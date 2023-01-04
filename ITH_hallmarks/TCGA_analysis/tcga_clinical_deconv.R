library(data.table)
library(magrittr)
library(ggplot2)
library(stringdist)
library(stringr)
library(survival)
library(readxl)
library(plyr)
library(matkot)

cancer_types <- c('BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'ESCA_AC', 'ESCA_ESCC', 'HNSC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PRAD', 'READ', 'STAD',
    'THCA', 'UCEC')

study_contrib <- readRDS('../data/study_contribution_per_MP.RDS')
names(study_contrib) <- gsub('  ', ' ', names(study_contrib))
names(study_contrib)[41] <- 'MP41 Unassigned'

study_tcga_map <- fread('../data/study_tcga_map.csv', na.strings = '')
cancer_types <- cancer_types[cancer_types %in% study_tcga_map$cancer_type]
setkey(study_tcga_map, cancer_type)

# Data from Liu et al. 2018 (https://doi.org/10.1016/j.cell.2018.02.052):
liu <- as.data.table(
    read_xlsx(
        '~/TCGA_data/Liu_1-s2.0-S0092867418302290-mmc1.xlsx',
        sheet = 1,
        na = c('[Unknown]', '[Not Evaluated]', '[Not Applicable]', '[Not Available]', '[Discrepancy]', '#N/A')
    )
)[, .(patient_id = gsub('-', '\\.', bcr_patient_barcode), cancer_type = type, OS = OS.time, OSc = OS, PFI = PFI.time, PFIc = PFI,
    stage = ajcc_pathologic_tumor_stage, grade = histological_grade, ther1 = str_to_title(treatment_outcome_first_course))]
setkey(liu, patient_id)

firehose <- fread('~/TCGA_data/clinical.csv', na.strings = '', key = 'patient_id')[, .(patient_id = patient_id,
    cancer_type = cancer_type, ln_n = as.character(number_of_lymphnodes_positive_by_he), ther2 = str_to_title(primary_therapy_outcome_success),
    ther3 = str_to_title(followup_treatment_success), t_stage = pathologic_t, n_stage = pathologic_n, m_stage = pathologic_m)]

clin <- merge(liu, firehose, by = c('patient_id', 'cancer_type'), all = TRUE)

meta_all <- lapply(
    cancer_types,
    function(ct) fread(paste0('~/TCGA_data/', ct, '/Cells.csv'))[, .(patient_id = patient_id, cancer_type = ct)]
) %>% rbindlist %>% unique
meta_all <- meta_all[patient_id %in% names(table(patient_id))[table(patient_id) == 1]]
setkey(meta_all, patient_id)
clin <- clin[patient_id %in% meta_all$patient_id]
clin[, cancer_type := do.call(`[`, list(meta_all, patient_id))$cancer_type] # Convert from TCGA label where relevant, e.g. 'ESCA' --> 'ESCA_AC'

clin_test <- lapply(cancer_types, function(ct) {
    
    cat(ct, '\n')
    
    if('Scores_deconv.csv' %in% dir(paste0('~/TCGA_data/', ct))) {
        scores <- fread(paste0('~/TCGA_data/', ct, '/Scores_deconv.csv'), key = 'sample_id')
    } else return(NULL)
    
    # Read in tumour metadata file and extract one sample ID per patient ID, according to order of priority of sample types:
    meta <- fread(paste0('~/TCGA_data/', ct, '/Cells.csv'), key = 'sample_type')[
        sample_type != 'normal',
        .SD[c('primary', 'primary_additional', 'recurrent', 'metastatic', 'metastatic_additional')][!is.na(sample_id)][1],
        keyby = patient_id
    ][clin[cancer_type == ct, patient_id]][!is.na(sample_id)]
    setkey(meta, sample_id)
    
    # Subset scores table with the above IDs:
    scores <- scores[sample_id %in% meta$sample_id]
    scores <- scores[, c(do.call(`[`, list(meta, sample_id))[, .(patient_id, purity)], .(meta_program = meta_program, score = score))]
    
    # Subset clinical data table with the above IDs:
    clin_ct <- clin[patient_id %in% scores$patient_id]
    
    # Define variable categories:
    
    # M stage:
    clin_ct[, m_stage := ifelse(grepl('m0', m_stage), 'Low', ifelse(grepl('m1', m_stage), 'High', NA))]
    
    # Therapy resistance:
    for(v in c('ther1', 'ther2', 'ther3')) clin_ct[, (v) := ifelse(grepl('^C|^No ', get(v)), 'Low', ifelse(is.na(get(v)), NA, 'High'))]
    clin_ct[, ther1 := {if(is.na(ther1)) {if(is.na(ther2)) NA else ther2} else if(is.na(ther2)) ther1 else if(ther1 == ther2) ther1 else NA},
        by = patient_id]
    clin_ct[, ther := 'a']
    clin_ct[, ther := if(any(c(ther1, ther3) == 'High', na.rm = TRUE)) 'High' else if(all(is.na(c(ther1, ther3)))) NA else 'Low', by = patient_id]
    clin_ct[, c('ther1', 'ther2', 'ther3') := NULL]
    
    # Lymphatic spread:
    clin_ct[!is.na(ln_n), ln_n := ifelse(ln_n == '0', 'Low', 'High')]
    clin_ct[, n_stage := ifelse(grepl('n0', n_stage), 'Low', ifelse(grepl('n1|n2|n3', n_stage), 'High', NA))]
    clin_ct[, ln := 'a']
    clin_ct[, ln := if(any(c(ln_n, n_stage) == 'High', na.rm = TRUE)) 'High' else if(all(is.na(c(ln_n, n_stage)))) NA else 'Low', by = patient_id]
    clin_ct[, c('ln_n', 'n_stage') := NULL]
    
    # In cases where we have intermediate levels, define boundary as wherever the cumulative sum of sample numbers passes 10:
    
    ints <- list()
    
    # Grade:
    clin_ct[, grade := mapvalues(grade, c('G1', 'G4', 'GX', 'GB', 'High Grade', 'Low Grade'), c('Low', 'High', NA, 'Low', 'High', 'Low'),
        warn_missing = FALSE)]
    ints$grade <- c('G2', 'G3')[c('G2', 'G3') %in% clin_ct$grade]
    
    # Stage:
    clin_ct[, stage := mapvalues(gsub('[ABC]$', '', stage), c('Stage 0', 'Stage I', 'Stage IV', 'IS', 'I/II NOS', 'Stage X'),
        c('Low', 'Low', 'High', 'Low', NA, NA), warn_missing = FALSE)]
    ints$stage <- c('Stage II', 'Stage III')[c('Stage II', 'Stage III') %in% clin_ct$stage]
    
    # T stage:
    clin_ct[, t_stage := mapvalues(str_extract(t_stage, 't[0-4]|tis'), c('tis', 't0', 't1', 't4'), c('Low', 'Low', 'Low', 'High'),
        warn_missing = FALSE)]
    ints$t_stage <- c('t2', 't3')[c('t2', 't3') %in% clin_ct$t_stage]
    
    for(v in names(ints)) {
        if(length(ints[[v]]) > 0) {
            sums <- clin_ct[!is.na(get(v)), .(N = .N), keyby = .(kv = get(v))][c('Low', ints[[v]], 'High')][!is.na(N), setNames(cumsum(N), kv)]
            if(!any(sums[names(sums) != 'High'] >= 10)) {
                clin_ct[, (v) := NA]
            } else {
                clin_ct[get(v) %in% names(sums)[1:min(which(sums >= 10))], (v) := 'Low']
                clin_ct[get(v) %in% names(sums)[(min(which(sums >= 10)) + 1):length(sums)], (v) := 'High']
            }
        }
    }
    
    scores <- scores[, cbind(.SD, do.call(`[`, list(clin_ct, patient_id))[, -c('patient_id', 'cancer_type')])]
    scores[, cancer_type := ct]
    setcolorder(scores, c('patient_id', 'cancer_type'))
    
    surv <- lapply(unique(scores$meta_program), function(mp) {
        lapply(c('OS', 'PFI'), function(v) {
            sdata <- scores[meta_program == mp & !is.na(get(v))]
            if(nrow(sdata) < 20) return(NULL)
            coxmod <- coxph(Surv(get(v), get(paste0(v, 'c'))) ~ score, data = sdata)
            data.table(cancer_type = ct, var_name = v, meta_program = mp, eff = coxmod$coefficients,# hr = exp(coxmod$coefficients),
                pval = summary(coxmod)$logtest['pvalue'])
        }) %>% rbindlist
    }) %>% rbindlist
    
    categ <- lapply(unique(scores$meta_program), function(mp) {
        lapply(c('stage', 'grade', 't_stage', 'm_stage', 'ln', 'ther'), function(v) {
            sdata <- scores[meta_program == mp & !is.na(get(v))]
            if(nrow(sdata) < 20) return(NULL) else if(sdata[, sum(get(v) == 'Low') < 10 | sum(get(v) == 'High') < 10]) return(NULL)
            data.table(cancer_type = ct, var_name = v, meta_program = mp,
                eff = sdata[, .SD[get(v) == 'High', mean(score)] - .SD[get(v) == 'Low', mean(score)]],
                pval = t.test(score ~ get(v), sdata)$p.value)
        }) %>% rbindlist
    }) %>% rbindlist
    
    return(list(res = rbind(surv, categ), scores = scores, clin = clin_ct))
    
})

clin_test_res <- rbindlist(lapply(clin_test, `[[`, 'res'))
clin_test_data <- rbindlist(lapply(clin_test, `[[`, 'clin'))
clin_test_scores <- rbindlist(lapply(clin_test, `[[`, 'scores'))

clin_test_res[, pval_adj := p.adjust(pval, method = 'BH')]
setkey(clin_test_res, cancer_type, meta_program)

fwrite(clin_test_res, '../data/clin_test_res_deconv.csv')
fwrite(clin_test_data, '../data/clin_test_data_deconv.csv')
fwrite(clin_test_scores, '../data/clin_test_scores_deconv.csv')
