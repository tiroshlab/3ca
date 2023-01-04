library(data.table)
library(readxl)
library(magrittr)
library(stringr)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(cowplot)
library(randomcoloR)
library(matkot)

cancer_types <- c('ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA_AC', 'ESCA_ESCC', 'GBM_IDH-WT', 'HNSC', 'KICH', 'KIRC', 'KIRP',
    'LAML', 'LGG_astro', 'LGG_IDH-WT', 'LGG_oligo', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM_primary',
    'SKCM_metastatic', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM')

cancer_genes <- as.data.table(read_xlsx('~/pan_cancer/cancer5000.xlsx', skip = 1, n_max = 260))[`Cancer5000-S (219)` == 1, gene]

hgnc_complete_set <- fread('../data/hgnc_complete_set.txt', key = 'symbol')[
    !(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])
][!grepl('unplaced', location), chrarm := str_extract(location, '[0-9]*[XY]*[pq]')][!is.na(chrarm)][
    order(as.numeric(mapvalues(gsub('[pq]', '', chrarm), c('X', 'Y'), c(23, 24))), str_extract(chrarm, '[pq]'))
]





# Identification of hypermutated subtypes in COAD, READ, STAD and UCEC:

# n_muts <- lapply(cancer_types, function(ct) {
#     if(!('Mutation_data.csv' %in% dir(paste0('~/TCGA_data/', ct)))) return(NULL)
#     mut_data <- fread(paste0('~/TCGA_data/', ct, '/Mutation_data.csv'))
#     mut_data[, .(cancer_type = ct, n_mut = .N), by = patient_id]
# }) %>% rbindlist
# n_muts <- n_muts[!(patient_id %in% names(table(patient_id))[table(patient_id) > 1])]
# n_muts[, patient_id := factor(patient_id, levels = patient_id[order(n_mut)])]
# n_muts[, subtype := NA][, subtype := as.character(subtype)]
# for(ct in c('COAD', 'READ', 'STAD', 'UCEC')) {
#     subtypes <- fread(paste0('~/TCGA_data/', ct, '/subtypes.csv'), na.strings = '', key = 'patient_id')
#     n_muts[cancer_type == ct, subtype := do.call(`[`, list(subtypes, as.character(patient_id)))$subtype]
# }
# 
# # A threshold of 500 mutations works well for COAD, READ and STAD, while 400 is better for UCEC.
# ggplot(n_muts[cancer_type == 'COAD'], aes(x = patient_id, y = n_mut, colour = subtype)) + geom_point() + scale_y_continuous(limits = c(0, 5000)) +
#     theme_test() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + geom_hline(yintercept = 500)
# ggplot(n_muts[cancer_type == 'READ'], aes(x = patient_id, y = n_mut, colour = subtype)) + geom_point() + scale_y_continuous(limits = c(0, 5000)) +
#     theme_test() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + geom_hline(yintercept = 500)
# ggplot(n_muts[cancer_type == 'STAD'], aes(x = patient_id, y = n_mut, colour = subtype)) + geom_point() + scale_y_continuous(limits = c(0, 5000)) +
#     theme_test() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + geom_hline(yintercept = 500)
# ggplot(n_muts[cancer_type == 'UCEC'], aes(x = patient_id, y = n_mut, colour = subtype)) + geom_point() + scale_y_continuous(limits = c(0, 5000)) +
#     theme_test() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + geom_hline(yintercept = 400)





# Mutations for cancer types without hypermutated subtypes:
test_mutations <- lapply(cancer_types[!(cancer_types %in% c('COAD', 'READ', 'STAD', 'UCEC'))], function(ct) {
    
    cat(ct, '\n')
    
    if(!('Scores.csv' %in% dir(paste0('~/TCGA_data/', ct)) & 'Mutation_data.csv' %in% dir(paste0('~/TCGA_data/', ct)))) return(NULL)
    
    meta <- fread(paste0('~/TCGA_data/', ct, '/Cells.csv'), key = 'sample_id')
    scores <- fread(paste0('~/TCGA_data/', ct, '/Scores.csv'))
    mut_data <- fread(paste0('~/TCGA_data/', ct, '/Mutation_data.csv'))
    
    mut_data <- mut_data[gene %in% cancer_genes]
    setkey(mut_data, gene, Variant_Classification)
    mut_data <- rbind(
        mut_data[Variant_Classification %in% c('Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site', 'Nonstop_Mutation',
            'Translation_Start_Site')],
        mut_data[
            mut_data[
                Variant_Classification %in% c('Missense_Mutation', 'In_Frame_Del', 'In_Frame_Ins'),
                .(n = length(unique(patient_id))),
                by = .(gene, Variant_Classification)
            ][n >= 2, -'n']
        ]
    )
    
    # Make a version of the scores data having only patient IDs, for easier searching:
    pscores <- copy(scores)[, patient_id := do.call(`[`, args = list(meta, sample_id))$patient_id][, -'sample_id']
    setcolorder(pscores, 'patient_id')
    setkey(pscores, patient_id)
    
    # Make patient IDs match between pscores and mut_data, removing any patient IDs which occur more than once in pscores:
    pscores <- pscores[patient_id %in% mut_data$patient_id]
    mut_data <- mut_data[patient_id %in% pscores$patient_id]
    pids_tab <- pscores[meta_program == meta_program[1], table(patient_id)]
    pids <- names(pids_tab)[pids_tab == 1]
    pscores <- pscores[patient_id %in% pids]
    mut_data <- mut_data[patient_id %in% pids]
    
    if(ct %in% c('BRCA', 'THCA')) {
        subtype_var <- switch(ct, BRCA = 'subtype_pam50', THCA = 'molecular_subtype')
        subtypes <- fread(paste0('~/TCGA_data/', ct, '/subtypes.csv'), na.strings = '')
        subtypes <- subtypes[!is.na(get(subtype_var)) & !is.na(patient_id) & patient_id %in% pids]
        mut_data <- mut_data[patient_id %in% subtypes$patient_id]
        ids_list <- slapply(unique(subtypes[[subtype_var]]), function(sbtp) subtypes[get(subtype_var) == sbtp, patient_id])
        names(ids_list) <- paste0(ct, '_', gsub(' ', '-', mapvalues(names(ids_list), c('positive', 'negative'), c('HPV-pos', 'HPV-neg'))))
    } else {
        ids_list <- setNames(list(pids), ct)
    }
    
    out <- lapply(names(ids_list), function(subtype) { # <subtype> could just be <ct>.
        
        ids_vec <- ids_list[[subtype]]
        
        # Get genes mutated in at least 10 samples and at least 10% of samples, and not mutated in at least 10 samples, then take top 20:
        genes <- mut_data[patient_id %in% ids_vec, .(N = length(unique(patient_id))), by = gene]
        genes <- genes[N >= 10 & N >= length(ids_vec)/10 & N <= length(ids_vec) - 10][order(-N)][1:min(.N, 20), gene]
        mut_data <- mut_data[gene %in% genes]
        
        if(nrow(mut_data) == 0) return(NULL)
        
        # t test for each gene/mutation type combination comparing scores between tumours with and without that mutation:
        mut_data[patient_id %in% ids_vec, {
            ids <- patient_id
            pscores[, categ := (patient_id %in% ids)]
            pscores[
                patient_id %in% ids_vec,
                .(cancer_type = subtype, event_type = 'mutation', n_event = sum(categ), n_patient = length(unique(patient_id)),
                    pval = t.test(score ~ categ)$p.value, diff = .SD[, .(m = mean(score)), keyby = categ][, diff(m)]),
                by = meta_program
            ]
        }, by = .(event = gene)]
        
    })
    
    if(all(sapply(out, is.null))) return(NULL) else return(rbindlist(out))
    
}) %>% rbindlist
setcolorder(test_mutations, c('cancer_type', 'event', 'event_type', 'n_event', 'n_patient'))

# Mutations for cancer types with hypermutated subtypes:
test_mutations_hm <- lapply(c('COAD', 'READ', 'STAD', 'UCEC'), function(ct) {
    
    meta <- fread(paste0('~/TCGA_data/', ct, '/Cells.csv'), key = 'sample_id')
    scores <- fread(paste0('~/TCGA_data/', ct, '/Scores.csv'))
    mut_data <- fread(paste0('~/TCGA_data/', ct, '/Mutation_data.csv'))
    
    # Count mutations:
    hm <- mut_data[, .(N = .N), keyby = patient_id][, hm := ifelse(N > ifelse(ct == 'UCEC', 400, 500), paste0(ct, '_MSI'), paste0(ct, '_non-MSI'))]
    mut_data[, hm := do.call(`[`, list(hm, patient_id))$hm]
    
    mut_data <- mut_data[gene %in% cancer_genes]
    setkey(mut_data, gene, Variant_Classification)
    mut_data <- rbind(
        mut_data[Variant_Classification %in% c('Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site', 'Nonstop_Mutation',
            'Translation_Start_Site')],
        mut_data[
            mut_data[
                Variant_Classification %in% c('Missense_Mutation', 'In_Frame_Del', 'In_Frame_Ins'),
                .(n = length(unique(patient_id))),
                by = .(gene, Variant_Classification)
            ][n >= 2, -'n']
        ]
    )
    
    # Make a version of the scores data having only patient IDs, for easier searching:
    pscores <- copy(scores)[, patient_id := do.call(`[`, args = list(meta, sample_id))$patient_id][, -'sample_id']
    setcolorder(pscores, 'patient_id')
    setkey(pscores, patient_id)
    
    # Make patient IDs match between pscores and mut_data, removing any patient IDs which occur more than once in pscores:
    pscores <- pscores[patient_id %in% mut_data$patient_id]
    mut_data <- mut_data[patient_id %in% pscores$patient_id]
    pids_tab <- pscores[meta_program == meta_program[1], table(patient_id)]
    pids <- names(pids_tab)[pids_tab == 1]
    pscores <- pscores[patient_id %in% pids]
    mut_data <- mut_data[patient_id %in% pids]
    hm <- hm[patient_id %in% pids]
    
    # Get genes mutated in at least 10 samples and at least 10% of samples, and not mutated in at least 10 samples, then take top 20:
    genes <- mut_data[,
        .(N = length(unique(patient_id))),
        by = .(msi = hm, gene = gene)
    ][, .(g = .SD[N >= 10 & N >= length(pids)/10 & N <= length(pids) - 10][order(-N)][1:min(.N, 20), gene]), by = msi]
    mut_data <- mut_data[, .SD[gene %in% genes[msi == hm, g]], by = hm]
    
    if(nrow(mut_data) == 0) return(NULL)
    
    # t test for each gene/mutation type combination comparing scores between tumours with and without that mutation:
    
    setnames(mut_data, 'hm', 'msi') # Avoids conflict between hm the variable and hm the data table
    mut_data[, {
        ids <- patient_id
        hm_ids <- do.call(function(x) hm[hm == x, patient_id], list(cancer_type))
        pscores[hm_ids, categ := (patient_id %in% ids)]
        pscores[
            hm_ids,
            .(event_type = 'mutation', n_event = sum(categ), n_patient = length(hm_ids),
                pval = t.test(score ~ categ)$p.value, diff = .SD[, .(m = mean(score)), keyby = categ][, diff(m)]),
            by = meta_program
        ]
    }, by = .(cancer_type = msi, event = gene)]
    
}) %>% rbindlist
setcolorder(test_mutations_hm, c('cancer_type', 'event', 'event_type', 'n_event', 'n_patient'))

# CNAs for cancer types without hypermutated subtypes:
test_cna <- lapply(cancer_types[!(cancer_types %in% c('COAD', 'READ', 'STAD', 'UCEC'))], function(ct) {
    
    cat(ct, '\n')
    
    if(!('Scores.csv' %in% dir(paste0('~/TCGA_data/', ct)) & 'CNV_GISTIC_data.csv' %in% dir(paste0('~/TCGA_data/', ct)))) return(NULL)
    meta <- fread(paste0('~/TCGA_data/', ct, '/Cells.csv'), key = 'sample_id')
    scores <- fread(paste0('~/TCGA_data/', ct, '/Scores.csv'))
    cna_data <- fread(paste0('~/TCGA_data/', ct, '/CNV_GISTIC_data.csv'))[, set_rownames(as.matrix(.SD), V1), .SDcols = -'V1']
    
    # Change column names to patient IDs (because they are initially sample IDs which do not match those in the expression data):
    pids <- apply(str_split_fixed(colnames(cna_data), '\\.', 4)[, 1:3], 1, paste, collapse = '.')
    colnames(cna_data) <- pids
    pids_tab <- table(pids)
    cna_data <- cna_data[, pids %in% meta[sample_id %in% scores$sample_id, patient_id] & !(pids %in% names(pids_tab)[pids_tab > 1])]
    
    # Make a version of the scores data having only patient IDs, for easier searching:
    pscores <- copy(scores)[, patient_id := do.call(`[`, args = list(meta, sample_id))$patient_id][patient_id %in% colnames(cna_data), -'sample_id']
    setcolorder(pscores, 'patient_id')
    setkey(pscores, patient_id)
    
    # Remove repeated patient IDs in pscores and make sure pscores and cna_data have the same patient IDs:
    pids_final_tab <- pscores[meta_program == meta_program[1], table(patient_id)]
    pids_final <- names(pids_final_tab)[pids_final_tab == 1]
    pscores <- pscores[patient_id %in% pids_final]
    cna_data <- cna_data[, colnames(cna_data) %in% pids_final]
    
    cna_data_sub <- cna_data[rownames(cna_data) %in% cancer_genes, ]
    
    if(ct %in% c('BRCA', 'THCA')) {
        subtype_var <- switch(ct, BRCA = 'subtype_pam50', THCA = 'molecular_subtype')
        subtypes <- fread(paste0('~/TCGA_data/', ct, '/subtypes.csv'), na.strings = '')
        subtypes <- subtypes[!is.na(get(subtype_var)) & !is.na(patient_id) & patient_id %in% pids]
        ids_list <- slapply(unique(subtypes[[subtype_var]]), function(sbtp) subtypes[get(subtype_var) == sbtp, patient_id])
        names(ids_list) <- paste0(ct, '_', gsub(' ', '-', mapvalues(names(ids_list), c('positive', 'negative'), c('HPV-pos', 'HPV-neg'))))
    } else {
        ids_list <- setNames(list(pids), ct)
    }
    
    test_gene <- lapply(names(ids_list), function(subtype) { # <subtype> could just be <ct>.
        
        ids_vec <- ids_list[[subtype]][ids_list[[subtype]] %in% colnames(cna_data)]
        
        # Get genes having at least 10 samples and at least 10% of samples with +/-2 CNA scores, and at least 10 samples without, then take top 20:
        genes <- lapply(
            rownames(cna_data_sub),
            function(g) list(gene = g, amp = sum(cna_data_sub[g, ids_vec] == 2), del = sum(cna_data_sub[g, ids_vec] == -2))
        ) %>% rbindlist
        genes <- melt(genes, id.var = 'gene', variable.name = 'type', value.name = 'n', variable.factor = FALSE)
        genes[, keep := ifelse(n >= 10 & n >= length(ids_vec)/10 & length(ids_vec) - n >= 10, TRUE, FALSE)]
        genes <- genes[keep == TRUE, .(gene = .SD[order(-n)][1:min(.N, 20), gene]), by = type]
        
        if(genes[, .(N = .N), keyby = type][c('amp', 'del'), all(is.na(N))]) {
            out <- NULL
        } else {
            # Run t tests for each of the above genes comparing MP scores of samples with and without the corresponding CNVs:
            out <- genes[, {
                gvec <- cna_data_sub[gene, ids_vec]
                pscores[
                    ids_vec,
                    categ := switch(type, amp = (patient_id %in% names(gvec)[gvec == 2]), del = (patient_id %in% names(gvec)[gvec == -2]))
                ]
                pscores[
                    ids_vec,
                    .(cancer_type = subtype, event_type = 'focal CNA', n_event = sum(categ), n_patient = length(ids_vec),
                        pval = t.test(score ~ categ)$p.value, diff = .SD[, .(m = mean(score)), keyby = categ][, diff(m)]),
                    by = meta_program
                ]
            }, by = .(type, gene)]
            out[, c('event', 'gene', 'type') := .(paste(gene, type), NULL, NULL)]
            setcolorder(out, c('cancer_type', 'event', 'event_type', 'n_event', 'n_patient'))
        }
        return(out)
        
    }) %>% rbindlist
    
    test_arm <- lapply(names(ids_list), function(subtype) { # <subtype> could just be <ct>.
        
        ids_vec <- ids_list[[subtype]][ids_list[[subtype]] %in% colnames(cna_data)]
        
        events <- hgnc_complete_set[
            locus_group == 'protein-coding gene' & symbol %in% rownames(cna_data), # We may lose some arms with this subsetting, e.g. chrY (p & q).
            {cm <- colMeans(cna_data[symbol, ids_vec, drop = FALSE]); data.table(patient_id = names(cm), ave = cm, n_gene = .N)},
            by = chrarm
        ][, event := ifelse(n_gene >= 100, ifelse(ave > 0.9, 'amp', ifelse(ave < -0.9, 'del', 'none')), as.character(NA))] # abs > 0.9 defines event
        
        events <- events[!is.na(event)][
            order(
                as.numeric(mapvalues(gsub('[pq]', '', chrarm), c('X', 'Y'), c(23, 24), warn_missing = FALSE)),
                str_extract(chrarm, '[pq]'),
                patient_id
            )
        ]
        
        out <- lapply(c('amp', 'del'), function(et) events[,
            if(sum(event == et) >= 10 & sum(event == et) >= length(ids_vec)/10 & sum(event != et) >= 10) {
                pscores[
                    ids_vec,
                    .(cancer_type = subtype, type = et, event_type = 'arm CNA', n_event = sum(event == et), n_patient = length(ids_vec),
                        pval = t.test(score ~ (event == et))$p.value, diff = .SD[, .(m = mean(score)), keyby = (event == et)][, diff(m)]),
                    by = meta_program
                ]
            } else NULL,
            by = chrarm
        ])
        out <- rbindlist(out[sapply(out, function(x) nrow(x) > 0)])
        out[, c('event', 'chrarm', 'type') := .(paste(chrarm, type), NULL, NULL)]
        setcolorder(out, c('cancer_type', 'event', 'event_type', 'n_event', 'n_patient'))
        return(out)
        
    }) %>% rbindlist
    
    return(rbind(test_gene, test_arm))
    
}) %>% rbindlist

# CNAs for cancer types with hypermutated subtypes:
test_cna_hm <- lapply(c('COAD', 'READ', 'STAD', 'UCEC'), function(ct) {
    
    cat(ct, '\n')
    
    meta <- fread(paste0('~/TCGA_data/', ct, '/Cells.csv'), key = 'sample_id')
    scores <- fread(paste0('~/TCGA_data/', ct, '/Scores.csv'))
    cna_data <- fread(paste0('~/TCGA_data/', ct, '/CNV_GISTIC_data.csv'))[, set_rownames(as.matrix(.SD), V1), .SDcols = -'V1']
    mut_data <- fread(paste0('~/TCGA_data/', ct, '/Mutation_data.csv'))
    
    # Count mutations:
    hm <- mut_data[, .(N = .N), keyby = patient_id][, hm := ifelse(N > ifelse(ct == 'UCEC', 400, 500), paste0(ct, '_MSI'), paste0(ct, '_non-MSI'))]
    
    # Change column names to patient IDs (because they are initially sample IDs which do not match those in the expression data):
    pids <- apply(str_split_fixed(colnames(cna_data), '\\.', 4)[, 1:3], 1, paste, collapse = '.')
    colnames(cna_data) <- pids
    pids_tab <- table(pids)
    cna_data <- cna_data[, pids %in% meta[sample_id %in% scores$sample_id, patient_id] & !(pids %in% names(pids_tab)[pids_tab > 1])]
    
    # Make a version of the scores data having only patient IDs, for easier searching:
    pscores <- copy(scores)[, patient_id := do.call(`[`, args = list(meta, sample_id))$patient_id][patient_id %in% colnames(cna_data), -'sample_id']
    setcolorder(pscores, 'patient_id')
    setkey(pscores, patient_id)
    
    # Remove repeated patient IDs in pscores and make sure pscores and cna_data have the same patient IDs:
    pids_final_tab <- pscores[meta_program == meta_program[1], table(patient_id)]
    pids_final <- names(pids_final_tab)[pids_final_tab == 1]
    pscores <- pscores[patient_id %in% pids_final]
    cna_data <- cna_data[, colnames(cna_data) %in% pids_final]
    
    # Subset patient IDs that have MSI annotations:
    cna_data <- cna_data[, colnames(cna_data) %in% hm$patient_id]
    hm <- hm[patient_id %in% colnames(cna_data)]
    pscores <- pscores[patient_id %in% colnames(cna_data)]
    
    cna_data_sub <- cna_data[rownames(cna_data) %in% cancer_genes, ]
    
    # Get genes having at least 10 samples and at least 10% of samples with +/-2 CNA scores, and at least 10 samples without, then take top 20:
    genes <- hm[, {
        genes_hm <- lapply(
            rownames(cna_data_sub),
            function(g) list(gene = g, amp = sum(cna_data_sub[g, patient_id] == 2), del = sum(cna_data_sub[g, patient_id] == -2))
        ) %>% rbindlist
        genes_hm <- melt(genes_hm, id.var = 'gene', variable.name = 'type', value.name = 'n', variable.factor = FALSE)
        genes_hm[, keep := ifelse(n >= 10 & n >= ncol(cna_data_sub)/10 & ncol(cna_data_sub) - n >= 10, TRUE, FALSE)]
        genes_hm[keep == TRUE, .(gene = .SD[order(-n)][1:min(.N, 20), gene]), by = type]
    }, by = .(msi = hm)]
    
    if(
        genes[, .(N = .N), keyby = .(msi, type)][CJ(c('COAD_MSI', 'COAD_non-MSI'), c('amp', 'del'))][
            !is.na(msi), # For some reason you get NA in the final msi column unless you add this !is.na step
            .(all_na = all(is.na(N))),
            by = msi
        ][, all(all_na)]
    ) {
        test_gene <- NULL
    } else {
        # Run t tests for each of the above genes comparing MP scores of samples with and without the corresponding CNVs:
        test_gene <- genes[, {
            ids <- do.call(function(x) hm[hm == x, patient_id], list(cancer_type))
            gvec <- cna_data_sub[gene, ]
            pscores[ids, categ := switch(type, amp = (patient_id %in% names(gvec)[gvec == 2]), del = (patient_id %in% names(gvec)[gvec == -2]))]
            pscores[
                ids,
                .(event_type = 'focal CNA', n_event = sum(categ), n_patient = length(ids),
                    pval = t.test(score ~ categ)$p.value, diff = .SD[, .(m = mean(score)), keyby = categ][, diff(m)]),
                by = meta_program
            ]
        }, by = .(cancer_type = msi, type = type, gene = gene)]
        test_gene[, c('event', 'gene', 'type') := .(paste(gene, type), NULL, NULL)]
        setcolorder(test_gene, c('cancer_type', 'event', 'event_type', 'n_event', 'n_patient'))
    }
    
    events <- hgnc_complete_set[
        locus_group == 'protein-coding gene' & symbol %in% rownames(cna_data), # We may lose some arms with this subsetting, e.g. chrY (p & q).
        {cm <- colMeans(cna_data[symbol, , drop = FALSE]); data.table(patient_id = names(cm), ave = cm, n_gene = .N)},
        by = chrarm
    ][, event := ifelse(n_gene >= 100, ifelse(ave > 0.9, 'amp', ifelse(ave < -0.9, 'del', 'none')), as.character(NA))] # abs > 0.98 defines event
    
    events <- events[!is.na(event)][
        order(
            as.numeric(mapvalues(gsub('[pq]', '', chrarm), c('X', 'Y'), c(23, 24), warn_missing = FALSE)),
            str_extract(chrarm, '[pq]'),
            patient_id
        )
    ][, msi := do.call(`[`, list(hm, patient_id))$hm]
    
    test_arm <- lapply(c('amp', 'del'), function(et) events[,
        if(sum(event == et) >= 10 & sum(event == et) >= .N/10 & sum(event != et) >= 10) { # .N = number of unique patient IDs in this MSI subtype
            do.call(`[`, list(pscores, patient_id))[,
                .(type = et, event_type = 'arm CNA', n_event = sum(event == et), n_patient = .N, # Same here, since it's per meta-program
                    pval = t.test(score ~ (event == et))$p.value, diff = .SD[, .(m = mean(score)), keyby = (event == et)][, diff(m)]),
                by = meta_program
            ]
        } else NULL,
        by = .(cancer_type = msi, chrarm = chrarm)
    ])
    test_arm <- rbindlist(test_arm[sapply(test_arm, function(x) nrow(x) > 0)])
    test_arm[, c('event', 'chrarm', 'type') := .(paste(chrarm, type), NULL, NULL)]
    setcolorder(test_arm, c('cancer_type', 'event', 'event_type', 'n_event', 'n_patient'))
    
    return(rbind(test_gene, test_arm))
    
}) %>% rbindlist

test_all <- rbind(test_mutations, test_mutations_hm, test_cna, test_cna_hm)

# Adjust p values using threshold of 0.4 (though for whole-arm CNAs we'll show only those passing a threshold of 0.6 in the figure):
test_all[abs(diff) > 0.4, pval_adj := p.adjust(pval, method = 'BH'), by = event_type]





# Summary plots:

pdata <- test_all[pval_adj < 0.05 & abs(diff) > ifelse(event_type == 'arm CNA', 0.6, 0.4)]

pdata1_pos <- pdata[event_type != 'arm CNA' & diff > 0, .(f = switch((.N > 1) + 1, cancer_type, 'Multiple')), by = .(event, meta_program)]
event_ord_pos <- pdata1_pos[, .(N = .N, ev = gsub(' amp| del', '', event)), by = event]
event_ord_pos[, N_gene := do.call(`[`, list(event_ord_pos[, .(N = mean(N)), keyby = ev], ev))$N]
pdata1_pos[, event := factor(event, levels = event_ord_pos[order(-N_gene, event), event])]
pdata1_pos[, meta_program := factor(meta_program, levels = .SD[, .(N = .N), by = meta_program][order(-N), meta_program])]
setkey(pdata1_pos, event, meta_program)

pdata1_neg <- pdata[event_type != 'arm CNA' & diff < 0, .(f = switch((.N > 1) + 1, cancer_type, 'Multiple')), by = .(event, meta_program)]
event_ord_neg <- pdata1_neg[, .(N = .N, ev = gsub(' amp| del', '', event)), by = event]
event_ord_neg[, N_gene := do.call(`[`, list(event_ord_neg[, .(N = mean(N)), keyby = ev], ev))$N]
pdata1_neg[, event := factor(event, levels = event_ord_neg[order(-N_gene, event), event])]
pdata1_neg[, meta_program := factor(meta_program, levels = .SD[, .(N = .N), by = meta_program][order(-N), meta_program])]
setkey(pdata1_neg, event, meta_program)

pdata_arm_pos <- pdata[event_type == 'arm CNA' & diff > 0, .(f = switch((.N > 1) + 1, cancer_type, 'Multiple')), by = .(event, meta_program)]
pdata_arm_pos[, meta_program := factor(meta_program, levels = .SD[, .(N = .N), by = meta_program][order(-N), meta_program])]
setkey(pdata_arm_pos, event, meta_program)
ev_pos <- unique(pdata_arm_pos$event)
pdata_arm_pos <- pdata_arm_pos[CJ(ev_pos, unique(meta_program))]
pdata_arm_pos[, event := factor(event, levels = ev_pos[order(as.numeric(str_extract(ev_pos, '[0-9]+')), str_extract(ev_pos, 'amp|del'))])]

pdata_arm_neg <- pdata[event_type == 'arm CNA' & diff < 0, .(f = switch((.N > 1) + 1, cancer_type, 'Multiple')), by = .(event, meta_program)]
pdata_arm_neg[, meta_program := factor(meta_program, levels = .SD[, .(N = .N), by = meta_program][order(-N), meta_program])]
setkey(pdata_arm_neg, event, meta_program)
ev_neg <- unique(pdata_arm_neg$event)
pdata_arm_neg <- pdata_arm_neg[CJ(ev_neg, unique(meta_program))]
pdata_arm_neg[, event := factor(event, levels = ev_neg[order(as.numeric(str_extract(ev_neg, '[0-9]+')), str_extract(ev_neg, 'amp|del'))])]

ldata <- rbind(pdata1_pos, pdata1_neg, pdata_arm_pos, pdata_arm_neg)[!is.na(f), .(categ = unique(f))]
set.seed(3559)
ldata[, col := distinctColorPalette(.N)]
ldata[, c('r', 'g', 'b') := lapply(c(1, 3, 5), function(i) substr(gsub('#', '', col), i, i+1))]
ldata[, c('r', 'g', 'b') := lapply(.SD, function(x) {
    digits <- mapvalues(str_split_fixed(x, '', 2), LETTERS[1:6], 10:15)
    as.numeric(digits[, 1])*16 + as.numeric(digits[, 2])
}), .SDcols = c('r', 'g', 'b')]
ldata[, dist_red_sq := (255 - r)^2 + g^2 + b^2]
col_mult <- ldata[categ == 'Multiple', col]
col_min <- ldata[dist_red_sq == min(dist_red_sq), col]
ldata[categ == 'Multiple', col := col_min]
ldata[categ == ldata[dist_red_sq == min(dist_red_sq), categ], col := col_mult]
setkey(ldata, categ)
colvec <- ldata[c('Multiple', sort(categ[categ != 'Multiple'])), setNames(col, categ)]

summary1_pos <- ggplot(pdata1_pos[CJ(unique(event), unique(meta_program))], aes(x = meta_program, y = event, fill = f)) +
    geom_tile(colour = 'grey90') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(na.value = 'white', values = colvec) +
    theme(panel.border = element_rect(size = 0.5, fill = NA), axis.text.x = element_text(angle = 55, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = 'Cancer type', title = 'Mutations and focal CNAs', subtitle = 'Positive effects')

summary1_neg <- ggplot(pdata1_neg[CJ(unique(event), unique(meta_program))], aes(x = meta_program, y = event, fill = f)) +
    geom_tile(colour = 'grey90') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(na.value = 'white', values = colvec) +
    theme(panel.border = element_rect(size = 0.5, fill = NA), axis.text.x = element_text(angle = 55, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = 'Cancer type', title = 'Mutations and focal CNAs', subtitle = 'Negative effects')

summary_arm_pos <- ggplot(pdata_arm_pos, aes(x = meta_program, y = event, fill = f)) +
    geom_tile(colour = 'grey90') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(na.value = 'white', values = colvec) +
    theme(panel.border = element_rect(size = 0.5, fill = NA), axis.text.x = element_text(angle = 55, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = 'Cancer type', title = 'Chromosome arms', subtitle = 'Positive effects')

summary_arm_neg <- ggplot(pdata_arm_neg, aes(x = meta_program, y = event, fill = f)) +
    geom_tile(colour = 'grey90') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(na.value = 'white', values = colvec) +
    theme(panel.border = element_rect(size = 0.5, fill = NA), axis.text.x = element_text(angle = 55, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = 'Cancer type', title = 'Chromosome arms', subtitle = 'Negative effects')

leg <- get_legend(ggplot(ldata) + geom_tile(aes(x = 0, y = 0, fill = categ)) + scale_fill_manual(name = 'Cancer type', values = colvec,
            guide = guide_legend(nrow = 5, ncol = 5, keywidth = unit(0.5, 'cm'), keyheight = unit(0.4, 'cm'))))
summary1_pos_grob <- ggplotGrob(summary1_pos + theme(legend.position = 'none'))
summary1_neg_grob <- ggplotGrob(summary1_neg + theme(legend.position = 'none'))
summary_arm_pos_grob <- ggplotGrob(summary_arm_pos + theme(legend.position = 'none'))
summary_arm_neg_grob <- ggplotGrob(summary_arm_neg + theme(legend.position = 'none'))

summary1_pos_widths <- c(0.2, 0, 0, 2.3, 0.5*pdata1_pos[, length(unique(meta_program))], 0, 0, 0, 0.2)
summary1_pos_heights <- c(0.2, 0, 0.5, 0.4, 0, 0, 0.3*pdata1_pos[, length(unique(event))], 3.5, 0, 0, 0, 0.2)
summary1_pos_grob$widths <- unit(summary1_pos_widths, 'cm')
summary1_pos_grob$heights <- unit(summary1_pos_heights, 'cm')

summary1_neg_widths <- c(0.2, 0, 0, 2.3, 0.5*pdata1_neg[, length(unique(meta_program))], 0, 0, 0, 0.2)
summary1_neg_heights <- c(0.2, 0, 0.5, 0.4, 0, 0, 0.3*pdata1_neg[, length(unique(event))], 3.5, 0, 0, 0, 0.2)
summary1_neg_grob$widths <- unit(summary1_neg_widths, 'cm')
summary1_neg_grob$heights <- unit(summary1_neg_heights, 'cm')

summary_arm_pos_widths <- c(1.2, 0, 0, 1.5, 0.5*pdata_arm_pos[, length(unique(meta_program))], 0, 0, 0, 0.2)
summary_arm_pos_heights <- c(0.2, 0, 0.5, 0.4, 0, 0, 0.3*pdata_arm_pos[, length(unique(event))], 3.5, 0, 0, 0, 0.2)
summary_arm_pos_grob$widths <- unit(summary_arm_pos_widths, 'cm')
summary_arm_pos_grob$heights <- unit(summary_arm_pos_heights, 'cm')

summary_arm_neg_widths <- c(0.2, 0, 0, 1.5, 0.5*pdata_arm_neg[, length(unique(meta_program))], 0, 0, 0, 0.2)
summary_arm_neg_heights <- c(0.2, 0, 0.5, 0.4, 0, 0, 0.3*pdata_arm_neg[, length(unique(event))], 3.5, 0, 0, 0, 0.2)
summary_arm_neg_grob$widths <- unit(summary_arm_neg_widths, 'cm')
summary_arm_neg_grob$heights <- unit(summary_arm_neg_heights, 'cm')

leg$heights <- unit(c(0.5, 0, 3, 0, 0), 'cm')

col1_height <- unit(sum(summary1_pos_heights), 'cm')
col2_height <- unit(max(sum(summary1_neg_heights), sum(summary_arm_pos_heights), sum(summary_arm_neg_heights)), 'cm')
leg_height <- sum(leg$heights)

summary1_pos_heights[12] <- col2_height + leg_height - unit(sum(summary1_pos_heights[1:11]), 'cm')
summary1_neg_heights[12] <- col2_height - unit(sum(summary1_neg_heights[1:11]), 'cm')
summary_arm_pos_heights[12] <- col2_height - unit(sum(summary_arm_pos_heights[1:11]), 'cm')
summary_arm_neg_heights[12] <- col2_height - unit(sum(summary_arm_neg_heights[1:11]), 'cm')

summary1_pos_grob$heights <- summary1_pos_heights
summary1_neg_grob$heights <- summary1_neg_heights
summary_arm_pos_grob$heights <- summary_arm_pos_heights
summary_arm_neg_grob$heights <- summary_arm_neg_heights

pdf(
    '../data/mut_summary.pdf',
    width = sum(summary1_pos_widths, summary1_neg_widths, summary_arm_pos_widths, summary_arm_neg_widths)/2.54,
    height = sum(summary1_pos_heights)/2.54
)
plot_grid(
    summary1_pos_grob,
    plot_grid(
        plot_grid(
            summary1_neg_grob,
            summary_arm_pos_grob,
            summary_arm_neg_grob,
            nrow = 1,
            ncol = 3,
            rel_widths = c(sum(summary1_neg_widths), sum(summary_arm_pos_widths), sum(summary_arm_neg_widths))
        ),
        leg,
        nrow = 2,
        ncol = 1,
        rel_heights = c(sum(summary1_neg_heights), leg_height)
    ),
    nrow = 1,
    ncol = 2,
    rel_widths = c(sum(summary1_pos_widths), sum(c(summary1_neg_widths, summary_arm_pos_widths, summary_arm_neg_widths)))
)
dev.off()

fwrite(
    test_all[pval_adj < 0.05, .(`Cancer type` = cancer_type, Event = event, `Event type` = event_type, `#patients with event` = n_event,
        `Total #patients` = n_patient, `Meta-program` = meta_program, `Difference in scores` = diff, `p value` = pval,
        `Adjusted p value` = pval_adj)],
    '../data/supp_table_mut.csv'
)





# CASP8 and Interferon/MHC-II in HNSC; also PDAC-classical and NSD1:

# TCGA HNSC paper suggests that CASP8 is mainly found in the basal (and maybe mesenchymal) subtype, so check if the association still holds when
# restricting to this (/these) subtype(s):
meta_hnsc <- fread('~/TCGA_data/HNSC/Cells.csv', key = 'sample_id', na.strings = '')
scores_hnsc <- fread('~/TCGA_data/HNSC/Scores.csv')[grep('MP17 Interferon|PDAC-classical', meta_program)]
mut_data_hnsc <- fread('~/TCGA_data/HNSC/Mutation_data.csv')
mut_data_hnsc <- mut_data_hnsc[gene %in% cancer_genes]
setkey(mut_data_hnsc, gene, Variant_Classification)
mut_data_hnsc <- rbind(
    mut_data_hnsc[Variant_Classification %in% c('Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site', 'Nonstop_Mutation',
        'Translation_Start_Site')],
    mut_data_hnsc[
        mut_data_hnsc[
            Variant_Classification %in% c('Missense_Mutation', 'In_Frame_Del', 'In_Frame_Ins'),
            .(n = length(unique(patient_id))),
            by = .(gene, Variant_Classification)
        ][n >= 2, -'n']
    ]
)
pscores_hnsc <- copy(scores_hnsc)[, patient_id := do.call(`[`, args = list(meta_hnsc, sample_id))$patient_id][, -'sample_id']
setcolorder(pscores_hnsc, 'patient_id')
setkey(pscores_hnsc, patient_id)
pscores_hnsc <- pscores_hnsc[patient_id %in% mut_data_hnsc$patient_id]
mut_data_hnsc <- mut_data_hnsc[patient_id %in% pscores_hnsc$patient_id]
pids_tab <- pscores_hnsc[meta_program == meta_program[1], table(patient_id)]
pids <- names(pids_tab)[pids_tab == 1]
pscores_hnsc <- pscores_hnsc[patient_id %in% pids]
mut_data_hnsc <- mut_data_hnsc[patient_id %in% pids]

pscores_hnsc <- dcast(pscores_hnsc[, .(patient_id, meta_program, score)], patient_id ~ meta_program, value.var = 'score')
genes <- c('CASP8', 'NSD1')
pscores_hnsc[, (genes) := lapply(genes, function(g) patient_id %in% mut_data_hnsc[gene == g, patient_id])]

# See how scores for this MP match with the trinity of CASP8-mut, HRAS-mut and TP53-WT:
pscores_hnsc[,
    categ3 := patient_id %in% mut_data_hnsc[gene == 'CASP8', patient_id] &
        patient_id %in% mut_data_hnsc[gene == 'HRAS', patient_id] &
        !(patient_id %in% mut_data_hnsc[gene == 'TP53', patient_id]),
    by = patient_id
]

bdata_hnsc_ifn <- melt(
    pscores_hnsc[, .(patient_id, `MP17 Interferon/MHC-II (I)`, CASP8, categ3)],
    id.vars = c('patient_id', 'MP17 Interferon/MHC-II (I)'),
    variable.name = 'event',
    value.name = 'event_pres'
)[, event := mapvalues(event, c('CASP8', 'categ3'), c('CASP8-mut', 'CASP8-mut, HRAS-mut,\nTP53-wt'))]
box_hnsc_ifn <- ggplot(bdata_hnsc_ifn) +
    geom_boxplot(aes(x = event_pres, y = `MP17 Interferon/MHC-II (I)`), outlier.shape = NA) +
    geom_point(aes(x = event_pres, y = `MP17 Interferon/MHC-II (I)`), position = position_jitter(width = 0.3), shape = 21, fill = 'black',
        stroke = 0.3, alpha = 0.3) +
    geom_text(
        data = data.table(
            x = c(1.5, 1.5),
            y = c(3.1, 3.1),
            label = pscores_hnsc[, paste('p =', c(
                signif(t.test(`MP17 Interferon/MHC-II (I)` ~ CASP8)$p.value, 2),
                signif(t.test(`MP17 Interferon/MHC-II (I)` ~ categ3)$p.value, 2)
            ))],
            event = c('CASP8-mut', 'CASP8-mut, HRAS-mut,\nTP53-wt')
        ),
        aes(x = x, y = y, label = label)
    ) +
    geom_segment(
        data = data.table(x = c(1, 1, 2, 1, 1, 2), xend = c(2, 1, 2, 2, 1, 2), y = rep(2.8, 6), yend = c(2.8, 2.7, 2.7, 2.8, 2.7, 2.7),
            event = c('CASP8-mut', 'CASP8-mut, HRAS-mut,\nTP53-wt')),
        aes(x = x, xend = xend, y = y, yend = yend)
    ) +
    facet_wrap(vars(event), strip.position = 'bottom') +
    theme_test() +
    theme(strip.placement = 'outside', strip.background = element_rect(fill = NA, colour = NA), strip.text = element_text(size = 11, vjust = 1)) +
    labs(x = NULL, y = 'Interferon/MHC-II (I) score', title = 'Interferon/MHC-II (I) and mutations in HNSC')

# mRNA subtypes:
subtypes_hnsc <- fread('~/TCGA_data/HNSC/subtypes.csv', na.strings = '', key = 'patient_id')
pscores_hnsc[, subtype := do.call(`[`, list(subtypes_hnsc, patient_id))$subtype]
pscores_hnsc[, subtype_classical := subtype == 'Classical']
# pscores_hnsc[, t.test(subtype_classical ~ pdac_high)] # No enrichment of the Classical subtype in the PDAC-high population.

clin_hnsc <- fread('~/TCGA_data/HNSC/Clinical.csv', key = 'patient_id')
pscores_hnsc[, site := do.call(`[`, list(clin_hnsc, patient_id))$anatomic_neoplasm_subdivision]
pscores_hnsc[, larynx := (site == 'larynx')]
# pscores_hnsc[, t.test(larynx ~ pdac_high)] # Larynx is enriched in the PDAC-high population.

pscores_hnsc[, larynx_nice := factor(ifelse(larynx, 'Larynx', 'Non-larynx'), levels = c('Non-larynx', 'Larynx'))]
pscores_hnsc[, nsd1_nice := factor(ifelse(NSD1, 'NSD1-mut', 'NSD1-wt'))]
bdata_hnsc_pdac_larynx <- rbind(
    pscores_hnsc[, .(id = patient_id, score = `MP30 PDAC-classical`, nsd1 = NSD1, categ = ifelse(larynx, 'Larynx', 'Non-larynx'))],
    pscores_hnsc[, .(id = patient_id, score = `MP30 PDAC-classical`, nsd1 = NSD1, categ = 'All tumors')]
)
bdata_hnsc_pdac_larynx[, nsd1 := factor(ifelse(nsd1, 'NSD1-mut', 'NSD1-wt'), levels = c('NSD1-wt', 'NSD1-mut'))]
bdata_hnsc_pdac_larynx[, categ := factor(categ, levels = c('All tumors', 'Non-larynx', 'Larynx'))]
box_hnsc_pdac_larynx <- ggplot(bdata_hnsc_pdac_larynx) +
    geom_boxplot(aes(x = categ, y = score, fill = nsd1, alpha = nsd1), outlier.shape = NA, width = 0.8,
        position = position_dodge2(padding = 0.1)) +
    geom_point(aes(x = categ, y = score, fill = nsd1, alpha = nsd1), shape = 21, stroke = 0.3,
        position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.55)) +
    geom_text(
        data = data.table(
            x = 1:3,
            y = rep(2.9, 3),
            label = paste('p =', signif(c(
                pscores_hnsc[, t.test(`MP30 PDAC-classical` ~ NSD1)$p.value],
                pscores_hnsc[larynx == FALSE, t.test(`MP30 PDAC-classical` ~ NSD1)$p.value],
                pscores_hnsc[larynx == TRUE, t.test(`MP30 PDAC-classical` ~ NSD1)$p.value]
            ), 2))
        ),
        aes(x = x, y = y, label = label)
    ) +
    geom_segment(
        data = data.table(
            x = c(0.8, 0.8, 1.2, 1.8, 1.8, 2.2, 2.8, 2.8, 3.2),
            xend = c(1.2, 0.8, 1.2, 2.2, 1.8, 2.2, 3.2, 2.8, 3.2),
            y = rep(2.7, 9),
            yend = rep(c(2.7, 2.65, 2.65), 3)
        ),
        aes(x = x, xend = xend, y = y, yend = yend)
    ) +
    scale_fill_manual(values = setNames(brewer.pal(9, 'Set1')[c(9, 6)], c('NSD1-wt', 'NSD1-mut'))) +
    scale_alpha_manual(values = c('NSD1-wt' = 0.25, 'NSD1-mut' = 0.5)) +
    theme_test() +
    theme(axis.text.x = element_text(size = 11), legend.text = element_text(size = 11)) +
    labs(title = 'PDAC-classical and NSD1 mutation in HNSC', x = NULL, y = 'PDAC-classical score', fill = NULL, alpha = NULL)

# Isolating PDAC-classical and mucinous components:

meta_luad <- fread('~/TCGA_data/LUAD/Cells.csv', key = 'sample_id', na.strings = '')
clin_luad <- fread('~/TCGA_data/LUAD/Clinical.csv', key = 'patient_id')
scores_luad <- fread('~/TCGA_data/LUAD/Scores.csv')
expmat_luad <- fread('~/TCGA_data/LUAD/Exp_data_TPM.csv')[, set_rownames(as.matrix(.SD), V1), .SDcols = -'V1']

meta_pdac <- fread('~/TCGA_data/PAAD/Cells.csv', key = 'sample_id', na.strings = '')
clin_pdac <- fread('~/TCGA_data/PAAD/Clinical.csv', key = 'patient_id')
scores_pdac <- fread('~/TCGA_data/PAAD/Scores.csv')
expmat_pdac <- fread('~/TCGA_data/PAAD/Exp_data_TPM.csv')[, set_rownames(as.matrix(.SD), V1), .SDcols = -'V1']

genes <- do.call(intersect, lapply(list(expmat_luad, expmat_pdac), function(x) {
    xgenes <- apply(x, 1, function(y) c(mean(y), var(y)))
    colnames(xgenes[, xgenes[1, ] > 3 & xgenes[2, ] > 0.2])
}))

muc_data <- copy(meta_luad[, .(sample_id, patient_id)])
muc_data[, muc := do.call(`[`, list(clin_luad, patient_id))[, grepl('^mucinous| mucinous', histological_type)]]
expmat_luad_sub <- t(apply(expmat_luad[genes, ], 1, function(x) x - mean(x)))
muc_diff <- rowMeans(expmat_luad_sub[, muc_data[muc == TRUE, sample_id]]) - rowMeans(expmat_luad_sub[, muc_data[muc == FALSE, sample_id]])

pdac_data <- scores_pdac[meta_program == 'MP30 PDAC-classical']
pdac_cor <- apply(expmat_pdac[genes, pdac_data$sample_id], 1, function(x) cor(x, pdac_data$score))
# all(names(pdac_cor) == names(muc_diff))

gene_scores <- data.table(gene = names(pdac_cor), pdac_cor = pdac_cor, muc_diff = muc_diff)

# Try differential expression between Classical and Basal PDAC tumours:
clas_data <- copy(meta_pdac[, .(sample_id, patient_id, subtype)])[, clas := (subtype == 'Classical')]
expmat_pdac_sub <- t(apply(expmat_pdac[genes, ], 1, function(x) x - mean(x)))
clas_diff <- rowMeans(expmat_pdac_sub[, clas_data[clas == TRUE, sample_id]]) - rowMeans(expmat_pdac_sub[, clas_data[clas == FALSE, sample_id]])
# all(gene_scores$gene == names(clas_diff))
gene_scores[, clas_diff := clas_diff]

scatter_spec1 <- ggplot(gene_scores) +
    geom_point(aes(x = pdac_cor, y = muc_diff), alpha = 0.5) +
    geom_hline(yintercept = c(0.55, 1.05), linetype = 'dashed', colour = 'tomato') +
    geom_vline(xintercept = c(0.4, 0.6), linetype = 'dashed', colour = 'tomato') +
    theme_test() +
    labs(x = 'Correlation with PDAC-classical score in PDAC', y = 'Differential expression in Mucinous LUAD',
        title = 'Classical PDAC and Mucinous LUAD specificity per gene')
pdac_spec1 <- gene_scores[pdac_cor > 0.6 & muc_diff < 0.55, gene]
joint1 <- gene_scores[pdac_cor > 0.6 & muc_diff > 1.05, gene]
muc_spec1 <- gene_scores[pdac_cor < 0.4 & muc_diff > 1.05, gene]

scatter_spec2 <- ggplot(gene_scores) +
    geom_point(aes(x = clas_diff, y = muc_diff), alpha = 0.5) +
    geom_hline(yintercept = c(0.7, 1.05), linetype = 'dashed', colour = 'tomato') +
    geom_vline(xintercept = c(0.35, 0.6), linetype = 'dashed', colour = 'tomato') +
    theme_test() +
    labs(x = 'Differential expression in Classical PDAC', y = 'Differential expression in Mucinous LUAD',
        title = 'Classical PDAC and Mucinous LUAD specificity per gene')
pdac_spec2 <- gene_scores[clas_diff > 0.6 & muc_diff < 0.7, gene]
joint2 <- gene_scores[clas_diff > 0.6 & muc_diff > 1.05, gene]
muc_spec2 <- gene_scores[clas_diff < 0.35 & muc_diff > 1.05, gene]

expmat_hnsc <- fread('~/TCGA_data/HNSC/Exp_data_TPM.csv')[, set_rownames(as.matrix(.SD), V1), .SDcols = -'V1']
expmat_hnsc <- expmat_hnsc[, colnames(expmat_hnsc) %in% scores_hnsc$sample_id]

setkey(meta_hnsc, patient_id)
expmat_hnsc <- expmat_hnsc[, meta_hnsc[sample_id %in% colnames(expmat_hnsc)][pscores_hnsc$patient_id, sample_id]]
setkey(meta_hnsc, sample_id)
colnames(expmat_hnsc) <- meta_hnsc[colnames(expmat_hnsc), patient_id]
# all(colnames(expmat_hnsc) == pscores_hnsc$patient_id) # TRUE

set.seed(1347)
pscores_hnsc[, c('pdac_spec1', 'joint1', 'muc_spec1') := lapply(
    list(pdac_spec1, joint1, muc_spec1),
    function(x) sig_score(expmat_hnsc, x, nbin = 50, n = 50)
)]
pscores_hnsc[, c('pdac_spec2', 'joint2', 'muc_spec2') := lapply(
    list(pdac_spec2, joint2, muc_spec2),
    function(x) sig_score(expmat_hnsc, x, nbin = 50, n = 50)
)]
bdata_spec <- rbind(
    melt(
        pscores_hnsc[, .(patient_id, NSD1, pdac_spec1, joint1, muc_spec1)],
        id.vars = c('patient_id', 'NSD1'),
        variable.name = 'geneset',
        value.name = 'score'
    ),
    melt(
        pscores_hnsc[, .(patient_id, NSD1, pdac_spec2, joint2, muc_spec2)],
        id.vars = c('patient_id', 'NSD1'),
        variable.name = 'geneset',
        value.name = 'score'
    )
)
bdata_spec[, NSD1 := factor(ifelse(NSD1, 'NSD1-mut', 'NSD1-wt'), levels = c('NSD1-wt', 'NSD1-mut'))]
bdata_spec[, categ := str_extract(geneset, '[12]')]
bdata_spec[, geneset := factor(
    mapvalues(gsub('[12]', '', geneset), c('pdac_spec', 'joint', 'muc_spec'), c('Classical PDAC', 'Shared', 'Mucinous LUAD')),
    levels = c('Mucinous LUAD', 'Classical PDAC', 'Shared')
)]
box_spec1 <- ggplot(bdata_spec[categ == 1]) +
    geom_boxplot(aes(x = geneset, y = score, fill = NSD1, alpha = NSD1), outlier.shape = NA, width = 0.8,
        position = position_dodge2(padding = 0.1)) +
    geom_point(aes(x = geneset, y = score, fill = NSD1, alpha = NSD1), shape = 21, stroke = 0.3,
        position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.55)) +
    geom_text(
        data = data.table(
            x = 1:3,
            y = rep(3.3, 3),
            label = pscores_hnsc[, paste('p =', signif(c(
                t.test(muc_spec1 ~ NSD1)$p.value,
                t.test(pdac_spec1 ~ NSD1)$p.value,
                t.test(joint1 ~ NSD1)$p.value
            ), 2))]
        ),
        aes(x = x, y = y, label = label)
    ) +
    geom_segment(
        data = data.table(
            x = c(0.8, 0.8, 1.2, 1.8, 1.8, 2.2, 2.8, 2.8, 3.2),
            xend = c(1.2, 0.8, 1.2, 2.2, 1.8, 2.2, 3.2, 2.8, 3.2),
            y = rep(3.1, 9),
            yend = rep(c(3.1, 3.05, 3.05), 3)
        ),
        aes(x = x, xend = xend, y = y, yend = yend)
    ) +
    scale_fill_manual(values = setNames(brewer.pal(9, 'Set1')[c(9, 6)], c('NSD1-wt', 'NSD1-mut'))) +
    scale_alpha_manual(values = c('NSD1-wt' = 0.25, 'NSD1-mut' = 0.5)) +
    theme_test() +
    theme(axis.text.x = element_text(size = 11, angle = 40, hjust = 1), legend.text = element_text(size = 11)) +
    labs(x = 'Gene set', y = 'Score', fill = NULL, alpha = NULL, title = 'Mucinous LUAD and Classical PDAC signatures in HNSC')
box_spec2 <- ggplot(bdata_spec[categ == 2]) +
    geom_boxplot(aes(x = geneset, y = score, fill = NSD1, alpha = NSD1), outlier.shape = NA, width = 0.8,
        position = position_dodge2(padding = 0.1)) +
    geom_point(aes(x = geneset, y = score, fill = NSD1, alpha = NSD1), shape = 21, stroke = 0.3,
        position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.55)) +
    geom_text(
        data = data.table(
            x = 1:3,
            y = rep(3.7, 3),
            label = pscores_hnsc[, paste('p =', signif(c(
                t.test(muc_spec2 ~ NSD1)$p.value,
                t.test(pdac_spec2 ~ NSD1)$p.value,
                t.test(joint2 ~ NSD1)$p.value
            ), 2))]
        ),
        aes(x = x, y = y, label = label)
    ) +
    geom_segment(
        data = data.table(
            x = c(0.8, 0.8, 1.2, 1.8, 1.8, 2.2, 2.8, 2.8, 3.2),
            xend = c(1.2, 0.8, 1.2, 2.2, 1.8, 2.2, 3.2, 2.8, 3.2),
            y = rep(3.5, 9),
            yend = rep(c(3.5, 3.45, 3.45), 3)
        ),
        aes(x = x, xend = xend, y = y, yend = yend)
    ) +
    scale_fill_manual(values = setNames(brewer.pal(9, 'Set1')[c(9, 6)], c('NSD1-wt', 'NSD1-mut'))) +
    scale_alpha_manual(values = c('NSD1-wt' = 0.25, 'NSD1-mut' = 0.5)) +
    theme_test() +
    theme(axis.text.x = element_text(size = 11, angle = 40, hjust = 1), legend.text = element_text(size = 11)) +
    labs(x = 'Gene set', y = 'Gene set score', fill = NULL, alpha = NULL, title = 'Mucinous LUAD and Classical PDAC signatures in HNSC')

pdf('../data/mut_hnsc_pdac.pdf', width = 32/2.54, height = 12/2.54)
plot_grid(
    box_hnsc_pdac_larynx + theme(legend.position = 'none', plot.margin = unit(c(5.5, 20, 5.5, 5.5), 'pt')),
    box_spec2 + theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 5.5, 20), 'pt')),
    get_legend(box_spec2),
    nrow = 1,
    ncol = 3,
    rel_widths = c(14, 14, 4)
)
dev.off()





# IFN response in LUAD:

meta_luad <- fread('~/TCGA_data/LUAD/Cells.csv', key = 'sample_id', na.strings = '')
scores_luad <- fread('~/TCGA_data/LUAD/Scores.csv')[grep('Alveolar|PDAC-classical|MP17|EMT-I[I]?$', meta_program)]
mut_data_luad <- fread('~/TCGA_data/LUAD/Mutation_data.csv')
mut_data_luad <- mut_data_luad[gene %in% cancer_genes]
setkey(mut_data_luad, gene, Variant_Classification)
mut_data_luad <- rbind(
    mut_data_luad[Variant_Classification %in% c('Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site', 'Nonstop_Mutation',
        'Translation_Start_Site')],
    mut_data_luad[
        mut_data_luad[
            Variant_Classification %in% c('Missense_Mutation', 'In_Frame_Del', 'In_Frame_Ins'),
            .(n = length(unique(patient_id))),
            by = .(gene, Variant_Classification)
        ][n >= 2, -'n']
    ]
)
pscores_luad <- copy(scores_luad)[, patient_id := do.call(`[`, args = list(meta_luad, sample_id))$patient_id][, -'sample_id']
setcolorder(pscores_luad, 'patient_id')
setkey(pscores_luad, patient_id)
pscores_luad <- pscores_luad[patient_id %in% mut_data_luad$patient_id]
mut_data_luad <- mut_data_luad[patient_id %in% pscores_luad$patient_id]
pids_tab <- pscores_luad[meta_program == meta_program[1], table(patient_id)]
pids <- names(pids_tab)[pids_tab == 1]
pscores_luad <- pscores_luad[patient_id %in% pids]
mut_data_luad <- mut_data_luad[patient_id %in% pids]

clin_luad <- fread('~/TCGA_data/LUAD/Clinical.csv', key = 'patient_id')
pscores_luad <- dcast(pscores_luad[, .(patient_id, meta_program, score)], patient_id ~ meta_program, value.var = 'score')
pscores_luad[, mucinous := grepl('^mucinous| mucinous', do.call(`[`, .(clin_luad, patient_id))$histological_type)]
genes <- c('TP53', 'XIRP2', 'KRAS', 'STK11', 'KEAP1', 'NFE2L2')
pscores_luad[, (genes) := lapply(genes, function(g) patient_id %in% mut_data_luad[gene == g, patient_id])]

pscores_luad_ifnmut <- melt(
    pscores_luad[, .(patient_id, `MP17 Interferon/MHC-II (I)`, STK11, KEAP1, NFE2L2)],
    id.vars = c('patient_id', 'MP17 Interferon/MHC-II (I)'),
    variable.name = 'gene',
    value.name = 'mut'
)[, mut := factor(ifelse(mut, 'Mutated', 'Wild-type'), levels = c('Wild-type', 'Mutated'))]
box_luad_ifn <- ggplot(pscores_luad_ifnmut) +
    geom_boxplot(aes(x = gene, y = `MP17 Interferon/MHC-II (I)`, fill = mut), width = 0.8, position = position_dodge2(padding = 0.1),
        outlier.shape = NA, alpha = 0.25) +
    geom_point(aes(x = gene, y = `MP17 Interferon/MHC-II (I)`, fill = mut), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.55),
        shape = 21, stroke = 0.3, alpha = 0.75) +
    geom_segment(
        aes(x = x, xend = xend, y = y, yend = yend),
        data = data.table(
            x = c(0.8, 0.8, 1.2, 1.8, 1.8, 2.2, 2.8, 2.8, 3.2),
            xend = c(1.2, 0.8, 1.2, 2.2, 1.8, 2.2, 3.2, 2.8, 3.2),
            y = rep(2.2, 9),
            yend = rep(c(2.2, 2.15, 2.15), 3)
        )
    ) +
    geom_text(
        aes(x = x, y = y, label = label),
        data = data.table(
            x = 1:3,
            y = 2.4,
            label = paste('p =', pscores_luad[, signif(c(
                t.test(`MP17 Interferon/MHC-II (I)` ~ STK11)$p.value,
                t.test(`MP17 Interferon/MHC-II (I)` ~ KEAP1)$p.value,
                t.test(`MP17 Interferon/MHC-II (I)` ~ NFE2L2)$p.value
            ), 2)])
        )
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(name = NULL, values = setNames(ggplot_colours(4)[c(2, 4)], c('Wild-type', 'Mutated'))) +
    theme_test() +
    theme(legend.justification = c(0, 1), axis.text.x = element_text(size = 12)) +
    labs(title = 'Interferon/MHC-II (I) and mutations in LUAD', x = NULL, y = 'Interferon/MHC-II (I) score')

# Interferon-related figures:
box_hnsc_ifn_grob <- ggplotGrob(box_hnsc_ifn + theme(plot.margin = unit(c(5.5, 15, 5.5, 5.5), 'pt')))
box_luad_ifn_grob <- ggplotGrob(box_luad_ifn + theme(axis.text.x = element_text(margin = margin(t = 5)), plot.margin = unit(c(5.5, 5.5, 5.5, 15), 'pt')))
box_hnsc_ifn_grob$heights <- unit(c(0.3, 0, 0.7, 0, 0, 0, 9, 0.5, 0.2, 0.8, 0, 0, 0, 0.3), 'cm')
box_luad_ifn_grob$heights <- unit(c(0.3, 0, 0.7, 0, 0, 0, 9, 0.8, 0, 0, 0, 1), 'cm')

pdf('../data/mut_ifn.pdf', width = 28/2.54, height = 11.8/2.54)
plot_grid(box_hnsc_ifn_grob, box_luad_ifn_grob, nrow = 1, ncol = 2, rel_widths = c(0.4, 0.6))
dev.off()





# Glutathione in KIRC:

meta_kirc <- fread('~/TCGA_data/KIRC/Cells.csv', key = 'sample_id', na.strings = '')
scores_kirc <- fread('~/TCGA_data/KIRC/Scores.csv')[grep('Glutathione|Cell Cy', meta_program)]
cna_data_kirc <- fread('~/TCGA_data/KIRC/CNV_GISTIC_data.csv')[, set_rownames(as.matrix(.SD), V1), .SDcols = -'V1']
pids_cna <- apply(str_split_fixed(colnames(cna_data_kirc), '\\.', 4)[, 1:3], 1, paste, collapse = '.')
colnames(cna_data_kirc) <- pids_cna
pids_cna_tab <- table(pids_cna)
cna_data_kirc <- cna_data_kirc[,
    pids_cna %in% meta_kirc[sample_id %in% scores_kirc$sample_id, patient_id] & !(pids_cna %in% names(pids_cna_tab)[pids_cna_tab > 1])
]
mut_data_kirc <- fread('~/TCGA_data/KIRC/Mutation_data.csv')
# mut_data_kirc <- mut_data_kirc[gene %in% cancer_genes]
setkey(mut_data_kirc, gene, Variant_Classification)
mut_data_kirc <- rbind(
    mut_data_kirc[Variant_Classification %in% c('Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site', 'Nonstop_Mutation',
        'Translation_Start_Site')],
    mut_data_kirc[
        mut_data_kirc[
            Variant_Classification %in% c('Missense_Mutation', 'In_Frame_Del', 'In_Frame_Ins'),
            .(n = length(unique(patient_id))),
            by = .(gene, Variant_Classification)
        ][n >= 2, -'n']
    ]
)
pscores_kirc <- copy(scores_kirc)[, patient_id := do.call(`[`, args = list(meta_kirc, sample_id))$patient_id]
pids_tab <- pscores_kirc[meta_program == meta_program[1], table(patient_id)]
pids <- names(pids_tab)[pids_tab == 1]
pscores_kirc <- pscores_kirc[patient_id %in% pids]
mut_data_kirc <- mut_data_kirc[patient_id %in% pids]
cna_data_kirc <- cna_data_kirc[, colnames(cna_data_kirc) %in% pids]
pscores_kirc <- dcast(pscores_kirc[, .(patient_id, meta_program, score)], patient_id ~ meta_program, value.var = 'score')
expmat_kirc <- fread('~/TCGA_data/KIRC/Exp_data_TPM.csv')[, set_rownames(as.matrix(.SD), V1), .SDcols = -'V1']
expmat_kirc <- expmat_kirc[, colnames(expmat_kirc) %in% scores_kirc$sample_id]
colnames(expmat_kirc) <- meta_kirc[colnames(expmat_kirc), patient_id]
expmat_kirc <- expmat_kirc[, colnames(expmat_kirc) %in% pscores_kirc$patient_id]
pscores_kirc[, glut_low := (`MP38 Glutathione` < -1)]
# pscores_kirc[, plot(sort(`MP38 Glutathione`))]
# abline(h = -1)

# Density barplot:
pscores_kirc[, glut_score_cut := cut(`MP38 Glutathione`, breaks = seq(-4.5, 2.1, by = 0.2))]
pscores_kirc[, glut_score_cut_mid := as.numeric(gsub('\\(', '', strsplit(as.character(glut_score_cut), ',')[[1]][1])) + 0.1, by = patient_id]
density_kirc_glut <- ggplot(pscores_kirc, aes(x = glut_score_cut_mid)) +
    geom_bar(width = 0.15) +
    theme_bw() +
    geom_vline(xintercept = -1, colour = 'red', linetype = 'dashed') +
    scale_x_continuous(limits = c(-4.6, 2.2), expand = c(0, 0)) +
    labs(x = 'Score', y = 'Count', title = 'Glutathione score distribution in KIRC')

# Mutation enrichments:
pscores_kirc[, c('VHL mut/del', 'PBRM1 mut/del') := lapply(c('VHL', 'PBRM1'), # Variables encompassing VHL/PBRM1 mutations and deletions
    function(g) patient_id %in% mut_data_kirc[gene == g, patient_id] | patient_id %in% colnames(cna_data_kirc)[cna_data_kirc[g, ] == -2])]
pscores_kirc[, 'CDKN2A del' := patient_id %in% colnames(cna_data_kirc)[cna_data_kirc['CDKN2A', ] == -2]]
bdata <- melt(
    pscores_kirc[, .(patient_id, glut_low, `VHL mut/del`, `PBRM1 mut/del`, `CDKN2A del`)],#), `MUC17 mut`, `PKHD1 mut`)],
    id.vars = c('patient_id', 'glut_low'),
    variable.name = 'event',
    value.name = 'event_present'
)[, .(perc = 100*sum(event_present)/.N), by = .(event, glut_low)]
bdata[, glut_categ := factor(ifelse(glut_low, 'Glutathione-low', 'Glutathione-high'), levels = c('Glutathione-low', 'Glutathione-high'))]
bar_kirc_enrich <- ggplot(bdata) +
    geom_col(aes(x = event, y = perc, group = glut_categ, fill = glut_categ), position = position_dodge2(), width = 0.7) +
    geom_segment(
        data = data.table(x = c(1:3 - 0.175, 1:3 - 0.175, 1:3 + 0.175), xend = c(1:3 + 0.175, 1:3 - 0.175, 1:3 + 0.175),
            y = c(45, 45, 15, 45, 45, 15, 45, 45, 15), yend = c(45, 45, 15, 44.5, 44.5, 14.5, 44.5, 44.5, 14.5)),
        aes(x = x, xend = xend, y = y, yend = yend)
    ) +
    geom_text(
        data = data.table(
            x = 1:3,
            y = c(46, 46, 16),
            label = pscores_kirc[,
                .(pval_adj = p.adjust(sapply(.SD, function(x) fisher.test(x, glut_low)$p.value), method = 'BH')),
                .SDcols = c('VHL mut/del', 'PBRM1 mut/del', 'CDKN2A del')#, 'MUC17 mut', 'PKHD1 mut')
            ][, ifelse(pval_adj < 0.01, ifelse(pval_adj < 0.001, '***', '**'), '*')]
        ),
        aes(x = x, y = y, label = label),
        size = 5
    ) +
    scale_fill_manual(name = NULL, values = setNames(ggplot_colours(2), c('Glutathione-low', 'Glutathione-high'))) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 11)) +
    labs(x = NULL, y = 'Percentage of tumours with alteration', title = 'Glutathione and genomic\nalterations in KIRC')

# Boxplot of cell cycle genes between glut-low and -high groups:
pscores_kirc[, cc_ave := rowMeans(.SD), .SDcols = c('MP1 Cell Cycle - G2/M', 'MP2 Cell Cycle - G1/S')]
pscores_kirc[, glut_categ := factor(ifelse(glut_low, 'Glutathione-low', 'Glutathione-high'), levels = c('Glutathione-low', 'Glutathione-high'))]
box_kirc_glut_cc <- ggplot(pscores_kirc) +
    geom_boxplot(aes(x = glut_categ, y = cc_ave, fill = glut_categ), width = 0.9, position = position_dodge2(), outlier.shape = NA, alpha = 0.25) +
    geom_point(aes(x = glut_categ, y = cc_ave, fill = glut_categ), position = position_jitterdodge(jitter.width = 1.4),
        shape = 21, stroke = 0.3, colour = 'black', alpha = 0.75) +
    geom_segment(data = data.table(x = c(1, 1, 2), xend = c(2, 1, 2), y = c(3, 3, 3), yend = c(3, 2.9, 2.9)),
        aes(x = x, xend = xend, y = y, yend = yend)) +
    annotate('text', x = 1.5, y = 3.2, label = pscores_kirc[, paste('p =', signif(t.test(cc_ave ~ glut_low)$p.value, 2))]) +
    theme_test() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = NULL, y = 'Average of G1/S and G2/M scores', fill = 'Category', title = 'Glutathione and cell cycle in KIRC')

# K-M plot:
clin_kirc <- fread('../data/clin_test_data.csv', na.strings = '')[cancer_type == 'KIRC']
sdata_kirc <- merge(clin_kirc[, .(patient_id, OS, OSc)], pscores_kirc[, .(patient_id, glut_low)], by = 'patient_id')
survmod <- survival::survfit(survival::Surv(OS, OSc) ~ glut_low, data = sdata_kirc)
km_kirc <- survminer::ggsurvplot(survmod, conf.int = TRUE, size = 0.5)$plot +
    annotate('text', x = 2000, y = 0.1, label = paste('p =', signif(survminer::surv_pvalue(survmod)$pval, 2))) +
    scale_colour_manual(
        name = 'Category',
        labels = c('glut_low=FALSE' = 'Glutathione-high', 'glut_low=TRUE' = 'Glutathione-low'),
        values = setNames(ggplot_colours(2), c('glut_low=TRUE', 'glut_low=FALSE'))
    ) +
    scale_fill_manual(
        name = 'Category',
        labels = c('glut_low=FALSE' = 'Glutathione-high', 'glut_low=TRUE' = 'Glutathione-low'),
        values = setNames(ggplot_colours(2), c('glut_low=TRUE', 'glut_low=FALSE'))
    ) +
    theme_test() +
    labs(x = 'Time (days)', title = 'Glutathione and survival in KIRC')

density_kirc_glut_grob <- ggplotGrob(density_kirc_glut + theme(plot.margin = unit(c(5.5, 10, 5.5, 5.5), 'pt')))
bar_kirc_enrich_grob <- ggplotGrob(bar_kirc_enrich + theme(legend.position = 'none', plot.margin = unit(c(5.5, 10, 5.5, 10), 'pt')))
km_kirc_grob <- ggplotGrob(km_kirc + theme(legend.position = 'none', plot.margin = unit(c(5.5, 10, 5.5, 10), 'pt')))
box_kirc_glut_cc_grob <- ggplotGrob(box_kirc_glut_cc + theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 5.5, 10), 'pt')))

density_kirc_glut_grob$heights <- unit(c(0.2, 0, 0.7, 0, 0, 0, 9, 0.5, 0.6, 0, 0, 0.2), 'cm')
bar_kirc_enrich_grob$heights <- unit(c(0.2, 0, 1.2, 0, 0, 0, 7.5, 2.1, 0, 0, 0, 0.2), 'cm')
km_kirc_grob$heights <- unit(c(0.2, 0, 0.7, 0, 0, 0, 8, 0.5, 0.6, 0, 0, 1.2), 'cm')
box_kirc_glut_cc_grob$heights <- unit(c(0.2, 0, 0.7, 0, 0, 0, 9.8, 0, 0, 0, 0, 0.5), 'cm')

pdf('../data/mut_kirc_glut.pdf', width = 42/2.54, height = sum(density_kirc_glut_grob$heights)/2.54)
plot_grid(density_kirc_glut_grob, bar_kirc_enrich_grob, km_kirc_grob, box_kirc_glut_cc_grob, get_legend(bar_kirc_enrich), nrow = 1, ncol = 5,
    rel_widths = c(12, 7, 11, 8, 4))
dev.off()
