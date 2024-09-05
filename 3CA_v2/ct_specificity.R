library(ggplot2)
library(ggtext)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(randomcoloR)
library(data.table)
library(magrittr)
library(plyr)
library(stringr)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'), encoding = 'UTF-8')

hgnc_complete_set <- fread('../data/hgnc_complete_set_2023-04-13.txt', key = 'symbol')
hgnc_complete_set <- hgnc_complete_set[!(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])]
alias_table <- make_alias_table(hgnc_complete_set)

cts <- c('Malignant', 'Epithelial', 'Fibroblast', 'T_cell', 'Macrophage', 'Endothelial', 'B_cell')





# Pseudobulks generated in pseudobulks.R:

pbulks_all <- slapply(gsub('.csv$', '', dir('../data/pseudobulks')), function(abbr) {
    cat(abbr, '\n')
    fread(paste0('..data/pseudobulks/', abbr, '.csv'))
})

pbulks <- slapply(cts, function(ct) {
    cat(ct, '\n')
    pbulks_ct_list <- lapply(pbulks_all, function(dt) dt[, c('gene', names(dt)[grep(gsub('_', '', ct), names(dt))]), with = FALSE])
    Reduce(function(x, y) merge(x, y, all = TRUE), pbulks_ct_list)
})

# The following generates a table containing study, cancer type, sample and cancer type associated to each pseudobulk ID:
samples_all <- paths_table[, {
    abbr <- paste(gsub(' | et al. ', '', unique(study)), gsub(' |/', '-', unique(cancer_type)), sep = '_')
    if(abbr %in% names(pbulks_all)) {
        cat(abbr, '\n')
        ids <- names(pbulks_all[[abbr]][, -'gene'])
        ids_smpl <- gsub(paste0(abbr, '_[0-9]_|_', paste(gsub('_', '', cts), collapse = '|_')), '', ids)
        samples <- fread(
            paste0('~/../shared/pan_cancer_datasets/', directory[1], '/samples.csv'),
            colClasses = c(sample = 'character', cancer_type = 'character', disease = 'character')
        )
        if(abbr == 'Puram2017_Head-and-Neck') { # Split HNSCC samples into HPV pos and neg
            samples <- samples[sample %in% c('5', '6', '16', '17', '18', '20', '22', '25', '26', '28')][, cancer_type := 'HNSCC (HPV-neg.)']
        } else if(abbr == 'Cillo2020_Head-and-Neck') {
            samples[, cancer_type := ifelse(p16_status == 'p16+', 'HNSCC (HPV-pos.)', 'HNSCC (HPV-neg.)')]
        } else if(abbr == 'Kürten2021_Head-and-Neck') {
            samples[, cancer_type := ifelse(additional_tumor_characterisics == 'HPV Positive', 'HNSCC (HPV-pos.)', 'HNSCC (HPV-neg.)')]
        } else if(abbr == 'UnpublishedHNSCCHTAN_Head-and-Neck') {
            samples[, cancer_type := ifelse(hpv == 'HPV+', 'HNSCC (HPV-pos.)', 'HNSCC (HPV-neg.)')]
        }
        if(abbr == 'Lee2020_Colorectal') { # Split CRC and Gastric into MSS and MSI
            samples[, cancer_type := ifelse(grepl('mss', genetic_hormonal_features), 'CRC (MSS)', 'CRC (MSI)')]
        } else if(abbr == 'Pelka2021_Colorectal') {
            samples[, cancer_type := ifelse(MSIStatus == 'MSS', 'CRC (MSS)', 'CRC (MSI)')]
        } else if(abbr == 'Zhang2018_Colorectal') {
            samples[, cancer_type := ifelse(msi_status == 'MSS', 'CRC (MSS)', 'CRC (MSI)')]
        } else if(abbr == 'Kumar2022_Gastric') {
            samples[, cancer_type := ifelse(is.na(subtype) | subtype == '', NA, ifelse(subtype == 'MSI', 'Gastric (MSI)', 'Gastric (MSS)'))]
        } else if(abbr == 'Zhang2019_Gastric') {
            samples[, cancer_type := NA] # There's only 1 true gastric cancer sample in this dataset, and it's not annotated with MSI status.
        }
        if('cancer_type' %in% names(samples)) {
            samples[, .(id = ids[ids_smpl == unique(sample)]), by = .(sample, disease = cancer_type)]
        } else if('disease' %in% names(samples)) samples[, .(id = ids[ids_smpl == unique(sample)]), by = .(sample, disease)]
    }
}, by = .(study, cancer_type)]
setkey(samples_all, id)




cor_data <- lapply(cts, function(ct) {
    
    cat(ct, '\n')
    
    # Get pseudobulks:
    pbulks_ct <- pbulks[[ct]][, set_rownames(as.matrix(.SD), gene), .SDcols = -'gene']
    pbulks_ct <- pbulks_ct[, colnames(pbulks_ct) %in% samples_all$id]
    
    # Log transform:
    pbulks_ct <- log2(pbulks_ct + 1)
    
    # Filter genes:
    pbulks_ct <- pbulks_ct[!(rownames(pbulks_ct) %in% hgnc_complete_set[grepl('RP[SL]', symbol) | grepl('[Rr]ibosom', gene_group), symbol]), ]
    pbulks_ct <- pbulks_ct[rownames(pbulks_ct) %in% names(sort(apply(pbulks_ct, 1, var, na.rm = TRUE), decreasing = TRUE)[1:5000]), ]
        
    # Compute correlation matrix:
    cor_ct <- cor(pbulks_ct, use = 'pairwise.complete.obs')
    
    # Prepare correlation table:
    cor_ct <- as.data.table(cor_ct, keep.rownames = 'sample_id1')
    cor_ct <- melt(cor_ct, id.vars = 'sample_id1', variable.name = 'sample_id2', value.name = 'corr', variable.factor = FALSE)
    cor_ct[, c('study1', 'disease1') := samples_all[sample_id1, .(paste(study, cancer_type, sep = ' - '), disease)]]
    cor_ct[, c('study2', 'disease2') := samples_all[sample_id2, .(paste(study, cancer_type, sep = ' - '), disease)]]
    cor_ct[, cell_type := ct]
    setcolorder(cor_ct, 'cell_type')
    
    return(cor_ct)
    
}) %>% rbindlist





cor_data <- cor_data[!is.na(corr) & !is.na(disease1) & !is.na(disease2) &
        !(disease1 %in% c('Normal', 'Premalignant', '')) & !(disease2 %in% c('Normal', 'Premalignant', '')) &
        !grepl('Other/Models', study1) & !grepl('Other/Models', study2)]

for(v in c('disease1', 'disease2')) {
    cor_data[, (v) := mapvalues(
        gsub(' Cancer', '', get(v)),
        c('Acute Myeloid Leukemia', 'Chronic Lymphocytic Leukemia', 'Chronic Myeloid Leukemia', 'Clear Cell Renal Cell Carcinoma', 'Colorectal',
            'Cutaneous Basal Cell Carcinoma', 'Cutaneous Squamous Cell Carcinoma', 'Cutaneous T Cell Lymphoma', 'Diffuse Large B Cell Lymphoma',
            'Follicular Lymphoma', 'Glioblastoma', 'Head and Neck Squamous Cell Carcinoma', 'Hepatocellular Carcinoma', 'Lung Adenocarcinoma',
            'Lung Squamous Cell Carcinoma', 'Multiple Myeloma', 'Nasopharyngeal Carcinoma', 'Neuroendocrine Tumor',
            'Pancreatic Ductal Adenocarcinoma', 'Small Cell Lung'),
        c('AML', 'CLL', 'CML', 'ccRCC', 'CRC', 'Skin BCC', 'Skin SCC', 'CTCL', 'DLBCL', 'FL', 'GBM', 'HNSCC', 'HCC', 'Lung Adeno.', 'Lung Squamous',
            'MM', 'Nasopharyngeal', 'NET', 'PDAC', 'SCLC'),
        warn_missing = FALSE
    )]
}

mdata <- cor_data[, {
    dlist <- .SD[,
        .(n_sample = length(unique(sample_id1)), n_study = length(unique(study1))),
        by = disease1
    ][n_sample >= 10 & n_study >= 2, disease1] # Changing to 3 removes the worst outliers but reduces statistical power
    .SD[
        disease1 %in% dlist & disease2 %in% dlist,
        .(mean_corr = mean(corr), n_sample = length(unique(sample_id1))),
        by = .(disease1, disease2, study1, study2)
    ]
}, by = cell_type]

mdata_merge <- Reduce(merge, list(
    mdata[disease1 == disease2 & study1 == study2, .(cell_type, disease = disease1, study = study1, a = mean_corr, n_a = n_sample)],
    mdata[
        disease1 == disease2 & study1 != study2,
        .(b = sum(mean_corr*n_sample)/sum(n_sample), n_b = sum(n_sample)),
        by = .(cell_type, disease = disease1, study = study1)
    ],
    mdata[
        disease1 != disease2 & study1 != study2,
        .(c = sum(mean_corr*n_sample)/sum(n_sample), n_c = sum(n_sample)),
        by = .(cell_type, disease = disease1, study = study1)
    ]
))
mdata_merge[, c('n', 'se_a', 'se_b', 'se_c', 'o_a', 'o_b', 'o_c') := c(.(.N, sd(a)/.N, sd(b)/(.N - 1), sd(c)/.N), lapply(.SD, function(x) {
    x > quantile(x, 0.75) + 1.5*IQR(x) | x < quantile(x, 0.25) - 1.5*IQR(x)
})), by = .(cell_type, disease), .SDcols = c('a', 'b', 'c')]
mdata_merge[disease == 'HNSCC (HPV-neg.)', disease := 'HNSCC (HPV<sup>\u2013</sup>)']
mdata_merge[disease == 'HNSCC (HPV-pos.)', disease := 'HNSCC (HPV<sup>+</sup>)']
mdata_merge[, cell_type := gsub('_', ' ', cell_type)]

spec_p_data <- mdata_merge[,
    .(spec_p = 1 - sum(a*n_a)/sum(n_a)),
    by = .(cell_type, disease)
][, cell_type := factor(cell_type, levels = .SD[, .(m = median(spec_p)), by = cell_type][order(m), cell_type])]
spec_d_data <- mdata_merge[
    se_b < 0.1 & se_c < 0.1,
    .(spec_d = sum(b*n_b)/sum(n_b) - sum(c*n_c)/sum(n_c)),
    by = .(cell_type, disease)
][, cell_type := factor(cell_type, levels = .SD[, .(m = median(spec_d)), by = cell_type][order(m), cell_type])]

test_mat_d <- dcast(spec_d_data[, .(disease, cell_type, spec_d)], disease ~ cell_type)
test_mat_d <- test_mat_d[, set_rownames(as.matrix(.SD), disease), .SDcols = -'disease']
pmat_d <- apply(test_mat_d, 2, function(x) apply(test_mat_d, 2, function(y) t.test(x = x, y = y, paired = TRUE)$p.value))
pdt_d <- pmat_d
for(i in 1:nrow(pmat_d)) for(j in 1:ncol(pmat_d)) if(i > j) pdt_d[i, j] <- NaN
pdt_d <- melt(as.data.table(pdt_d, keep.rownames = 'ct1'), id.var = 'ct1', variable.name = 'ct2', value.name = 'pval')
pdt_d <- pdt_d[is.finite(pval)]
pdt_d[, pval_adj := p.adjust(pval, 'BH')]

test_mat_p <- dcast(spec_p_data[, .(disease, cell_type, spec_p)], disease ~ cell_type)
test_mat_p <- test_mat_p[, set_rownames(as.matrix(.SD), disease), .SDcols = -'disease']
pmat_p <- apply(test_mat_p, 2, function(x) apply(test_mat_p, 2, function(y) t.test(x = x, y = y, paired = TRUE)$p.value))
pdt_p <- pmat_p
for(i in 1:nrow(pmat_p)) for(j in 1:ncol(pmat_p)) if(i > j) pdt_p[i, j] <- NaN
pdt_p <- melt(as.data.table(pdt_p, keep.rownames = 'ct1'), id.var = 'ct1', variable.name = 'ct2', value.name = 'pval')
pdt_p <- pdt_p[is.finite(pval)]
pdt_p[, pval_adj := p.adjust(pval, 'BH')]
seg_data <- pdt_p[pval_adj < 0.05, .(
    pval = pval_adj,
    x = spec_p_data[cell_type == ct1, unique(as.numeric(cell_type))],
    xend = spec_p_data[cell_type == ct2, unique(as.numeric(cell_type))]
), by = .(ct1, ct2)]
seg_data[, l := ifelse(pval < 0.001, '***', ifelse(pval < 0.01, '**', ifelse(pval < 0.05, '*', '')))]
seg_data[, c('y', 'yend') := .(0.25 + .I/66, 0.25 + .I/66)]

ggplot(spec_d_data, aes(x = cell_type, y = spec_d)) +
    geom_hline(yintercept = 0, linetype = 'dotted', colour = 'black') +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = disease), shape = 21, size = 3, position = position_jitter(w = 0.2)) +
    scale_fill_manual(values = distinctColorPalette(spec_d_data[, length(unique(disease))])) +
    geom_segment(
        aes(x = x, xend = xend, y = y, yend = yend),
        data = data.table(
            x = c(1, 2.8, 2.8, 7, 3.2, 3.2, 6),
            xend = c(5, 2.8, 7, 7, 3.2, 6, 6),
            y = c(0.21, 0.21, 0.38, 0.38, 0.21, 0.33, 0.33),
            yend = c(0.21, 0.38, 0.38, 0.37, 0.33, 0.33, 0.32)
        )
    ) +
    annotate(geom = 'text', x = c(5, 4.5), y = c(0.39, 0.34), label = c('**', '*'), size = 7) +
    theme_test() +
    theme(
        axis.text.x = element_text(angle = 35, hjust = 1, size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        legend.text = element_markdown(size = 14),
        legend.title = element_text(size = 16),
        legend.position = 'bottom',
        legend.box.margin = margin(t = 30),
        legend.justification = 'right'
    ) +
    guides(fill = guide_legend(nrow = 7, ncol = 3, title.position = 'top')) +
    labs(x = 'Cell type', y = 'Cancer type specificity', fill = 'Cancer type') +
    ylim(spec_d_data[, min(spec_d)], 0.4)
ggplot(spec_p_data, aes(x = cell_type, y = spec_p)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = disease), shape = 21, size = 3, position = position_jitter(w = 0.2)) +
    scale_fill_manual(values = distinctColorPalette(spec_p_data[, length(unique(disease))])) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend), data = seg_data) +
    geom_segment(aes(x = x, xend = x, y = y, yend = y - 0.002), data = seg_data) +
    geom_segment(aes(x = xend, xend = xend, y = y, yend = y - 0.002), data = seg_data) +
    geom_text(aes(x = (x + xend)/2, y = y + 0.002, label = l), data = seg_data, size = 7) +
    theme_test() +
    theme(
        axis.text.x = element_text(angle = 35, hjust = 1, size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        legend.text = element_markdown(size = 14),
        legend.title = element_text(size = 16)
    ) +
    labs(x = 'Cell type', y = 'Patient specificity', fill = 'Cancer type')
