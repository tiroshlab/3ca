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

cancer_genes <- as.data.table(read_xlsx('/home/labs/tirosh/tyler/pan_cancer/cancer5000.xlsx', skip = 1, n_max = 260))[`Cancer5000-S (219)` == 1, gene]

hgnc_complete_set <- fread('../data/hgnc_complete_set.txt', key = 'symbol')[
    !(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])
][!grepl('unplaced', location), chrarm := str_extract(location, '[0-9]*[XY]*[pq]')][!is.na(chrarm)][
    order(as.numeric(mapvalues(gsub('[pq]', '', chrarm), c('X', 'Y'), c(23, 24))), str_extract(chrarm, '[pq]'))
]





meta_hnsc <- fread('/home/labs/tirosh/tyler/TCGA_data/HNSC/Cells.csv', key = 'sample_id', na.strings = '')
scores_hnsc <- fread('/home/labs/tirosh/tyler/TCGA_data/HNSC/Scores.csv')[grep('MP17 Interferon|PDAC-classical', meta_program)]
mut_data_hnsc <- fread('/home/labs/tirosh/tyler/TCGA_data/HNSC/Mutation_data.csv')
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

subtypes_hnsc <- fread('../../TCGA_data/HNSC/subtypes.csv', na.strings = '', key = 'patient_id')
pscores_hnsc[, subtype := do.call(`[`, list(subtypes_hnsc, patient_id))$subtype]
pscores_hnsc[, subtype_classical := subtype == 'Classical']
# pscores_hnsc[, t.test(subtype_classical ~ pdac_high)] # No enrichment of the Classical subtype in the PDAC-high population.

clin_hnsc <- fread('../../TCGA_data/HNSC/Clinical.csv', key = 'patient_id')
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

# Just the "All tumors" boxes:
box_hnsc <- ggplot(bdata_hnsc_pdac_larynx[categ == 'All tumors']) +
    geom_boxplot(aes(x = nsd1, y = score), outlier.shape = NA, width = 0.8) +
    geom_point(aes(x = nsd1, y = score), fill = 'grey80', size = 2, shape = 21, stroke = 0.3, position = position_jitter(width = 0.3)) +
    geom_text(
        data = data.table(x = 1.5, y = 2.9, label = paste('p =', signif(c(pscores_hnsc[, t.test(`MP30 PDAC-classical` ~ NSD1)$p.value]), 2))),
        aes(x = x, y = y, label = label)
    ) +
    geom_segment(
        data = data.table(x = c(1, 1, 2), xend = c(2, 1, 2), y = rep(2.7, 3), yend = c(2.7, 2.65, 2.65)),
        aes(x = x, xend = xend, y = y, yend = yend)
    ) +
    theme_test() +
    theme(axis.text.x = element_text(size = 11), plot.title = element_text(hjust = 0.5)) +
    labs(title = 'PDAC-classical and\nNSD1 mutation in HNSC', x = NULL, y = 'PDAC-classical score')





library(seriation)

scores_patient_mut_all <- fread('../data/scores_patient_mut.csv')

exp_data <- slapply(c('LUAD', 'PAAD'), function(ct) {
    expmat <- fread(paste0('/home/labs/tirosh/tyler/TCGA_data/', ct, '/Exp_data_TPM.csv'))[, set_rownames(as.matrix(.SD), V1), .SDcols = -'V1']
    meta <- fread(paste0('/home/labs/tirosh/tyler/TCGA_data/', ct, '/Cells.csv'), na.strings = '')[
        sample_type != 'normal' & patient_id %in% scores_patient_mut_all[cancer_type == ct, patient_id]
    ]
    expmat <- expmat[, meta$sample_id]
    colnames(expmat) <- meta$patient_id
    return(list(expmat = expmat, meta = meta))
})
top_genes <- slapply(
    names(exp_data),
    function(ct) head(names(sort(apply(exp_data[[ct]]$expmat, 1, function(x) median(abs(x - median(x)))), decreasing = TRUE)), 2500)
)
cormat <- slapply(names(exp_data), function(ct) cor(t(apply(exp_data[[ct]]$expmat[top_genes[[ct]], ], 1, function(x) x - median(x)))))
clust <- slapply(names(exp_data), function(ct) seriate(dist(cormat[[ct]]), method = 'HC_average')[[1]])
pdac_clust <- labels(as.dendrogram(clust$LUAD)[[1]][[2]][[1]])

luad_clinical <- fread('/home/labs/tirosh/tyler/TCGA_data/LUAD/Clinical.csv', key = 'patient_id')
hist_data <- copy(scores_patient_mut_all[cancer_type == 'LUAD' & grepl('PDAC Classical', meta_program)])[,
    c('pdac_high', 'hist_type') := .((score > 1.5), do.call(function(x) luad_clinical[x, histological_type], list(x = patient_id)))
][, mucinous := grepl('^mucinous| mucinous', hist_type)] # This regex ensures "nonmucinous" is not counted
setkey(hist_data, patient_id)

top_genes_int <- intersect(top_genes$LUAD, top_genes$PAAD)
cor_luad_paad <- cor(exp_data$LUAD$expmat[top_genes_int, ], exp_data$PAAD$expmat[top_genes_int, ])
bxplt_data <- data.table(patient_id = rownames(cor_luad_paad), mean_cor = rowMeans(cor_luad_paad))[, pdac_clust := (patient_id %in% pdac_clust)]
setkey(scores_patient_mut_all, patient_id)
bxplt_data[, pdac_high := do.call(function(x) scores_patient_mut_all[grepl('PDAC Classical', meta_program)][x, score > 1.5], list(x = patient_id))]
bxplt_data[, mucinous := hist_data$mucinous] # Checked that this gives TRUE: all(bxplt_data$patient_id == hist_data$patient_id)

box_luad <- ggplot(bxplt_data, aes(x = mapvalues(pdac_high, c(TRUE, FALSE), c('PDAC score\n\u2265 1.5', 'PDAC score\n< 1.5')), y = mean_cor)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(
        aes(fill = factor(mapvalues(mucinous, c(TRUE, FALSE), c('Mucinous', 'Non-mucinous')), levels = c('Non-mucinous', 'Mucinous'))),
        size = 2, shape = 21, stroke = 0.3, position = position_jitter(width = 0.3)
    ) +
    scale_fill_manual(values = c('Non-mucinous' = 'grey80', 'Mucinous' = 'deeppink2')) +
    theme_test() +
    theme(axis.text.x = element_text(size = 11)) +
    labs(x = NULL, y = 'Mean correlation with PDAC samples', title = NULL, fill = 'LUAD subtype')





box_hnsc_grob <- ggplotGrob(box_hnsc)
box_luad_grob <- ggplotGrob(box_luad)
widths_hnsc <- c(0.3, 0, 0.5, 0.4, 8, 0, 0, 0, 0.3)
widths_luad <- c(0.3, 0, 0.5, 0.7, 8, 0, 0, 0.3, 3, 0, 0.3)
heights_hnsc <- c(0.3, 0, 1, 0, 0, 0, 10, 0.8, 0, 0, 0, 0.3)
heights_luad <- c(0.3, 0, 1, 0, 0, 0, 10, 0.8, 0, 0, 0, 0.3)
box_hnsc_grob$widths <- unit(widths_hnsc, 'cm')
box_luad_grob$widths <- unit(widths_luad, 'cm')
box_hnsc_grob$heights <- unit(heights_hnsc, 'cm')
box_luad_grob$heights <- unit(heights_luad, 'cm')

cairo_pdf('../data/pdac_classical_main.pdf', height = sum(heights_luad)/2.54, width = sum(c(widths_luad, widths_hnsc))/2.54)
plot_grid(box_luad_grob, box_hnsc_grob, nrow = 1, ncol = 2, rel_widths = c(sum(widths_luad), sum(widths_hnsc)))
dev.off()
