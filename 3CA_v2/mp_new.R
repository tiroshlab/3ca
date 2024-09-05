library(data.table)
library(magrittr)
library(Matrix)
library(readxl)
library(stringr)
library(ggplot2)
library(scales)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', encoding = 'UTF-8', key = c('study', 'cancer_type'))

hgnc_complete_set <- fread('../data/hgnc_complete_set_2023-04-13.txt', key = 'symbol')
hgnc_complete_set <- hgnc_complete_set[!(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])]
alias_table <- make_alias_table(hgnc_complete_set)

ds <- as.data.table(read_xlsx('../3ca_manuscript/revision/dataset_summary.xlsx'))
setkey(ds, Study, Category)





custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))





d <- '~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF'





nmf_intersect_ccle <- readRDS(paste0(d, '/nmf_intersect_ccle_New_Batch.RDS')) # Jaccard similarity
Cluster_list <- readRDS(paste0(d, '/Cluster_list_Final_New_Batch.RDS'))
MP_list <- readRDS(paste0(d, '/MP_list_Final_New_Batch.RDS'))

nmf_ord <- rbindlist(lapply(names(Cluster_list), function(x) data.table(mp_name = x, nmf = Cluster_list[[x]])))
mp_htmp_data <- as.data.table(nmf_intersect_ccle[nmf_ord$nmf, nmf_ord$nmf], keep.rownames = 'mp1')
mp_htmp_data <- melt(mp_htmp_data, id.vars = 'mp1', variable.name = 'mp2', value.name = 'sim')
mp_htmp_data[, c('mp1', 'mp2') := .(factor(mp1, levels = nmf_ord$nmf), factor(mp2, levels = rev(nmf_ord$nmf)))]
mp_htmp_rect <- nmf_ord[, .N, by = mp_name][, .(mp_name, N, cumN = cumsum(N))][, c('mn', 'mx') := .(cumN - N + 0.5, cumN + 0.5)]

ggplot(mp_htmp_data) +
    geom_raster(aes(x = mp1, y = mp2, fill = 100*sim/(100 - sim))) +
    scale_fill_gradientn(name = 'Similarity (Jaccard index)', colours = custom_magma, limits = c(2, 25), oob = squish) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_rect(aes(xmin = mn, xmax = mx, ymin = mp_htmp_rect[, max(cumN) + 1 - mx], ymax = mp_htmp_rect[, max(cumN) + 1 - mn]),
        data = mp_htmp_rect, color = "red", linetype = "dashed", fill = NA, linewidth = 0.5) +
    theme(axis.ticks = element_blank(), panel.border = element_rect(fill = NA, colour = 'black'), panel.background = element_blank(),
        axis.line = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 16), legend.title = element_text(size = 16),
        legend.text = element_text(size = 14), legend.justification = 'right', legend.direction = 'horizontal', legend.position = 'bottom') +
    guides(fill = guide_colourbar(barheight = unit(0.5, 'cm'), barwidth = unit(4, 'cm'), frame.colour = 'black', ticks.colour = 'black')) +
    labs(x = paste0('NMF Programs (n=', mp_htmp_rect[, max(cumN)], ')'), y = 'NMF Programs')





Cancer_MPs <- as.data.table(read_xlsx(paste0(d, '/MP_by_cell_type/New/Cancer.xlsx')))
intersect_MPs_Jaccard <- apply(Cancer_MPs, 2, function(x) apply(MP, 2, function(y) length(intersect(x,y))/length(unique(c(x,y))))) 
intersect_MPs_Jaccard_meltI <- reshape2::melt(intersect_MPs_Jaccard) 

mp_int_data <- as.data.table(intersect_MPs_Jaccard_meltI)
setnames(mp_int_data, c('Var1', 'Var2', 'value'), c('mp_new', 'mp_old', 'sim'))
mp_int_data[, mp_new_n := as.numeric(gsub('MP', '', str_extract(mp_new, 'MP[0-9]+')))]
mp_int_data[, mp_old_n := as.numeric(gsub('MP', '', str_extract(mp_old, 'MP[0-9]+')))]
mp_int_data[, mp_new_name := gsub('Cylce', 'Cycle', gsub('^ ', '', gsub('MP[0-9]+ ', '', mp_new)))]
mp_int_data[, mp_old_name := gsub('EMT ', 'EMT-', gsub('Cylce', 'Cycle', gsub('^ ', '', gsub('MP[0-9]+ ', '', mp_old))))]
mp_int_data[, mp_new_name := factor(mp_new_name, levels = unique(mp_new_name[order(mp_new_n)]))]
mp_int_data[, mp_old_name := factor(mp_old_name, levels = unique(mp_old_name[order(mp_old_n)]))]

mp_match <- mp_int_data[,
    .(mp_old_name = if(max(sim) > 0.3 & all(sim[sim < max(sim)] < 0.2)) as.character(mp_old_name)[sim == max(sim)] else as.character(NA)),
    by = .(mp_new_name, mp_new_n)
]
mp_int_data[, mp_new_name := factor(mp_new_name, levels = mp_match[order(is.na(mp_old_name)), mp_new_name])]
mp_int_data[, mp_old_name := factor(mp_old_name, levels = mp_match[!is.na(mp_old_name), c(mp_old_name, 'PDAC-related')])]

ggplot(data = mp_int_data, aes(x = mp_new_name, y = mp_old_name, fill = 100*sim)) +
    geom_tile(colour = 'grey90') +
    scale_fill_gradientn(name = 'Similarity (Jaccard index)', colours = custom_magma, limits = c(2, 25), oob = squish) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(panel.border = element_rect(fill = NA), axis.ticks.length = unit(1.5, 'pt'),
        axis.text.x = element_text(size = 10.5, angle = 55, hjust = 1), axis.text.y = element_text(size = 10.5),
        axis.title = element_text(size = 16), legend.title = element_text(size = 16), legend.text = element_text(size = 14),
        legend.justification = 'bottom') +
    guides(fill = guide_colourbar(barheight = unit(4, 'cm'), barwidth = unit(0.5, 'cm'))) +
    labs(x = 'New meta-programs', y = 'Previous meta-programs')





nuc_seq_prop <- fread('../3ca_manuscript/revision/nuc_seq_porportion.csv')

ggplot(nuc_seq_prop, aes(x = total_samples_num - nuc_seq_num, y = nuc_seq_num, label = MP)) +
    geom_point(shape = 24, size = 3, stroke = 1, fill = 'plum1') +
    scale_x_continuous(limits = c(-50, 1200), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-3, 60), expand = c(0, 0)) +
    geom_text_repel(data = nuc_seq_prop[MP == 'Unassigned 2'], nudge_x = 20, nudge_y = 4, size = 5) +
    geom_text_repel(data = nuc_seq_prop[MP == 'Cilia 2'], nudge_x = 20, nudge_y = 3, size = 5) +
    geom_text_repel(data = nuc_seq_prop[MP == 'Cell Cycle - G1/S'], nudge_x = 20, nudge_y = 5, size = 5) +
    geom_text_repel(data = nuc_seq_prop[MP == 'Stress 1'], nudge_x = 20, nudge_y = 3, size = 5) +
    geom_text_repel(data = nuc_seq_prop[MP == 'Cell cycle single-nucleus'], nudge_x = 200, nudge_y = -4, size = 5) +
    geom_text_repel(data = nuc_seq_prop[MP == 'Cell Cycle - G2/M'], nudge_y = 4, size = 5) +
    geom_text_repel(data = nuc_seq_prop[MP == 'Glioma single-nucleus'], nudge_x = 260, nudge_y = 1, size = 5) +
    theme_half_open() +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
    labs(x = 'Number of single-cell samples', y = 'Number of single-nucleus samples')





rep_smpl <- data.table(
    stdy = c('Yue et al. 2020', 'Hwang et al. 2022'),
    smpl = c('91610', 'U2'),
    tech = c('single-cell', 'single-nucleus')
)
rep_smpl_paths <- slapply(rep_smpl$smpl, function(s) {
    sdt <- rep_smpl[smpl == s]
    paths_table[sdt$stdy, if(s %in% fread(cells)$sample) .SD[, .(genes, expmat, study, tech = sdt$tech)], by = cells]
})
set.seed(3320)
pdata_cilia <- slapply(names(rep_smpl_paths), function(s) {
    p <- rep_smpl_paths[[s]]
    cat(p$study, '\n')
    cells <- suppressWarnings(fread(p$cells, na.strings = '', colClasses = c(cell_name = 'character', sample = 'character')))
    if(cells[, .N > length(unique(cell_name))]) cells[, cell_name := paste(cell_name, .I, sep = '_')]
    genes <- fread(p$genes, header = FALSE)$V1
    expmat <- readMM(p$expmat)
    expmat <- expmat[genes != '', ]
    genes <- genes[genes != '']
    expmat <- expmat[genes %in% names(table(genes))[table(genes) == 1], ]
    genes <- genes[genes %in% names(table(genes))[table(genes) == 1]]
    genes <- update_symbols_fast(genes, alias_table)
    rownames(expmat) <- genes
    colnames(expmat) <- cells$cell_name
    cells <- cells[col_nnz(expmat) >= 1000 & cell_type == 'Malignant' & sample == s]
    cells <- cells[, if(.N > 350) .SD[sample(1:.N, 300, replace = FALSE)] else .SD]
    expmat <- expmat[, cells$cell_name]
    expmat <- round(log_transform(median(colSums(expmat))*to_frac(expmat)), 4)
    mp1 <- MP_list$`Cilia 1`; mp1 <- mp1[mp1 %in% genes]
    mp2 <- MP_list$`Cilia 2`; mp2 <- mp2[mp2 %in% genes]
    if(p$tech == 'single-cell') mp_score <- mp1 else mp_score <- mp2
    scores <- cells[, .(cell_name, score = sig_score(expmat[, cell_name], mp_score, nbin = 100, n = 50))]
    setkey(scores, cell_name)
    expmat_sub <- expmat[union(mp1, mp2), ]
    expmat_cent <- as.matrix(expmat_sub - rowMeans(expmat_sub))
    score_corr <- scores[, .(symbol = mp_score, corr = cor(t(expmat_cent[mp_score, cell_name]), score)[, 1])]
    setkey(score_corr, symbol)
    out <- as.data.table(expmat_cent, keep.rownames = 'gene')
    out <- melt(out, id.vars = 'gene', variable.name = 'cell', value.name = 'exp_level', variable.factor = FALSE)
    setkey(cells, cell_name)
    out[, c('study', 'tech', 'sample_name') := c(p[, .(study, tech)], list(s))]
    out[, c('score', 'score_corr') := .(scores[cell, score], score_corr[gene, corr])]
    setcolorder(out, c('study', 'tech', 'sample_name', 'cell', 'gene'))
    return(out)
}) %>% rbindlist

pdata_cilia_final <- copy(pdata_cilia)
pdata_cilia_final <- pdata_cilia_final[gene %in% pdata_cilia_final[, .(n = length(unique(sample_name))), by = gene][n == 2, gene]]
pdata_cilia_final <- pdata_cilia_final[gene %in% pdata_cilia_final[!is.na(score_corr), unique(gene)]]
pdata_cilia_final[, gene_cat := ifelse(
    gene %in% intersect(MP_list$`Cilia 1`, MP_list$`Cilia 2`),
    'Shared',
    ifelse(gene %in% MP_list$`Cilia 1`, 'Single-cell signature', 'Single-nucleus signature'))
]
pdata_cilia_final[, gene := factor(
    gene,
    levels = unique(.SD[!is.na(score_corr), .(gene, score_corr)])[, .(score_corr = mean(score_corr)), by = gene][order(score_corr), gene]
)]
pdata_cilia_final[, cell := factor(cell, levels = unique(.SD)[order(-score), cell]), .SDcols = c('cell', 'score')]
pdata_cilia_final[, gene_cat := factor(gene_cat, levels = c('Single-cell signature', 'Shared', 'Single-nucleus signature'))]
pdata_cilia_final[, sample_name_ext := paste0(gsub(' .+$', '', study), ' - ', sample_name, '\n(', tech, ')')]
pdata_cilia_final[, sample_name_ext := factor(sample_name_ext, levels = c('Yue - 91610\n(single-cell)', 'Hwang - U2\n(single-nucleus)'))]
pdata_cilia_final <- pdata_cilia_final[ # Filter genes
    gene_cat == 'Shared' |
        gene %in% unique(pdata_cilia_final[!is.na(score_corr) & gene_cat != 'Shared', .(gene, score_corr)])[score_corr > 0.2, gene]
]

ggplot(pdata_cilia_final) +
    geom_raster(aes(x = cell, y = gene, fill = exp_level)) +
    facet_grid(cols = vars(sample_name_ext), rows = vars(gene_cat), scales = 'free', space = 'free_y') +
    scale_fill_gradientn(colours = rev(brewer.pal(9, 'RdBu')), limits = c(-3, 3), breaks = c(-3, 0, 3), oob = squish) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, margin = margin(r = 7)), axis.ticks.x = element_blank(), legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, margin = margin(b = 7)), panel.border = element_rect(fill = NA),
        plot.title = element_text(size = 16), strip.text = element_text(size = 14), strip.background = element_rect(fill = NA, colour = NA),
        strip.clip = 'off') +
    guides(fill = guide_colourbar(barheight = unit(0.5, 'cm'), barwidth = unit(3, 'cm'), frame.colour = 'black', ticks.colour = 'black',
        title.position = 'right')) +
    labs(x = 'Cells', y = 'Gene', fill = 'Relative expression level (log2 ratio)', title = 'Cilia')





rep_smpl <- data.table(
    stdy = c('Griffiths et al. 2021', 'Pal et al. 2021'),
    smpl = c('P21_S', 'ER_positive_0167'),
    tech = c('single-nucleus', 'single-cell')
)
rep_smpl_paths <- slapply(rep_smpl$smpl, function(s) {
    sdt <- rep_smpl[smpl == s]
    paths_table[sdt$stdy, if(s %in% fread(cells)$sample) .SD[, .(genes, expmat, study, tech = sdt$tech)], by = cells]
})
set.seed(3400)
pdata_cc <- slapply(names(rep_smpl_paths), function(s) {
    p <- rep_smpl_paths[[s]]
    cat(p$study, '\n')
    cells <- suppressWarnings(fread(p$cells, na.strings = '', colClasses = c(cell_name = 'character', sample = 'character')))
    if(cells[, .N > length(unique(cell_name))]) cells[, cell_name := paste(cell_name, .I, sep = '_')]
    genes <- fread(p$genes, header = FALSE)$V1
    expmat <- readMM(p$expmat)
    expmat <- expmat[genes != '', ]
    genes <- genes[genes != '']
    expmat <- expmat[genes %in% names(table(genes))[table(genes) == 1], ]
    genes <- genes[genes %in% names(table(genes))[table(genes) == 1]]
    genes <- update_symbols_fast(genes, alias_table)
    rownames(expmat) <- genes
    colnames(expmat) <- cells$cell_name
    cells <- cells[col_nnz(expmat) >= 1000 & cell_type == 'Malignant' & sample == s]
    cells <- cells[, if(.N > 350) .SD[sample(1:.N, 300, replace = FALSE)] else .SD]
    expmat <- expmat[, cells$cell_name]
    expmat <- round(log_transform(median(colSums(expmat))*to_frac(expmat)), 4)
    mp1 <- union(MP_list$`Cell Cycle - G2/M`, MP_list$`Cell Cycle - G1/S`); mp1 <- mp1[mp1 %in% genes]
    mp2 <- MP_list$`Cell cycle single-nucleus`; mp2 <- mp2[mp2 %in% genes]
    if(p$tech == 'single-cell') mp_score <- mp1 else mp_score <- mp2
    scores <- cells[, .(cell_name, score = sig_score(expmat[, cell_name], mp_score, nbin = 100, n = 50))]
    setkey(scores, cell_name)
    expmat_sub <- expmat[union(mp1, mp2), ]
    expmat_cent <- as.matrix(expmat_sub - rowMeans(expmat_sub))
    score_corr <- scores[, .(symbol = mp_score, corr = cor(t(expmat_cent[mp_score, cell_name]), score)[, 1])]
    setkey(score_corr, symbol)
    out <- as.data.table(expmat_cent, keep.rownames = 'gene')
    out <- melt(out, id.vars = 'gene', variable.name = 'cell', value.name = 'exp_level', variable.factor = FALSE)
    setkey(cells, cell_name)
    out[, c('study', 'tech', 'sample_name') := c(p[, .(study, tech)], list(s))]
    out[, c('score', 'score_corr') := .(scores[cell, score], score_corr[gene, corr])]
    setcolorder(out, c('study', 'tech', 'sample_name', 'cell', 'gene'))
    return(out)
}) %>% rbindlist

pdata_cc_final <- copy(pdata_cc)
pdata_cc_final <- pdata_cc_final[gene %in% pdata_cc_final[, .(n = length(unique(sample_name))), by = gene][n == 2, gene]]
pdata_cc_final <- pdata_cc_final[gene %in% pdata_cc_final[!is.na(score_corr), unique(gene)]]
pdata_cc_final[, gene_cat := ifelse(
    gene %in% intersect(union(MP_list$`Cell Cycle - G2/M`, MP_list$`Cell Cycle - G1/S`), MP_list$`Cell cycle single-nucleus`),
    'Shared',
    ifelse(gene %in% MP_list$`Cell cycle single-nucleus`, 'Single-nucleus signature', 'Single-cell signature'))
]
pdata_cc_final[, gene := factor(
    gene,
    levels = unique(.SD[!is.na(score_corr), .(gene, score_corr)])[, .(score_corr = mean(score_corr)), by = gene][order(score_corr), gene]
)]
pdata_cc_final[, cell := factor(cell, levels = unique(.SD)[order(-score), cell]), .SDcols = c('cell', 'score')]
pdata_cc_final[, gene_cat := factor(gene_cat, levels = c('Single-cell signature', 'Shared', 'Single-nucleus signature'))]
pdata_cc_final[, sample_name_ext := paste0(gsub(' .+$', '', study), ' - ', sample_name, '\n(', tech, ')')]
pdata_cc_final[, sample_name_ext := factor(
    sample_name_ext,
    levels = c('Pal - ER_positive_0167\n(single-cell)', 'Griffiths - P21_S\n(single-nucleus)')
)]
pdata_cc_final <- pdata_cc_final[ # Filter genes
    gene_cat == 'Shared' | gene %in% unique(pdata_cc_final[!is.na(score_corr) & gene_cat != 'Shared', .(gene, score_corr, gene_cat)])[
        score_corr > ifelse(gene_cat == 'Single-cell signature', 0.55, 0.45),
        gene
    ]
]

ggplot(pdata_cc_final) +
    geom_raster(aes(x = cell, y = gene, fill = exp_level)) +
    facet_grid(cols = vars(sample_name_ext), rows = vars(gene_cat), scales = 'free', space = 'free_y') +
    scale_fill_gradientn(colours = rev(brewer.pal(9, 'RdBu')), limits = c(-3, 3), breaks = c(-3, 0, 3), oob = squish) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, margin = margin(r = 7)), axis.ticks.x = element_blank(), legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, margin = margin(b = 7)), panel.border = element_rect(fill = NA),
        plot.title = element_text(size = 16), strip.text = element_text(size = 14), strip.background = element_rect(fill = NA, colour = NA),
        strip.clip = 'off') +
    guides(fill = guide_colourbar(barheight = unit(0.5, 'cm'), barwidth = unit(3, 'cm'), frame.colour = 'black', ticks.colour = 'black',
        title.position = 'right')) +
    labs(x = 'Cells', y = 'Gene', fill = 'Relative expression level (log2 ratio)', title = 'Cell cycle')





p <- paths_table[grep('Pelka', study), if('C161_T' %in% fread(cells)$sample) .SD, by = cells]
cells <- suppressWarnings(fread(p$cells, na.strings = '', colClasses = c(cell_name = 'character', sample = 'character')))
genes <- fread(p$genes, header = FALSE)$V1
expmat <- readMM(p$expmat)
expmat <- expmat[genes != '', ]
genes <- genes[genes != '']
expmat <- expmat[genes %in% names(table(genes))[table(genes) == 1], ]
genes <- genes[genes %in% names(table(genes))[table(genes) == 1]]
genes <- update_symbols_fast(genes, alias_table)
rownames(expmat) <- genes
colnames(expmat) <- cells$cell_name
cells <- cells[sample == 'C161_T' & cell_type == 'Malignant' & col_nnz(expmat) >= 1000][sample(1:.N, 500)] # Downsample to 500 cells
expmat <- expmat[, cells$cell_name]
expmat <- round(log_transform(median(colSums(expmat))*to_frac(expmat)), 4)
mp46 <- MP[`MP46 CRC stemness` %in% genes, `MP46 CRC stemness`]
set.seed(5528)
scores_mp46 <- sig_score(expmat, mp46, nbin = 100, n = 100)
expmat_sub <- expmat[mp46, ]
expmat_cent <- as.matrix(expmat_sub - rowMeans(expmat_sub))
score_corr <- cor(t(expmat_cent), scores_mp46)[, 1]
gene_cond <- rowMeans(expmat_sub) > 0.1 & score_corr > 0.1
pdata_mp46 <- as.data.table(expmat_cent[gene_cond, ], keep.rownames = 'gene')
pdata_mp46 <- melt(pdata_mp46, id.vars = 'gene', variable.name = 'cell', value.name = 'exp_level')
pdata_mp46[, cell := factor(cell, levels = names(sort(scores_mp46, decreasing = TRUE)))]
pdata_mp46[, gene := factor(gene, levels = names(sort(score_corr[gene_cond])))]

ggplot(pdata_mp46) +
    geom_raster(aes(x = cell, y = gene, fill = exp_level)) +
    scale_fill_gradientn(colours = rev(brewer.pal(9, 'RdBu')), limits = c(-3, 3), breaks = c(-3, 0, 3), oob = squish) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, margin = margin(r = 7)),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, margin = margin(b = 7)),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(size = 16)
    ) +
    guides(fill = guide_colourbar(barheight = unit(3, 'cm'), barwidth = unit(0.5, 'cm'), frame.colour = 'black', ticks.colour = 'black')) +
    labs(x = 'Cells', y = 'Gene', fill = 'Relative\nexpression\nlevel\n(log2 ratio)',
        title = 'CRC stemness - Pelka - C161_T')
