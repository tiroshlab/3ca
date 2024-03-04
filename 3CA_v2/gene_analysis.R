library(data.table)
library(magrittr)
library(ggplot2)
library(plyr)
library(stringr)
library(Matrix)
library(RColorBrewer)
library(scales)
library(cowplot)
library(randomcoloR)
library(ggrepel)
library(seriation)
library(ggtext)
library(ggpubr)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'), encoding = 'UTF-8')

hgnc_complete_set <- fread('../data/hgnc_complete_set_2023-04-13.txt', key = 'symbol')
hgnc_complete_set <- hgnc_complete_set[!(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])]
hgnc_complete_set <- hgnc_complete_set[locus_group == 'protein-coding gene']
alias_table <- make_alias_table(hgnc_complete_set)





diseases <- lapply(transpose(as.list(unique(paths_table[cancer_type != 'Other/Models', .(study, cancer_type)]))), function(r) {
    cat(r, '\n')
    samples_path <- paste0('~/../shared/pan_cancer_datasets/', paths_table[as.list(r), directory[1]], '/samples.csv')
    if(file.exists(samples_path)) {
        samples <- fread(samples_path)[!is.na(sample)]
        samples[, sample := as.character(sample)]
        samples[cancer_type == '', cancer_type := NA]
        setkey(samples, sample)
    } else {warning("Samples file doesn't exist for ", r[1], ", ", r[2]); return(NULL)}
    if('cancer_type' %in% names(samples)) {
        return(samples[, .(study = r[1], cancer_type = r[2], sample = sample, disease = cancer_type)])
    } else {warning("Samples file for ", r[1], ", ", r[2], " doesn't contain cancer_type column"); return(NULL)}
}) %>% rbindlist

to_exclude <- list(c('Chen et al. 2020', 'Head and Neck'), c('Sun et al. 2021', 'Liver/Biliary'), c('Li et al. 2017', 'Colorectal')) %>%
    transpose %>% as.data.table %>% setNames(c('study', 'cancer_type'))

setkey(diseases, study, cancer_type)
diseases_to_include <- diseases[!to_exclude][
    !is.na(disease) & !(disease %in% c('Normal', 'Premalignant')),
    .SD[, .(n_sample = .N, n_study = nrow(unique(.SD))), by = disease, .SDcols = c('study', 'cancer_type')][
        (n_study >= 2 & n_sample >= 10) | (n_sample >= 20),
        disease
    ]
]
to_include <- diseases[!to_exclude][disease %in% diseases_to_include, .(study, cancer_type)] %>% unique

cell_types <- readRDS('../data/gene_plots_cell_types.rds')
cell_types <- cell_types[!(cell_types %in% c('Pericyte', 'Oligodendrocyte'))] # Remove pericytes & oligodendrocytes - too cancer-type-specific

data_genes <- lapply(transpose(as.list(to_include)), function(r) {
    
    cat(r, '\n')
    if(!('gene_ave.csv' %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])))) return(NULL)
    
    gene_ave <- fread(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/gene_ave.csv'), colClasses = c(sample = 'character'))
    gene_ave <- gene_ave[n_cell >= 10 & cell_type %in% cell_types & symbol %in% hgnc_complete_set$symbol]
    if(gene_ave[, unique(.SD)[, sum(n_cell)], .SDcols = c('cell_type', 'sample', 'group', 'n_cell')] < 100) return(NULL)
    
    samples_path <- paste0('~/../shared/pan_cancer_datasets/', paths_table[as.list(r), directory[1]], '/samples.csv')
    if(file.exists(samples_path)) {
        samples <- fread(samples_path, colClasses = c(sample = 'character'), na.strings = '')[!is.na(sample) & cancer_type %in% diseases_to_include]
        setkey(samples, sample)
    }
    
    # Split HNSCC samples into HPV pos and neg:
    if(all(r == c('Puram et al. 2017', 'Head and Neck'))) {
        samples <- samples[sample %in% c('5', '6', '16', '17', '18', '20', '22', '25', '26', '28')][, cancer_type := 'HNSCC (HPV-neg.)']
    } else if(all(r == c('Cillo et al. 2020', 'Head and Neck'))) {
        samples[, cancer_type := ifelse(p16_status == 'p16+', 'HNSCC (HPV-pos.)', 'HNSCC (HPV-neg.)')]
    } else if(all(r == c('Kürten et al. 2021', 'Head and Neck'))) {
        samples[, cancer_type := ifelse(additional_tumor_characterisics == 'HPV Positive', 'HNSCC (HPV-pos.)', 'HNSCC (HPV-neg.)')]
    } else if(all(r == c('Unpublished HNSCC HTAN', 'Head and Neck'))) {
        samples[, cancer_type := ifelse(hpv == 'HPV+', 'HNSCC (HPV-pos.)', 'HNSCC (HPV-neg.)')]
    }
    
    # Split CRC and Gastric into MSS and MSI:
    if(all(r == c('Lee et al. 2020', 'Colorectal'))) {
        samples[, cancer_type := ifelse(grepl('mss', genetic_hormonal_features), 'CRC (MSS)', 'CRC (MSI)')]
    } else if(all(r == c('Pelka et al. 2021', 'Colorectal'))) {
        samples[, cancer_type := ifelse(MSIStatus == 'MSS', 'CRC (MSS)', 'CRC (MSI)')]
    } else if(all(r == c('Zhang et al. 2018', 'Colorectal'))) {
        samples[, cancer_type := ifelse(msi_status == 'MSS', 'CRC (MSS)', 'CRC (MSI)')]
    } else if(all(r == c('Kumar et al. 2022', 'Gastric'))) {
        samples[, cancer_type := ifelse(is.na(subtype) | subtype == '', NA, ifelse(subtype == 'MSI', 'Gastric (MSI)', 'Gastric (MSS)'))]
    } else if(all(r == c('Zhang et al. 2019', 'Gastric'))) {
        samples[, cancer_type := NA] # There's only 1 true gastric cancer sample in this dataset, and it's not annotated with MSI status.
    }
    
    samples <- samples[!is.na(cancer_type)]; if(nrow(samples) == 0) return(NULL)
    
    gene_ave <- gene_ave[sample %in% samples$sample] # Could put n >= 50 here, but we might want to take averages per study first
    if(gene_ave[, unique(.SD)[, sum(n_cell)], .SDcols = c('cell_type', 'sample', 'group', 'n_cell')] < 100) return(NULL)
    gene_ave[, c('disease', 'tech') := do.call(`[`, list(samples, sample))[, .(cancer_type, technology)]]
    gene_ave[, c('study', 'cancer_type') := as.list(r)]
    
    return(gene_ave)
    
}) %>% rbindlist

setnames(data_genes, 'n_cell', 'n')

data_genes[
    cell_type %in% c('Macrophage', 'Myeloid', 'Monocyte'),
    c('cell_type', 'n', 'ave', 'prop_pos') := .('Macrophage', sum(n), sum(ave*n)/sum(n), sum(prop_pos*n)/sum(n)),
    by = .(disease, study, tech, sample, symbol)
]
data_genes <- unique(data_genes)

data_genes <- data_genes[!(cancer_type == 'Brain' & cell_type == 'Fibroblast')]

# Retain genes that have at least one value in at least 20 cancer types:
data_genes <- data_genes[symbol %in% data_genes[, .(n = length(unique(disease))), by = symbol][n >= 20, symbol]]

# Edit disease names:
data_genes[, disease := mapvalues(
    gsub(' Cancer', '', disease),
    c('Acute Myeloid Leukemia', 'Chronic Myeloid Leukemia', 'Clear Cell Renal Cell Carcinoma', 'Colorectal', 'Cutaneous Basal Cell Carcinoma',
        'Cutaneous Squamous Cell Carcinoma', 'Diffuse Large B Cell Lymphoma', 'Glioblastoma', 'Hepatocellular Carcinoma', 'Lung Adenocarcinoma',
        'Lung Squamous Cell Carcinoma', 'Multiple Myeloma', 'Neuroendocrine Tumor', 'Pancreatic Ductal Adenocarcinoma', 'Small Cell Lung'),
    c('AML', 'CML', 'ccRCC', 'CRC', 'Skin BCC', 'Skin SCC', 'DLBCL', 'GBM', 'HCC', 'Lung Adeno.', 'Lung Squamous', 'MM', 'NET', 'PDAC', 'SCLC'),
    warn_missing = FALSE
)]





ave_data <- data_genes[,
    .(ave = sum(ave*n)/sum(n), prop_pos = sum(prop_pos*n)/sum(n), n_sample = .N),
    by = .(symbol, cell_type, disease, study, tech) # Take means across samples within each study, disease and tech
][,
    .(m = mean(ave), p = mean(prop_pos), n_sample = sum(n_sample)),
    by = .(symbol, cell_type, disease, study) # First take mean across datasets of the same study and disease with different tech values
][,
    .(m = weighted.mean(m, n_sample), p = weighted.mean(p, n_sample)),
    by = .(symbol, cell_type, disease) # Weighted mean across studies of the same disease
]

ave_data[, cell_type := gsub('_', ' ', cell_type)]
ave_data[, cell_type := factor(cell_type, levels = c('Malignant', 'Macrophage', 'T cell', 'NK cell', 'B cell', 'Plasma',
    'Dendritic', 'Mast', 'Fibroblast', 'Endothelial', 'Epithelial'))]





ss_data <- ave_data[, .(m = mean(m), p = mean(p), v = var(m), med = median(m), iqr = IQR(m)), by = .(symbol, cell_type)]
ss_data[, spec := ifelse(med == max(med) & med > 0, 1/(m[order(-m)][2] + 1), as.numeric(NA)), by = symbol]
ss_data[, sens := ifelse(med == max(med) & med > 0, med, as.numeric(NA)), by = symbol]
ss_data[!is.na(sens), sens := sens/max(sens)]

# Facetted plot for all cell types:

labs_manual <- data.table(
    symbol = c(
        'MCM7', 'NME1', # Malignant
        'C1QC', 'C1QA', 'AIF1', 'LYZ', # Macrophage
        'CD3D', 'CD3G', 'IL32', 'CD2', # T cell
        'KLRF1', 'GNLY', 'NKG7', 'KLRD1', # NK cell
        'MS4A1', 'BANK1', 'CD79A', 'CD79B', # B cell
        'JCHAIN', 'MZB1', 'TNFRSF17', 'SPAG4', # Plasma
        'LGALS2', 'FLT3', # Dendritic
        'TPSB2', 'CPA3', 'MS4A2', 'HDC', # Mast
        'COL1A2', 'DCN', 'LUM', 'PDGFRB', # Fibroblast
        'RAMP2', 'VWF', 'IFI27', 'CLEC14A', # Endothelial
        'KRT7', 'SLPI', 'KRT19', 'PIGR' # Epithelial
    ),
    cell_type = c(rep('Malignant', 2), rep('Macrophage', 4), rep('T cell', 4), rep('NK cell', 4), rep('B cell', 4), rep('Plasma', 4),
        rep('Dendritic', 2), rep('Mast', 4), rep('Fibroblast', 4), rep('Endothelial', 4), rep('Epithelial', 4))
)
setkey(ss_data, symbol, cell_type)
ss_data[labs_manual, lab := symbol]
ss_data[, cell_type := factor(cell_type, levels = c('Malignant', 'Macrophage', 'T cell', 'NK cell', 'B cell', 'Plasma',
    'Dendritic', 'Mast', 'Fibroblast', 'Endothelial', 'Epithelial'))]

ggplot(ss_data[!is.na(sens)]) +
    geom_point(aes(x = spec, y = sens), colour = 'grey20', alpha = 0.5) +
    facet_wrap(vars(cell_type), nrow = 4, ncol = 3) +
    labs(x = 'Specificity', y = 'Sensitivity') +
    geom_text_repel(aes(x = spec, y = sens, label = lab), size = 4.5, nudge_x = 0.05, nudge_y = 0.1) +
    scale_x_continuous(breaks = c(0.2, 0.6, 1), labels = c('0.2', '0.6', '1')) +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
    theme_test() +
    theme(
        axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, margin = margin(t = 10)),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = 'aliceblue')
    )

# Dot plot of marker genes:

marker_data <- ss_data[symbol %in% labs_manual$symbol]
ct_ord <- c('Mast', 'Fibroblast', 'Endothelial', 'NK cell', 'T cell', 'B cell', 'Plasma', 'Macrophage', 'Dendritic', 'Epithelial', 'Malignant')
marker_data[, cell_type := factor(cell_type, levels = ct_ord)]
marker_data[, symbol := factor(symbol, levels = marker_data[!is.na(lab)][order(-cell_type, sens), symbol])]

ggplot(marker_data, aes(x = cell_type, y = symbol)) +
    geom_point(aes(fill = m, size = 100*p), shape = 21, colour = 'black') +
    scale_fill_gradientn(colours = brewer.pal(9, 'YlOrRd'), limits = c(0, 8), oob = squish) +
    scale_radius(range = c(1, 7), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
        labels = c('0' = '0', '25' = '25', '50' = '50', '75' = '75', '100' = '100')) +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1, size = 16),
        axis.text.y = element_markdown(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.grid = element_line(linewidth = 0.4, colour = 'grey92'),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = 'black'),
        legend.key = element_rect(fill = NA)
    ) +
    guides(fill = guide_colourbar(frame.colour = 'black', ticks.colour = 'black', order = 1), size = guide_legend(order = 2)) +
    labs(x = 'Cell type', y = 'Gene', fill = 'Mean', size = '% expressing\ncells')





# Markers for malignant cells within each cancer type:

ss_data_m <- ave_data[
    !(disease %in% c('Lung Squamous', 'HCC', 'CML')), # Remove carcinomas without normal epithelial cells, and CML, which only has malignant cells
    if(.N >= 5 && 'Malignant' %in% cell_type) {
        mm <- m[cell_type == 'Malignant'] # Define sens and spec only if this gene's expression in highest in malignant cells:
        if(mm > 0 & mm == max(m)) .(spec = 1/(max(m[cell_type != 'Malignant']) + 1), sens = mm)
    },
    by = .(symbol, disease)
]
ss_data_m[, sens := sens/max(sens)]

# Facetted plot for all cancer types:

labs_manual_m <- data.table(
    symbol = c(
        'EIF4A1', 'CAT', 'SPINK2', 'NPR3', # AML
        'TRPS1', 'ANKRD30A', 'ESR1', 'BMPR1B', # Breast
        'NNMT', 'FABP7', 'CP', 'PNCK', # ccRCC
        'AREG', 'GDF15', 'DPEP1', 'KRT23', # CRC (MSS)
        'HMGA1', 'IRAG2', 'GTSF1', 'RASSF6', # DLBCL
        'KRT5', 'CXCL14', 'TP63', 'KLK5', # HNSCC (HPV-neg.)
        'KRT17', 'CDKN2A', 'CALML5', 'LY6K', # HNSCC (HPV-pos.)
        'NAPSA', 'CEACAM6', 'ABCC3', 'CEACAM5', # Lung Adeno.
        'PMEL', 'TYR', 'MIA', 'CDH19', # Melanoma
        'GNAS', 'NME2', 'IRAK1', 'RRAGD', # MM
        'STMN2', 'TUBB2B', 'ELAVL4', 'MARCHF11', # Neuroblastoma
        'CD9', 'S100A1', 'CAPS', 'RAB25', # Ovarian
        'CD24', 'AHNAK2', 'ANKRD36B', 'MUC16', # PDAC
        'MESP1', 'MBOAT2', 'PCSK1N', 'SERPINB11', # Prostate
        'STMN1', 'HES6', 'CHGA', 'SCG3', # SCLC
        'BCAM', 'MXRA8', 'WT1', 'CRB2' # Wilms Tumor
    ),
    disease = lapply(ss_data_m[, sort(unique(disease))], rep, 4) %>% unlist
)
setkey(ss_data_m, symbol, disease)
ss_data_m[labs_manual_m, lab := symbol]
ss_data_m[, disease := factor(disease, levels = sort(unique(disease)))]

ggplot(ss_data_m) +
    geom_point(aes(x = spec, y = sens), colour = 'grey20', alpha = 0.5) +
    facet_wrap(vars(disease), nrow = 4, ncol = 4) +
    labs(x = 'Specificity', y = 'Sensitivity') +
    geom_text_repel(aes(x = spec, y = sens, label = lab), size = 4.5, nudge_x = 0.05, nudge_y = 0.1) +
    scale_x_continuous(breaks = c(0.2, 0.6, 1), labels = c('0.2', '0.6', '1')) +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
    theme_test() +
    theme(
        axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, margin = margin(t = 10)),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = 'aliceblue')
    )

# Dot plot for malignant cell markers:

genes_dot_m <- data.table(
    disease = ss_data_m[, sort(unique(disease))],
    symbol = c('SPINK2', 'ESR1', 'FABP7', 'DPEP1', 'GTSF1', 'TP63', 'CDKN2A', 'CEACAM5', 'PMEL', 'IRAK1', 'ELAVL4', 'S100A1', 'CD24', 'PCSK1N',
        'CHGA', 'WT1')
)
mdata <- copy(ave_data)
setkey(mdata, disease, symbol)
mdata <- mdata[genes_dot_m]
mdata[, disease := factor(disease, levels = rev(genes_dot_m$disease))]
mdata[, disease_num := as.numeric(disease)]
setkey(genes_dot_m, disease)

ggplot(mdata, aes(x = cell_type, y = disease_num)) +
    geom_point(aes(fill = m, size = 100*p), shape = 21, colour = 'black') +
    scale_fill_gradientn(colours = brewer.pal(9, 'YlOrRd'), limits = c(0, 8), oob = squish) +
    scale_radius(range = c(1, 7), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
        labels = c('0' = '0', '25' = '25', '50' = '50', '75' = '75', '100' = '100')) +
    scale_x_discrete(expand = c(0, 0.8)) +
    scale_y_continuous(name = 'Cancer type', breaks = mdata[, unique(disease_num)], labels = mdata[, unique(disease)],
        sec.axis = dup_axis(name = 'Gene', labels = genes_dot_m[levels(mdata$disease), rev(symbol)]), # No idea why we need to use rev() on symbol!
        expand = c(0, 0.8)) +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y.left = element_text(size = 16, margin = margin(r = 10)),
        axis.title.y.right = element_text(size = 16, margin = margin (l = 10)),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = 'black'),
        panel.grid.major = element_line(linewidth = 0.4, colour = 'grey92'),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
    ) +
    guides(fill = guide_colourbar(frame.colour = 'black', ticks.colour = 'black', order = 1), size = guide_legend(order = 2)) +
    labs(x = 'Cell type', fill = 'Mean', size = '% expressing\ncells')





# Cancer-type-specific genes:

ss_data_v <- ave_data[cell_type == 'Malignant']
ss_data_v[, spec := ifelse(m == max(m) & m > 0, 1/(m[order(-m)][2] + 1), as.numeric(NA)), by = symbol]
ss_data_v[, sens := ifelse(m == max(m) & m > 0, m, as.numeric(NA)), by = symbol]
ss_data_v[!is.na(sens), sens := sens/max(sens)]

# Facetted plot for all cancer types:

labs_manual_v <- data.table(
    symbol = c(
        'PCNP', 'CFD', 'SPN', 'CLEC12A', # AML
        'TRPS1', 'ANKRD30A', 'ESR1', 'PRLR', # Breast
        'CRYAB', 'NDUFA4L2', 'ANGPTL4', 'SLC3A1', # ccRCC
        'TESC', 'CDA', 'CXCL6', 'NPSR1', # Cholangiocarcinoma
        'H3-5', 'MPL', 'PPIAL4G', 'CCDC27', # CML
        'AREG', 'CXCL3', 'DPEP1', 'CDX1', # CRC (MSS)
        'CD79A', 'CD79B', 'NCF1', 'TCL1A', # DLBCL
        'PLTP', 'ELN', 'C9orf24', 'VWA3A', # Ependymoma
        'PTPRZ1', 'MT3', 'PMP2', 'OLIG1', # GBM
        'ALB', 'RBP4', 'APOA2', 'APOC3', # HCC
        'KRT14', 'KRT6B', 'CLCA2', 'LGALS7', # HNSCC (HPV-neg.)
        'S100A9', 'CSTA', 'CALML5', 'FDCSP', # HNSCC (HPV-pos.)
        'NAPSA', 'SFTPB', 'SCGB3A2', 'SCGB1A1', # Lung Adeno.
        'AKR1C2', 'ADH7', 'PLCXD2', 'CHST7', # Lung Squamous
        'SDCBP', 'MLANA', 'PMEL', 'TYR', # Melanoma
        'JCHAIN', 'TENT5C', 'SLAMF7', 'BMP6', # MM
        'STMN2', 'NPY', 'HAND2', 'PHOX2A', # Neuroblastoma
        'RHEX', 'CMTM7', 'CLDN6', 'SOX17', # Ovarian
        'MMP7', 'GCNT3', 'SYT8', 'CLDN18', # PDAC
        'PDLIM5', 'KLK2', 'KLK3', 'KLK4', # Prostate
        'HES6', 'ASCL1', 'INSM1', 'GRP', # SCLC
        'CCN2', 'PODXL', 'TCF21', 'NPHS2' # Wilms Tumor
    ),
    disease = lapply(ss_data_v[, sort(unique(disease))], rep, 4) %>% unlist
)
setkey(ss_data_v, symbol, disease)
ss_data_v[labs_manual_v, lab := symbol]
ss_data_v[, disease := factor(disease, levels = sort(unique(disease)))]

ggplot(ss_data_v) +
    geom_point(aes(x = spec, y = sens), colour = 'grey20', alpha = 0.5) +
    facet_wrap(vars(disease), nrow = 6, ncol = 4) +
    labs(x = 'Specificity', y = 'Sensitivity') +
    geom_text_repel(aes(x = spec, y = sens, label = lab), size = 4.5, nudge_x = 0.05, nudge_y = 0.1) +
    scale_x_continuous(breaks = c(0.2, 0.6, 1), labels = c('0.2', '0.6', '1')) +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
    theme_test() +
    theme(
        axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, margin = margin(t = 10)),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = 'aliceblue')
    )

# Dot plot for cancer-type-specific genes:

genes_dot_v <- data.table(
    disease = sort(unique(ss_data_v$disease)),
    symbol = c('CLEC12A', 'ANKRD30A', 'ANGPTL4', 'CXCL6', 'MPL', 'CDX1', 'TCL1A', 'C9orf24', 'PMP2', 'APOA2', 'KRT14', 'FDCSP', 'SFTPB', 'AKR1C2',
        'PMEL', 'TENT5C', 'HAND2', 'SOX17', 'CLDN18', 'KLK3', 'GRP', 'TCF21')
)
vdata <- ss_data_v[symbol %in% genes_dot_v$symbol]
vdata[, c('symbol', 'disease') := .(factor(symbol, levels = genes_dot_v$symbol), factor(disease, levels = rev(genes_dot_v$disease)))]

ggplot(vdata, aes(x = symbol, y = disease)) +
    geom_point(aes(fill = m, size = 100*p), shape = 21, colour = 'black') +
    scale_fill_gradientn(colours = brewer.pal(9, 'YlOrRd'), limits = c(0, 8), oob = squish) +
    scale_radius(range = c(1, 7), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
        labels = c('0' = '0', '25' = '25', '50' = '50', '75' = '75', '100' = '100')) +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1, size = 14),
        axis.text.y = element_markdown(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.grid = element_line(linewidth = 0.4, colour = 'grey92'),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = 'black'),
        legend.key = element_rect(fill = NA)
    ) +
    guides(fill = guide_colourbar(frame.colour = 'black', order = 1), size = guide_legend(order = 2)) +
    labs(x = 'Gene', y = 'Cancer type', fill = 'Mean', size = '% expressing\ncells')





# Cancer-type-specific genes for other cell types:

ss_data_vct <- copy(ave_data)
ss_data_vct[, spec_inv := ifelse(m == max(m) & m > 0, m[order(-m)][2], as.numeric(NA)), by = .(cell_type, symbol)]
ss_data_vct[!is.na(spec_inv), spec := 1/(spec_inv + 1)]
ss_data_vct[, sens := ifelse(m == max(m) & m > 0, m, as.numeric(NA)), by = .(cell_type, symbol)]

cell_types_thresh <- lapply(transpose(as.list(unique(paths_table[, .(study, cancer_type)]))), function(r) {
    cat(r, '\n')
    cts <- lapply(paths_table[as.list(r), cells], function(p) {
        cells <- suppressWarnings(fread(p, na.strings = ''))
        if(!('cell_type' %in% names(cells))) return(NULL)
        cells <- cells[!is.na(cell_type) & cell_type != 'Unassigned']
        out <- cells[, .(N = .N/nrow(cells)), by = cell_type][N >= 0.01, cell_type] # Cell types that constitute at least 1% of this dataset
        # Check if 'Epithelial' cells are definitely normal epithelial cells:
        if(
            'Epithelial' %in% out && (
                ('malignant' %in% names(cells) && cells[cell_type == 'Epithelial' & malignant == 'no', .N/nrow(cells)] < 0.01) |
                    !('Malignant' %in% cells$cell_type)
            )
        ) {out <- out[out != 'Epithelial']}
        if(length(out) > 0) return(out)
    })
    cts <- cts[!sapply(cts, is.null)]
    if(length(cts) > 0) return(Reduce(intersect, cts))
})
cell_types_ssn <- unlist(cell_types_thresh[!sapply(cell_types_thresh, is.null)])
cell_types_ssn[cell_types_ssn %in% c('Myeloid', 'Monocyte')] <- 'Macrophage'
cell_types_ssn <- table(cell_types_ssn)
cell_types_ssn <- names(cell_types_ssn)[cell_types_ssn >= 30 & names(cell_types_ssn) != ''] # Cell types that pass 1% threshold in >= 30 datasets

ss_data_vct <- ss_data_vct[cell_type %in% gsub('_', ' ', cell_types_ssn)]

ssn_scheme <- ggplot(ss_data_vct[!is.na(sens) & cell_type == 'Malignant' & disease == 'HCC'], aes(x = spec_inv, y = sens)) +
    geom_point(colour = 'grey50', alpha = 0.5, size = 2) +
    geom_abline(aes(slope = 1, intercept = 0), colour = 'red', linewidth = 1) +
    annotate(geom = 'text', label = 'y = x', x = 2.5, y = 2.5, hjust = -0.1, vjust = 1.1, colour = 'red', size = 5) +
    geom_abline(
        aes(slope = 1, intercept = i, colour = c), data = data.table(i = c(1, 1.5, 2, 3))[, c := factor(as.character(i), levels = rev(i))],
        linewidth = 0.8,
        linetype = 'longdash'
    ) +
    scale_colour_manual(values = viridis::viridis(5)[c(1, 3:5)]) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 10), breaks = c(0, 3, 6, 9)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 12.1), breaks = c(0:3, 6, 9, 12)) +
    theme(
        axis.line.y.left = element_line(linewidth = 0.5),
        axis.line.x.bottom = element_line(linewidth = 0.5),
        panel.background = element_rect(colour = NA, fill = NA),
        axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        plot.title = element_text(size = 16),
        legend.text = element_markdown(size = 14),
        legend.title = element_text(size = 16),
        legend.position = c(0.9, 0.25),
        legend.key = element_rect(fill = NA)
    ) +
    labs(x = 'Highest mean expression in malignant cells\namong non-HCC cancer types', y = 'Mean expression in malignant HCC cells',
        colour = 'Threshold', title = 'Specificity of gene expression to malignant cells in HCC')

ssn <- ss_data_vct[
    !is.na(sens),
    setNames(lapply(seq(1, 3.5, 0.5), function(x) sum(sens > spec_inv + x)), paste0('n', seq(1, 3.5, 0.5))),
    by = .(disease = disease, cell_type = as.character(cell_type))
]
ssn_q <- ssn[, slapply(.SD, function(x) quantile(x, c(0.25, 0.5, 0.75))), by = cell_type, .SDcols = paste0('n', seq(1, 3.5, 0.5))]
ssn_q[, q := rep(c(0.25, 0.5, 0.75), .N/3)]
ssn_q <- melt(ssn_q, id.vars = c('cell_type', 'q'), variable.name = 'thresh', value.name = 'n')[, thresh := gsub('n', '', thresh)]
ssn_q <- dcast(ssn_q, cell_type + thresh ~ q)
ssn_q[, cell_type := factor(cell_type, levels = ssn_q[, .(m = mean(`0.5`)), by = cell_type][order(m), cell_type])]
ssn_q[, cell_type_num := as.numeric(cell_type)]
ssn_plot <- ggplot(ssn_q[thresh %in% c('1', '1.5', '2', '3')]) +
    geom_point(aes(x = cell_type, y = `0.5`, colour = thresh), size = 2) +
    geom_line(aes(x = cell_type, y = `0.5`, colour = thresh, group = thresh), linewidth = 0.8) +
    scale_colour_manual(values = viridis::viridis(5)[c(5:3, 1)]) +
    theme(
        axis.line.y.left = element_line(linewidth = 0.5),
        axis.line.x.bottom = element_line(linewidth = 0.5),
        panel.background = element_rect(colour = NA, fill = NA),
        axis.text.x = element_text(size = 14, angle = 35, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        plot.title = element_text(size = 16),
        legend.text = element_markdown(size = 14),
        legend.title = element_text(size = 16),
        legend.key = element_rect(fill = NA)
    ) +
    labs(x = 'Cell type', y = 'Number of cancer-type-specific genes', colour = 'Threshold')

ssn_b <- melt(ssn, id.vars = c('disease', 'cell_type'), variable.name = 'thresh', value.name = 'n')[, thresh := gsub('n', '', thresh)]
ssn_b[, cell_type := factor(cell_type, levels = levels(ssn_q$cell_type))]
ssn_b[disease == 'HNSCC (HPV-neg.)', disease := 'HNSCC (HPV<sup>\u2013</sup>)']
ssn_b[disease == 'HNSCC (HPV-pos.)', disease := 'HNSCC (HPV<sup>+</sup>)']
ssn_box <- ggplot(ssn_b[thresh %in% c('1', '1.5', '2', '3')], aes(x = cell_type, y = log10(n + 1))) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(fill = disease), shape = 21, position = position_jitter(width = 0.2)) +
    scale_fill_manual(values = distinctColorPalette(length(unique(ssn_b$disease)))) +
    facet_wrap(facets = vars(thresh), nrow = 2, ncol = 2) +
    theme_test() +
    theme(
        axis.text.x = element_text(size = 14, angle = 35, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_markdown(size = 16, margin = margin(r = 10)),
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 16),
        legend.text = element_markdown(size = 14),
        legend.title = element_text(size = 16)
    ) +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    labs(x = 'Cell type', y = 'log<sub>10</sub>(1 + number of cancer-type-specific genes)', fill = 'Cancer type',
        title = 'Number of cancer-type-specific genes above different thresholds')
