library(data.table)
library(ggplot2)
library(magrittr)
library(Matrix)
library(stringr)
library(plyr)
library(cowplot)
library(RColorBrewer)
library(scales)
library(readxl)
library(ggpubr)
library(ggtext)
library(matkot)

try(library(randomcoloR), silent = TRUE)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'))

hgnc_complete_set <- fread('../data/hgnc_complete_set_2023-04-13.txt', key = 'symbol')
hgnc_complete_set <- hgnc_complete_set[!(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])]
alias_table <- make_alias_table(hgnc_complete_set)

cancer_genes <- as.data.table(read_xlsx('../cancer5000.xlsx', sheet = 3, skip = 1))
cancer_genes <- cancer_genes[combined...32 < 0.05, sort(update_symbols_fast(gene, alias_table))]





# Choose studies to include:
diseases <- lapply(transpose(as.list(unique(paths_table[cancer_type != 'Other/Models', .(study, cancer_type)]))), function(r) {
    cat(r, '\n')
    samples_path <- paste0('~/../shared/pan_cancer_datasets/', paths_table[as.list(r), directory[1]], '/samples.csv')
    if(file.exists(samples_path)) {
        samples <- fread(samples_path)[!is.na(sample)] # Setting na.strings = '' doesn't seem to work here...
        samples[, sample := as.character(sample)]
        samples[cancer_type == '', cancer_type := NA]
        setkey(samples, sample)
    } else {
        warning("Samples file doesn't exist for ", r[1], ", ", r[2])
        return(NULL)
    }
    if('cancer_type' %in% names(samples)) {
        return(samples[, .(study = r[1], cancer_type = r[2], sample = sample, disease = cancer_type)])
    } else {
        warning("Samples file for ", r[1], ", ", r[2], " doesn't contain cancer_type column")
        return(NULL)
    }
}) %>% rbindlist

to_exclude <- list(c('Chen et al. 2020', 'Head and Neck'), c('Sun et al. 2021', 'Liver/Biliary')) %>%
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





# Make table of cell type proportions per cancer type and study, and per sample:
set.seed(1852) # Set seed because we downsample in cases where cells of some type have extremely high counts in certain samples.
cc_prop <- lapply(transpose(as.list(to_include)), function(r) {
    
    cat(r, '\n')
    if(!('data_cc.RDS' %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])))) return(NULL)
    
    plot_data <- readRDS(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/data_cc.RDS'))
    nullcond <- sapply(plot_data, function(x) ifelse(is.null(x), TRUE, all(sapply(x[names(x) != 'path'], is.null))))
    if(all(nullcond)) return(NULL)
    
    paths <- apply(paths_table[as.list(r), .(cells, genes, expmat)], 1, as.list, simplify = FALSE)
    
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
        samples[, cancer_type := NA]
    }
    
    samples <- samples[!is.na(cancer_type)]; if(nrow(samples) == 0) return(NULL)
    
    rout <- lapply(which(!nullcond), function(i) {
        
        dt <- plot_data[[i]]$ccdata
        dt <- dt[sample %in% samples$sample][, c('disease', 'tech') := do.call(`[`, list(samples, sample))[, .(cancer_type, technology)]]
        if(!(nrow(dt) >= 100 & 'cell_type' %in% names(dt))) return(NULL)
        dt <- dt[, .(cell_name, cell_type, sample, disease, tech, phase)]
        if(any(c('Myeloid', 'Monocyte') %in% dt$cell_type)) dt[cell_type %in% c('Myeloid', 'Monocyte'), cell_type := 'Macrophage']
        
        # Add complexity column to dt:
        cells <- suppressWarnings(fread(paths[[i]]$cells, na.strings = '', colClasses = c(cell_name = 'character', sample = 'character')))
        setkey(cells, cell_name)
        if('complexity' %in% names(cells)) dt[, n_gene := do.call(`[`, list(cells, cell_name, 'complexity'))] else dt[, n_gene := as.numeric(NA)]
        return(dt)
        
    }) %>% rbindlist
    
    if(is.null(rout) || nrow(rout) == 0) return(NULL)
    
    # Cell cycle proportions across all cells in a dataset (per disease and tech):
    rout_all <- rout[, # Require at least 100 cells altogether, across samples
        if(.N >= 100) .(sample = 'all', n_sample = length(unique(sample)), n_sample_thresh = .SD[, .N, by = sample][N >= 10, .N], n_cell = .N,
            n_g1s = sum(phase == 'G1/S'), n_g2m = sum(phase == 'G2/M'), n_int = sum(phase == 'Intermediate'), n_gene = mean(n_gene)),
        by = .(cell_type, disease, tech)
    ]
    
    # As above but imposing bounds on the number of cells per sample, so that samples with very many cells don't skew the average:
    rout_all_b <- rout[,
        if(.N >= 100) { # Require at least 100 cells altogether, across samples
            # Here, for a given cell type, we downsample cells in samples with exceptionally many (or few) cells of that type.
            bounds <- .SD[, .(N = .N), by = sample][, quantile(N, c(0.75, 0.25)) + c(1.5, -1.5)*floor(IQR(N))]
            .SD[, if(.N >= bounds[2]) {if(.N <= bounds[1]) .SD else .SD[sample(1:.N, bounds[1])]}, by = sample][,
                if(.N >= 100) .(sample = 'all_b', n_sample = length(unique(sample)), n_sample_thresh = .SD[, .N, by = sample][N >= 10, .N],
                    n_cell = .N, n_g1s = sum(phase == 'G1/S'), n_g2m = sum(phase == 'G2/M'), n_int = sum(phase == 'Intermediate'),
                    n_gene = mean(n_gene))
            ]
        },
        by = .(cell_type, disease, tech)
    ]
    
    # Cell cycle proportions per sample:
    rout_smpl <- rout[, # Require at least 100 cells per sample
        if(.N >= 100) .(n_sample = NA, n_sample_thresh = NA, n_cell = .N, n_g1s = sum(phase == 'G1/S'), n_g2m = sum(phase == 'G2/M'),
            n_int = sum(phase == 'Intermediate'), n_gene = mean(n_gene)),
        by = .(cell_type, disease, tech, sample)
    ]
    
    # Check nrow in the following, because if one table has no rows then it will have only 4 columns, when it should have 10:
    rout <- rbindlist(list(rout_all, rout_all_b, rout_smpl)[c(nrow(rout_all) > 0, nrow(rout_all_b) > 0, nrow(rout_smpl) > 0)])
    
    if(!is.null(rout) && nrow(rout) > 0) {
        rout <- unique(rout)
        rout[, c('cancer_type', 'study') := .(r[2], r[1])]
        setcolorder(rout, c('cancer_type', 'study'))
        return(rout)
    }
    
}) %>% rbindlist

cc_prop <- cc_prop[!(cancer_type == 'Brain' & cell_type == 'Fibroblast')]

cc_prop[, disease := mapvalues(
    gsub(' Cancer', '', disease),
    c('Acute Myeloid Leukemia', 'Chronic Myeloid Leukemia', 'Clear Cell Renal Cell Carcinoma', 'Colorectal', 'Cutaneous Basal Cell Carcinoma',
        'Cutaneous Squamous Cell Carcinoma', 'Diffuse Large B Cell Lymphoma', 'Glioblastoma', 'Hepatocellular Carcinoma', 'Lung Adenocarcinoma',
        'Lung Squamous Cell Carcinoma', 'Multiple Myeloma', 'Neuroendocrine Tumor', 'Pancreatic Ductal Adenocarcinoma', 'Small Cell Lung'),
    c('AML', 'CML', 'ccRCC', 'CRC', 'Skin BCC', 'Skin SCC', 'DLBCL', 'GBM', 'HCC', 'Lung Adeno.', 'Lung Squamous', 'MM', 'NET', 'PDAC', 'SCLC'),
    warn_missing = FALSE
)]





# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Summary figures
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# Subset cell types for which we computed cell cycle proportions in at least 10 studies:
common_cts <- cc_prop[sample == 'all_b', .(n = nrow(unique(.SD))), by = cell_type, .SDcols = c('study', 'disease')][n >= 10, cell_type]





# Barplot of cell cycle proportions per cell type across all studies and cancer types:

# Average cell cycle proportions across datasets with different tech within the same study, weighting by number of samples in each:
pdata1 <- cc_prop[
    sample == 'all_b' & cell_type %in% common_cts,
    .(prop = weighted.mean((n_g1s + n_g2m + n_int)/n_cell, n_sample_thresh)),
    by = .(cell_type, study, disease)
]

pdata1 <- pdata1[, .(prop = mean(prop), dev = sd(prop)/sqrt(.N)), by = cell_type]
pdata1[, cell_type := gsub('_', ' ', cell_type)]
pdata1[, cell_type := factor(cell_type, levels = cell_type[order(-prop)])]
sdata1 <- pdata1[, .(
    x = which(levels(pdata1$cell_type) == cell_type) + c(0, -0.1, -0.1),
    xend = which(levels(pdata1$cell_type) == cell_type) + c(0, 0.1, 0.1),
    y = prop + dev*c(-1, -1, 1),
    yend = prop + dev*c(1, -1, 1)
), by = cell_type]
cc_bar_global <- ggplot() +
    geom_col(aes(x = cell_type, y = 100*prop), width = 0.8, colour = 'black', linewidth = 0.5, fill = 'darkseagreen2', data = pdata1) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend), data = sdata1[, .(x = x, xend = xend, y = 100*y, yend = 100*yend)], linewidth = 0.5) +
    scale_y_continuous(expand = c(0, 0.2)) +
    theme_pubclean() +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1, size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, margin = margin(r = 12)),
        axis.ticks.length.x = unit(0, 'pt'),
        plot.title = element_text(size = 18)
    ) +
    labs(x = NULL, y = '% cycling cells', title = 'Mean percentage of cycling cells per cell type\nacross all studies')





# Barplots of cell cycle proportions per study, grouped by cancer type, including points showing weighted mean across studies within each cancer
# type. We take the weighted mean across tech within each study:

cc_bar_ct <- slapply(common_cts, function(ct) {
    diseases_ct <- cc_prop[
        sample == 'all_b' & cell_type == ct,
        .(n_sample_thresh = sum(n_sample_thresh), n_study = nrow(unique(.SD))),
        by = disease,
        .SDcols = c('study', 'cancer_type')
    ][(n_study >= 2 & n_sample_thresh >= 10) | (n_sample_thresh >= 20), disease]
    pdata_ct <- cc_prop[
        sample == 'all_b' & cell_type == ct & disease %in% diseases_ct,
        .(
            prop = weighted.mean((n_g1s + n_g2m + n_int)/n_cell, n_sample_thresh), # Weighted mean of cell cycle prop across tech within each study
            n_sample_thresh = sum(n_sample_thresh) # Keep number of samples to use for weighted mean across studies, and for bar annotations.
        ),
        by = .(disease, study)
    ]
    pdata_ct <- pdata_ct[, # Remove outliers
        .(study = study, prop = prop, n_sample_thresh = n_sample_thresh, outlier = prop > quantile(prop, 0.75) + 1.5*IQR(prop)),
        by = disease
    ][outlier == FALSE, -'outlier']
    pdata_ct[disease == 'HNSCC (HPV-neg.)', disease := 'HNSCC (HPV<sup>\u2013</sup>)']
    pdata_ct[disease == 'HNSCC (HPV-pos.)', disease := 'HNSCC (HPV<sup>+</sup>)']
    pdata_ct[, disease := factor(disease, levels = .SD[, .(m = weighted.mean(prop, n_sample_thresh)), by = disease][order(m), disease])]
    pdata_ct[, study := paste(study, '-', disease)]
    pdata_ct[, study := factor(study, levels = .SD[order(prop), study])]
    pdata_ct[, n_sample_categ := ifelse(n_sample_thresh <= 5, '\u2264 5', ifelse(n_sample_thresh <= 20, '6 - 20', '> 20'))]
    pdata_ct[, n_sample_categ := factor(n_sample_categ, levels = c('\u2264 5', '6 - 20', '> 20'))]
    title_ct <- mapvalues(
        ct,
        c('B_cell', 'Dendritic', 'Endothelial', 'Epithelial', 'Fibroblast', 'Macrophage', 'Malignant', 'Mast', 'NK_cell', 'Pericyte', 'Plasma',
            'T_cell'),
        c('B cells', 'Dendritic cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Macrophages',
            'Mean percentage of cycling malignant cells per cancer type', 'Mast cells', 'NK cells', 'Pericytes', 'Plasma cells', 'T cells'),
        warn_missing = FALSE
    )
    bar_ct <- ggplot(pdata_ct) +
        geom_col(aes(x = as.numeric(disease), y = 100*prop, group = study, fill = n_sample_categ), colour = 'black', linewidth = 0.3,
            position = position_dodge2(preserve = 'single')) +
        scale_fill_manual(values = c('\u2264 5' = 'white', '6 - 20' = 'grey', '> 20' = 'grey25')) +
        scale_x_continuous(
            limits = pdata_ct[, c(0.4, length(levels(disease)) + 0.6)],
            expand = c(0, 0),
            breaks = pdata_ct[, (1:length(levels(disease)))],
            labels = pdata_ct[, setNames(levels(disease), 1:length(unique(disease)))],
            sec.axis = sec_axis(
                ~.,
                breaks = pdata_ct[, (1:length(levels(disease)))],
                labels = pdata_ct[, .(n_sample_thresh = sum(n_sample_thresh)), keyby = disease][pdata_ct[, levels(disease)], n_sample_thresh],
                name = 'Number of samples per cancer type'
            )
        ) +
        geom_point(
            aes(x = x, y = perc),
            data = pdata_ct[, .(perc = 100*weighted.mean(prop, n_sample_thresh)), by = disease][order(perc)][, x := .I],
            size = 4, stroke = 1.5, shape = 3, colour = 'deeppink'
        ) +
        scale_y_continuous(limits = pdata_ct[, c(0, max(10, 100*max(prop)))], breaks = seq(0, 100, if(ct == 'Malignant') 10 else 5),
            expand = c(0, 1)) +
        theme_bw() +
        theme(
            axis.text.x.bottom = element_markdown(angle = 45, hjust = 1, vjust = 1, size = 16),
            axis.text.x.top = element_text(size = 12),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_text(size = 16, margin = margin(r = 10)),
            axis.title.x.top = element_text(margin = margin(b = 8), size = 14),
            axis.ticks.x = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(margin = margin(b = 20), size = 18),
            legend.spacing.y = unit(7, 'pt'),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14),
            panel.border = element_rect(fill = NA, colour = NA)
        ) +
        guides(fill = guide_legend(byrow = TRUE, keyheight = unit(6, 'pt'), keywidth = unit(8, 'pt'))) +
        labs(x = NULL, y = '% cycling cells', fill = 'Number of\nsamples per\nstudy', title = title_ct)
    return(list(plot = bar_ct, data = pdata_ct))
})





# Example cell cycle scatter plot:
data_cc_example <- readRDS('../data/study_plots/Brain/Venteicher et al. 2017/data_cc.RDS')
data_cc_example <- data_cc_example[[1]]$ccdata[cell_type == 'Malignant']
data_cc_example[, phase := factor(phase, levels = c('Not cycling', 'G1/S', 'Intermediate', 'G2/M'))]
cc_example <- ggplot() +
    annotate(geom = 'text', size = 6, x = 2.3, y = 3.1,
        label = paste(percent(data_cc_example[, nrow(.SD[phase != 'Not cycling'])/.N], accuracy = 0.1), 'cycling cells')) +
    geom_point(aes(x = g1s_score, y = g2m_score, colour = phase, size = phase), alpha = 0.5, data = data_cc_example) +
    scale_colour_manual(values = setNames(c('grey', brewer.pal(3, 'Dark2')), c('Not cycling', 'G1/S', 'Intermediate', 'G2/M'))) +
    scale_size_manual(values = c(`Not cycling` = 1, `G1/S` = 2, Intermediate = 2, `G2/M` = 2)) +
    guides(
        colour = guide_legend(title = 'Phase', label.position = 'left', override.aes = list(size = c(2, 3.5, 3.5, 3.5))),
        size = guide_legend(title = 'Phase')
    ) +
    theme_pubr() +
    theme(
        axis.title.x = element_text(size = 16, margin = margin(t = 10)),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        axis.text = element_text(size = 14),
        legend.position = 'right',
        legend.direction = 'vertical',
        legend.justification = 'right',
        legend.key = element_rect(fill = NA, colour = NA),
        legend.key.size = unit(16, 'pt'),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, hjust = 1),
        plot.title = element_text(size = 18, margin = margin(b = 15))
    ) +
    labs(x = 'G1/S score', y = 'G2/M score', title = 'Malignant cells in Venteicher et al. 2017')





# Heatmap of average cell cycle proportions per cell type and cancer type (weighted mean across tech within study, then across studies):

pdata_htmp <- cc_prop[
    sample == 'all_b' & cell_type %in% common_cts,
    .(prop = weighted.mean((n_g1s + n_g2m + n_int)/n_cell, n_sample_thresh), n_sample_thresh = sum(n_sample_thresh)),
    by = .(cell_type, disease, study)
]

# For each cell type, apply thresholds for number of samples and studies having data for that cell type:
setkey(pdata_htmp, cell_type, disease)
diseases_htmp <- cc_prop[
    sample == 'all_b' & cell_type %in% common_cts,
    .(n_sample_thresh = sum(n_sample_thresh), n_study = nrow(unique(.SD))),
    keyby = .(cell_type, disease),
    .SDcols = c('study', 'cancer_type')
][(n_study >= 2 & n_sample_thresh >= 10) | (n_sample_thresh >= 20), .(cell_type, disease)]
pdata_htmp <- pdata_htmp[diseases_htmp]

pdata_htmp <- pdata_htmp[, .(prop = weighted.mean(prop, n_sample_thresh)), keyby = .(cell_type, disease)]
pdata_htmp <- pdata_htmp[CJ(unique(cell_type), unique(disease))]

# So that no row or column has fewer than 3 non-NA values:
pdata_htmp <- pdata_htmp[cell_type %in% pdata_htmp[, .(N = sum(!is.na(prop))), by = cell_type][N >= 3, cell_type]]

pdata_htmp[, cell_type := gsub('_', ' ', cell_type)]
pdata_htmp[disease == 'HNSCC (HPV-neg.)', disease := 'HNSCC (HPV<sup>\u2013</sup>)']
pdata_htmp[disease == 'HNSCC (HPV-pos.)', disease := 'HNSCC (HPV<sup>+</sup>)']
pdata_htmp <- pdata_htmp[!(disease %in% pdata_htmp[cell_type == 'Malignant'][is.na(prop), disease])] # Remove diseases lacking malignant cell data
pdata_htmp[, cell_type := factor(cell_type, levels = .SD[, .(m = mean(prop, na.rm = TRUE)), by = cell_type][order(m), cell_type])]
pdata_htmp[, disease := factor(disease, levels = .SD[cell_type == 'Malignant', disease[order(prop)]])]
pdata_htmp[, na := ifelse(is.na(prop), 'na', 'non_na')]
cc_prop_htmp <- ggplot(pdata_htmp) +
    geom_tile(aes(x = disease, y = cell_type, fill = 100*prop), colour = 'grey95') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(colours = rev(brewer.pal(11, 'Spectral')), limits = c(0, 40), breaks = seq(0, 40, 10), oob = squish, na.value = 'grey95') +
    theme(
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        legend.title = element_text(size = 16, margin = margin(b = 7)),
        legend.text = element_text(size = 14),
        legend.key.height = unit(30, 'pt'),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(size = 18)
    ) +
    guides(fill = guide_colourbar(frame.colour = 'black', ticks.colour = 'black')) +
    labs(x = 'Cancer type', y = 'Cell type', fill = '% cycling\ncells', title = 'Mean percentage of cycling cells')





# Phase bias at study level:

# Per study:

phase_bar_ct <- slapply(common_cts, function(ct) {
    diseases_ct <- cc_prop[
        sample == 'all_b' & cell_type == ct & n_g1s + n_g2m + n_int >= 50, # Minimum 50 cycling cells per study
        .(n_sample_thresh = sum(n_sample_thresh), n_study = nrow(unique(.SD))),
        by = disease,
        .SDcols = c('study', 'cancer_type')
    ][(n_study >= 2 & n_sample_thresh >= 10) | (n_sample_thresh >= 20), disease]
    pdata_ct <- cc_prop[
        sample == 'all_b' & cell_type == ct & disease %in% diseases_ct & n_g1s + n_g2m + n_int >= 50,
        .(
            prop = weighted.mean((n_g1s - n_g2m)/(n_g1s + n_g2m + n_int), n_sample_thresh), # Weighted mean of across tech within each study
            n_sample_thresh = sum(n_sample_thresh) # Keep number of samples to use for weighted mean across studies, and for bar annotations.
        ),
        by = .(disease, study)
    ]
    pdata_ct <- pdata_ct[, # Remove outliers
        .(study = study, prop = prop, n_sample_thresh = n_sample_thresh,
            outlier = prop > quantile(prop, 0.75) + 1.5*IQR(prop) | prop < quantile(prop, 0.25) - 1.5*IQR(prop)),
        by = disease
    ][outlier == FALSE, -'outlier']
    pdata_ct[disease == 'HNSCC (HPV-neg.)', disease := 'HNSCC (HPV<sup>\u2013</sup>)']
    pdata_ct[disease == 'HNSCC (HPV-pos.)', disease := 'HNSCC (HPV<sup>+</sup>)']
    pdata_ct[, disease := factor(disease, levels = .SD[, .(m = weighted.mean(prop, n_sample_thresh)), by = disease][order(m), disease])]
    pdata_ct[, study := paste(study, '-', disease)]
    pdata_ct[, study := factor(study, levels = .SD[order(prop), study])]
    pdata_ct[, n_sample_categ := ifelse(n_sample_thresh <= 5, '\u2264 5', ifelse(n_sample_thresh <= 20, '6 - 20', '> 20'))]
    pdata_ct[, n_sample_categ := factor(n_sample_categ, levels = c('\u2264 5', '6 - 20', '> 20'))]
    title_ct <- mapvalues(
        ct,
        c('B_cell', 'Dendritic', 'Endothelial', 'Epithelial', 'Fibroblast', 'Macrophage', 'Malignant', 'Mast', 'NK_cell', 'Pericyte', 'Plasma',
            'T_cell'),
        c('B cells', 'Dendritic cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Macrophages',
            'Mean phase bias of malignant cells per cancer type', 'Mast cells', 'NK cells', 'Pericytes', 'Plasma cells', 'T cells'),
        warn_missing = FALSE
    )
    bar_ct <- ggplot(pdata_ct) +
        geom_col(aes(x = as.numeric(disease), y = prop, group = study, fill = n_sample_categ), colour = 'black', linewidth = 0.3,
            position = position_dodge2(preserve = 'single')) +
        scale_fill_manual(values = c('\u2264 5' = 'white', '6 - 20' = 'grey', '> 20' = 'grey25')) +
        scale_x_continuous(
            limits = pdata_ct[, c(0.4, length(levels(disease)) + 0.6)],
            expand = c(0, 0),
            breaks = pdata_ct[, (1:length(levels(disease)))],
            labels = pdata_ct[, setNames(levels(disease), 1:length(unique(disease)))],
            sec.axis = sec_axis(
                ~.,
                breaks = pdata_ct[, (1:length(levels(disease)))],
                labels = pdata_ct[, .(n_sample_thresh = sum(n_sample_thresh)), keyby = disease][pdata_ct[, levels(disease)], n_sample_thresh],
                name = 'Number of samples per cancer type'
            )
        ) +
        geom_point(aes(x = x, y = prop), data = pdata_ct[, .(prop = weighted.mean(prop, n_sample_thresh)), by = disease][order(prop)][, x := .I],
            size = 4, stroke = 1.5, shape = 3, colour = 'deeppink') +
        scale_y_continuous(breaks = seq(-1, 1, if(ct %in% c('Malignant', 'Macrophage', 'T_cell', 'Epithelial')) 0.2 else 0.1)) +
        theme_bw() +
        theme(
            axis.text.x.bottom = element_markdown(angle = 45, hjust = 1, size = 16),
            axis.text.x.top = element_text(size = 12),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_markdown(size = 16, margin = margin(r = 10)),
            axis.title.x.top = element_text(margin = margin(b = 8), size = 14),
            axis.ticks.x = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(margin = margin(b = 20), size = 18),
            legend.spacing.y = unit(7, 'pt'),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14),
            panel.border = element_rect(fill = NA, colour = NA)
        ) +
        guides(fill = guide_legend(byrow = TRUE, keyheight = unit(6, 'pt'), keywidth = unit(8, 'pt'))) +
        labs(x = NULL, y = '(n<sub>G1/S</sub> - n<sub>G2/M</sub>) / n<sub>cycling</sub>', fill = 'Number of\nsamples per\nstudy', title = title_ct)
    return(list(plot = bar_ct, data = pdata_ct))
})





# Per cell type and cancer type:

# Because several cell types have many missing data points, show only Malignant cells, T cells and Macrophages:

pdata_phase_htmp <- cc_prop[
    sample == 'all_b' & cell_type %in% c('Malignant', 'T_cell', 'Macrophage') & n_g1s + n_g2m + n_int >= 50, # Minimum 50 cycling cells per study
    .(prop = weighted.mean((n_g1s - n_g2m)/(n_g1s + n_g2m + n_int), n_sample_thresh), n_sample_thresh = sum(n_sample_thresh)),
    by = .(cell_type, disease, study)
]
setkey(pdata_phase_htmp, cell_type, disease)
diseases_htmp <- cc_prop[
    sample == 'all_b' & cell_type %in% c('Malignant', 'T_cell', 'Macrophage') & n_g1s + n_g2m + n_int >= 50,
    .(n_sample_thresh = sum(n_sample_thresh), n_study = nrow(unique(.SD))),
    keyby = .(cell_type, disease),
    .SDcols = c('study', 'cancer_type')
][(n_study >= 2 & n_sample_thresh >= 10) | (n_sample_thresh >= 20), .(cell_type, disease)]
pdata_phase_htmp <- pdata_phase_htmp[diseases_htmp]
pdata_phase_htmp <- pdata_phase_htmp[, .(prop = weighted.mean(prop, n_sample_thresh)), keyby = .(cell_type, disease)]
pdata_phase_htmp <- pdata_phase_htmp[CJ(unique(cell_type), unique(disease))]
pdata_phase_htmp[, cell_type := gsub('_', ' ', cell_type)]
pdata_phase_htmp[, cell_type := factor(cell_type, levels = .SD[, .(m = mean(prop, na.rm = TRUE)), by = cell_type][order(-m), cell_type])]
pdata_phase_htmp[disease == 'HNSCC (HPV-neg.)', disease := 'HNSCC (HPV<sup>\u2013</sup>)']
pdata_phase_htmp[disease == 'HNSCC (HPV-pos.)', disease := 'HNSCC (HPV<sup>+</sup>)']
pdata_phase_htmp[, disease := factor(disease, levels = .SD[cell_type == 'Malignant', disease[order(prop)]])]
phase_htmp <- ggplot(pdata_phase_htmp) +
    geom_tile(aes(x = disease, y = cell_type, fill = prop), colour = 'grey95') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(colours = rev(brewer.pal(9, 'RdBu')), limits = c(-0.1, 0.5),
        values = rescale(c(seq(-0.1, 0.05, 0.05), seq(0.1, 0.5, 0.1))), breaks = c(-0.1, 0.1, 0.3, 0.5), oob = squish, na.value = 'grey95') +
    theme(
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.direction = 'horizontal',
        legend.key.width = unit(30, 'pt'),
        legend.box.margin = margin(l = 15),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(size = 18)
    ) +
    guides(fill = guide_colourbar(frame.colour = 'black', ticks.colour = 'black', title.position = 'top', title.hjust = 0.5)) +
    labs(x = 'Cancer type', y = 'Cell type', fill = expression(frac(n[G1/S]-n[G2/M], n[cycling])), title = 'Mean phase bias')





# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Follow-ups
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# Correlation of cell cycle between cell types across samples:

cts <- c('Malignant', 'Macrophage', 'T_cell', 'NK_cell', 'B_cell', 'Plasma', 'Dendritic', 'Mast', 'Fibroblast', 'Pericyte', 'Endothelial',
    'Epithelial')
cordata <- cc_prop[
    !grepl('^all', sample) & cell_type %in% cts,
    .(disease = disease, study = study, sample = sample, cell_type = cell_type, prop = (n_g1s + n_g2m + n_int)/n_cell)
]
setkey(cordata, study, sample)
cordata[, prop_cent := prop - mean(prop), by = .(disease, study, cell_type)] # Centre proportions per study and cell type
cordata <- cordata[, {
    smpls <- unique(.SD[, .(study, sample)])
    lapply(1:(length(cts) - 1), function(i) lapply((i + 1):length(cts), function(j) {
        v1 <- .SD[cell_type == cts[i]][smpls, prop_cent]; v2 <- .SD[cell_type == cts[j]][smpls, prop_cent]
        cond <- !is.na(v1) & !is.na(v2)
        if(sum(cond) >= 10) return(
            data.table(ct1 = cts[i], ct2 = cts[j])[,
                c('corr', 'pval') := cor.test(v1[cond], v2[cond], method = 'spearman')[c('estimate', 'p.value')]
            ]
        )
    })) %>% unlist(recursive = FALSE) %>% rbindlist
}, by = disease]
cordata[, pval_adj := p.adjust(pval, method = 'BH')]

# Triangular faceted plot where each panel is a bar plot with one bar per cancer type (and coloured by cancer type):

# Filter cell types to keep those with the highest number of correlation values:
cts_filt <- cordata[, sapply(common_cts, function(ct) sum(ct1 == ct) + sum(ct2 == ct))]
cts_filt <- names(cts_filt)[cts_filt >= 20]
cts_filt <- cts_filt[order(sapply(cts_filt, function(x) which(x == cts)))]

pdata <- cordata[ct1 %in% cts_filt & ct2 %in% cts_filt]

pdata[disease == 'HNSCC (HPV-neg.)', disease := 'HNSCC (HPV^(\u2013))']
pdata[disease == 'HNSCC (HPV-pos.)', disease := 'HNSCC (HPV^(+))']
set.seed(6881)
disease_colours <- unique(pdata[, .(disease)])[, setNames(distinctColorPalette(.N), disease)]
corr_plots <- lapply(1:(length(cts_filt) - 1), function(i) lapply(length(cts_filt):(i + 1), function(j) {
    pdata_ij <- pdata[ct1 == cts_filt[i] & ct2 == cts_filt[j]]
    pdata_ij[, disease := factor(disease, levels = disease[order(corr)])]
    out <- ggplot(pdata_ij) +
        geom_col(aes(x = disease, y = corr, fill = disease)) +
        geom_text(aes(x = disease, y = corr + ifelse(corr > 0, 0.05, -0.1), label = '*'), data = pdata_ij[pval_adj < 0.05], size = 5) +
        scale_fill_manual(values = disease_colours) +
        scale_x_discrete(expand = c(0.05, (pdata[, .N, by = .(ct1, ct2)][, max(N)] - pdata_ij[, .N])/2)) +
        scale_y_continuous(sec.axis = dup_axis(name = 'Correlation'), limits = c(-0.5, 1), breaks = c(-0.4, 0, 0.4, 0.8)) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 16, hjust = 0.5),
            axis.title.y.left = element_text(size = 16, margin = margin(r = 7)),
            axis.title.y.right = element_text(size = 14, margin = margin(l = 15)),
            axis.ticks.y.left = element_blank(),
            axis.ticks.length.y.left = unit(0, 'pt'),
            axis.text.y.left = element_blank(),
            axis.text.y.right = element_text(size = 12),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.length.x = unit(0, 'pt'),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank()
        )
    if(i == 1) out <- out + labs(title = gsub('_', ' ', cts_filt[j]))
    if(j == length(cts_filt)) out <- out + labs(y = gsub('_', ' ', cts_filt[i])) else out <- out + labs(y = NULL)
    if(j != i + 1) out <- out + theme(axis.text.y.right = element_blank(), axis.title.y.right = element_blank(),
        axis.ticks.y.right = element_blank())
    return(out)
}))
corr_plots <- lapply(
    corr_plots,
    function(x) if(length(x) < length(cts_filt) - 1) {
        c(x, replicate(length(cts_filt) - 1 - length(x), ggplot() + theme_void(), simplify = FALSE))
    } else x
)

corr_grobs <- lapply(corr_plots, function(x) lapply(x, function(y) ggplotGrob(y + theme(legend.position = 'none'))))
for(i in 1:length(corr_grobs)) {
    # Heights common to all plots:
    for(j in 1:6) corr_grobs[[i]][[j]]$heights[c(1:2, 4:12)] <- unit(c(0.1, 0, 0, 0, 0, 3.5, 0, 0, 0, 0, 0.1), 'cm')
    # Allowing title space only for plots on the first row:
    if(i == 1) {
        for(j in 1:6) corr_grobs[[i]][[j]]$heights[3] <- unit(0.7, 'cm')
    } else {
        for(j in 1:6) corr_grobs[[i]][[j]]$heights[3] <- unit(0, 'cm')
    }
    for(j in 1:6) {
        # Widths common to all plots (even the blank ones):
        corr_grobs[[i]][[j]]$widths[c(1:2, 4:5, 8:9)] <- unit(c(0.1, 0, 0, 3.5, 0, 0.1), 'cm')
        # Allowing y axis title space only for leftmost column:
        if(j == 1) corr_grobs[[i]][[j]]$widths[3] <- unit(0.7, 'cm') else corr_grobs[[i]][[j]]$widths[3] <- unit(0, 'cm')
        # Allowing right y axis text/title space only for rightmost plot in each row:
        if(j == 7 - i) corr_grobs[[i]][[j]]$widths[6:7] <- unit(c(0.7, 0.8), 'cm') else corr_grobs[[i]][[j]]$widths[6:7] <- unit(c(0, 0), 'cm')
    }
}

hwidths <- sapply(corr_grobs, function(x) sum(as.numeric(x[[1]]$heights)))
rwidths <- lapply(corr_grobs, function(x) sapply(x, function(y) sum(as.numeric(y$widths))))

corr_grobs[[5]][[5]] <- get_legend(
    ggplot(pdata) +
        geom_col(aes(x = disease, y = corr, fill = disease)) +
        scale_fill_manual(values = disease_colours) +
        theme(legend.position = c(0.7, 0.5), legend.title = element_text(size = 16), legend.text = element_markdown(size = 14)) +
        guides(fill = guide_legend(nrow = ceiling(length(disease_colours)/2), byrow = TRUE)) +
        labs(fill = 'Cancer type')
)

plot_grid(
    plotlist = lapply(1:length(corr_grobs), function(i) plot_grid(plotlist = corr_grobs[[i]], nrow = 1, ncol = 6, rel_widths = rwidths[[i]])),
    nrow = 6, ncol = 1, rel_heights = hwidths
)

# Global correlations, across cancer types:

cordata_global <- cc_prop[
    !grepl('^all', sample) & cell_type %in% cts,
    .(disease = disease, study = study, sample = sample, cell_type = cell_type, prop = (n_g1s + n_g2m + n_int)/n_cell)
]
setkey(cordata_global, disease, study, sample)
cordata_global[, prop_cent := prop - mean(prop), by = .(disease, study, cell_type)] # Centre proportions per study and cell type
cordata_global <- cordata_global[, {
    smpls <- unique(.SD[, .(disease, study, sample)])
    lapply(1:(length(cts) - 1), function(i) lapply((i + 1):length(cts), function(j) {
        v1 <- .SD[cell_type == cts[i]][smpls, prop_cent]; v2 <- .SD[cell_type == cts[j]][smpls, prop_cent]
        cond <- !is.na(v1) & !is.na(v2)
        if(sum(cond) >= 10) return(
            data.table(ct1 = cts[i], ct2 = cts[j])[,
                c('corr', 'pval') := cor.test(v1[cond], v2[cond], method = 'spearman')[c('estimate', 'p.value')]
            ]
        )
    })) %>% unlist(recursive = FALSE) %>% rbindlist
}]
cordata_global[, pval_adj := p.adjust(pval, method = 'BH')]

pdata_global <- copy(cordata_global)
setkey(pdata_global, ct1, ct2)
pdata_global <- pdata_global[rbindlist(lapply(1:(length(cts) - 1), function(i) data.table(ct1 = cts[i], ct2 = cts[(i + 1):length(cts)])))]
pdata_global[, ct1 := factor(gsub('_', ' ', ct1), levels = rev(gsub('_', ' ', cts)))]
pdata_global[, ct2 := factor(gsub('_', ' ', ct2), levels = rev(gsub('_', ' ', cts)))]

corr_plot_global <- ggplot(pdata_global) +
    geom_tile(aes(x = ct2, y = ct1, fill = corr), colour = 'grey90', linewidth = 0.3) +
    geom_text(
        aes(x = ct2, y = ct1, label = lab),
        data = pdata_global[,
            .(ct1 = ct1, ct2 = ct2, lab = ifelse(pval_adj < 0.001, '***', ifelse(pval_adj < 0.01, '**', ifelse(pval_adj < 0.05, '*', ''))))
        ],
        size = 6,
        vjust = 0.75 # Pushes asterisks downwards, towards centre of tile
    ) +
    geom_segment(
        aes(x = x, xend = xend, y = y, yend = yend),
        data = data.table(x = c(0:10, 1:11) + 0.5, xend = c(1:11, 1:11) + 0.5, y = c(0:10, 0:10) + 0.5, yend = c(0:10, 1:11) + 0.5)
    ) +
    scale_x_discrete(expand = c(0, 0), position = 'top') +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(
        colours = rev(brewer.pal(9, 'RdBu'))[2:9],
        values = rescale(c(seq(-0.3, 0, 0.1), seq(0.15, 0.6, 0.15)), c(0, 1)),
        limits = c(-0.3, 0.6),
        oob = squish,
        breaks = c(-0.3, 0, 0.2, 0.4, 0.6)
    ) +
    guides(fill = guide_colourbar(frame.colour = 'black', ticks.colour = 'black')) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = -30, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = c(0.8, 0.3),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.key.height = unit(30, 'pt'),
        plot.title = element_text(size = 18, margin = margin(b = 15)),
        plot.title.position = 'plot'
    ) +
    labs(x = NULL, y = NULL, fill = 'Correlation', title = 'Correlation of proliferation rates between cell types')





# Association between phase bias and mutations:

cancer_types <- c('ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA_AC', 'ESCA_ESCC', 'GBM_IDH-WT', 'HNSC', 'KICH', 'KIRC', 'KIRP',
    'LAML', 'LGG_astro', 'LGG_IDH-WT', 'LGG_oligo', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM_primary',
    'SKCM_metastatic', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM')

cc_sig <- readRDS('../data/cc_sigs_consensus.rds')

set.seed(1284)
cc_scores <- lapply(cancer_types, function(ct) {
    cat(ct, '\n')
    expmat <- fread(paste0('~/TCGA_data/', ct, '/Exp_data_TPM.csv'))[, set_rownames(as.matrix(.SD), V1), .SDcols = -'V1']
    meta <- fread(paste0('~/TCGA_data/', ct, '/Cells.csv'))[sample_type != 'normal']
    expmat <- expmat[, meta$sample_id]
    cc_sig_ct <- slapply(cc_sig, function(x) x[x %in% rownames(expmat)])
    set.seed(140)
    scores <- lapply(
        names(cc_sig_ct),
        function(p) meta[, .(cancer_type = ct, sample_id = sample_id, phase = p, score = sig_score(expmat, cc_sig_ct[[p]], nbin = 50, n = 50))]
    ) %>% rbindlist
}) %>% rbindlist

cc_diff_mut <- lapply(cancer_types, function(ct) {
    cat(ct, '\n')
    if(!('Mutation_data.csv' %in% dir(paste0('~/TCGA_data/', ct)))) return(NULL)
    
    meta <- fread(paste0('~/TCGA_data/', ct, '/Cells.csv'), key = 'sample_id')
    scores <- cc_scores[cancer_type == ct]
    mut_data <- fread(paste0('~/TCGA_data/', ct, '/Mutation_data.csv'))#, key = 'gene')
    
    # Define MSI and MSS subtypes of COAD, READ, STAD and UCEC by counting mutations:
    if(ct %in% c('COAD', 'READ', 'STAD', 'UCEC')) {
        ct_alt <- mut_data[, .(N = .N), keyby = patient_id]
        ct_alt[, ct_alt := ifelse(N > ifelse(ct == 'UCEC', 400, 500), paste0(ct, '_MSI'), paste0(ct, '_MSS'))]
    }
    
    # Separate HPV-pos and -neg HNSC tumours:
    if(ct == 'HNSC') {
        subtypes_hnsc <- fread('~/TCGA_data/HNSC/subtypes.csv', key = 'patient_id', na.strings = '')
        ct_alt <- subtypes_hnsc[unique(mut_data$patient_id), .(patient_id, hpv_status)]
        ct_alt[!is.na(hpv_status), ct_alt := paste0('HNSC_HPV-', str_extract(hpv_status, '^...'))]
        setkey(ct_alt, patient_id)
    }
    
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
    setkey(mut_data, gene)
    
    # Make a version of the scores data having only patient IDs, for easier searching:
    pscores <- copy(scores)[, patient_id := do.call(`[`, args = list(meta, sample_id))$patient_id][, -'sample_id']
    setcolorder(pscores, 'patient_id')
    setkey(pscores, patient_id)
    
    # Make patient IDs match between pscores and mut_data, removing any patient IDs which occur more than once in pscores:
    pscores <- pscores[patient_id %in% mut_data$patient_id]
    mut_data <- mut_data[patient_id %in% pscores$patient_id]
    pids_tab <- pscores[phase == phase[1], table(patient_id)]
    pids <- names(pids_tab)[pids_tab == 1]
    pscores <- pscores[patient_id %in% pids]
    mut_data <- mut_data[patient_id %in% pids]
    
    if(ct %in% c('COAD', 'READ', 'STAD', 'UCEC', 'HNSC')) {
        pscores[, cancer_type := do.call(`[`, list(ct_alt, patient_id))$ct_alt]
    } else pscores[, cancer_type := ct]
    
    out <- pscores[, .(diff = score[phase == 'g1s'] - score[phase == 'g2m']), by = .(patient_id, cancer_type)]
    cgs <- cancer_genes[cancer_genes %in% mut_data$gene]
    out[, (cgs) := lapply(cgs, function(x) patient_id %in% mut_data[x, patient_id])]
    
    return(out[!is.na(cancer_type)])
    
}) %>% rbindlist(fill = TRUE)

# Remember "diff" is G1/S score minus G2/M score, so positive means more G1/S and negative means more G2/M.

# p53 and Rb boxplots for paper, for main and supp:

ct_p53_main <- cc_diff_mut[!is.na(TP53) & !is.na(diff), .(sum_mut = sum(TP53), N = .N), by = cancer_type]
ct_p53_main <- ct_p53_main[sum_mut >= 10 & N - sum_mut >= 10 & N >= 100 & sum_mut/N >= 0.05, cancer_type]
pdata_p53 <- cc_diff_mut[cancer_type %in% ct_p53_main, .(id = patient_id, cancer_type = cancer_type, diff = diff, gene = TP53)]
pdata_p53[, gene := ifelse(gene, '*TP53*-mut', '*TP53*-WT')]
pdata_p53[grep('_', cancer_type), cancer_type := gsub('_', ' - ', cancer_type)]
pdata_p53_pval <- pdata_p53[, {tout <- t.test(diff ~ gene); .(mean_diff = diff(tout$estimate), pval = tout$p.value)}, by = cancer_type]
pdata_p53_pval[, pval_adj := p.adjust(pval, 'BH')]
pdata_p53_pval[, sig := ifelse(pval_adj < 0.001, '***', ifelse(pval_adj < 0.01, '**', ifelse(pval_adj < 0.05, '*', 'ns')))]
ct_p53_main <- pdata_p53_pval[sig != 'ns', cancer_type]
pdata_p53_main <- pdata_p53[cancer_type %in% ct_p53_main]; pdata_p53_pval_main <- pdata_p53_pval[cancer_type %in% ct_p53_main]
pdata_p53_main[, cancer_type := factor(
    cancer_type,
    levels = pdata_p53_pval_main[order(rowMeans(data.table(r1 = order(order(-mean_diff)), r2 = order(order(pval))))), cancer_type]
)]
plot_p53_main <- ggplot(pdata_p53_main) +
    geom_boxplot(aes(x = cancer_type, y = diff, fill = gene), colour = 'grey50', size = 0.6, outlier.shape = NA, width = 0.6, alpha = 0.25,
        position = position_dodge(width = 0.7)) +
    geom_point(aes(x = cancer_type, y = diff, colour = gene), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
        size = 0.7, alpha = 0.65) +
    geom_segment(
        aes(x = x, xend = xend, y = y, yend = yend),
        data = lapply(
            1:nrow(pdata_p53_pval_main),
            function(i) data.table(x = c(-0.175, -0.175, 0.175) + i, xend = c(0.175, -0.175, 0.175) + i, y = rep(1.35, 3), yend = c(1.35, 1.33, 1.33))
        ) %>% rbindlist
    ) +
    geom_text(aes(x = cancer_type, y = 1.4, label = sig), data = pdata_p53_pval_main, size = 6) +
    theme_pubr() +
    theme(
        axis.text.x = element_text(angle = 35, hjust = 1, size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        plot.title = element_markdown(size = 18, margin = margin(b = 15)),
        legend.text = element_markdown(size = 16, margin = margin(l = -30)),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = 'right',
        legend.spacing.x = unit(40, 'pt'),
        legend.key.size = unit(30, 'pt'),
        legend.box.spacing = unit(0, 'pt')
    ) +
    labs(x = NULL, y = '[G1/S score] \u2013 [G2/M score]', colour = NULL, fill = NULL, title = 'Phase bias and *TP53* mutations in TCGA samples')

pdata_p53[, cancer_type := factor(
    cancer_type,
    levels = pdata_p53_pval[order(rowMeans(data.table(r1 = order(order(-mean_diff)), r2 = order(order(pval))))), cancer_type]
)]
pdata_p53_pval[, vj := ifelse(sig == 'ns', 0, 0.4)]
plot_p53_supp <- ggplot(pdata_p53) +
    geom_boxplot(aes(x = cancer_type, y = diff, fill = gene), colour = 'grey50', size = 0.6, outlier.shape = NA, width = 0.6, alpha = 0.25,
        position = position_dodge(width = 0.7)) +
    geom_point(aes(x = cancer_type, y = diff, colour = gene), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
        size = 0.7, alpha = 0.65) +
    geom_segment(
        aes(x = x, xend = xend, y = y, yend = yend),
        data = lapply(
            1:nrow(pdata_p53_pval),
            function(i) data.table(x = c(-0.175, -0.175, 0.175) + i, xend = c(0.175, -0.175, 0.175) + i, y = rep(1.35, 3), yend = c(1.35, 1.33, 1.33))
        ) %>% rbindlist
    ) +
    geom_text(aes(x = cancer_type, y = 1.4, label = sig, vjust = vj), data = pdata_p53_pval, size = 6) +
    theme_pubr() +
    theme(
        axis.text.x = element_text(angle = 35, hjust = 1, size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        plot.title = element_markdown(size = 18, margin = margin(b = 15)),
        legend.text = element_markdown(size = 16, margin = margin(l = -30)),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = 'right',
        legend.spacing.x = unit(40, 'pt'),
        legend.key.size = unit(30, 'pt'),
        legend.box.spacing = unit(0, 'pt')
    ) +
    labs(x = NULL, y = '[G1/S score] \u2013 [G2/M score]', colour = NULL, fill = NULL, title = 'Phase bias and *TP53* mutations in TCGA samples')

ct_rb <- cc_diff_mut[!is.na(RB1) & !is.na(diff), .(sum_mut = sum(RB1), N = .N), by = cancer_type]
ct_rb <- ct_rb[sum_mut >= 10 & N - sum_mut >= 10 & N >= 100 & sum_mut/N >= 0.05, cancer_type]
pdata_rb <- cc_diff_mut[cancer_type %in% ct_rb, .(id = patient_id, cancer_type = cancer_type, diff = diff, gene = RB1)]
pdata_rb[, gene := ifelse(gene, '*RB1*-mut', '*RB1*-WT')]
pdata_rb[grep('_', cancer_type), cancer_type := gsub('_', ' - ', cancer_type)]
pdata_rb_pval <- pdata_rb[, {tout <- t.test(diff ~ gene); .(mean_diff = diff(tout$estimate), pval = tout$p.value)}, by = cancer_type]
pdata_rb_pval[, pval_adj := p.adjust(pval, 'BH')]
pdata_rb_pval[, sig := ifelse(pval_adj < 0.001, '***', ifelse(pval_adj < 0.01, '**', ifelse(pval_adj < 0.05, '*', 'ns')))]
pdata_rb[, cancer_type := factor(
    cancer_type,
    levels = pdata_rb_pval[order(rowMeans(data.table(r1 = order(order(mean_diff)), r2 = order(order(pval))))), cancer_type]
)]
pdata_rb_pval[, vj := ifelse(sig == 'ns', 0, 0.4)]
plot_rb <- ggplot(pdata_rb) +
    geom_boxplot(aes(x = cancer_type, y = diff, fill = gene), colour = 'grey50', size = 0.6, outlier.shape = NA, width = 0.6, alpha = 0.25,
        position = position_dodge(width = 0.7)) +
    geom_point(aes(x = cancer_type, y = diff, colour = gene), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
        size = 0.7, alpha = 0.65) +
    geom_segment(
        aes(x = x, xend = xend, y = y, yend = yend),
        data = lapply(
            1:nrow(pdata_rb_pval),
            function(i) data.table(x = c(-0.175, -0.175, 0.175) + i, xend = c(0.175, -0.175, 0.175) + i, y = rep(1.35, 3), yend = c(1.35, 1.33, 1.33))
        ) %>% rbindlist
    ) +
    geom_text(aes(x = cancer_type, y = 1.4, label = sig, vjust = vj), data = pdata_rb_pval, size = 6) +
    theme_pubr() +
    theme(
        axis.text.x = element_text(angle = 35, hjust = 1, size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        plot.title = element_markdown(size = 18, margin = margin(b = 15)),
        legend.text = element_markdown(size = 16, margin = margin(l = -30)),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.justification = 'right',
        legend.spacing.x = unit(40, 'pt'),
        legend.key.size = unit(30, 'pt'),
        legend.box.spacing = unit(0, 'pt')
    ) +
    labs(x = NULL, y = '[G1/S score] \u2013 [G2/M score]', colour = NULL, fill = NULL, title = 'Phase bias and *RB1* mutations in TCGA samples')





# Dot plot showing consistency of effect on phase bias for many cancer genes:

test_genes <- names(cc_diff_mut[, -c('patient_id', 'cancer_type', 'diff')])
test_data <- cc_diff_mut[,
    lapply(test_genes, function(g) {
        if(.SD[!is.na(get(g)) & !is.na(diff), sum(get(g)) >= 10 & .N - sum(get(g)) >= 10 & .N >= 100 & sum(get(g))/.N >= 0.05]) {
            gt <- t.test(diff ~ get(g))
            .(gene = g, mean_diff = diff(gt$estimate), pval = gt$p.value)
        }
    }) %>% rbindlist,
    by = cancer_type
]
test_data[abs(mean_diff) > 0.1, pval_adj := p.adjust(pval, 'BH')]

dot_data <- test_data[gene %in% test_data[, .(N = .N, n_sig = sum(pval_adj[!is.na(pval_adj)] < 0.05)), by = gene][N > 1 & n_sig > 0, gene]]
dot_data[pval < 1e-05, pval := 1e-05] # Squish high significance values down to limit of -log10(p) = 5
setkey(dot_data, cancer_type, gene)
dot_genes_ord <- dot_data[dot_data[, CJ(unique(cancer_type), unique(gene))]]
dot_genes_ord[is.na(mean_diff), mean_diff := 0]
dot_genes_ord <- dot_genes_ord[, .(m = mean(mean_diff)), by = gene][order(m), gene]
dot_data[, gene := factor(gene, levels = dot_genes_ord)]
dot_data[grep('_', cancer_type), cancer_type := gsub('_', ' - ', cancer_type)]
dot_test_genes <- ggplot(dot_data) +
    geom_point(aes(x = gene, y = cancer_type, fill = mean_diff, size = -log10(pval)), shape = 21) +
    scale_radius(limits = c(0, 5), breaks = c(0, 2.5, 5), range = c(1, 5.5)) +
    scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(9, 'RdBu')), limits = c(-0.3, 0.3), breaks = c(-0.3, 0, 0.3), oob = squish) +
    scale_x_discrete(expand = c(0, 0.7)) +
    scale_y_discrete(expand = c(0, 0.7)) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
    ) +
    guides(fill = guide_colourbar(barwidth = unit(15, 'pt'), frame.colour = 'black', ticks.colour = 'black')) +
    labs(x = 'Gene', y = 'Cancer type', fill = expression(paste(Delta, '(Score)')), size = 'Significance\n(-log10(p))')





# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Effect of sequencing technology: re-make figures using only 10x data
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# Global summary bar plot:
pdata1_10x <- cc_prop[
    tech == '10x' & sample == 'all_b' & cell_type %in% common_cts,
    .(prop = (n_g1s + n_g2m + n_int)/n_cell),
    by = .(cell_type, study, disease)
]
pdata1_10x <- pdata1_10x[, .(prop = mean(prop), dev = sd(prop)/sqrt(.N)), by = cell_type]
pdata1_10x[, cell_type := gsub('_', ' ', cell_type)]
pdata1_10x[, cell_type := factor(cell_type, levels = cell_type[order(-prop)])]
sdata1_10x <- pdata1_10x[, .(
    x = which(levels(pdata1_10x$cell_type) == cell_type) + c(0, -0.1, -0.1),
    xend = which(levels(pdata1_10x$cell_type) == cell_type) + c(0, 0.1, 0.1),
    y = prop + dev*c(-1, -1, 1),
    yend = prop + dev*c(1, -1, 1)
), by = cell_type]
cc_bar_global_10x <- ggplot() +
    geom_col(aes(x = cell_type, y = 100*prop), width = 0.8, colour = 'black', linewidth = 0.5, fill = 'darkseagreen2', data = pdata1_10x) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend), data = sdata1_10x[, .(x = x, xend = xend, y = 100*y, yend = 100*yend)],
        linewidth = 0.5) +
    scale_y_continuous(expand = c(0, 0.2)) +
    theme_pubclean() +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1, size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, margin = margin(r = 12)),
        axis.ticks.length.x = unit(0, 'pt'),
        plot.title = element_text(size = 18)
    ) +
    labs(x = NULL, y = '% cycling cells', title = 'Mean percentage of cycling cells per cell type\nacross all 10x datasets')

cc_bar_ct_10x <- slapply(common_cts, function(ct) {
    diseases_ct <- cc_prop[
        tech == '10x' & sample == 'all_b' & cell_type == ct,
        .(n_sample_thresh = sum(n_sample_thresh), n_study = nrow(unique(.SD))),
        by = disease,
        .SDcols = c('study', 'cancer_type')
    ][(n_study >= 2 & n_sample_thresh >= 10) | (n_sample_thresh >= 20), disease]
    pdata_ct <- cc_prop[
        tech == '10x' & sample == 'all_b' & cell_type == ct & disease %in% diseases_ct,
        .(prop = (n_g1s + n_g2m + n_int)/n_cell, n_sample_thresh = n_sample_thresh),
        by = .(disease, study)
    ]
    pdata_ct <- pdata_ct[, # Remove outliers
        .(study = study, prop = prop, n_sample_thresh = n_sample_thresh, outlier = prop > quantile(prop, 0.75) + 1.5*IQR(prop)),
        by = disease
    ][outlier == FALSE, -'outlier']
    pdata_ct[disease == 'HNSCC (HPV-neg.)', disease := 'HNSCC (HPV<sup>\u2013</sup>)']
    pdata_ct[disease == 'HNSCC (HPV-pos.)', disease := 'HNSCC (HPV<sup>+</sup>)']
    pdata_ct[, disease := factor(disease, levels = .SD[, .(m = weighted.mean(prop, n_sample_thresh)), by = disease][order(m), disease])]
    pdata_ct[, study := paste(study, '-', disease)]
    pdata_ct[, study := factor(study, levels = .SD[order(prop), study])]
    pdata_ct[, n_sample_categ := ifelse(n_sample_thresh <= 5, '\u2264 5', ifelse(n_sample_thresh <= 20, '6 - 20', '> 20'))]
    pdata_ct[, n_sample_categ := factor(n_sample_categ, levels = c('\u2264 5', '6 - 20', '> 20'))]
    title_ct <- mapvalues(
        ct,
        c('B_cell', 'Dendritic', 'Endothelial', 'Epithelial', 'Fibroblast', 'Macrophage', 'Malignant', 'Mast', 'NK_cell', 'Pericyte', 'Plasma',
            'T_cell'),
        c('B cells', 'Dendritic cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Macrophages',
            'Mean percentage of cycling malignant cells per cancer type', 'Mast cells', 'NK cells', 'Pericytes', 'Plasma cells', 'T cells'),
        warn_missing = FALSE
    )
    bar_ct <- ggplot(pdata_ct) +
        geom_col(aes(x = as.numeric(disease), y = 100*prop, group = study, fill = n_sample_categ), colour = 'black', linewidth = 0.3,
            position = position_dodge2(preserve = 'single')) +
        scale_fill_manual(values = c('\u2264 5' = 'white', '6 - 20' = 'grey', '> 20' = 'grey25')) +
        scale_x_continuous(
            limits = pdata_ct[, c(0.4, length(levels(disease)) + 0.6)],
            expand = c(0, 0),
            breaks = pdata_ct[, (1:length(levels(disease)))],
            labels = pdata_ct[, setNames(levels(disease), 1:length(unique(disease)))],
            sec.axis = sec_axis(
                ~.,
                breaks = pdata_ct[, (1:length(levels(disease)))],
                labels = pdata_ct[, .(n_sample_thresh = sum(n_sample_thresh)), keyby = disease][pdata_ct[, levels(disease)], n_sample_thresh],
                name = 'Number of samples per cancer type'
            )
        ) +
        geom_point(
            aes(x = x, y = perc),
            data = pdata_ct[, .(perc = 100*weighted.mean(prop, n_sample_thresh)), by = disease][order(perc)][, x := .I],
            size = 4, stroke = 1.5, shape = 3, colour = 'deeppink'
        ) +
        scale_y_continuous(limits = pdata_ct[, c(0, max(10, 100*max(prop)))], breaks = seq(0, 100, if(ct == 'Malignant') 10 else 5),
            expand = c(0, 1)) +
        theme_bw() +
        theme(
            axis.text.x.bottom = element_markdown(angle = 45, hjust = 1, vjust = 1, size = 16),
            axis.text.x.top = element_text(size = 12),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_text(size = 16, margin = margin(r = 10)),
            axis.title.x.top = element_text(margin = margin(b = 8), size = 14),
            axis.ticks.x = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(margin = margin(b = 20), size = 18),
            legend.spacing.y = unit(7, 'pt'),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14),
            panel.border = element_rect(fill = NA, colour = NA)
        ) +
        guides(fill = guide_legend(byrow = TRUE, keyheight = unit(6, 'pt'), keywidth = unit(8, 'pt'))) +
        labs(x = NULL, y = '% cycling cells', fill = 'Number of\nsamples per\nstudy', title = title_ct)
    return(list(plot = bar_ct, data = pdata_ct))
})

# Phase bias at study level:
phase_bar_ct_10x <- slapply(common_cts, function(ct) {
    diseases_ct <- cc_prop[
        tech == '10x' & sample == 'all_b' & cell_type == ct & n_g1s + n_g2m + n_int >= 50, # Minimum 50 cycling cells per study
        .(n_sample_thresh = sum(n_sample_thresh), n_study = nrow(unique(.SD))),
        by = disease,
        .SDcols = c('study', 'cancer_type')
    ][(n_study >= 2 & n_sample_thresh >= 10) | (n_sample_thresh >= 20), disease]
    pdata_ct <- cc_prop[
        tech == '10x' & sample == 'all_b' & cell_type == ct & disease %in% diseases_ct & n_g1s + n_g2m + n_int >= 50,
        .(prop = weighted.mean((n_g1s - n_g2m)/(n_g1s + n_g2m + n_int), n_sample_thresh), n_sample_thresh = n_sample_thresh),
        by = .(disease, study)
    ]
    pdata_ct <- pdata_ct[, # Remove outliers
        .(study = study, prop = prop, n_sample_thresh = n_sample_thresh,
            outlier = prop > quantile(prop, 0.75) + 1.5*IQR(prop) | prop < quantile(prop, 0.25) - 1.5*IQR(prop)),
        by = disease
    ][outlier == FALSE, -'outlier']
    pdata_ct[disease == 'HNSCC (HPV-neg.)', disease := 'HNSCC (HPV<sup>\u2013</sup>)']
    pdata_ct[disease == 'HNSCC (HPV-pos.)', disease := 'HNSCC (HPV<sup>+</sup>)']
    pdata_ct[, disease := factor(disease, levels = .SD[, .(m = weighted.mean(prop, n_sample_thresh)), by = disease][order(m), disease])]
    pdata_ct[, study := paste(study, '-', disease)]
    pdata_ct[, study := factor(study, levels = .SD[order(prop), study])]
    pdata_ct[, n_sample_categ := ifelse(n_sample_thresh <= 5, '\u2264 5', ifelse(n_sample_thresh <= 20, '6 - 20', '> 20'))]
    pdata_ct[, n_sample_categ := factor(n_sample_categ, levels = c('\u2264 5', '6 - 20', '> 20'))]
    title_ct <- mapvalues(
        ct,
        c('B_cell', 'Dendritic', 'Endothelial', 'Epithelial', 'Fibroblast', 'Macrophage', 'Malignant', 'Mast', 'NK_cell', 'Pericyte', 'Plasma',
            'T_cell'),
        c('B cells', 'Dendritic cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Macrophages',
            'Mean phase bias of malignant cells per cancer type', 'Mast cells', 'NK cells', 'Pericytes', 'Plasma cells', 'T cells'),
        warn_missing = FALSE
    )
    bar_ct <- ggplot(pdata_ct) +
        geom_col(aes(x = as.numeric(disease), y = prop, group = study, fill = n_sample_categ), colour = 'black', linewidth = 0.3,
            position = position_dodge2(preserve = 'single')) +
        scale_fill_manual(values = c('\u2264 5' = 'white', '6 - 20' = 'grey', '> 20' = 'grey25')) +
        scale_x_continuous(
            limits = pdata_ct[, c(0.4, length(levels(disease)) + 0.6)],
            expand = c(0, 0),
            breaks = pdata_ct[, (1:length(levels(disease)))],
            labels = pdata_ct[, setNames(levels(disease), 1:length(unique(disease)))],
            sec.axis = sec_axis(
                ~.,
                breaks = pdata_ct[, (1:length(levels(disease)))],
                labels = pdata_ct[, .(n_sample_thresh = sum(n_sample_thresh)), keyby = disease][pdata_ct[, levels(disease)], n_sample_thresh],
                name = 'Number of samples per cancer type'
            )
        ) +
        geom_point(aes(x = x, y = prop), data = pdata_ct[, .(prop = weighted.mean(prop, n_sample_thresh)), by = disease][order(prop)][, x := .I],
            size = 4, stroke = 1.5, shape = 3, colour = 'deeppink') +
        scale_y_continuous(breaks = seq(-1, 1, if(ct %in% c('Malignant', 'Macrophage', 'T_cell', 'Epithelial')) 0.2 else 0.1)) +
        theme_bw() +
        theme(
            axis.text.x.bottom = element_markdown(angle = 45, hjust = 1, size = 16),
            axis.text.x.top = element_text(size = 12),
            axis.text.y = element_text(size = 14),
            axis.title.y = element_markdown(size = 16, margin = margin(r = 10)),
            axis.title.x.top = element_text(margin = margin(b = 8), size = 14),
            axis.ticks.x = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(margin = margin(b = 20), size = 18),
            legend.spacing.y = unit(7, 'pt'),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14),
            panel.border = element_rect(fill = NA, colour = NA)
        ) +
        guides(fill = guide_legend(byrow = TRUE, keyheight = unit(6, 'pt'), keywidth = unit(8, 'pt'))) +
        labs(x = NULL, y = '(n<sub>G1/S</sub> - n<sub>G2/M</sub>) / n<sub>cycling</sub>', fill = 'Number of\nsamples per\nstudy', title = title_ct)
    return(list(plot = bar_ct, data = pdata_ct))
})

# Global correlation of cell cycle between cell types across samples:
cts <- c('Malignant', 'Macrophage', 'T_cell', 'NK_cell', 'B_cell', 'Plasma', 'Dendritic', 'Mast', 'Fibroblast', 'Pericyte', 'Endothelial',
    'Epithelial')
cordata_10x <- cc_prop[
    tech == '10x' & !grepl('^all', sample) & cell_type %in% cts,
    .(disease = disease, study = study, sample = sample, cell_type = cell_type, prop = (n_g1s + n_g2m + n_int)/n_cell)
]
setkey(cordata_10x, disease, study, sample)
cordata_10x[, prop_cent := prop - mean(prop), by = .(disease, study, cell_type)] # Centre proportions per study and cell type
cordata_10x <- cordata_10x[, {
    smpls <- unique(.SD[, .(disease, study, sample)])
    lapply(1:(length(cts) - 1), function(i) lapply((i + 1):length(cts), function(j) {
        v1 <- .SD[cell_type == cts[i]][smpls, prop_cent]; v2 <- .SD[cell_type == cts[j]][smpls, prop_cent]
        cond <- !is.na(v1) & !is.na(v2)
        if(sum(cond) >= 10) return(
            data.table(ct1 = cts[i], ct2 = cts[j])[,
                c('corr', 'pval') := cor.test(v1[cond], v2[cond], method = 'spearman')[c('estimate', 'p.value')]
            ]
        )
    })) %>% unlist(recursive = FALSE) %>% rbindlist
}]
cordata_10x[, pval_adj := p.adjust(pval, method = 'BH')]
pdata_10x <- copy(cordata_10x)
setkey(pdata_10x, ct1, ct2)
pdata_10x <- pdata_10x[rbindlist(lapply(1:(length(cts) - 1), function(i) data.table(ct1 = cts[i], ct2 = cts[(i + 1):length(cts)])))]
pdata_10x[, ct1 := factor(gsub('_', ' ', ct1), levels = rev(gsub('_', ' ', cts)))]
pdata_10x[, ct2 := factor(gsub('_', ' ', ct2), levels = rev(gsub('_', ' ', cts)))]
corr_plot_10x <- ggplot(pdata_10x) +
    geom_tile(aes(x = ct2, y = ct1, fill = corr), colour = 'grey90', linewidth = 0.3) +
    geom_text(
        aes(x = ct2, y = ct1, label = lab),
        data = pdata_10x[,
            .(ct1 = ct1, ct2 = ct2, lab = ifelse(pval_adj < 0.001, '***', ifelse(pval_adj < 0.01, '**', ifelse(pval_adj < 0.05, '*', ''))))
        ],
        size = 6,
        vjust = 0.75 # Pushes asterisks downwards, towards centre of tile
    ) +
    geom_segment(
        aes(x = x, xend = xend, y = y, yend = yend),
        data = data.table(x = c(0:10, 1:11) + 0.5, xend = c(1:11, 1:11) + 0.5, y = c(0:10, 0:10) + 0.5, yend = c(0:10, 1:11) + 0.5)
    ) +
    scale_x_discrete(expand = c(0, 0), position = 'top') +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(
        colours = rev(brewer.pal(9, 'RdBu'))[2:9],
        values = rescale(c(seq(-0.3, 0, 0.1), seq(0.15, 0.6, 0.15)), c(0, 1)),
        limits = c(-0.3, 0.6),
        oob = squish,
        breaks = c(-0.3, 0, 0.2, 0.4, 0.6)
    ) +
    guides(fill = guide_colourbar(frame.colour = 'black', ticks.colour = 'black')) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = -30, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = c(0.8, 0.3),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.key.height = unit(30, 'pt'),
        plot.title = element_text(size = 18, margin = margin(b = 15)),
        plot.title.position = 'plot'
    ) +
    labs(x = NULL, y = NULL, fill = 'Correlation', title = 'Correlation of proliferation rates between cell types')





# ----------------------------------------------------------------------------------------------------------------------------------------------------
# Effect of complexity and number of cells
# ----------------------------------------------------------------------------------------------------------------------------------------------------

tdata2 <- cc_prop[sample == 'all_b' & tech %in% c('10x', 'SmartSeq2')]
tdata2_sub <- tdata2[, .(n_ss2 = sum(tech == 'SmartSeq2'), n_10x = sum(tech == '10x')), by = .(cell_type, disease)][n_ss2 > 0 & n_10x > 0]
tdata2_sub <- tdata2_sub[cell_type %in% tdata2_sub[, .(n_ss2 = sum(n_ss2), n_10x = sum(n_10x)), by = cell_type][n_ss2 >= 5 & n_10x >= 5, cell_type]]
setkey(tdata2, cell_type, disease)
tdata2 <- tdata2[tdata2_sub[, .(cell_type, disease)]]
tdata2[, prop := (n_g1s + n_g2m + n_int)/n_cell]
tdata2[n_g1s + n_g2m + n_int >= 50, phase_bias := (n_g1s - n_g2m)/(n_g1s + n_g2m + n_int)]
tdata2[, cell_type := gsub('_', ' ', cell_type)]

# Get candidate outliers by looking at residuals:
tdata2[!is.na(prop) & !is.na(n_gene), o_gene := {
    r1 <- lm(prop ~ n_gene)$residuals; r2 <- lm(n_gene ~ prop)$residuals
    r1 > quantile(r1, 0.75) + 1.5*IQR(r1) | r1 < quantile(r1, 0.25) - 1.5*IQR(r1) |
        r2 > quantile(r2, 0.75) + 1.5*IQR(r2) | r2 < quantile(r2, 0.25) - 1.5*IQR(r2)
}, by = .(tech, cell_type)]
tdata2[!is.na(prop) & !is.na(n_cell), o_cell := {
    r1 <- lm(prop ~ n_cell)$residuals; r2 <- lm(n_cell ~ prop)$residuals
    r1 > quantile(r1, 0.75) + 1.5*IQR(r1) | r1 < quantile(r1, 0.25) - 1.5*IQR(r1) |
        r2 > quantile(r2, 0.75) + 1.5*IQR(r2) | r2 < quantile(r2, 0.25) - 1.5*IQR(r2)
}, by = .(tech, cell_type)]

# Manually choose final list of outliers:
tdata2[, o_gene := FALSE]
tdata2[tech == '10x' & cell_type == 'Macrophage' & prop > 0.14, o_gene := TRUE]
tdata2[tech == '10x' & cell_type == 'T cell' & n_gene > 3000, o_gene := TRUE]
tdata2[tech == 'SmartSeq2' & cell_type == 'T cell' & prop > 0.5, o_gene := TRUE]
tdata2[, o_cell := FALSE]
tdata2[tech == '10x' & cell_type == 'T cell' & prop > 0.24, o_cell := TRUE]
tdata2[tech == 'SmartSeq2' & cell_type == 'Macrophage' & prop > 0.13, o_cell := TRUE]
tdata2[tech == 'SmartSeq2' & cell_type == 'T cell' & prop > 0.5, o_cell := TRUE]

rpos <- unique(tdata2[, .(tech, cell_type)])[, .(v = c('n_gene', 'n_cell')), keyby = .(tech, cell_type)]
rpos[v == 'n_gene', c('x', 'y') := .(c(1, 1, 1, 2, 2, 2, 1, 1), rep(2, 8))]
rpos[v == 'n_cell', c('x', 'y') := .(c(1, 1, 1, 2, 1, 2, 1, 2), rep(2, 8))]
setkey(rpos, tech, cell_type, v)

# Proportion of cycling cells against complexity:
scatter_gene <- slapply(c('10x', 'SmartSeq2'), function(tch) {
    slapply(tdata2[, sort(unique(cell_type))], function(ct) {
        dt <- tdata2[tech == tch & cell_type == ct & !is.na(n_gene)][, c('n_gene', 'prop') := .(n_gene/1000, 100*prop)]
        xlims <- dt[, {
            l1 <- if(ceiling(min(n_gene)) >= floor(max(n_gene))) floor(min(n_gene)) else min(n_gene)
            c(l1, max(max(n_gene), l1 + 1))
        }]
        ylims <- c(0, dt[, max(10, max(prop))])
        out <- ggplot(dt, aes(x = n_gene, y = prop)) +
            geom_point() +
            scale_x_continuous(breaks = seq(1, tdata2[!is.na(n_gene), floor(max(n_gene))], 1), limits = xlims) +
            scale_y_continuous(breaks = 100*seq(0, tdata2[, floor(max(10*prop))/10], 0.1), limits = ylims) +
            theme(
                panel.grid.minor = element_blank(),
                panel.grid.major = element_line(linewidth = 0.3, colour = 'grey85'),
                panel.border = element_rect(fill = NA),
                panel.background = element_rect(fill = 'white'),
                strip.background = element_rect(fill = NA),
                strip.text.x = element_text(size = 16, margin = margin(l = 0, r = 0, t = 2, b = 8)),
                strip.text.y = element_text(size = 16, margin = margin(l = 8, r = 2, t = 0, b = 0)),
                axis.text = element_text(size = 14),
                axis.title = element_blank()
            )
        rp <- rpos[.(tch, ct, 'n_gene')]
        if(any(dt$o_gene)) {
            out <- out + geom_smooth(data = dt[o_gene == FALSE], method = 'lm', colour = 'blue', se = FALSE)
            out <- out + geom_smooth(data = dt, method = 'lm', colour = 'red', se = FALSE)
            out <- out + annotate('text', label = dt[o_gene == FALSE, paste('r =', signif(cor(prop, n_gene), 2))], size = 5, colour = 'blue',
                x = xlims[rp$x], y = ylims[rp$y], hjust = rp$x - 1, vjust = rp$y - 1)
            out <- out + annotate('text', label = dt[, paste('r =', signif(cor(prop, n_gene), 2))], size = 5, colour = 'red',
                x = xlims[rp$x], y = ylims[rp$y] - (rp$y - 1.5)*diff(ylims/5), hjust = rp$x - 1, vjust = rp$y - 1)
        } else {
            out <- out + geom_smooth(data = dt, method = 'lm', colour = 'purple', se = FALSE)
            out <- out + annotate('text', label = dt[, paste('r =', signif(cor(prop, n_gene), 2))], size = 5, colour = 'purple',
                x = xlims[rp$x], y = ylims[rp$y], hjust = rp$x - 1, vjust = rp$y - 1)
        }
        if(tch == '10x' & ct != 'T cell') out <- out + facet_grid(cols = vars(cell_type))
        if(tch == '10x' & ct == 'T cell') out <- out + facet_grid(cols = vars(cell_type), rows = vars(tech))
        if(tch == 'SmartSeq2' & ct == 'T cell') out <- out + facet_grid(rows = vars(tech))
        return(out)
    })
})
scatter_gene_grob <- lapply(unlist(scatter_gene, recursive = FALSE), ggplotGrob)
lgrob_gene <- length(scatter_gene_grob)
for(i in 1:lgrob_gene) {
    l <- length(scatter_gene_grob[[i]]$widths)
    if(i %in% c(1, lgrob_gene/2 + 1)) {
        scatter_gene_grob[[i]]$widths[1] <- unit(12, 'mm') # Left margin
    } else {
        scatter_gene_grob[[i]]$widths[1] <- unit(2, 'mm') # Left margin
    }
    scatter_gene_grob[[i]]$widths[l] <- unit(2, 'mm') # Right margin
}
for(i in 1:lgrob_gene) scatter_gene_grob[[i]]$widths[4] <- unit(7, 'mm') # Y axis text space
for(i in 1:lgrob_gene) scatter_gene_grob[[i]]$widths[5] <- unit(60, 'mm') # Plot area
for(i in (1:2)*lgrob_gene/2) scatter_gene_grob[[i]]$widths[6] <- unit(7.5, 'mm') # Facet strip space
for(i in 1:lgrob_gene) {
    l <- length(scatter_gene_grob[[i]]$heights)
    scatter_gene_grob[[i]]$heights[1] <- unit(2, 'mm') # Top margin
    if(i %in% (lgrob_gene/2 + 1):lgrob_gene) {
        scatter_gene_grob[[i]]$heights[l] <- unit(10, 'mm') # Bottom margin
    } else {
        scatter_gene_grob[[i]]$heights[l] <- unit(2, 'mm') # Bottom margin
    }
}
for(i in 1:(lgrob_gene/2)) scatter_gene_grob[[i]]$heights[7] <- unit(7.5, 'mm') # Facet strip space
for(i in 1:(lgrob_gene/2)) scatter_gene_grob[[i]]$heights[8:9] <- unit(c(60, 5), 'mm') # Plot area and X axis text space
for(i in (lgrob_gene/2 + 1):lgrob_gene) scatter_gene_grob[[i]]$heights[7:8] <- unit(c(60, 5), 'mm') # Plot area and X axis text space

# Proportion of cycling cells against number of captured cells:
scatter_cell <- slapply(c('10x', 'SmartSeq2'), function(tch) {
    slapply(tdata2[, sort(unique(cell_type))], function(ct) {
        dt <- tdata2[tech == tch & cell_type == ct][, c('n_cell', 'prop') := .(log10(n_cell), 100*prop)]
        xlims <- c(dt[, if(ceiling(min(n_cell)) == floor(max(n_cell))) floor(min(n_cell)) else min(n_cell)], dt[, max(n_cell)])
        ylims <- c(0, dt[, max(10, max(prop))])
        out <- ggplot(dt, aes(x = n_cell, y = prop)) +
            geom_point() +
            scale_x_continuous(breaks = seq(1, tdata2[, floor(max(n_cell))], 1), limits = xlims) +
            scale_y_continuous(breaks = 100*seq(0, tdata2[, floor(max(10*prop))/10], 0.1), limits = ylims) +
            theme(
                panel.grid.minor = element_blank(),
                panel.grid.major = element_line(linewidth = 0.3, colour = 'grey85'),
                panel.border = element_rect(fill = NA),
                panel.background = element_rect(fill = 'white'),
                strip.background = element_rect(fill = NA),
                strip.text.x = element_text(size = 16, margin = margin(l = 0, r = 0, t = 2, b = 8)),
                strip.text.y = element_text(size = 16, margin = margin(l = 8, r = 2, t = 0, b = 0)),
                axis.text = element_text(size = 14),
                axis.title = element_blank()
            )
        rp <- rpos[.(tch, ct, 'n_cell')]
        if(any(dt$o_cell)) {
            out <- out + geom_smooth(data = dt[o_cell == FALSE], method = 'lm', colour = 'blue', se = FALSE)
            out <- out + geom_smooth(data = dt, method = 'lm', colour = 'red', se = FALSE)
            out <- out + annotate('text', label = dt[o_cell == FALSE, paste('r =', signif(cor(prop, n_cell), 2))], size = 5, colour = 'blue',
                x = xlims[rp$x], y = ylims[rp$y], hjust = rp$x - 1, vjust = rp$y - 1)
            out <- out + annotate('text', label = dt[, paste('r =', signif(cor(prop, n_cell), 2))], size = 5, colour = 'red',
                x = xlims[rp$x], y = ylims[rp$y] - (rp$y - 1.5)*diff(ylims/5), hjust = rp$x - 1, vjust = rp$y - 1)
        } else {
            out <- out + geom_smooth(data = dt, method = 'lm', colour = 'purple', se = FALSE)
            out <- out + annotate('text', label = dt[, paste('r =', signif(cor(prop, n_cell), 2))], size = 5, colour = 'purple',
                x = xlims[rp$x], y = ylims[rp$y], hjust = rp$x - 1, vjust = rp$y - 1)
        }
        if(tch == '10x' & ct != 'T cell') out <- out + facet_grid(cols = vars(cell_type))
        if(tch == '10x' & ct == 'T cell') out <- out + facet_grid(cols = vars(cell_type), rows = vars(tech))
        if(tch == 'SmartSeq2' & ct == 'T cell') out <- out + facet_grid(rows = vars(tech))
        return(out)
    })
})
scatter_cell_grob <- lapply(unlist(scatter_cell, recursive = FALSE), ggplotGrob)
lgrob_cell <- length(scatter_cell_grob)
for(i in 1:lgrob_cell) {
    l <- length(scatter_cell_grob[[i]]$widths)
    if(i %in% c(1, lgrob_cell/2 + 1)) {
        scatter_cell_grob[[i]]$widths[1] <- unit(12, 'mm') # Left margin
    } else {
        scatter_cell_grob[[i]]$widths[1] <- unit(2, 'mm') # Left margin
    }
    scatter_cell_grob[[i]]$widths[l] <- unit(2, 'mm') # Right margin
}
for(i in 1:lgrob_cell) scatter_cell_grob[[i]]$widths[4] <- unit(7, 'mm') # Y axis text space
for(i in 1:lgrob_cell) scatter_cell_grob[[i]]$widths[5] <- unit(60, 'mm') # Plot area
for(i in (1:2)*lgrob_cell/2) scatter_cell_grob[[i]]$widths[6] <- unit(7.5, 'mm') # Facet strip space
for(i in 1:lgrob_cell) {
    l <- length(scatter_cell_grob[[i]]$heights)
    scatter_cell_grob[[i]]$heights[1] <- unit(2, 'mm') # Top margin
    if(i %in% (lgrob_cell/2 + 1):lgrob_cell) {
        scatter_cell_grob[[i]]$heights[l] <- unit(10, 'mm') # Bottom margin
    } else {
        scatter_cell_grob[[i]]$heights[l] <- unit(2, 'mm') # Bottom margin
    }
}
for(i in 1:(lgrob_cell/2)) scatter_cell_grob[[i]]$heights[7] <- unit(7.5, 'mm') # Facet strip space
for(i in 1:(lgrob_cell/2)) scatter_cell_grob[[i]]$heights[8:9] <- unit(c(60, 5), 'mm') # Plot area and X axis text space
for(i in (lgrob_cell/2 + 1):lgrob_cell) scatter_cell_grob[[i]]$heights[7:8] <- unit(c(60, 5), 'mm') # Plot area and X axis text space





# Phase bias:

tcts_phase <- tdata2[!is.na(phase_bias), .N, by = .(cell_type, tech)][N >= 5, .(n_tech = length(unique(tech))), by = cell_type]
tcts_phase <- tcts_phase[n_tech == 2, cell_type]
tdata2_phase <- tdata2[!is.na(phase_bias) & cell_type %in% tcts_phase]

# Get candidate outliers by looking at residuals:
tdata2_phase[!is.na(phase_bias) & !is.na(n_gene), o_gene := {
    r1 <- lm(phase_bias ~ n_gene)$residuals; r2 <- lm(n_gene ~ phase_bias)$residuals
    r1 > quantile(r1, 0.75) + 1.5*IQR(r1) | r1 < quantile(r1, 0.25) - 1.5*IQR(r1) |
        r2 > quantile(r2, 0.75) + 1.5*IQR(r2) | r2 < quantile(r2, 0.25) - 1.5*IQR(r2)
}, by = .(tech, cell_type)]
tdata2_phase[!is.na(phase_bias) & !is.na(n_cell), o_cell := {
    r1 <- lm(phase_bias ~ n_cell)$residuals; r2 <- lm(n_cell ~ phase_bias)$residuals
    r1 > quantile(r1, 0.75) + 1.5*IQR(r1) | r1 < quantile(r1, 0.25) - 1.5*IQR(r1) |
        r2 > quantile(r2, 0.75) + 1.5*IQR(r2) | r2 < quantile(r2, 0.25) - 1.5*IQR(r2)
}, by = .(tech, cell_type)]

# Manually choose final list of outliers:
tdata2_phase[, o_gene := FALSE]
tdata2_phase[tech == '10x' & cell_type == 'T cell' & (phase_bias < -0.16 | n_gene > 3000), o_gene := TRUE]
tdata2_phase[tech == 'SmartSeq2' & cell_type == 'Malignant' & n_gene < 4000, o_gene := TRUE]
tdata2_phase[tech == 'SmartSeq2' & cell_type == 'T cell' & phase_bias < -0.45, o_gene := TRUE]
tdata2_phase[, o_cell := FALSE]
tdata2_phase[tech == '10x' & cell_type == 'T cell' & phase_bias < -0.16, o_cell := TRUE]
tdata2_phase[tech == 'SmartSeq2' & cell_type == 'T cell' & phase_bias < -0.45, o_cell := TRUE]

rpos_phase <- unique(tdata2_phase[, .(tech, cell_type)])[, .(v = c('n_gene', 'n_cell')), keyby = .(tech, cell_type)]
rpos_phase[v == 'n_gene', c('x', 'y') := .(c(2, 2, 2, 1), c(1, 1, 2, 1))]
rpos_phase[v == 'n_cell', c('x', 'y') := .(c(1, 1, 1, 1), c(2, 1, 2, 2))]
setkey(rpos_phase, tech, cell_type, v)

# Phase bias against complexity:
scatter_gene_phase <- slapply(c('10x', 'SmartSeq2'), function(tch) {
    slapply(tdata2_phase[, sort(unique(cell_type))], function(ct) {
        dt <- tdata2_phase[tech == tch & cell_type == ct & !is.na(n_gene)][, n_gene := n_gene/1000]
        xlims <- c(dt[, if(ceiling(min(n_gene)) == floor(max(n_gene))) floor(min(n_gene)) else min(n_gene)], dt[, max(n_gene)])
        ylims <- dt[, c(min(0, min(phase_bias)), max(0.2, max(phase_bias)))]
        out <- ggplot(dt, aes(x = n_gene, y = phase_bias)) +
            geom_point() +
            scale_x_continuous(breaks = seq(1, tdata2_phase[!is.na(n_gene), floor(max(n_gene))], 1), limits = xlims) +
            scale_y_continuous(breaks = tdata2_phase[, seq(floor(min(10*phase_bias)/10), ceiling(max(10*phase_bias)/10), 0.2)], limits = ylims) +
            theme(
                panel.grid.minor = element_blank(),
                panel.grid.major = element_line(linewidth = 0.3, colour = 'grey85'),
                panel.border = element_rect(fill = NA),
                panel.background = element_rect(fill = 'white'),
                strip.background = element_rect(fill = NA),
                strip.text.x = element_text(size = 16, margin = margin(l = 0, r = 0, t = 2, b = 8)),
                strip.text.y = element_text(size = 16, margin = margin(l = 8, r = 2, t = 0, b = 0)),
                axis.text = element_text(size = 14),
                axis.title = element_blank()
            )
        rp <- rpos_phase[.(tch, ct, 'n_gene')]
        if(any(dt$o_gene)) {
            out <- out + geom_smooth(data = dt[o_gene == FALSE], method = 'lm', colour = 'blue', se = FALSE)
            out <- out + geom_smooth(data = dt, method = 'lm', colour = 'red', se = FALSE)
            out <- out + annotate('text', label = dt[o_gene == FALSE, paste('r =', signif(cor(phase_bias, n_gene), 2))], size = 5, colour = 'blue',
                x = xlims[rp$x], y = ylims[rp$y], hjust = rp$x - 1, vjust = rp$y - 1)
            out <- out + annotate('text', label = dt[, paste('r =', signif(cor(phase_bias, n_gene), 2))], size = 5, colour = 'red',
                x = xlims[rp$x], y = ylims[rp$y] - (rp$y - 1.5)*diff(ylims/5), hjust = rp$x - 1, vjust = rp$y - 1)
        } else {
            out <- out + geom_smooth(data = dt, method = 'lm', colour = 'purple', se = FALSE)
            out <- out + annotate('text', label = dt[, paste('r =', signif(cor(phase_bias, n_gene), 2))], size = 5, colour = 'purple',
                x = xlims[rp$x], y = ylims[rp$y], hjust = rp$x - 1, vjust = rp$y - 1)
        }
        if(tch == '10x' & ct != 'T cell') out <- out + facet_grid(cols = vars(cell_type))
        if(tch == '10x' & ct == 'T cell') out <- out + facet_grid(cols = vars(cell_type), rows = vars(tech))
        if(tch == 'SmartSeq2' & ct == 'T cell') out <- out + facet_grid(rows = vars(tech))
        return(out)
    })
})
scatter_gene_phase_grob <- lapply(unlist(scatter_gene_phase, recursive = FALSE), ggplotGrob)
lgrob_gene_phase <- length(scatter_gene_phase_grob)
for(i in 1:lgrob_gene_phase) {
    l <- length(scatter_gene_phase_grob[[i]]$widths)
    if(i %in% c(1, lgrob_gene_phase/2 + 1)) {
        scatter_gene_phase_grob[[i]]$widths[1] <- unit(12, 'mm') # Left margin
    } else {
        scatter_gene_phase_grob[[i]]$widths[1] <- unit(2, 'mm') # Left margin
    }
    scatter_gene_phase_grob[[i]]$widths[l] <- unit(2, 'mm') # Right margin
}
for(i in 1:lgrob_gene_phase) scatter_gene_phase_grob[[i]]$widths[4] <- unit(9, 'mm') # Y axis text space
for(i in 1:lgrob_gene_phase) scatter_gene_phase_grob[[i]]$widths[5] <- unit(60, 'mm') # Plot area
for(i in (1:2)*lgrob_gene_phase/2) scatter_gene_phase_grob[[i]]$widths[6] <- unit(7.5, 'mm') # Facet strip space
for(i in 1:lgrob_gene_phase) {
    l <- length(scatter_gene_phase_grob[[i]]$heights)
    scatter_gene_phase_grob[[i]]$heights[1] <- unit(2, 'mm') # Top margin
    if(i %in% (lgrob_gene_phase/2 + 1):lgrob_gene_phase) {
        scatter_gene_phase_grob[[i]]$heights[l] <- unit(10, 'mm') # Bottom margin
    } else {
        scatter_gene_phase_grob[[i]]$heights[l] <- unit(2, 'mm') # Bottom margin
    }
}
for(i in 1:(lgrob_gene_phase/2)) scatter_gene_phase_grob[[i]]$heights[7] <- unit(7.5, 'mm') # Facet strip space
for(i in 1:(lgrob_gene_phase/2)) scatter_gene_phase_grob[[i]]$heights[8:9] <- unit(c(60, 5), 'mm') # Plot area and X axis text space
for(i in (lgrob_gene_phase/2 + 1):lgrob_gene_phase) scatter_gene_phase_grob[[i]]$heights[7:8] <- unit(c(60, 5), 'mm') # Plot area and X axis text

# Phase bias against number of captured cells:
scatter_cell_phase <- slapply(c('10x', 'SmartSeq2'), function(tch) {
    slapply(tdata2_phase[, sort(unique(cell_type))], function(ct) {
        dt <- tdata2_phase[tech == tch & cell_type == ct][, n_cell := log10(n_cell)]
        xlims <- c(dt[, if(ceiling(min(n_cell)) == floor(max(n_cell))) floor(min(n_cell)) else min(n_cell)], dt[, max(n_cell)])
        ylims <- dt[, c(min(0, min(phase_bias)), max(0.2, max(phase_bias)))]
        out <- ggplot(dt, aes(x = n_cell, y = phase_bias)) +
            geom_point() +
            scale_x_continuous(breaks = seq(1, tdata2_phase[!is.na(n_cell), floor(max(n_cell))], 1), limits = xlims) +
            scale_y_continuous(breaks = tdata2_phase[, seq(floor(min(10*phase_bias)/10), ceiling(max(10*phase_bias)/10), 0.2)], limits = ylims) +
            theme(
                panel.grid.minor = element_blank(),
                panel.grid.major = element_line(linewidth = 0.3, colour = 'grey85'),
                panel.border = element_rect(fill = NA),
                panel.background = element_rect(fill = 'white'),
                strip.background = element_rect(fill = NA),
                strip.text.x = element_text(size = 16, margin = margin(l = 0, r = 0, t = 2, b = 8)),
                strip.text.y = element_text(size = 16, margin = margin(l = 8, r = 2, t = 0, b = 0)),
                axis.text = element_text(size = 14),
                axis.title = element_blank()
            )
        rp <- rpos_phase[.(tch, ct, 'n_cell')]
        if(any(dt$o_cell)) {
            out <- out + geom_smooth(data = dt[o_cell == FALSE], method = 'lm', colour = 'blue', se = FALSE)
            out <- out + geom_smooth(data = dt, method = 'lm', colour = 'red', se = FALSE)
            out <- out + annotate('text', label = dt[o_cell == FALSE, paste('r =', signif(cor(phase_bias, n_cell), 2))], size = 5, colour = 'blue',
                x = xlims[rp$x], y = ylims[rp$y], hjust = rp$x - 1, vjust = rp$y - 1)
            out <- out + annotate('text', label = dt[, paste('r =', signif(cor(phase_bias, n_cell), 2))], size = 5, colour = 'red',
                x = xlims[rp$x], y = ylims[rp$y] - (rp$y - 1.5)*diff(ylims/5), hjust = rp$x - 1, vjust = rp$y - 1)
        } else {
            out <- out + geom_smooth(data = dt, method = 'lm', colour = 'purple', se = FALSE)
            out <- out + annotate('text', label = dt[, paste('r =', signif(cor(phase_bias, n_cell), 2))], size = 5, colour = 'purple',
                x = xlims[rp$x], y = ylims[rp$y], hjust = rp$x - 1, vjust = rp$y - 1)
        }
        if(tch == '10x' & ct != 'T cell') out <- out + facet_grid(cols = vars(cell_type))
        if(tch == '10x' & ct == 'T cell') out <- out + facet_grid(cols = vars(cell_type), rows = vars(tech))
        if(tch == 'SmartSeq2' & ct == 'T cell') out <- out + facet_grid(rows = vars(tech))
        return(out)
    })
})
scatter_cell_phase_grob <- lapply(unlist(scatter_cell_phase, recursive = FALSE), ggplotGrob)
lgrob_cell_phase <- length(scatter_cell_phase_grob)
for(i in 1:lgrob_cell_phase) {
    l <- length(scatter_cell_phase_grob[[i]]$widths)
    if(i %in% c(1, lgrob_cell_phase/2 + 1)) {
        scatter_cell_phase_grob[[i]]$widths[1] <- unit(12, 'mm') # Left margin
    } else {
        scatter_cell_phase_grob[[i]]$widths[1] <- unit(2, 'mm') # Left margin
    }
    scatter_cell_phase_grob[[i]]$widths[l] <- unit(2, 'mm') # Right margin
}
for(i in 1:lgrob_cell_phase) scatter_cell_phase_grob[[i]]$widths[4] <- unit(9, 'mm') # Y axis text space
for(i in 1:lgrob_cell_phase) scatter_cell_phase_grob[[i]]$widths[5] <- unit(60, 'mm') # Plot area
for(i in (1:2)*lgrob_cell_phase/2) scatter_cell_phase_grob[[i]]$widths[6] <- unit(7.5, 'mm') # Facet strip space
for(i in 1:lgrob_cell_phase) {
    l <- length(scatter_cell_phase_grob[[i]]$heights)
    scatter_cell_phase_grob[[i]]$heights[1] <- unit(2, 'mm') # Top margin
    if(i %in% (lgrob_cell_phase/2 + 1):lgrob_cell_phase) {
        scatter_cell_phase_grob[[i]]$heights[l] <- unit(10, 'mm') # Bottom margin
    } else {
        scatter_cell_phase_grob[[i]]$heights[l] <- unit(2, 'mm') # Bottom margin
    }
}
for(i in 1:(lgrob_cell_phase/2)) scatter_cell_phase_grob[[i]]$heights[7] <- unit(7.5, 'mm') # Facet strip space
for(i in 1:(lgrob_cell_phase/2)) scatter_cell_phase_grob[[i]]$heights[8:9] <- unit(c(60, 5), 'mm') # Plot area and X axis text space
for(i in (lgrob_cell_phase/2 + 1):lgrob_cell_phase) scatter_cell_phase_grob[[i]]$heights[7:8] <- unit(c(60, 5), 'mm') # Plot area and X axis text





w_prop <- sapply(scatter_gene_grob[1:(lgrob_gene/2)], function(x) sum(as.numeric(x$widths)))
w_phase <- sapply(scatter_gene_phase_grob[1:(lgrob_gene_phase/2)], function(x) sum(as.numeric(x$widths)))
h_prop <- sapply(scatter_gene_grob[c(1, lgrob_gene/2 + 1)], function(x) sum(as.numeric(x$heights)))
h_phase <- sapply(scatter_gene_phase_grob[c(1, lgrob_gene_phase/2 + 1)], function(x) sum(as.numeric(x$heights)))

corr_prop_gene <- t.test(tdata2[!is.na(n_gene), .(corr = cor(prop, n_gene)), by = .(tech, cell_type)]$corr)
corr_prop_gene_o <- t.test(tdata2[!is.na(n_gene) & o_gene == FALSE, .(corr = cor(prop, n_gene)), by = .(tech, cell_type)]$corr)
corr_prop_cell <- t.test(tdata2[!is.na(n_cell), .(corr = cor(prop, log10(n_cell))), by = .(tech, cell_type)]$corr)
corr_prop_cell_o <- t.test(tdata2[!is.na(n_cell) & o_cell == FALSE, .(corr = cor(prop, log10(n_cell))), by = .(tech, cell_type)]$corr)
corr_phase_gene <- t.test(tdata2_phase[!is.na(n_gene), .(corr = cor(phase_bias, n_gene)), by = .(tech, cell_type)]$corr)
corr_phase_gene_o <- t.test(tdata2_phase[!is.na(n_gene) & o_gene == FALSE, .(corr = cor(phase_bias, n_gene)), by = .(tech, cell_type)]$corr)
corr_phase_cell <- t.test(tdata2_phase[!is.na(n_cell), .(corr = cor(phase_bias, log10(n_cell))), by = .(tech, cell_type)]$corr)
corr_phase_cell_o <- t.test(tdata2_phase[!is.na(n_cell) & o_cell == FALSE, .(corr = cor(phase_bias, log10(n_cell))), by = .(tech, cell_type)]$corr)
title_gene <- ggplot() + theme_void() + theme(plot.title = element_text(size = 16), plot.subtitle = element_markdown(size = 14)) + labs(
    title = 'Proportion of cycling cells vs. number of detected genes',
    subtitle = paste0('<span style="color:blue">Mean correlation: r = ', signif(corr_prop_gene_o$estimate, 2), ', p = ',
        signif(corr_prop_gene_o$p.value, 2), '</span> <span style="color:red">(with outliers: r = ', signif(corr_prop_gene$estimate, 2), ', p = ',
        signif(corr_prop_gene$p.value, 2), ')</span>')
)
title_cell <- ggplot() + theme_void() + theme(plot.title = element_text(size = 16), plot.subtitle = element_markdown(size = 14)) + labs(
    title = 'Proportion of cycling cells vs. number of captured cells',
    subtitle = paste0('<span style="color:blue">Mean correlation: r = ', signif(corr_prop_cell_o$estimate, 2), ', p = ',
        signif(corr_prop_cell_o$p.value, 2), '</span> <span style="color:red">(with outliers: r = ', signif(corr_prop_cell$estimate, 2), ', p = ',
        signif(corr_prop_cell$p.value, 2), ')</span>')
)
title_gene_phase <- ggplot() + theme_void() + theme(plot.title = element_text(size = 16), plot.subtitle = element_markdown(size = 14)) + labs(
    title = 'Phase bias vs. number of detected genes',
    subtitle = paste0('<span style="color:blue">Mean corr.: r = ', signif(corr_phase_gene_o$estimate, 2), ', p = ',
        signif(corr_phase_gene_o$p.value, 2), '</span> <span style="color:red">(+outliers: r = ', signif(corr_phase_gene$estimate, 2), ', p = ',
        signif(corr_phase_gene$p.value, 2), ')</span>')
)
title_cell_phase <- ggplot() + theme_void() + theme(plot.title = element_text(size = 16), plot.subtitle = element_markdown(size = 14)) + labs(
    title = 'Phase bias vs. number of captured cells',
    subtitle = paste0('<span style="color:blue">Mean corr.: r = ', signif(corr_phase_cell_o$estimate, 2), ', p = ',
        signif(corr_phase_cell_o$p.value, 2), '</span> <span style="color:red">(+outliers: r = ', signif(corr_phase_cell$estimate, 2), ', p = ',
        signif(corr_phase_cell$p.value, 2), ')</span>')
)

plot_grid(
    plot_grid(
        plot_grid(
            lab_plot('a'),
            title_gene,
            plot_grid(plotlist = scatter_gene_grob, nrow = 2, ncol = lgrob_gene/2, rel_widths = w_prop, rel_heights = h_prop) +
                draw_label('(Number of detected genes)/1000', x = 0.51, y = 0.03, size = 16) +
                draw_label('% cycling cells', x = 0.01, y = 0.51, size = 16, angle = 90),
            ggplot() + theme_void(),
            lab_plot('b'),
            title_cell,
            plot_grid(plotlist = scatter_cell_grob, nrow = 2, ncol = lgrob_cell/2, rel_widths = w_prop, rel_heights = h_prop) +
                draw_label('log10(Number of cells)', x = 0.51, y = 0.03, size = 16) +
                draw_label('% cycling cells', x = 0.01, y = 0.51, size = 16, angle = 90),
            nrow = 7, ncol = 1, rel_heights = c(20, 20, sum(h_prop), 20, 20, 20, sum(h_prop))
        ),
        ggplot() + theme_void(),
        nrow = 1, ncol = 2, rel_widths = c(sum(w_prop), 2*sum(w_phase) + 20 - sum(w_prop))
    ),
    ggplot() + theme_void(),
    plot_grid(
        lab_plot('c'),
        ggplot() + theme_void(),
        lab_plot('d'),
        title_gene_phase,
        ggplot() + theme_void(),
        title_cell_phase,
        plot_grid(plotlist = scatter_gene_phase_grob, nrow = 2, ncol = lgrob_gene_phase/2, rel_widths = w_phase, rel_heights = h_phase) +
            draw_label('(Number of detected genes)/1000', x = 0.51, y = 0.03, size = 16) +
            draw_label(expression((n['G1/S']~-~'n'['G2/M'])~'/'~n['cycling']), x = 0.02, y = 0.51, size = 16, angle = 90),
        ggplot() + theme_void(),
        plot_grid(plotlist = scatter_cell_phase_grob, nrow = 2, ncol = lgrob_cell_phase/2, rel_widths = w_phase, rel_heights = h_phase) +
            draw_label('log10(Number of cells)', x = 0.51, y = 0.03, size = 16) +
            draw_label(expression((n['G1/S']~-~'n'['G2/M'])~'/'~n['cycling']), x = 0.02, y = 0.51, size = 16, angle = 90),
        nrow = 3, ncol = 3, rel_widths = c(sum(w_phase), 20, sum(w_phase)), rel_heights = c(20, 20, sum(h_phase))
    ),
    nrow = 3, ncol = 1, rel_heights = c(2*sum(h_prop) + 100, 20, 40 + sum(h_phase))
)
