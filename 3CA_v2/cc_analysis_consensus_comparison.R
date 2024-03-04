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





set.seed(1852)
cc_prop <- lapply(c('data_cc_consensus.rds', 'data_cc.RDS'), function(fn) rbindlist(lapply(transpose(as.list(to_include)), function(r) {
    
    cat(r, '\n')
    
    # Only difference is to change 'data_cc_consensus.rds' to 'data_cc.RDS' in the following two lines:
    if(!(fn %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])))) return(NULL)
    plot_data <- readRDS(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/', fn))
    
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
        samples[, cancer_type := NA] # There's only 1 true gastric cancer sample in this dataset, and it's not annotated with MSI status.
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
        if(.N >= 100) .(sample = 'all', n_sample = length(unique(sample)), n_cell = .N, n_g1s = sum(phase == 'G1/S'),
            n_g2m = sum(phase == 'G2/M'), n_int = sum(phase == 'Intermediate'), n_gene = mean(n_gene)),
        by = .(cell_type, disease, tech)
    ]
    
    # As above but imposing bounds on the number of cells per sample, so that samples with very many cells don't skew the average:
    rout_all_b <- rout[,
        if(.N >= 100) { # Require at least 100 cells altogether, across samples
            bounds <- .SD[, .(N = .N), by = sample][, quantile(N, c(0.75, 0.25)) + c(1.5, -1.5)*floor(IQR(N))]
            .SD[, if(.N >= bounds[2]) {if(.N <= bounds[1]) .SD else .SD[sample(1:.N, bounds[1])]}, by = sample][,
                if(.N >= 100) .(sample = 'all_b', n_sample = length(unique(sample)), n_cell = .N, n_g1s = sum(phase == 'G1/S'),
                    n_g2m = sum(phase == 'G2/M'), n_int = sum(phase == 'Intermediate'), n_gene = mean(n_gene))
            ]
        },
        by = .(cell_type, disease, tech)
    ]
    
    # Cell cycle proportions per sample:
    rout_smpl <- rout[, # Require at least 100 cells per sample
        if(.N >= 100) .(n_sample = NA, n_cell = .N, n_g1s = sum(phase == 'G1/S'), n_g2m = sum(phase == 'G2/M'),
            n_int = sum(phase == 'Intermediate'), n_gene = mean(n_gene)),
        by = .(cell_type, disease, tech, sample)
    ]
    
    rout <- rbindlist(list(rout_all, rout_all_b, rout_smpl)[c(nrow(rout_all) > 0, nrow(rout_all_b) > 0, nrow(rout_smpl) > 0)])
    
    if(!is.null(rout) && nrow(rout) > 0) {
        rout <- unique(rout)
        rout[, c('cancer_type', 'study') := .(r[2], r[1])]
        setcolorder(rout, c('cancer_type', 'study'))
        return(rout)
    }
    
}))) %>% setNames(c('cons', 'bes'))

for(x in c('cons', 'bes')) {
    cc_prop[[x]] <- cc_prop[[x]][!(cancer_type == 'Brain' & cell_type == 'Fibroblast')][, disease := mapvalues(
        gsub(' Cancer', '', disease),
        c('Acute Myeloid Leukemia', 'Chronic Myeloid Leukemia', 'Clear Cell Renal Cell Carcinoma', 'Colorectal', 'Cutaneous Basal Cell Carcinoma',
            'Cutaneous Squamous Cell Carcinoma', 'Diffuse Large B Cell Lymphoma', 'Glioblastoma', 'Hepatocellular Carcinoma', 'Lung Adenocarcinoma',
            'Lung Squamous Cell Carcinoma', 'Multiple Myeloma', 'Neuroendocrine Tumor', 'Pancreatic Ductal Adenocarcinoma', 'Small Cell Lung'),
        c('AML', 'CML', 'ccRCC', 'CRC', 'Skin BCC', 'Skin SCC', 'DLBCL', 'GBM', 'HCC', 'Lung Adeno.', 'Lung Squamous', 'MM', 'NET', 'PDAC', 'SCLC'),
        warn_missing = FALSE
    )]
}

common_cts <- cc_prop$cons[sample == 'all_b', .(n = nrow(unique(.SD))), by = cell_type, .SDcols = c('study', 'disease')][n >= 10, cell_type]





cor_data <- do.call(merge, lapply(names(cc_prop), function(x) rbind(
    cc_prop[[x]][
        sample == 'all_b' & cell_type %in% common_cts,
        setNames(.('Cell cycle proportion', weighted.mean((n_g1s + n_g2m + n_int)/n_cell, n_sample)), c('v', paste0('val_', x))),
        by = .(cell_type, study, disease)
    ],
    cc_prop[[x]][
        sample == 'all_b' & cell_type %in% common_cts & n_g1s + n_g2m + n_int >= 50,
        setNames(.('Phase bias', weighted.mean((n_g1s - n_g2m)/(n_g1s + n_g2m + n_int), n_sample)), c('v', paste0('val_', x))),
        by = .(cell_type, study, disease)
    ]
)))

ggplot(cor_data, aes(x = val_cons, y = val_bes)) +
    geom_point(size = 3, shape = 21, fill = 'lightblue') +
    facet_wrap(vars(v), nrow = 1, ncol = 2, scales = 'free') +
    geom_text(
        aes(x = x, y = y, label = l),
        cor_data[,
            .(x = min(val_cons) + 0.18*diff(range(val_cons)), y = min(val_bes) + 0.9*diff(range(val_bes)),
                l = paste('r =', signif(cor(val_bes, val_cons), 2))),
            by = v
        ],
        size = 8
    ) +
    geom_abline(slope = 1, intercept = 0, colour = 'tomato', linewidth = 1) +
    theme_bw() +
    theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = 'grey95'),
        aspect.ratio = 1
    ) +
    labs(x = 'Consensus signature', y = 'Bespoke signature', title = 'Cell cycle quantification with bespoke and consensus gene signatures')
