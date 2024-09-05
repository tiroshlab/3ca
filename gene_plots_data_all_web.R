library(data.table)
library(magrittr)
library(stringr)
library(Matrix)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'), encoding = 'UTF-8')

cell_types <- readRDS('../data/gene_plots_cell_types.rds')

hgnc_complete_set <- fread('../data/hgnc_complete_set.txt', key = 'symbol')
hgnc_complete_set <- hgnc_complete_set[!(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])]

# The following averages over studies for a given cancer type:

to_include <- unique(paths_table[
    cancer_type != 'Other/Models' &
        !grepl('Unpublished', study) &
        !(study == 'Chen et al. 2020' & cancer_type == 'Head and Neck') &
        !(study == 'Sun et al. 2021' & cancer_type == 'Liver/Biliary'),
    .(study, cancer_type)
])

gene_ave <- lapply(transpose(as.list(to_include)), function(r) {
    
    cat(r, '\n')
    if(!('gene_ave.csv' %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])))) return(NULL)
    
    rout <- fread(
        paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/gene_ave.csv'),
        colClasses = c(cell_type = 'character', symbol = 'character'),
        key = c('cell_type', 'symbol')
    )[cell_type %in% cell_types]
    rout[, c('study', 'cancer_type') := as.list(r)]
    
    samples_path <- paste0('/home/labs/tirosh/shared/pan_cancer_datasets/', paths_table[as.list(r), directory[1]], '/samples.csv')
    samples <- fread(samples_path, colClasses = c(sample = 'character'), na.strings = '')
    samples <- samples[!is.na(sample) & !is.na(cancer_type) & !(cancer_type %in% c('Normal', 'Premalignant'))]
    if(nrow(samples) == 0) return(NULL)
    
    rout <- rout[sample %in% samples$sample] # This excludes the sample == 'all' category
    
    if(all(r == c('Jerby-Arnon et al. 2021', 'Sarcoma'))) { # This dataset is unusual because it has the same sample names in 10x and SS2 datasets
        rout[, tech := group_name]
    } else {
        setkey(samples, sample)
        rout[, tech := do.call(`[`, list(samples, sample))$technology]
    }
    
    setcolorder(rout, c('cancer_type', 'study', 'group', 'group_name', 'tech'))
    
    return(rout)
    
}) %>% rbindlist

gene_ave[
    cell_type %in% c('Macrophage', 'Myeloid', 'Monocyte'),
    c('cell_type', 'n_cell', 'ave', 'prop_pos') := .('Macrophage', sum(n_cell), sum(ave*n_cell)/sum(n_cell), sum(prop_pos*n_cell)/sum(n_cell)),
    by = .(cancer_type, study, tech, sample, symbol)
]
gene_ave <- unique(gene_ave)

gene_ave <- gene_ave[!(cancer_type == 'Brain' & cell_type == 'Fibroblast')]

gene_ave <- gene_ave[symbol %in% hgnc_complete_set$symbol]





# Retain genes that have at least one value in all but at most 3 cancer types:
gene_ave_study <- gene_ave[symbol %in% gene_ave[, .(n = length(unique(cancer_type))), by = symbol][n >= max(n) - 3, symbol]]

gene_ave_study <- gene_ave_study[,
    if(sum(n_cell) >= 10) .(ave = sum(ave*n_cell)/sum(n_cell), prop_pos = sum(prop_pos*n_cell)/sum(n_cell), n_cell = sum(n_cell), n_sample = .N,
        n_sample_thresh = sum(n_cell >= 10)),
    by = .(symbol, cell_type, cancer_type, study, tech) # Mean across samples for each study, cancer type and tech (require >=10 cells in each case)
][,
    .(ave = mean(ave), prop_pos = mean(prop_pos), n_cell = sum(n_cell), n_sample = sum(n_sample), n_sample_thresh = sum(n_sample_thresh)),
    by = .(symbol, cell_type, cancer_type, study) # Mean across datasets of the same study and cancer type but different tech
]

gene_ave_all <- gene_ave_study[,
    .(study = 'all', ave = weighted.mean(ave, n_sample_thresh + 1), prop_pos = weighted.mean(prop_pos, n_sample_thresh + 1), n_cell = sum(n_cell),
        n_sample = sum(n_sample), n_sample_thresh = sum(n_sample_thresh)),
    by = .(symbol, cell_type, cancer_type) # Weighted mean across studies of the same disease, weighted by (number of samples with >= 10 cells) + 1
]

gene_plots_data_all_web <- rbind(gene_ave_study, gene_ave_all, use.names = TRUE)

fwrite(gene_plots_data_all_web, '../data/gene_plots_data_all_web.csv')





unique_ct <- gene_plots_data_all_web[study == 'all', setNames(CJ(unique(cell_type), unique(cancer_type)), c('cell_type', 'cancer_type'))]
unique_ct[, study := 'all']
unique_study <- gene_plots_data_all_web[
    study != 'all',
    setNames(CJ(unique(cell_type), unique(paste(study, cancer_type, sep = ' - '))), c('cell_type', 'study'))
][, c('study', 'cancer_type') := as.data.table(str_split_fixed(study, ' - ', 2))]
unique_dt <- rbind(unique_ct, unique_study, use.names = TRUE)

# Save tables for individual genes, to speed up the writing of the Rmd files:
setkey(gene_plots_data_all_web, cell_type, cancer_type, study)
for(g in gene_plots_data_all_web[, sort(unique(symbol))]) {
    cat(g, '\n')
    # No gene names have underscores in them: sum(grepl('_', gene_plots_data_all_web[, sort(unique(symbol))]))
    # Deal with gene names with slashes in them by replacing with underscore:
    fwrite(
        gene_plots_data_all_web[symbol == g, -'symbol'][unique_dt],
        paste0('../data/gene_plots/gene_plots_data_all_web/', gsub('/', '_', g), '.csv')
    )
}
