library(data.table)
library(magrittr)
library(Matrix)
library(stringr)
library(plyr)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'), encoding = 'UTF-8')

hgnc_complete_set <- fread('../data/hgnc_complete_set.txt', key = 'symbol')
hgnc_complete_set <- hgnc_complete_set[!(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])]
alias_table <- make_alias_table(hgnc_complete_set)

set.seed(763)
data_genes_all <- lapply(transpose(as.list(unique(paths_table[, .(study, cancer_type)]))), function(r) {
    cat(r[1], '-', r[2], '\n')
    rdir <- paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])
    if('gene_ave.csv' %in% dir(rdir)) gene_ave <- fread(paste0(rdir, '/gene_ave.csv'))[!is.na(ave)] else return(NULL)
    if(nrow(gene_ave) == 0) return(NULL)
    gene_ave <- gene_ave[!is.na(ave), if(unique(.SD[, .(group, sample, n_cell)])[, sum(n_cell) >= 10]) .SD[,
        .(study = r[1], cancer_type = r[2], ave = sum(ave*n_cell)/sum(n_cell)),
        by = symbol
    ], by = .(group = cell_type)]
    setcolorder(gene_ave, c('study', 'cancer_type'))
    setkey(gene_ave, study, cancer_type, group, symbol)
    # Shuffle the genes before ranking so that the zero genes don't simply get ordered alphabetically:
    gene_ave[, ave_rank := .SD[sample(1:.N), order(order(ave))[order(symbol)]]/.N, .SDcols = c('symbol', 'ave'), by = group]
}) %>% rbindlist

# Bins for non-malignant cell types:

# Get the non-malignant cell types occuring in sufficiently many datasets:
cell_types_nm <- data_genes_all[
    !(group %in% c('', 'Unassigned', 'Malignant')) &
        !(study == 'Chen et al. 2020' & cancer_type == 'Head and Neck') &
        !(study == 'Sun et al. 2021' & cancer_type == 'Liver/Biliary'),
    unique(.SD),
    .SDcols = c('study', 'cancer_type', 'group')
][, .N, by = group][N >= 3, group]

genes_nm <- fread('../data/gene_plots_data_all_web.csv', select = 'symbol')$symbol %>% unique

bins_nm <- data_genes_all[
    symbol %in% genes_nm & group %in% cell_types_nm &
        !(study == 'Chen et al. 2020' & cancer_type == 'Head and Neck') &
        !(study == 'Sun et al. 2021' & cancer_type == 'Liver/Biliary'),
    .(ave = mean(ave[!is.na(ave)]), ave_rank = mean(ave_rank[!is.na(ave_rank)])),
    by = .(symbol, group)
][, cbind(.SD[order(ave_rank)], bin = cut(1:.N, 15, labels = FALSE)), by = group]

fwrite(bins_nm, '../data/bins_nm.csv')
