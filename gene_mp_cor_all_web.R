library(data.table)
library(magrittr)
library(stringr)
library(Matrix)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'), encoding = 'UTF-8')

hgnc_complete_set <- fread('../data/hgnc_complete_set.txt', key = 'symbol')
hgnc_complete_set <- hgnc_complete_set[!(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])]

genes <- fread('../data/gene_plots_data_all_web.csv', select = 'symbol')[symbol %in% hgnc_complete_set$symbol, unique(symbol)]

to_include <- unique(paths_table[
    cancer_type != 'Other/Models' &
        !grepl('Unpublished', study) &
        !(study == 'Chen et al. 2020' & cancer_type == 'Head and Neck') &
        !(study == 'Sun et al. 2021' & cancer_type == 'Liver/Biliary'),
    .(study, cancer_type)
])

gene_mp_cor <- lapply(transpose(as.list(to_include)), function(r) {
    
    cat(r, '\n')
    if(!('gene_mp_cor.csv' %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])))) return(NULL)
    
    rout <- fread(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/gene_mp_cor.csv'))
    rout <- rout[gene %in% genes]
        
    rout[, c('study', 'cancer_type') := .(r[1], r[2])]
    setcolorder(rout, c('cancer_type', 'study'))
    
}) %>% rbindlist

dtkey <- unique(gene_mp_cor[, .(gene, cell_type, cancer_type, meta_program)])[order(gene, cell_type, cancer_type, meta_program)]
gene_mp_cor_ave <- gene_mp_cor[
    !is.na(corr),
    if(sum(n_cell) >= 50 & sum(n_sample) >= 3) .(mean_corr = weighted.mean(corr, n_sample_thresh)),
    keyby = .(gene, cell_type, cancer_type, meta_program)
]
gene_mp_cor_ave <- gene_mp_cor_ave[dtkey]

fwrite(gene_mp_cor_ave, '../data/gene_mp_cor_all_web.csv')

# Save tables for individual genes, to speed up the writing of the Rmd files:
setkey(gene_mp_cor_ave, gene)
for(g in sort(genes[genes %in% gene_mp_cor_ave$gene])) {
    cat(g, '\n')
    fwrite(gene_mp_cor_ave[g, -'gene'], paste0('../data/gene_plots/gene_mp_cor/', g, '.csv'))
}
