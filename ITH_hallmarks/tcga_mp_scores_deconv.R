library(data.table)
library(magrittr)
library(stringr)
library(readxl)
library(matkot)

source('functions.R')

cancer_types <- c('BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'ESCA_AC', 'ESCA_ESCC', 'HNSC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'PRAD', 'READ', 'STAD',
    'THCA', 'UCEC')

hgnc_complete_set <- fread('../data/hgnc_complete_set.txt', key = 'ensembl_gene_id')[
    !(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])
]
alias_table <- make_alias_table(hgnc_complete_set)

study_contrib <- readRDS('../data/study_contribution_per_MP.RDS')
mps <- as.list(read_xlsx('../data/mps_malignant.xlsx')) # Also already ordered by MP number
mps <- slapply(mps, update_symbols_fast, alias_table)

# Make sure MP names are consistent:
names(mps) <- gsub('  ', ' ', names(mps))
names(study_contrib) <- gsub('  ', ' ', names(study_contrib))
names(study_contrib)[41] <- 'MP41 Unassigned'

study_tcga_map <- fread('../data/study_tcga_map.csv', na.strings = '')
cancer_types <- cancer_types[cancer_types %in% study_tcga_map$cancer_type]
setkey(study_tcga_map, cancer_type)

for(ct in cancer_types) {
    
    cat(ct, '\n')
    
    mps_ct <- mps[sapply(study_contrib, function(x) any(x %in% study_tcga_map[ct, study]))]
    
    meta <- fread(paste0('~/TCGA_data/', ct, '/Cells.csv'))[sample_type != 'normal']
    
    if(startsWith(ct, 'ESCA')) {
        expmat <- fread('~/pan_cancer/data/TCGA_deconvolved_scaled/ESCA_Epithelial.cells.txt')
    } else expmat <- fread(paste0('~/pan_cancer/data/TCGA_deconvolved_scaled/', ct, '_Epithelial.cells.txt'))
    expmat <- expmat[, set_rownames(as.matrix(.SD), Gene), .SDcols = -'Gene']
    meta <- meta[sample_id %in% colnames(expmat)]
    expmat <- expmat[, meta$sample_id]
    
    mps_ct <- slapply(mps_ct, function(x) x[x %in% rownames(expmat)])
    
    set.seed(140)
    scores <- lapply(
        names(mps_ct),
        function(mp) meta[, .(sample_id = sample_id, meta_program = mp, score = colMeans(expmat[mps_ct[[mp]], ]))]
    ) %>% rbindlist
    
    fwrite(scores, paste0('~/TCGA_data/', ct, '/Scores_deconv.csv'))
    
}
