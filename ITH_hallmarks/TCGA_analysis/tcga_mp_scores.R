library(data.table)
library(magrittr)
library(stringr)
library(readxl)
library(matkot)

source('functions.R')

cancer_types <- c('ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA_AC', 'ESCA_ESCC', 'GBM_IDH-WT', 'HNSC', 'KICH', 'KIRC', 'KIRP',
    'LAML', 'LGG_astro', 'LGG_IDH-WT', 'LGG_oligo', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM_primary',
    'SKCM_metastatic', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM')

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
    
    expmat <- fread(paste0('~/TCGA_data/', ct, '/Exp_data_TPM.csv'))[, set_rownames(as.matrix(.SD), V1), .SDcols = -'V1']
    meta <- fread(paste0('~/TCGA_data/', ct, '/Cells.csv'))[sample_type != 'normal']
    expmat <- expmat[, meta$sample_id]
    
    mps_ct <- slapply(mps_ct, function(x) x[x %in% rownames(expmat)])
    
    set.seed(140)
    scores <- lapply(
        names(mps_ct),
        function(mp) meta[, .(sample_id = sample_id, meta_program = mp, score = sig_score(expmat, mps_ct[[mp]], nbin = 50, n = 50))]
    ) %>% rbindlist
    
    fwrite(scores, paste0('~/TCGA_data/', ct, '/Scores.csv'))
    
}
