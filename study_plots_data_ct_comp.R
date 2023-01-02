# Command line arguments to Rscript should supply study name and cancer type (in that order).
r = commandArgs(trailingOnly = TRUE)

library(data.table)
library(ggplot2)
library(magrittr)
library(Matrix)
library(stringr)
library(plyr)
library(matkot)

source('functions.R')

study_folder_map <- fread('../data/study_folder_map.csv', encoding = 'UTF-8')
canonical_markers <- fread('../data/canonical_markers.csv')
paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'))

paths <- apply(paths_table[as.list(r), .(cells, genes, expmat)], 1, as.list, simplify = FALSE)

out <- lapply(paths, function(p) {
        
    if(endsWith(p$expmat, 'mtx')) {
        
        cells <- suppressWarnings(fread(p$cells, na.strings = ''))
        
        cells$cell_name <- as.character(cells$cell_name) # In case cell names are the same as row numbers
        
        if(cells[, .N > length(unique(cell_name))]) cells[, cell_name := paste(cell_name, .I, sep = '_')]
        
        genes <- fread(p$genes, header = FALSE)$V1
        expmat <- readMM(p$expmat)
        rownames(expmat) <- genes
        colnames(expmat) <- cells$cell_name
        
        # Remove low-complexity cells and normalise to log TPM/10:
        cells <- cells[col_nnz(expmat) >= 1000]
        expmat <- round(log_transform(1e+05*to_frac(expmat[, cells$cell_name])), 4)
        
        ct_comp_data_out <- ct_comp_data(expmat, cells, markers = canonical_markers)
        
        return(list(data = ct_comp_data_out, path = p))
        
    } else return(list(data = NULL, path = p))
    
})

if(!all(sapply(out, function(x) is.null(x$data)))) {
    
    ct <- gsub('/', '-', r[2])
    if(!(ct %in% dir('../data/study_plots'))) {dir.create(paste0('../data/study_plots/', ct))}
    if(!(r[1] %in% dir(paste0('../data/study_plots/', ct)))) {dir.create(paste0('../data/study_plots/', ct, '/', r[1]))}
    
    saveRDS(out, paste0('../data/study_plots/', ct, '/', r[1], '/data_ct_comp.RDS'))
    
}
