# The command line arguments should supply study name and cancer type (in that order).
r = commandArgs(trailingOnly = TRUE)





library(data.table)
library(ggplot2)
library(magrittr)
library(Matrix)
library(stringr)
library(plyr)
library(irlba)
library(uwot)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', encoding = 'UTF-8', key = c('study', 'cancer_type'))





paths <- apply(paths_table[as.list(r), .(cells, genes, expmat)], 1, as.list, simplify = FALSE)

out <- lapply(paths, function(p) {
    
    if(endsWith(p$expmat, 'mtx')) {
        
        cells <- suppressWarnings(fread(p$cells, na.strings = ''))
        
        cells$cell_name <- as.character(cells$cell_name)
        
        if(cells[, .N > length(unique(cell_name))]) cells[, cell_name := paste(cell_name, .I, sep = '_')]
        
        genes <- fread(p$genes, header = FALSE)$V1
        expmat <- readMM(p$expmat)
        rownames(expmat) <- genes
        colnames(expmat) <- cells$cell_name
        
        # Remove low-complexity cells and normalise to log TPM/10:
        cells <- cells[complexity >= 1000]
        expmat <- round(log_transform(1e+05*to_frac(expmat[, cells$cell_name])), 4)
        
        # Filtered gene list, after removing lowly expressed genes:
        if('sample' %in% names(cells)) {
            filtered_genes <- cells[
                sample %in% cells[, .(N = .N), by = sample][N >= 10, sample],
                .(gene = genes[rowMeans(log_transform(expmat[, cell_name], reverse = TRUE)) >= 1e+05*(2^4 - 1)/1e+06]),
                by = sample
            ]$gene %>% unique
        } else {
            filtered_genes <- genes[rowMeans(log_transform(expmat, reverse = TRUE)) >= 1e+05*(2^4 - 1)/1e+06]
        }
        
        ct_umap_data_out <- ct_umap_data(expmat[filtered_genes, ], cells)
        
        return(list(data = ct_umap_data_out, path = p))
        
    } else return(NULL)
    
})

if(!all(sapply(out, function(x) is.null(x$data)))) {
    
    ct <- gsub('/', '-', r[2])
    if(!(ct %in% dir('../data/study_plots'))) {dir.create(paste0('../data/study_plots/', ct))}
    if(!(r[1] %in% dir(paste0('../data/study_plots/', ct)))) {dir.create(paste0('../data/study_plots/', ct, '/', r[1]))}
    
    saveRDS(out, paste0('../data/study_plots/', ct, '/', r[1], '/data_ct_umap.RDS'))
    
}
