# The command line arguments should supply study name and cancer type (in that order).
r = commandArgs(trailingOnly = TRUE)





library(data.table)
library(magrittr)
library(plyr)
library(stringr)
library(Matrix)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'))
paths <- copy(paths_table[as.list(r)]); setkey(paths, group)

hgnc_complete_set <- fread('../data/hgnc_complete_set_2023-04-13.txt', key = 'symbol')
hgnc_complete_set <- hgnc_complete_set[!(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])]
alias_table <- make_alias_table(hgnc_complete_set)

genes_all <- fread('../data/gene_plots_data_all_web.csv', select = 'symbol')[symbol %in% hgnc_complete_set$symbol, unique(symbol)]
cts <- c('Malignant', 'B_cell', 'Endothelial', 'Epithelial', 'Fibroblast', 'Macrophage', 'T_cell')
rabbr <- paste(gsub(' | et al. ', '', r[1]), gsub(' |/', '-', r[2]), sep = '_')





pbulks <- lapply(paths$group, function(g) {
    
    cells <- suppressWarnings(fread(paths[g, cells], na.strings = '', colClasses = c(cell_name = 'character', sample = 'character')))
    if(!all(c('cell_type', 'sample') %in% names(cells)) | !endsWith(paths[g, expmat], 'mtx')) return(NULL)
    if(cells[, .N > length(unique(cell_name))]) cells[, cell_name := paste(cell_name, .I, sep = '_')]
    genes <- fread(paths[g, genes], header = FALSE)$V1
    expmat <- readMM(paths[g, expmat])
    expmat <- expmat[genes != '', ] # In case any genes are empty strings
    genes <- genes[genes != '']
    expmat <- expmat[genes %in% names(table(genes))[table(genes) == 1], ] # Removes repeated gene names
    genes <- genes[genes %in% names(table(genes))[table(genes) == 1]]
    genes <- update_symbols_fast(genes, alias_table) # Update gene symbols
    rownames(expmat) <- genes
    colnames(expmat) <- cells$cell_name
    
    genes <- genes[genes %in% genes_all & genes %in% names(table(genes))[table(genes) == 1]]
    
    cells <- cells[cell_type %in% cts & col_nnz(expmat) >= 1000] # Filter cell types and remove low-complexity cells
    setkey(cells, sample, cell_type)
    cells <- cells[cells[, .N, by = .(sample, cell_type)][N >= 30, -'N']]
    
    if(nrow(cells) < 30) return(NULL)
    expmat <- round(1e+05*to_frac(expmat[genes, cells$cell_name]), 4) # Normalise to TPM/10
    
    Reduce(merge, apply(unique(cells[, .(sample, cell_type)]), 1, function(sct) {
        out <- as.data.table(rowMeans(expmat[, cells[as.list(sct), cell_name]]), keep.rownames = TRUE)
        names(out) <- c('gene', paste(rabbr, g, sct[1], gsub('_', '', sct[2]), sep = '_'))
        return(out)
    }, simplify = FALSE))
    
})

pbulks <- Reduce(function(x, y) merge(x, y, by = 'gene', all = TRUE), pbulks[!sapply(pbulks, is.null)])

if(!is.null(pbulks)) fwrite(pbulks, paste0('../data/pseudobulks/', rabbr, '.csv'))
