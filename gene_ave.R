r = commandArgs(trailingOnly = TRUE)





library(data.table)
library(magrittr)
library(Matrix)
library(stringr)
library(plyr)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'))

hgnc_complete_set <- fread('../data/hgnc_complete_set_2023-04-13.txt', key = 'symbol')
hgnc_complete_set <- hgnc_complete_set[!(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])]
alias_table <- make_alias_table(hgnc_complete_set)

paths <- copy(paths_table[as.list(r)]); setkey(paths, group)

out <- rbindlist(lapply(paths$group, function(grp) {
    
    cells <- suppressWarnings(fread(paths[grp, cells], na.strings = '', colClasses = c(cell_name = 'character', sample = 'character')))
    
    if(!all(c('cell_type', 'sample') %in% names(cells)) | !endsWith(paths[grp, expmat], 'mtx')) {
        return(data.table(group = grp, group_name = paths[grp, group_name], cell_type = character(), sample = character(), symbol = character(),
            ave = numeric(), prop_pos = numeric(), n_cell = integer()))
    }
    
    genes <- fread(paths[grp, genes], header = FALSE)$V1
    expmat <- readMM(paths[grp, expmat])
    expmat <- expmat[genes %in% names(table(genes))[table(genes) == 1], ]
    genes <- genes[genes %in% names(table(genes))[table(genes) == 1]]
    genes <- update_symbols_fast(genes, alias_table)
    rownames(expmat) <- genes
    colnames(expmat) <- cells$cell_name
    
    # Remove low-complexity cells and normalise to log TPM/10:
    cells <- cells[col_nnz(expmat) >= 1000]
    if(nrow(cells) < 10) return(data.table(group = grp, group_name = paths[grp, group_name], cell_type = character(), sample = character(),
        symbol = character(), ave = numeric(), prop_pos = numeric(), n_cell = integer()))
    expmat <- round(log_transform(1e+05*to_frac(expmat[, cells$cell_name])), 4)
    
    pout <- cells[
        !is.na(cell_type) & !is.na(sample),
        .(symbol = genes, ave = rowMeans(expmat[, cell_name, drop = FALSE]), prop_pos = row_nnz(expmat[, cell_name, drop = FALSE])/.N, n_cell = .N),
        by = .(cell_type, sample)
    ]
    pout[, c('group', 'group_name') := .(grp, paths[grp, group_name])]
    setcolorder(pout, c('group', 'group_name'))
    
    return(pout)
    
}), use.names = TRUE)

if(nrow(out) > 0) {
    ct <- gsub('/', '-', r[2])
    if(!(ct %in% dir('../data/study_plots'))) {dir.create(paste0('../data/study_plots/', ct))}
    if(!(r[1] %in% dir(paste0('../data/study_plots/', ct)))) {dir.create(paste0('../data/study_plots/', ct, '/', r[1]))}
    fwrite(out, paste0('../data/study_plots/', ct, '/', r[1], '/gene_ave.csv'))
}
