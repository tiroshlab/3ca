# The command line arguments should supply study name and cancer type (in that order).
r = commandArgs(trailingOnly = TRUE)





library(data.table)
library(magrittr)
library(Matrix)
library(stringr)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'))

hgnc_complete_set <- fread('../data/hgnc_complete_set_2023-04-13.txt', key = 'symbol')
hgnc_complete_set <- hgnc_complete_set[!(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])]
alias_table <- make_alias_table(hgnc_complete_set)

paths <- copy(paths_table[as.list(r)]); setkey(paths, group)

rdir <- paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])

if('mp_scores_all.csv' %in% dir(rdir)) {
    scores <- fread(paste0(rdir, '/mp_scores_all.csv'), colClasses = c(cell_name = 'character'))
    if(nrow(scores) == 0 | !('score' %in% names(scores))) scores <- NULL
} else scores <- NULL

if(!is.null(scores)) {
    
    out <- lapply(paths$group, function(grp) {
        
        cells <- suppressWarnings(fread(paths[grp, cells], na.strings = '', colClasses = c(cell_name = 'character', sample = 'character')))
        
        if(
            scores[group == grp, .N == 0 | all(is.na(score))] |
                !all(c('cell_type', 'sample') %in% names(cells)) |
                !endsWith(paths[grp, expmat], 'mtx')
        ) {
            return(data.table(group = numeric(), group_name = character(), cell_type = character(), meta_program = character(),
                gene = character(), corr = numeric(), n_cell = integer(), n_sample = integer(), n_sample_thresh = integer()))
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
        if(nrow(cells) < 10) return(data.table(group = numeric(), group_name = character(), cell_type = character(), meta_program = character(),
            gene = character(), corr = numeric(), n_cell = integer(), n_sample = integer(), n_sample_thresh = integer()))
        expmat <- round(log_transform(1e+05*to_frac(expmat[, cells$cell_name])), 4)
        
        setkey(cells, cell_name)
        scores[group == grp, c('cell_type', 'sample') := do.call(`[`, list(cells, cell_name))[, .(cell_type, sample)]]
        setkey(scores, sample, cell_name)
        # In the following, I'm using a threshold of 10, in a few places, for minimum number of cells.
        grp_out <- scores[group == grp][, # Chaining so we can use the key even if the same cell/sample names are present in multiple groups
            if(length(unique(cell_name)) >= 10) {
                cat(unique(cell_type), '\n')
                ids <- unique(.SD[, .(sample, cell_name)])
                bounds <- ids[, .N, by = sample][, quantile(N, c(0.75, 0.25)) + c(1.5, -1.5)*floor(IQR(N))]
                ids <- ids[, if(.N >= bounds[2]) {if(.N <= bounds[1]) .(cell_name) else .(cell_name = sample(cell_name, bounds[1]))}, keyby = sample]
                ids_n <- ids[, .N, keyby = sample]
                expmat_sub <- expmat[, ids$cell_name]
                rm_all <- rowMeans(expmat_sub)
                # The following binds together matrices across samples, which might lead to memory issues for large datasets. But we need it to
                # compute the correlation across samples.
                expmat_sub <- t(Reduce(cbind, lapply(ids_n$sample, function(smpl) {
                    if(ids_n[smpl, N] < 10) {
                        apply(expmat_sub[, ids[smpl, cell_name], drop = FALSE], 2, function(x) x - rm_all)
                    } else {
                        rm_smpl <- rowMeans(expmat_sub[, ids[smpl, cell_name]])
                        apply(expmat_sub[, ids[smpl, cell_name]], 2, function(x) x - rm_smpl)
                    }
                })))
                setkey(ids, NULL)
                out <- .SD[ids, .(gene = rownames(expmat), corr = cor(score, expmat_sub)[1, ]), by = meta_program]
                out[, c('n_cell', 'n_sample', 'n_sample_thresh') := ids_n[, .(nrow(ids), .N, sum(N >= 10))]]
                out
            },
            by = cell_type
        ]
        
        grp_out[, c('group', 'group_name') := .(grp, paths[grp, group_name])]
        setcolorder(grp_out, c('group', 'group_name'))
        
        return(grp_out)
        
    }) %>% rbindlist(use.names = TRUE)
    
    if(nrow(out) > 0) fwrite(out, paste0(rdir, '/gene_mp_cor.csv'))
    
}
