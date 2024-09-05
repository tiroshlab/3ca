# The command line arguments should supply study name and cancer type (in that order).
r = commandArgs(trailingOnly = TRUE)





library(data.table)
library(magrittr)
library(readxl)
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

meta_programs <- lapply(
    c('Malignant', 'B cells', 'Endothelial', 'Epithelial', 'Fibroblasts', 'Macrophages', 'CD4 T cells', 'CD8 T cells'),
    function(ct) {
        dt <- as.data.table(read_xlsx('../data/meta_programs_2023-07-13.xlsx', sheet = ct))
        dt[, (names(dt)) := lapply(dt, update_symbols_fast, alias_table)]
        return(dt)
    }
)

# Change some MP names, especially for CD4 and CD8 T cells, so we can merge them into a single T cell category:
names(meta_programs[[1]])[c(12:15, 30, 32, 34, 39)] <- gsub('-', ' ', names(meta_programs[[1]])[c(12:15, 30, 32, 34, 39)])
names(meta_programs[[1]])[c(19, 35, 37, 41)] <- c('EpiSen', 'Hemato-related I', 'Hemato-related II', 'MP41 (Unassigned)')
names(meta_programs[[2]])[3] <- 'Cell Cycle'
names(meta_programs[[3]])[c(1, 6)] <- c('Notch signaling', 'Cell Cycle')
names(meta_programs[[4]])[5] <- 'Cell Cycle'
names(meta_programs[[5]])[c(6, 10, 15)] <- c('Cell Cycle', 'Metal response', 'Lipid metabolism')
names(meta_programs[[6]])[c(3, 9, 12)] <- c('Cell Cycle', 'Proteasomal degradation', 'Unfolded protein response')
names(meta_programs[[7]]) <- paste('CD4 -', names(meta_programs[[7]]))
names(meta_programs[[7]])[c(1, 3)] <- c('CD4 - Treg', 'CD4 - Cell Cycle')
names(meta_programs[[8]]) <- paste('CD8 -', names(meta_programs[[8]]))
names(meta_programs[[8]])[c(2, 10)] <- c('CD8 - Cell Cycle', 'CD8 - Heat shock')

# Merge CD4 and CD8 T cells:
meta_programs[[7]] <- cbind(meta_programs[[7]], meta_programs[[8]])
meta_programs <- meta_programs[1:7]

cts <- c('Malignant', 'B_cell', 'Endothelial', 'Epithelial', 'Fibroblast', 'Macrophage', 'T_cell')
names(meta_programs) <- cts





rdir <- paste0('../study_plots/', gsub('/', '-', r[2]), '/', r[1])

gene_ave <- fread(paste0(rdir, '/gene_ave.csv'))[!is.na(ave)]

if(nrow(gene_ave) > 0) {
    bins_nm <- fread('../data/bins_nm.csv')
    gene_ave_list <- list(gene_ave[cell_type == 'Malignant' & n_cell >= 10, .(group, symbol, sample, cell_type, ave)])
    if(unique(gene_ave[cell_type == 'Malignant', .(group, sample, n_cell)])[, sum(n_cell) >= 10]) gene_ave_list <- c(gene_ave_list, list(
        gene_ave[
            cell_type == 'Malignant',
            .(group = 'all', sample = 'all', cell_type = 'Malignant', ave = sum(ave*n_cell)/sum(n_cell)),
            by = symbol
        ]
    ))
    gene_ave_list <- c(gene_ave_list, list(gene_ave[,
        if(unique(.SD[, .(cell_type, n_cell)])[, sum(n_cell) >= 10]) .SD[, .(cell_type = 'all', ave = sum(ave*n_cell)/sum(n_cell)), by = symbol],
        by = .(group, sample)
    ]))
    if(unique(gene_ave[, .(sample, cell_type, n_cell)])[, sum(n_cell) >= 10]) gene_ave_list <- c(gene_ave_list, list(
        gene_ave[, .(group = 'all', sample = 'all', cell_type = 'all', ave = sum(ave*n_cell)/sum(n_cell)), by = symbol]
    ))
    gene_ave <- rbindlist(gene_ave_list, use.names = TRUE)
    set.seed(7327) # Shuffle the genes before ranking so that the zero genes don't simply get ordered alphabetically:
    gene_ave[sample(1:.N), ave_rank := order(order(ave))/.N, by = .(group, sample, cell_type)]
    bins_m <- gene_ave[
        cell_type == 'Malignant',
        cbind(.SD[order(ave_rank)], bin = cut(1:.N, 15, labels = FALSE)),
        by = .(group, sample),
        .SDcols = -'cell_type'
    ]
    bins_all <- gene_ave[
        cell_type == 'all',
        cbind(.SD[order(ave_rank)], bin = cut(1:.N, 15, labels = FALSE)),
        by = .(group, sample),
        .SDcols = -'cell_type'
    ]
}





if(nrow(gene_ave) > 0) out <- lapply(paths$group, function(g) {
    
    cat('Preparing data\n')
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
    cells <- cells[cell_type %in% cts & col_nnz(expmat) >= 1000] # Filter cell types and remove low-complexity cells
    if(nrow(cells) < 30) return(NULL)
    expmat <- round(log_transform(1e+05*to_frac(expmat[, cells$cell_name])), 4) # Normalise to log TPM/10
    
    if(sum(unique(bins_nm$symbol) %in% genes) < 1000) bins_nm[, symbol := str_to_title(symbol)] # In case this is mouse data
    if(sum(unique(bins_nm$symbol) %in% genes) < 1000) return(NULL)
    bins_nm <- bins_nm[symbol %in% genes, .(group, symbol, bin = as.character(bin))]
    bins_m <- bins_m[symbol %in% genes & group %in% c(g, 'all'), .(group = mapvalues(sample, 'all', 'Malignant'), symbol, bin = as.character(bin))]
    bins_all <- bins_all[symbol %in% genes & group %in% c(g, 'all'), .(group = sample, symbol, bin = as.character(bin))]
    
    # Consolidate bins tables (and add <group> column to <cells>, for easy matching with <bins>):
    cat('Consolidating bins tables\n')
    if(all(unique(cells$sample) %in% bins_m$group)) {
        bins <- rbind(bins_nm[group %in% cells$cell_type], bins_m[group != 'Malignant'])
        cells[cell_type %in% c(unique(bins_nm$group), 'Malignant'), group := ifelse(cell_type == 'Malignant', sample, cell_type)]
    } else { # Use "Malignant" group for samples with too few malignant cells:
        bins <- rbind(bins_nm[group %in% cells$cell_type], bins_m)
        cells[
            cell_type %in% c(unique(bins_nm$group), 'Malignant'),
            group := ifelse(cell_type != 'Malignant' | !(sample %in% bins_m$group), cell_type, sample)
        ]
    }
    if(!all(cells[cell_type != 'Malignant', unique(cell_type)] %in% bins_nm$group)) {
        bins <- rbind(bins, bins_all[group == 'all']) # Use "all" group for non-malignant cell types not listed in bins_nm
        cells[cell_type != 'Malignant' & !(cell_type %in% bins_nm$group), group := 'all']
    }
    
    setkey(bins, group, symbol)
    bins_copy <- copy(bins); setkey(bins_copy, group, bin)
    
    cat('Computing bin averages\n')
    set.seed(9432)
    bin_aves <- cells[, {
        ct <- unique(cell_type)
        cat(ct, '\n')
        ct_genes <- sort(unique(unlist(meta_programs[[ct]])))
        if(sum(ct_genes %in% genes) < 50) ct_genes <- str_to_title(ct_genes) # In case this is mouse data
        ct_genes <- ct_genes[ct_genes %in% genes]
        if(length(ct_genes) < 50) NULL else {
            .SD[, {
                grp <- unique(group)
                cat('\t', grp, '\n')
                sdout <- bins_copy[
                    .(grp, bins[.(grp, ct_genes), as.character(bin)]),
                    if(.N >= 50) {
                        if(.N <= 100) {
                            .(gene = ct_genes[.GRP], bin_genes = symbol)
                        } else .(gene = ct_genes[.GRP], bin_genes = sample(symbol, 100, replace = FALSE))
                    },
                    by = .EACHI
                ]
                expmat_sub <- expmat[unique(sdout$bin_genes), cell_name, drop = FALSE]
                sdout[, .(cell_name = cell_name, ave = colMeans(expmat_sub[bin_genes, , drop = FALSE])), by = gene]
            }, by = group]
        }
    }, by = cell_type]
    
    bin_ave_mat <- dcast(bin_aves[, .(gene, cell_name, ave)], gene ~ cell_name)[, set_rownames(as.matrix(.SD), gene), .SDcols = -'gene']
    
    cat('Computing scores\n')
    scores <- cells[, {
        ct <- unique(cell_type)
        cat(ct, '\n')
        c(
            .(cell_name = cell_name),
            slapply(names(meta_programs[[ct]]), function(mp_name) {
                cat('\t', mp_name, '\n')
                mp <- meta_programs[[ct]][[mp_name]]
                if(sum(mp %in% genes & mp %in% bin_aves$gene) < 10) mp <- str_to_title(mp)
                mp <- mp[mp %in% genes & mp %in% bin_aves$gene]
                if(length(mp) >= 30) {
                    xout <- expmat[mp, cell_name, drop = FALSE] - bin_ave_mat[mp, cell_name, drop = FALSE]
                    xout[xout > 3] <- 3; xout[xout < -3] <- 3 # Cap relative expression levels to reduce influence of extreme values.
                    colMeans(xout, na.rm = TRUE) # Vector of scores - values for some cells/groups may be NA/NaN.
                }
            })
        ) %>% as.data.table %>% melt(id.vars = 'cell_name', variable.name = 'meta_program', value.name = 'score', variable.factor = FALSE)
    }, by = cell_type]
    
    scores <- scores[, -'cell_type']
    scores[, group := g]
    setcolorder(scores, 'group')
    return(scores)
    
}) %>% rbindlist

if(nrow(gene_ave) > 0) fwrite(out, paste0(rdir, '/mp_scores_all.csv'))
