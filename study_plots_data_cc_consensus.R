# The command line arguments should supply study name and cancer type (in that order).
r = commandArgs(trailingOnly = TRUE)





library(data.table)
library(ggplot2)
library(magrittr)
library(Matrix)
library(stringr)
library(plyr)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'))

hgnc_complete_set <- fread('../data/hgnc_complete_set.txt', key = 'symbol')[
    !(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])
]
alias_table <- make_alias_table(hgnc_complete_set)

bins_nm <- fread('../data/bins_nm.csv')

# Study-specific bins for malignant cells and whole samples:
rdir <- paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])
gene_ave <- fread(paste0(rdir, '/gene_ave.csv'))
gene_ave_list <- list(gene_ave[cell_type == 'Malignant' & n_cell >= 10, .(group, symbol, sample, cell_type, ave)])
if(unique(gene_ave[cell_type == 'Malignant', .(group, sample, n_cell)])[, sum(n_cell) >= 10]) gene_ave_list <- c(gene_ave_list, list(
    gene_ave[cell_type == 'Malignant', .(group = 'all', sample = 'all', cell_type = 'Malignant', ave = sum(ave*n_cell)/sum(n_cell)), by = symbol]
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
gene_ave[, ave_rank := .SD[sample(1:.N), order(order(ave))[order(symbol)]]/.N, .SDcols = c('symbol', 'ave'), by = .(group, sample, cell_type)]
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

# Consensus cell cycle signatures:
genesets <- readRDS('../data/cc_sigs_consensus.rds')





paths <- apply(paths_table[as.list(r), .(cells, genes, expmat)], 1, as.list, simplify = FALSE)

out <- lapply(1:length(paths), function(paths_i) {
    
    p <- paths[[paths_i]]
    
    if(!endsWith(p$expmat, 'mtx')) return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
    
    # Prepare data:
    cat('Preparing data\n')
    cells <- suppressWarnings(fread(p$cells, na.strings = ''))
    cells$cell_name <- as.character(cells$cell_name) # In case cell names are the same as row numbers
    if('sample' %in% names(cells)) cells$sample <- as.character(cells$sample)
    if(cells[, .N > length(unique(cell_name))]) cells[, cell_name := paste(cell_name, .I, sep = '_')]
    genes <- fread(p$genes, header = FALSE)$V1
    expmat <- readMM(p$expmat)
    expmat <- expmat[genes != '', ]
    genes <- genes[genes != '']
    expmat <- expmat[genes %in% names(table(genes))[table(genes) == 1], ]
    genes <- genes[genes %in% names(table(genes))[table(genes) == 1]]
    genes <- update_symbols_fast(genes, alias_table) # Update gene symbols
    rownames(expmat) <- genes
    colnames(expmat) <- cells$cell_name
    cells <- cells[col_nnz(expmat) >= 1000] # Remove low-complexity cells
    if(nrow(cells) < 30) return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
    expmat <- round(log_transform(1e+05*to_frac(expmat[, cells$cell_name])), 4) # Normalise to log TPM/10
    
    cat('Obtaining cell cycle signatures\n')
    g1s <- genesets$g1s; g2m <- genesets$g2m
    
    # If too few genes are in expmat, it might be mouse data, so convert them (de-capitalising seems to work better than biomaRt!):
    if(sum(g1s %in% rownames(expmat)) < 10 | sum(g2m %in% rownames(expmat)) < 10) {
        g1s <- str_to_title(g1s)
        if(sum(g1s %in% rownames(expmat)) >= 30) {
            g2m <- str_to_title(g2m)
            if(sum(g2m %in% rownames(expmat)) >= 30) {
                g1s <- g1s[g1s %in% rownames(expmat)]
                g2m <- g2m[g2m %in% rownames(expmat)]
            } else return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
        } else return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
    } else {
        g1s <- g1s[g1s %in% rownames(expmat)]
        g2m <- g2m[g2m %in% rownames(expmat)]
    }
    
    bins_nm <- bins_nm[symbol %in% rownames(expmat)]
    bins_m <- bins_m[symbol %in% rownames(expmat) & group %in% c(paths_i, 'all'), -'group']; setnames(bins_m, 'sample', 'group')
    bins_all <- bins_all[symbol %in% rownames(expmat) & group %in% c(paths_i, 'all'), -'group']; setnames(bins_all, 'sample', 'group')
    binslist <- list(bins_nm, bins_m, bins_all)[sapply(list(bins_nm, bins_m, bins_all), function(x) nrow(x) > 0)]
    g1s <- g1s[g1s %in% Reduce(intersect, lapply(binslist, function(x) x$symbol))]
    g2m <- g2m[g2m %in% Reduce(intersect, lapply(binslist, function(x) x$symbol))]
    
    if(length(g1s) < 20 | length(g2m) < 20) return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
    
    if('cell_type' %in% names(cells)) {
        cells <- cells[!is.na(cell_type) & cell_type != 'Unassigned']
        cells <- cells[cell_type %in% cells[, .(N = .N), by = cell_type][N >= 30, cell_type]] # Only consider cell types with at least 30 cells
        expmat <- expmat[, cells$cell_name]
    }
    
    # Consolidate bins tables:
    cat('Consolidating bins tables\n')
    if('cell_type' %in% names(cells)) {
        if('sample' %in% names(cells)) {
            if(all(unique(cells$sample) %in% bins_m$group)) {
                bins <- rbind(bins_nm[group %in% cells$cell_type], bins_m[group != 'all'])
            } else {
                bins <- rbind(bins_nm[group %in% cells$cell_type], bins_m) # Use "Malignant" group for samples with too few malignant cells
            }
        } else {
            bins <- rbind(bins_nm[group %in% cells$cell_type], bins_m[group == 'all'])
        }
        if(!all(cells[cell_type != 'Malignant', unique(cell_type)] %in% bins_nm$group)) {
            bins <- rbind(bins, bins_all[group == 'all']) # Use "all" group for non-malignant cell types not listed in bins_nm
        }
    } else if('sample' %in% names(cells)) {
        bins <- bins_all # Use "all" group for samples with too few cells
    } else {
        bins <- bins_all[group == 'all']
    }
    setkey(bins, group, symbol)
    
    # Stop if any bin doesn't have enough genes to sample from:
    if(bins[, {
        bins_grp <- .SD[symbol %in% c(g1s, g2m), unique(bin)]
        .SD[bin %in% bins_grp, .(N = .N), by = bin]
    }, by = group][, any(N < 100)]) return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
    
    # Compute matrices of bin averages for the cell cycle signature genes, to be used for several computations:
    cat('Computing bin averages\n')
    bin_aves <- slapply(c('g1s', 'g2m'), function(x) t(sapply(get(x), function(g) bin_ave(g, bins, cells, expmat))))
    
    # Initial classification of cycling cells:
    cat('Initial classification of cycling cells\n')
    ccdata <- cc_data(cells, g1s, g2m, expmat, bin_aves)$ccdata
    
    ccdata[,
        c('phase_count', 'phase_score', 'phase_ave') := lapply(c('count', 'score', 'ave'), function(vt) switch(
            (get(paste0('g1s_', vt, '_padj')) < 0.05 | get(paste0('g2m_', vt, '_padj')) < 0.05) + 1,
            'Not cycling',
            'Cycling'
        )),
        by = cell_name
    ]
    
    ccdata[, phase_cons := switch(all(c(phase_count, phase_score) == 'Cycling') + 1, 'Not cycling', 'Cycling'), by = cell_name]
        
    # ggplot(ccdata[cell_type == 'Malignant'], aes(x = g1s_score, y = g2m_score)) + geom_point(aes(colour = phase_cons)) + theme_test()
    
    cat('Redefining cell cycle gene bins\n')
    ccbins <- data.table(grp = character(), symb = character(), r = numeric())
    if('cell_type' %in% names(cells)) {
        ct_change <- ccdata[, .(prop_cyc = sum(phase_cons == 'Cycling')/.N), by = cell_type][prop_cyc > 0.1, cell_type]
        if('sample' %in% names(cells)) {
            if('Malignant' %in% ct_change && ccdata[phase_cons == 'Not cycling' & cell_type == 'Malignant' & sample %in% bins$group, .N > 0]) {
                ccbins <- rbind(ccbins, ccdata[
                    phase_cons == 'Not cycling' & cell_type == 'Malignant' & sample %in% bins$group,
                    .(symb = bins[unique(sample), symbol],
                        r = order(order(rowMeans(expmat[bins[unique(sample), symbol], cell_name, drop = FALSE])))/nrow(expmat)),
                    by = .(grp = sample)
                ])
            }
            if('Malignant' %in% ct_change && ccdata[phase_cons == 'Not cycling' & cell_type == 'Malignant' & !(sample %in% bins$group), .N > 0]) {
                ccbins <- rbind(ccbins, ccdata[
                    phase_cons == 'Not cycling' & cell_type == 'Malignant' & !(sample %in% bins$group),
                    .(grp = 'Malignant', symb = bins['Malignant', symbol],
                        r = order(order(rowMeans(expmat[bins['Malignant', symbol], cell_name, drop = FALSE])))/nrow(expmat))
                ])
            }
            if(ccdata[phase_cons == 'Not cycling' & cell_type != 'Malignant' & cell_type %in% bins$group, .N > 0]) {
                ccbins <- rbind(ccbins, ccdata[
                    phase_cons == 'Not cycling' & cell_type != 'Malignant' & cell_type %in% bins$group & cell_type %in% ct_change,
                    .(symb = bins[unique(cell_type), symbol],
                        r = order(order(rowMeans(expmat[bins[unique(cell_type), symbol], cell_name, drop = FALSE])))/nrow(expmat)),
                    by = .(grp = cell_type)
                ])
            }
            if(ccdata[phase_cons == 'Not cycling' & cell_type != 'Malignant' & !(cell_type %in% bins$group), .N > 0]) {
                ccbins <- rbind(ccbins, ccdata[
                    phase_cons == 'Not cycling' & cell_type != 'Malignant' & !(cell_type %in% bins$group) & cell_type %in% ct_change,
                    .(grp = 'all', symb = bins['all', symbol],
                        r = order(order(rowMeans(expmat[bins['all', symbol], cell_name, drop = FALSE])))/nrow(expmat))
                ])
            }
            ccbins <- ccbins[symb %in% c(g1s, g2m)]
        } else {
            if(ccdata[phase_cons == 'Not cycling' & cell_type %in% bins$group, .N > 0]) {
                ccbins <- rbind(ccbins, ccdata[
                    phase_cons == 'Not cycling' & cell_type %in% bins$group & cell_type %in% ct_change,
                    .(symb = bins[unique(cell_type), symbol],
                        r = order(order(rowMeans(expmat[bins[unique(cell_type), symbol], cell_name, drop = FALSE])))/nrow(expmat)),
                    by = .(grp = cell_type)
                ])
            }
            if(ccdata[phase_cons == 'Not cycling' & !(cell_type %in% bins$group), .N > 0]) {
                ccbins <- rbind(ccbins, ccdata[
                    phase_cons == 'Not cycling' & !(cell_type %in% bins$group) & cell_type %in% ct_change,
                    .(grp = 'all', symb = bins['all', symbol],
                        r = order(order(rowMeans(expmat[bins['all', symbol], cell_name, drop = FALSE])))/nrow(expmat))
                ])
            }
            ccbins <- ccbins[symb %in% c(g1s, g2m)]
        }
    } else if('sample' %in% names(cells) & ccdata[, sum(phase_cons == 'Cycling')/.N > 0.1]) {
        if(ccdata[phase_cons == 'Not cycling' & sample %in% bins$group, .N > 0]) {
            ccbins <- rbind(ccbins, ccdata[
                phase_cons == 'Not cycling' & sample %in% bins$group,
                .(symb = bins[unique(sample), symbol],
                    r = order(order(rowMeans(expmat[bins[unique(sample), symbol], cell_name, drop = FALSE])))/nrow(expmat)),
                by = .(grp = sample)
            ])
        }
        if(ccdata[phase_cons == 'Not cycling' & !(sample %in% bins$group), .N > 0]) {
            ccbins <- rbind(ccbins, ccdata[
                phase_cons == 'Not cycling' & !(sample %in% bins$group),
                .(grp = 'all', symb = bins['all', symbol],
                    r = order(order(rowMeans(expmat[bins['all', symbol], cell_name, drop = FALSE])))/nrow(expmat))
            ])
        }
        ccbins <- ccbins[symb %in% c(g1s, g2m)]
    } else if(ccdata[, sum(phase_cons == 'Cycling')/.N > 0.1]) {
        ccbins <- rbind(ccbins, ccdata[
            phase_cons == 'Not cycling',
            .(grp = 'all', symb = bins$symbol, r = order(order(rowMeans(expmat[bins$symbol, cell_name, drop = FALSE])))/nrow(expmat))
        ][symb %in% c(g1s, g2m)])
    }
    
    if(nrow(ccbins) > 0) {
        
        setkey(ccbins, grp, symb)
        bin_range <- bins[, .(rr = range(ave_rank)), keyby = .(group, bin)]
        ccbins <- ccbins[, bin_range[grp, .(rm = r %between% rr), by = bin], by = .(grp, symb)][,
            .(bin = ifelse(any(rm), bin[rm], 1L)),
            keyby = .(grp, symb)
        ]
        bins[symbol %in% c(g1s, g2m), bin := switch((group %in% ccbins$grp) + 1, bin, ccbins[.(group, symbol), bin]), by = group]
        
        # Stop if any bin doesn't have enough genes to sample from:
        if(bins[, {
            bins_grp <- .SD[symbol %in% c(g1s, g2m), unique(bin)]
            .SD[bin %in% bins_grp, .(N = .N), by = bin]
        }, by = group][, any(N < 100)]) return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
        
        # Recompute matrices of bin averages for the cell cycle signature genes using the updated bins:
        cat('Recomputing bin averages\n')
        bin_aves <- slapply(c('g1s', 'g2m'), function(x) t(sapply(get(x), function(g) bin_ave(g, bins, cells, expmat))))
        
        # Final classification of cycling cells using cells initially classified as non-cycling to define null distribution:
        cat('Final classification of cycling cells\n')
        if('cell_type' %in% names(ccdata)) {
            null_ids <- ccdata[,
                .(ids = switch((sum(phase_cons == 'Cycling')/.N > 0.1) + 1, cell_name, cell_name[phase_cons == 'Not cycling'])),
                by = cell_type
            ]$ids
        } else {
            null_ids <- ccdata[, switch((sum(phase_cons == 'Cycling')/.N > 0.1) + 1, cell_name, cell_name[phase_cons == 'Not cycling'])]
        }
        ccdata <- cc_data(ccdata, g1s, g2m, expmat, bin_aves, null_ids)$ccdata
        
        ccdata[,
            c('phase_count', 'phase_score', 'phase_ave') := lapply(c('count', 'score', 'ave'), function(vt) switch(
                (get(paste0('g1s_', vt, '_padj')) < 0.05 | get(paste0('g2m_', vt, '_padj')) < 0.05) + 1,
                'Not cycling',
                'Cycling'
            )),
            by = cell_name
        ]
        
        ccdata[, phase_cons := switch(all(c(phase_count, phase_score) == 'Cycling') + 1, 'Not cycling', 'Cycling'), by = cell_name]
        
    }
    
    ccdata[,
        c('phase_g1s', 'phase_g2m') := .(all(c(g1s_score_padj, g1s_count_padj) < 0.05), all(c(g2m_score_padj, g2m_count_padj) < 0.05)),
        by = cell_name
    ]
    if('cell_type' %in% names(ccdata)) {
        ccdata[,
            c('lr_g1s', 'lr_g2m') := .(
                glm(phase_g1s ~ g1s_score, family = binomial(), .SD)$fitted.values > 0.5,
                glm(phase_g2m ~ g2m_score, family = binomial(), .SD)$fitted.values > 0.5
            ),
            by = cell_type
        ]
    } else {
        ccdata[,
            c('lr_g1s', 'lr_g2m') := .(
                glm(phase_g1s ~ g1s_score, family = binomial(), .SD)$fitted.values > 0.5,
                glm(phase_g2m ~ g2m_score, family = binomial(), .SD)$fitted.values > 0.5
            )
        ]
    }
    ccdata[, cycling := ifelse(lr_g1s | lr_g2m, 'Cycling', 'Not cycling')]
    ccdata[, c('phase_g1s', 'phase_g2m') := NULL]
    
    # Re-centre the scores relative to means of the non-cycling cells:
    if('cell_type' %in% names(ccdata)) {
        ccdata[,
            c('g1s_score', 'g2m_score') := lapply(.(g1s_score, g2m_score), function(x) x - mean(x[cycling == 'Not cycling'])),
            by = cell_type
        ]
    } else {
        ccdata[, c('g1s_score', 'g2m_score') := lapply(.(g1s_score, g2m_score), function(x) x - mean(x[cycling == 'Not cycling']))]
    }
    
    # ggplot(ccdata[cell_type == 'Malignant'], aes(x = g1s_score, y = g2m_score)) + geom_point(aes(colour = cycling)) + theme_test()
    
    # Define separate data for heatmap:
    cat('Constructing heatmap data\n')
    annot_cols <- c('cell_name', 'g1s_score', 'g2m_score', 'cycling')
    if('cell_type' %in% names(ccdata)) {
    # if(!is.null(ctvar)) {
        hdata <- ccdata[,
            merge(.SD, as.data.table(apply(expmat[c(g1s, g2m), cell_name], 1, function(x) x - mean(x)), keep.rownames = 'cell_name')),
            .SDcols = annot_cols,
            by = cell_type
        ][, -'cell_type']
    } else {
        hdata <- ccdata[,
            merge(.SD, as.data.table(apply(expmat[c(g1s, g2m), ], 1, function(x) x - mean(x)), keep.rownames = 'cell_name')),
            .SDcols = annot_cols
        ]
    }
    
    cat('Done!\n\n')
    return(list(ccdata = ccdata, hdata = hdata, g1s = g1s, g2m = g2m, path = p))
    
})

if(!all(sapply(out, function(x) is.null(x$ccdata)))) {
    ct <- gsub('/', '-', r[2])
    if(!(ct %in% dir('../data/study_plots'))) {dir.create(paste0('../data/study_plots/', ct))}
    if(!(r[1] %in% dir(paste0('../data/study_plots/', ct)))) {dir.create(paste0('../data/study_plots/', ct, '/', r[1]))}
    saveRDS(out, paste0('../data/study_plots/', ct, '/', r[1], '/data_cc_consensus.rds'))
}
