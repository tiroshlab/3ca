# The command line arguments should supply study name and cancer type (in that order).
r = commandArgs(trailingOnly = TRUE)





library(data.table)
library(magrittr)
library(Matrix)
library(caTools)
library(stringr)
library(plyr)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'))

hgnc_complete_set <- fread('../data/hgnc_complete_set_2023-04-13.txt', key = 'ensembl_gene_id')
hgnc_complete_set <- hgnc_complete_set[!(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])]
alias_table <- make_alias_table(hgnc_complete_set)

gene_positions <- fread('../data/gene_positions.csv')

# This is how gene_positions was defined:
# library(biomaRt)
# ensembl_mart <- useMart('ensembl')
# ensembl_data <- useDataset('hsapiens_gene_ensembl', mart = ensembl_mart)
# gene_positions <- as.data.table(getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position'), mart = ensembl_data))
# fwrite(gene_positions, '../data/gene_positions.csv')

gene_positions[, c('symbol', 'location') := do.call(`[`, list(hgnc_complete_set, ensembl_gene_id))[, .(symbol, location_sortable)]]
gene_positions <- gene_positions[!is.na(symbol) & chromosome_name %in% c(as.character(1:22), 'X', 'Y')]
gene_positions[, chromosome_name := mapvalues(chromosome_name, as.character(1:9), paste0(0, 1:9))]
setkey(gene_positions, symbol)





paths <- apply(paths_table[as.list(r), .(cells, genes, expmat)], 1, as.list, simplify = FALSE)
ct <- gsub('/', '-', r[2])
out <- lapply(paths, function(p) {
    
    if(!endsWith(p$expmat, 'mtx')) return(list(cna_data = NULL, ref_cells = NULL, ref_range = NULL, path = p))
    
    cat('Preparing data\n')
    cells <- suppressWarnings(fread(p$cells, na.strings = ''))
    if(!('cell_type' %in% names(cells)) || length(unique(cells$cell_type)) == 1 || !('Malignant' %in% cells$cell_type)) {
        return(list(cna_data = NULL, ref_cells = NULL, ref_range = NULL, path = p))
    }
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
    genes <- genes[genes %in% gene_positions$symbol]; expmat <- expmat[genes, ]
    colnames(expmat) <- cells$cell_name
    cells <- cells[col_nnz(expmat) >= 1000] # Remove low-complexity cells
    if(nrow(cells) < 30) return(list(cna_data = NULL, ref_cells = NULL, ref_range = NULL, path = p))
    expmat <- round(log_transform(1e+05*to_frac(expmat[, cells$cell_name])), 4) # Normalise to log TPM/10
    
    # Get reference cells, prioritising Fibroblasts, then Endothelial, then Macrophage, then whatever is most abundant.
    cat('Defining reference\n')
    ref_cts_global <- cells[cell_type != 'Malignant', .(N = .N), keyby = cell_type][N >= 100]
    ref_cts_global <- rbind(
        ref_cts_global[c('Fibroblast', 'Endothelial', 'Macrophage')],
        ref_cts_global[!(cell_type %in% c('Fibroblast', 'Endothelial', 'Macrophage'))][order(-N)]
    )[!is.na(N), cell_type[complete.cases(.SD)][1:min(2, .N)]]
    if(length(ref_cts_global) < 2) return(list(cna_data = NULL, ref_cells = NULL, ref_range = NULL, path = p))
    ref_cells <- cells[cell_type %in% ref_cts_global, .(cell_name = cell_name, cell_type = cell_type, ref_type = 'global')]
    ref_cells_smpl <- cells[
        cell_type != 'Malignant',
        {
            ref_cts_smpl <- .SD[, .(N = .N), keyby = cell_type][N >= 100]
            ref_cts_smpl <- ref_cts_smpl[c('Fibroblast', 'Endothelial', 'Macrophage'), cell_type[complete.cases(.SD)]]
            if(length(ref_cts_smpl) >= 2) .SD[cell_type %in% ref_cts_smpl, .(cell_name = cell_name, cell_type = cell_type)]
        },
        by = .(ref_type = sample)
    ]
    if(nrow(ref_cells_smpl) > 0) ref_cells <- rbind(ref_cells, ref_cells_smpl, use.names = TRUE)
    rm(ref_cells_smpl)

    cat('Computing CNA values\n')
    rms <- rowMeans(expmat)
    genes_top <- head(names(rms)[order(-rms)], 5000)
    expmat <- expmat[genes_top, ]; rms <- rms[genes_top]
    expmat <- expmat - rms # Centre rows/genes
    expmat[expmat > 3] <- 3; expmat[expmat < -3] <- -3 # Restrict range to limit influence of extreme values
    expdt <- as.data.table(as.matrix(expmat), keep.rownames = 'gene')
    expdt <- melt(expdt, id.vars = 'gene', variable.name = 'cell_name', value.name = 'value', variable.factor = FALSE)
    expdt[, c('chr', 'start_pos') := gene_positions[gene, .(chromosome_name, start_position)]]
    expdt <- expdt[order(chr, start_pos)]
    expdt[, cna := log2(runmean(2^value, k = 100)), by = .(chr, cell_name)] # Computing the CNA values
    expdt[, cna := cna - median(cna), keyby = cell_name]
    ref_range <- ref_cells[, do.call(`[`, list(expdt, cell_name))[, .(mean_cna = mean(cna)), by = gene], by = .(ref_type, cell_type)]
    ref_range <- ref_range[, .(mn = min(mean_cna), mx = max(mean_cna)), keyby = .(ref_type, gene)]
    setkey(cells, cell_name)
    expdt[, sample := do.call(`[`, list(cells, cell_name))[, as.character(sample)]]

    cat('Correcting CNA values w.r.t. reference\n')
    expdt[, # Correcting the CNA values w.r.t. reference cells
        cna := {
            rt <- unique(sample)
            if(!(rt %in% ref_range$ref_type)) rt <- 'global'
            out <- .SD[, .(cna, do.call(`[`, list(ref_range, .(rt, gene)))[, .(mn, mx)])]
            out[, ifelse(cna > mx + 0.1, cna - mx - 0.1, ifelse(cna < mn - 0.1, cna - mn + 0.1, 0))]
        },
        by = sample
    ]
    
    cat('Done!\n\n')
    return(list(cna_data = expdt, ref_cells = ref_cells, ref_range = ref_range, path = p))
    
})

if(!all(sapply(out, function(x) is.null(x$cna_data)))) {

    ct <- gsub('/', '-', r[2])
    if(!(ct %in% dir('../data/study_plots'))) {dir.create(paste0('../data/study_plots/', ct))}
    if(!(r[1] %in% dir(paste0('../data/study_plots/', ct)))) {dir.create(paste0('../data/study_plots/', ct, '/', r[1]))}

    saveRDS(out, paste0('../data/study_plots/', ct, '/', r[1], '/data_cna.rds'))

}
