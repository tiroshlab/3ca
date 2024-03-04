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





paths <- apply(paths_table[as.list(r), .(cells, genes, expmat)], 1, as.list, simplify = FALSE)

out <- lapply(paths, function(p) {
    
    if(!endsWith(p$expmat, 'mtx')) return(NULL)
    
    cat('Preparing data\n')
    cells <- suppressWarnings(fread(p$cells, na.strings = ''))
    if(!all(c('cell_type', 'sample') %in% names(cells))) return(NULL)
    cells$cell_name <- as.character(cells$cell_name)
    cells$sample <- as.character(cells$sample)
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
    if(nrow(cells) < 30) return(NULL)
    expmat <- round(log_transform(1e+05*to_frac(expmat[, cells$cell_name])), 4) # Normalise to log TPM/10
    
    cat('Assigning MP to each cell\n')
    set.seed(4387)
    scores <- cells[
        cell_type %in% cts,
        if(.N >= 100) c( # Cut-off of 100 cells of this cell type in this sample
            .(cell_name = cell_name),
            slapply(
                meta_programs[[unique(cell_type)]],
                function(mp) {
                    if(sum(mp %in% rownames(expmat)) < 10) mp <- str_to_title(mp)
                    if(sum(mp %in% rownames(expmat)) >= 30) sig_score(expmat[, cell_name], mp[mp %in% rownames(expmat)], nbin = 50, n = 50)
                }
            )
        ) %>% as.data.table %>% melt(id.vars = 'cell_name', variable.name = 'meta_program', value.name = 'score', variable.factor = FALSE),
        by = .(cell_type, sample)
    ]
    if(nrow(scores) == 0) return(NULL)
    # Assign MP to each cell (take the one with highest score):
    assignments <- scores[, .(max_score = max(score), meta_program = meta_program[which.max(score)]), by = .(cell_type, sample, cell_name)]
    assignments[, a_ct := ifelse(max_score < 1, as.character(NA), meta_program)]
    # In each cell type and sample, remove MPs that were assigned to less than 5% of those cells that were given an assignment (i.e. that had max
    # score >= 1):
    assignments[!is.na(a_ct), a_smpl := {a <- table(a_ct)/.N; ifelse(a_ct %in% names(a)[a < 0.05], NA, a_ct)}, by = .(cell_type, sample)]
    # Similar but more lenient criterion for a_ct, to get rid of spurious assignments:
    assignments[!is.na(a_ct), a_ct := {a <- table(a_ct)/.N; ifelse(a_ct %in% names(a)[a < 0.01], NA, a_ct)}, by = .(cell_type, sample)]
    setkey(assignments, cell_type, sample)
    
    cat('Constructing plot data\n')
    pie_data <- slapply(cts[cts %in% assignments$cell_type], function(ct) {
        pie_out <- assignments[ct, {N_ct <- .N; .SD[, .(prop = .N/N_ct), by = a_ct]}][order(-prop)]
        pie_out[is.na(a_ct), a_ct := 'Unassigned']
        pie_out[prop < 0.03, a_ct := 'Rare MPs (each <3%)']
        pie_out[a_ct == 'Rare MPs (each <3%)', prop := sum(prop)]
        pie_out <- unique(pie_out)
        pie_out[, a_ct := factor(a_ct, levels = c(a_ct[!(a_ct %in% c('Rare MPs (each <3%)', 'Unassigned'))], 'Rare MPs (each <3%)', 'Unassigned'))]
        pie_out <- pie_out[order(a_ct)]
        pie_out[, cumprop := c(0, cumsum(prop)[-.N])]
        return(pie_out)
    })
    bar_data <- slapply(cts[cts %in% assignments$cell_type], function(ct) {
        slapply(assignments[ct, unique(sample)], function(smpl) {
            bar_out <- assignments[.(ct, as.character(smpl))][, {N_smpl <- .N; .SD[, .(prop = .N/N_smpl), by = a_smpl]}][!is.na(a_smpl)][order(prop)]
            if(nrow(bar_out) == 0) return(NULL)
            bar_out[, a_smpl := factor(a_smpl, levels = a_smpl)]
            return(bar_out)
        })
    })
    bar_data <- slapply(bar_data, function(x) x[!sapply(x, is.null)])
    htmp_data <- slapply(cts[cts %in% assignments$cell_type], function(ct) {
        slapply(assignments[ct, unique(sample)], function(smpl) {
            if(!(smpl %in% names(bar_data[[ct]]))) return(NULL)
            mps <- bar_data[[ct]][[as.character(smpl)]][, levels(a_smpl)]
            g <- slapply(meta_programs[[ct]][, ..mps], function(mp) {
                if(sum(mp %in% rownames(expmat)) < 10) mp <- str_to_title(mp)
                mp[mp %in% rownames(expmat)]
            })
            asub <- assignments[!is.na(a_smpl)][.(ct, as.character(smpl))]
            if(nrow(asub) > 100) asub <- asub[, {n <- floor(100*.N/nrow(asub)); .SD[sample(1:.N, n, replace = FALSE)]}, by = a_smpl]
            htmp_out_list <- slapply(mps, function(mp_name) {
                mat <- expmat[g[[mp_name]], asub$cell_name, drop = FALSE]
                mat_sub <- mat[, asub[a_smpl == mp_name, cell_name], drop = FALSE]
                dt <- as.data.table(as.matrix(mat - rowMeans(mat)), keep.rownames = 'gene') # centre per gene and convert to data.table
                dt <- melt(dt, id.vars = 'gene', variable.name = 'cell_name', value.name = 'value')[, a_smpl := mp_name]
                return(list(data = dt, genes = names(sort(rowMeans(mat_sub))), cells = names(sort(colMeans(mat_sub)))))
            })
            htmp_out <- rbindlist(lapply(htmp_out_list, `[[`, 'data'))
            htmp_out[, c('gene_a', 'a_smpl') := .(paste(gene, gsub(' |-|/|\\(|\\)', '', a_smpl), sep = '_'), factor(a_smpl, levels = rev(mps)))]
            cells_ord <- unlist(lapply(htmp_out_list[rev(mps)], `[[`, 'cells'))
            genes_ord <- unlist(lapply(mps, function(mp_name) paste(htmp_out_list[[mp_name]]$genes, gsub(' |-|/|\\(|\\)', '', mp_name), sep = '_')))
            htmp_out[, c('gene_a', 'cell_name') := .(factor(gene_a, levels = genes_ord), factor(cell_name, levels = cells_ord))]
            htmp_out[, gene_num := as.numeric(gene_a)]
            return(htmp_out)
        })
    })
    htmp_data <- slapply(htmp_data, function(x) x[!sapply(x, is.null)])
    
    cat('Done!\n')
    return(list(scores = scores, assignments = assignments, pie_data = pie_data, bar_data = bar_data, heatmap_data = htmp_data, path = p))
    
})

if(!all(sapply(out, function(x) is.null(x$assignments)))) {
    
    ct <- gsub('/', '-', r[2])
    if(!(ct %in% dir('../data/study_plots'))) {dir.create(paste0('../data/study_plots/', ct))}
    if(!(r[1] %in% dir(paste0('../data/study_plots/', ct)))) {dir.create(paste0('../data/study_plots/', ct, '/', r[1]))}
    
    saveRDS(out, paste0('../data/study_plots/', ct, '/', r[1], '/data_dist.rds'))
    
}
