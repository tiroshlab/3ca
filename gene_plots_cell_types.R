library(data.table)
library(magrittr)
library(ggplot2)
library(stringr)
library(Matrix)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'), encoding = 'UTF-8')





# For each dataset, find cell types that constitute at least 1% of the data:
cell_types_thresh <- lapply(transpose(as.list(unique(paths_table[, .(study, cancer_type)]))), function(r) {
        
    cat(r, '\n')
    
    cts <- lapply(paths_table[as.list(r), cells], function(p) {
        cells <- suppressWarnings(fread(p, na.strings = ''))
        if(!('cell_type' %in% names(cells))) return(NULL)
        cells <- cells[!is.na(cell_type) & cell_type != 'Unassigned']
        out <- cells[, .(N = .N/nrow(cells)), by = cell_type][N >= 0.01, cell_type]
        # Check if 'Epithelial' cells are definitely normal epithelial cells:
        if(
            'Epithelial' %in% out && (
                ('malignant' %in% names(cells) && cells[cell_type == 'Epithelial' & malignant == 'no', .N/nrow(cells)] < 0.01) |
                    !('Malignant' %in% cells$cell_type)
            )
        ) {out <- out[out != 'Epithelial']}
        if(length(out) > 0) return(out)
    })
    
    cts <- cts[!sapply(cts, is.null)]
    
    if(length(cts) > 0) return(Reduce(intersect, cts))
    
})

# Take cell types that pass the 1% threshold in at least 5 datasets:
cell_types <- table(unlist(cell_types_thresh[!sapply(cell_types_thresh, is.null)]))
cell_types <- names(cell_types)[cell_types >= 5 & names(cell_types) != '']

saveRDS(cell_types, '../data/gene_plots_cell_types.rds')
