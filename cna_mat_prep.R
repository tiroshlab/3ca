# The command line arguments should supply study name and cancer type (in that order).
r = commandArgs(trailingOnly = TRUE)





library(data.table)
library(magrittr)
library(matkot)

paths_table <- fread('../data/paths_table.csv', encoding = 'UTF-8', key = c('study', 'cancer_type'))

d <- paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])
if('data_cna.rds' %in% dir(d)) {
    cna <- readRDS(paste0(d, '/data_cna.rds')) # Genes are already in order
    cna_cond <- sapply(cna, function(x) is.null(x) | is.null(x$cna_data))
    if(!all(cna_cond)) {
        for(i in (1:length(cna))[!cna_cond]) { # Indices where data is not NULL
            out <- cna[[i]]$cna_data[, .(gene, cell_name, cna = round(cna, 4))]
            out <- dcast(out, cell_name ~ gene)[, c('cell_name', out[, unique(gene)]), with = FALSE]
            if(sum(!cna_cond) > 1) {
                if(!('CNA matrix' %in% dir(d))) dir.create(paste0(d, '/CNA matrix'))
                suff <- paths_table[as.list(r)][i, if(group_name != '') group_name else paste0('group', i)]
                fwrite(out, paste0(d, '/CNA matrix/CNA_matrix_', suff, '.csv'))
            } else fwrite(out, paste0(d, '/CNA matrix.csv'))
        }
    } else cat('No CNA data\n')
} else cat('No CNA data\n')
