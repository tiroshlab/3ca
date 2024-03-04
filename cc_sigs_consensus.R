library(data.table)
library(ggplot2)
library(magrittr)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'))

sigs <- lapply(transpose(as.list(unique(paths_table[cancer_type != 'Other/Models', .(study, cancer_type)]))), function(r) {
    cat(r, '\n')
    if(!('data_cc.RDS' %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])))) return(NULL)
    plot_data <- readRDS(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/data_cc.RDS'))
    nullcond <- sapply(plot_data, function(x) ifelse(is.null(x), TRUE, all(sapply(x[names(x) != 'path'], is.null))))
    if(all(nullcond)) return(NULL)
    paths <- apply(paths_table[as.list(r), .(cells, genes, expmat)], 1, as.list, simplify = FALSE)
    rout <- lapply(which(!nullcond), function(i) list(r = r, group = i, g1s = plot_data[[i]]$g1s, g2m = plot_data[[i]]$g2m))
}) %>% unlist(recursive = FALSE)

g1s_tab <- table(unlist(lapply(sigs, `[[`, 'g1s')))
g2m_tab <- table(unlist(lapply(sigs, `[[`, 'g2m')))

sigs_cons <- list(g1s = names(g1s_tab)[order(-g1s_tab)][1:50], g2m = names(g2m_tab)[order(-g2m_tab)][1:50])

saveRDS(sigs_cons, '../data/cc_sigs_consensus.rds')
