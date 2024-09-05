library(data.table)
library(magrittr)
library(matkot)

paths_table <- fread('../data/paths_table.csv', encoding = 'UTF-8', key = c('study', 'cancer_type'))

for(r in transpose(as.list(unique(paths_table[, .(study, cancer_type)])))) tryCatch({
    cat(r, '-')
    d <- paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])
    if('data_dist.rds' %in% dir(d)) {
        data_dist <- readRDS(paste0(d, '/data_dist.rds'))
        data_dist_cond <- sapply(data_dist, function(x) is.null(x) | is.null(x$scores))
        if(!all(data_dist_cond)) {
            for(i in (1:length(data_dist))[!data_dist_cond]) { # Indices where data is not NULL
                out <- data_dist[[i]]$scores[, .(cell_name, sample, cell_type, meta_program, score = round(score, 4))]
                if(sum(!data_dist_cond) > 1) {
                    if(!('MP scores' %in% dir(d))) dir.create(paste0(d, '/MP scores'))
                    suff <- paths_table[as.list(r)][i, if(group_name != '') group_name else paste0('group', i)]
                    fwrite(out, paste0(d, '/MP scores/MP_scores_', suff, '.csv'))
                } else fwrite(out, paste0(d, '/MP scores.csv'))
            }
        } else cat('No MP score data')
    } else cat('No MP score data')
    cat('\n')
}, error = function(e) print('ERROR\n'))
