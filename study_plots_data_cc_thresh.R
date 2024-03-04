library(data.table)
library(magrittr)
library(matkot)

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'))

thresh <- lapply(transpose(as.list(unique(paths_table[, .(study, cancer_type)]))), function(r) {
    if(!('data_cc.RDS' %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])))) return(NULL)
    plot_data <- readRDS(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/data_cc.RDS'))
    nullcond <- sapply(plot_data, function(x) ifelse(is.null(x), TRUE, all(sapply(x[names(x) != 'path'], is.null))))
    if(all(nullcond)) return(NULL)
    lapply(which(!nullcond), function(i) {
        if('cell_type' %in% names(plot_data[[i]]$ccdata)) {
            plot_data[[i]]$ccdata[,
                .(g1s = min(g1s_score[lr_g1s]), g2m = min(g2m_score[lr_g2m]), study = r[1], cancer_type = r[2]),
                by = cell_type
            ]
        } else {
            plot_data[[i]]$ccdata[, .(cell_type = NA, g1s = min(g1s_score[lr_g1s]), g2m = min(g2m_score[lr_g2m]), study = r[1], cancer_type = r[2])]
        }
    }) %>% rbindlist
}) %>% rbindlist

# plot(thresh[is.finite(g1s), sort(g1s)]); abline(h = 0.7)
# plot(thresh[is.finite(g2m), sort(g2m)]); abline(h = 0.65)

for(r in transpose(as.list(unique(paths_table[, .(study, cancer_type)])))) {
    cat(r, '\n')
    if(!('data_cc.RDS' %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])))) next
    plot_data <- readRDS(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/data_cc.RDS'))
    nullcond <- sapply(plot_data, function(x) ifelse(is.null(x), TRUE, all(sapply(x[names(x) != 'path'], is.null))))
    if(all(nullcond)) next
    for(i in which(!nullcond)) {
        
        ccdata <- copy(plot_data[[i]]$ccdata)
        ccdata[, phase := 'Cycling']
        ccdata[g1s_score < 0.7 & g2m_score < 0.65, phase := 'Not cycling']
        
        # Define phases and intermediates by fold change threshold:
        ccdata[phase == 'Cycling' & g1s_score > g2m_score & g1s_score - 0.5 < 2*(g2m_score - 0.5), phase := 'Intermediate']
        ccdata[phase == 'Cycling' & g2m_score > g1s_score & g2m_score - 0.5 < 2*(g1s_score - 0.5), phase := 'Intermediate']
        ccdata[phase == 'Cycling', phase := c('G1/S', 'G2/M')[which.max(.(g1s_score, g2m_score))], by = cell_name]
        
        # ggplot(ccdata[cell_type == 'Malignant'], aes(x = g1s_score, y = g2m_score)) + geom_point(aes(colour = phase)) + theme_test()
        
        setkey(ccdata, cell_name)
        plot_data[[i]]$ccdata <- ccdata
        
        hdata <- copy(plot_data[[i]]$hdata)
        if('data.table' %in% class(hdata)) {
            setkey(hdata, cell_name)
            hdata[, phase := do.call(`[`, list(ccdata, cell_name))$phase]
            annot_cols <- c('cell_name', 'g1s_score', 'g2m_score', 'phase')
            hdata <- slapply(
                switch(
                    ('cell_type' %in% names(ccdata)) + 1,
                    list(ccdata$cell_name),
                    slapply(
                        ccdata[phase != 'Not cycling', .(N = .N), by = cell_type][N >= 50, cell_type],
                        function(ct) ccdata[cell_type == ct, cell_name]
                    )
                ),
                function(idvec) {
                    # Downsample cell IDs if one phase (probably non-cycling cells) takes up more than two thirds of the data:
                    cell_ids <- hdata[
                        idvec,
                        .(id = switch((.N > hdata[idvec][phase != p, 2*.N]) + 1, cell_name, sample(cell_name, hdata[idvec][phase != p, 2*.N]))),
                        by = .(p = phase)
                    ]$id
                    out_data <- copy(hdata)[cell_name %in% cell_ids, -'cycling']
                    # Take cumulative sum of numbers of cells for each phase, to be used in ordering the cells and plotting division
                    # lines in the heatmap:
                    setkey(out_data, phase)
                    const <- out_data[cell_name %in% cell_ids][c('G2/M', 'Intermediate', 'G1/S')][!is.na(cell_name), .(N = .N), by = phase]
                    const <- const[, setNames(c(0, cumsum(N)), c(phase, 'Not cycling'))]
                    # Get ordered cell ID list:
                    cell_ids <- out_data[c('G2/M', 'Intermediate', 'G1/S', 'Not cycling')][
                        !is.na(cell_name),
                        .(cell_name = cell_name, cell_order = order(-pmax(g1s_score, g2m_score)) + const[phase]),
                        by = phase
                    ][, cell_name[cell_order]]
                    # Melt heatmap data and order cells:
                    out_data <- melt(out_data, id.vars = annot_cols, variable.name = 'gene', value.name = 'exp_level')
                    out_data[, cell_name := factor(cell_name, levels = cell_ids)]
                    # Order the genes:
                    setkey(out_data, phase)
                    if(all(c('G1/S', 'G2/M') %in% out_data$phase)) {
                        out_data[, gene := factor(gene, levels = out_data[
                            c('G1/S', 'G2/M'),
                            .SD[gene %in% plot_data[[i]][[tolower(gsub('/', '', phase))]], .(ave_exp = mean(exp_level)), by = gene],
                            by = phase
                        ][, .(ordered_genes = gene[order(ave_exp)]), by = phase]$ordered_genes)]
                    } else { # Already checked there are enough cycling cells in idvec, so can assume that 'G1/S' or 'G2/M' is present
                        out_data[, gene := factor(gene, levels = out_data[phase %in% c('G1/S', 'G2/M'), .(ave_exp = mean(exp_level)), by = gene][,
                            Reduce(c, lapply(
                                c('G1/S', 'G2/M'),
                                function(p) .SD[gene %in% plot_data[[i]][[tolower(gsub('/', '', p))]]][order(ave_exp), gene]
                            ))
                        ])]
                    }
                    return(list(data = out_data, vline = const))
                }
            )
            plot_data[[i]]$hdata <- hdata
        } else plot_data[[i]]$hdata <- NULL
        
    }
    saveRDS(plot_data, paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/data_cc.RDS'))
}
