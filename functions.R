make_alias_table <- function(hgnc_complete_set) {
    out <- hgnc_complete_set[,
        rbind(
            data.table(symbol_alt = str_split(alias_symbol, '\\|')[[1]], type = 'alias'),
            data.table(symbol_alt = str_split(prev_symbol, '\\|')[[1]], type = 'prev')
        ),
        by = symbol
    ][symbol_alt == '', symbol_alt := NA]
    setkey(out, symbol_alt)
    return(out)
}

update_symbols_fast <- function(symbols, alias_table) {
    if(length(symbols) != length(unique(symbols))) warning('Not all symbols are unique!')
    setkey(alias_table, symbol_alt)
    symbols[grep('C[0-9]+ORF[0-9]+', symbols)] <- gsub('ORF', 'orf', symbols[grep('C[0-9]+ORF[0-9]+', symbols)])
    cond <- !(symbols %in% alias_table$symbol)
    # The following addition to <cond> ensures two different symbols don't map to the same "prev" symbol...
    cond <- cond & !(symbols %in% alias_table[type == 'prev'][symbols[cond]][symbol %in% names(table(symbol))[table(symbol) > 1], symbol_alt])
    # ...And the following makes sure we have all unique prev symbols:
    cond <- cond & !(symbols %in% alias_table[type == 'prev'][symbols[cond], names(table(symbol_alt))[table(symbol_alt) > 1]])
    if(any(cond)) {
        # alias_table[type == 'prev'][symbols[cond]]
        symbols[cond] <- alias_table[type == 'prev'][symbols[cond]][,
            .(newsymb = { # Match symbols that have at least one "prev" symbol that isn't already in symbols:
                cands <- symbol[!is.na(symbol) & !(symbol %in% symbols)]
                ifelse(length(cands) == 1, cands, symbol_alt) # Convert only those symbols that have exactly one such "prev" symbol
            }),
            by = symbol_alt
        ]$newsymb
    }
    cond <- !(symbols %in% alias_table$symbol)
    # As above, the following additions to <cond> ensure two different symbols don't map to the same "alias", and we have all unique aliases:
    cond <- cond & !(symbols %in% alias_table[type == 'alias'][symbols[cond]][symbol %in% names(table(symbol))[table(symbol) > 1], symbol_alt])
    cond <- cond & !(symbols %in% alias_table[type == 'alias'][symbols[cond], names(table(symbol_alt))[table(symbol_alt) > 1]])
    if(any(cond)) {
        symbols[cond] <- alias_table[type == 'alias'][symbols[cond]][,
            .(newsymb = { # Match symbols that have at least one alias that isn't already in symbols:
                cands <- symbol[!is.na(symbol) & !(symbol %in% symbols)]
                ifelse(length(cands) == 1, cands, symbol_alt) # Convert only those symbols that have exactly one such alias
            }),
            by = symbol_alt
        ]$newsymb
    }
    return(symbols)
}





dirr <- function(x) dir(x, recursive = TRUE)





ct_comp_data <- function(expmat, cells, markers, max_genes = 24) {

    if('cell_type' %in% names(cells)) {
        pie_data <- cells[, .(N = .N/nrow(cells)), by = .(cell_type = gsub('_', ' ', mapvalues(cell_type, NA, 'Unassigned', warn_missing = FALSE)))]
    } else if('cell_subtype' %in% names(cells)) {
        pie_data <- cells[, .(N = .N/nrow(cells)), by = .(cell_type = gsub('_', ' ', mapvalues(cell_subtype, NA, 'Unassigned', warn_missing = FALSE)))]
    } else if('malignant' %in% names(cells)) {
        pie_data <- cells[, .(N = .N/nrow(cells)), by = .(cell_type = mapvalues(malignant, c('yes', 'no'), c('Malignant', 'Non-malignant')))]
    } else {
        return(NULL)
    }

    pie_data[, cell_type := factor(cell_type, levels = cell_type[order(-N)])]
    pie_data <- pie_data[order(cell_type)][, cumN := c(0, cumsum(N)[-.N])]

    if('cell_type' %in% names(cells)) {
        
        markers <- markers[gene %in% rownames(expmat)]

        if(is.null(key(markers))) {
            setkey(markers, cell_type)
            warning("No key set for <markers>: setting key to 'cell_type'.")
        }

        immune_cts <- c('B_cell', 'Basophil', 'Dendritic', 'Macrophage', 'Mast', 'Monocyte', 'Neutrophil', 'NK_cell', 'Plasma', 'T_cell')

        cts <- cells[cell_type %in% markers$cell_type, .(N = .N), by = cell_type][N >= 10][order(-N), cell_type] # Cell types with >= 10 cells
        if(length(cts) == 0) return(list(pie_data = pie_data))
        if(any(immune_cts %in% cts)) {
            first_immune <- min(which(cts %in% immune_cts))
            if(first_immune == 1) {
                cts <- unique(c('Immune', cts))
            } else {
                cts <- unique(c(cts[1:(first_immune - 1)], 'Immune', cts[first_immune:length(cts)]))
            }
        }
        
        # Take top 2 genes for each cell type:
        genes <- markers[cts, .(g = gene[1:min(.N, 2)]), by = cell_type]

        marker_plot_data <- cells[
            cell_type %in% cts,
            do.call(
                function(x) c(
                    setNames(markers[cts][gene %in% genes$g], c('marker_ct', 'gene')),
                    .(ave = rowMeans(x), percent_pos = 100*row_nnz(x)/ncol(x))
                ),
                args = list(x = expmat[genes$g, cell_name, drop = FALSE])
            ),
            by = cell_type
        ]

        if('Immune' %in% cts) {
            marker_plot_data <- marker_plot_data[
                gene %in% c(
                    markers['Immune', gene],
                    marker_plot_data[marker_ct == cell_type, .(gene = gene, keep = percent_pos >= 50)][, gene[keep]]
                )
            ]
        } else {
            marker_plot_data <- marker_plot_data[
                gene %in% marker_plot_data[marker_ct == cell_type, .(gene = gene, keep = percent_pos >= 50)][, gene[keep]]
            ]
        }
        
        # If necessary, iteratively remove one gene per cell type until the total number is equal to <max_genes>:
        while(
            marker_plot_data[, length(unique(gene)) > max_genes] & !all(marker_plot_data[marker_ct == cell_type, .(N = .N), by = cell_type]$N == 1)
        ) {
            marker_plot_data <- marker_plot_data[
                !(gene %in% marker_plot_data[
                    marker_ct == cell_type & cell_type %in% marker_plot_data[marker_ct == cell_type, .(N = .N), by = cell_type][N > 1, cell_type],
                    .(gene = gene, r = order(order(apply(data.table(r1 = order(order(ave)), r2 = order(order(percent_pos))), 1, mean))))
                ][order(r), gene[1]])
            ]
        }

        marker_plot_data[,
            c('cell_type', 'gene') := .(factor(gsub('_', ' ', cell_type), levels = gsub('_', ' ', cts)), factor(gene, levels = unique(gene)))
        ]

        return(list(pie_data = pie_data, marker_plot_data = marker_plot_data))

    } else {
        return(list(pie_data = pie_data))
    }

}





ct_comp_plot <- function(pie_data, marker_plot_data = NULL, colours = NULL, legends_arrange = c('vertical', 'horizontal')) {
    
    legends_arrange <- match.arg(legends_arrange)

    pie <- ggplot(pie_data, aes(x = '', y = N, fill = cell_type)) +
        geom_bar(width = 1, stat = 'identity') +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_manual(
            name = 'Cell type',
            values = switch(
                is.null(colours) + 1,
                colours,
                switch(('randomcoloR' %in% .packages(TRUE)) + 1, ggplot_colours(nrow(pie_data)), distinctColorPalette(nrow(pie_data)))
            )
        ) +
        coord_polar('y', start = 0, direction = -1) +
        geom_text(data = pie_data[N >= 0.05], aes(x = 1.3, y = 1 - cumN - N/2, label = percent(N, accuracy = 0.1))) +
        theme(
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.background = element_rect(fill = NA)
        )

    out <- list(pie = pie)

    if(!is.null(marker_plot_data)) {

        marker_plot <- ggplot(marker_plot_data, aes(x = cell_type, y = gene)) +
            geom_point(aes(fill = ave, size = percent_pos), shape = 21, colour = 'black') +
            scale_fill_gradientn(
                colours = brewer.pal(9, 'YlOrRd'),
                limits = c(2, 8),
                breaks = c(2, 4, 6, 8),
                labels = c('2' = '\u2264 2', '4' = '4', '6' = '6', '8' = '\u2265 8'),
                oob = squish
            ) +
            scale_radius(range = c(2, 7), limits = c(50, 100), breaks = c(50, 75, 100), labels = c('50' = '50', '75' = '75', '100' = '100')) +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle = 55, hjust = 1),
                axis.title.x = element_text(margin = margin(t = 10)),
                axis.title.y = element_text(margin = margin(r = 10)),
                panel.grid = element_line(size = 0.4),
                strip.background = element_blank(),
                strip.text = element_blank(),
                legend.box = legends_arrange
            ) +
            guides(fill = guide_colourbar(frame.colour = 'black')) +
            labs(x = 'Cell type', y = 'Gene', fill = 'Mean', size = '% expressing\ncells')

        out <- c(out, list(marker_plot = marker_plot))

    }

    return(out)

}





ct_umap_data <- function(expmat, cells) {
    if(!any(c('cell_type', 'sample') %in% names(cells))) return(NULL)
    if(ncol(expmat) > 100) {
        exp_pca <- irlba(t(apply(expmat, 1, function(x) x - mean(x))), nv = 50)
        exp_umap <- umap(exp_pca$v, n_neighbors = max(min(100, floor(ncol(expmat)/100)), 15), spread = 5, min_dist = 0.1)
    } else {
        exp_umap <- umap(as.matrix(t(expmat)), n_neighbors = 15, spread = 5, min_dist = 0.1)
    }
    clust_data <- copy(cells)[, .SD, .SDcols = c('cell_name', c('cell_type', 'sample')[c('cell_type', 'sample') %in% names(cells)])]
    if('sample' %in% names(clust_data)) clust_data$sample <- as.character(clust_data$sample)
    if('cell_type' %in% names(clust_data)) clust_data[, cell_type := gsub('_', ' ', mapvalues(cell_type, NA, 'Unassigned', warn_missing = FALSE))]
    setkey(clust_data, cell_name)
    clust_data <- clust_data[colnames(expmat)]
    clust_data[, c('UMAP 1', 'UMAP 2') := as.data.table(exp_umap)]
    return(clust_data)
}





ct_umap_plot <- function(data, col, colours = NULL, title = NULL, legend_title = waiver()) {
    ggplot(data, aes(x = `UMAP 1`, y = `UMAP 2`, colour = get(col))) +
        geom_point(size = 0.5) +
        scale_colour_manual(
            name = legend_title,
            values = switch(
                is.null(colours) + 1,
                colours,
                switch(
                    ('randomcoloR' %in% .packages(TRUE)) + 1,
                    ggplot_colours(length(unique(data[[col]]))),
                    randomcoloR::distinctColorPalette(length(unique(data[[col]])))
                )
            )
        ) +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        theme_test() +
        labs(title = title)
}





cc_plot <- function(ccdata, hdata, g1s, g2m, min_pts = 30, ngene = 15) {

    if(nrow(ccdata) < min_pts) return(NULL)
    
    if(!('cell_type' %in% names(ccdata))) {
        
        scatter <- ggplot(ccdata, aes(x = g1s_score, y = g2m_score, colour = phase)) +
            scale_colour_manual(values = setNames(c('grey', brewer.pal(3, 'Dark2')), c('Not cycling', 'G1/S', 'Intermediate', 'G2/M'))) +
            geom_point(alpha = 0.5) +
            annotate(
                'text',
                x = ccdata[, min(g1s_score) + 0.85*diff(range(g1s_score))],
                y = ccdata[, min(g2m_score) + 0.95*diff(range(g2m_score))],
                label = paste(percent(ccdata[, nrow(.SD[phase != 'Not cycling'])/.N], accuracy = 0.1), 'cycling cells')
            ) +
            theme_test()
        
        htmp_pair <- slapply(
            c('g1s', 'g2m'),
            function(x) {
                ggplot(hdata[[1]]$data[gene %in% get(x)], aes(x = cell_name, y = gene, fill = exp_level)) +
                    geom_raster() +
                    scale_fill_gradientn(
                        colours = brewer.pal(9, 'Reds'),
                        limits = c(1, 5),
                        breaks = 1:5,
                        labels = c('1' = '\u2264 1', '2' = '2', '3' = '3', '4' = '4', '5' = '\u2265 5'),
                        oob = squish
                    ) +
                    scale_x_discrete(expand = c(0, 0)) +
                    scale_y_discrete(expand = c(0, 0)) +
                    geom_vline(xintercept = hdata[[1]]$vline[-1] + 0.5) +
                    theme(
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.ticks.length = unit(0, 'pt'),
                        axis.title = element_blank(),
                        panel.border = element_rect(fill = NA, size = 1),
                        legend.justification = c(0, 0.2),
                        plot.margin = unit(c(0, 5.5, 0, 0), 'pt')
                    ) +
                    labs(fill = 'Relative\nexpression\nlevel')
            }
        )
        
        annots <- slapply(
            c('g1s', 'g2m'),
            function(sig) {
                out <- tableGrob(
                    data.table(a = c(rev(tail(levels(hdata[[1]]$data$gene)[levels(hdata[[1]]$data$gene) %in% get(sig)], ngene)), '')),
                    theme = ttheme_minimal(
                        rowhead = list(padding = unit(c(0, 0), "pt"), fg_params = list(cex = 0)),
                        colhead = list(padding = unit(c(0, 0), "pt"), fg_params = list(cex = 0)),
                        core = list(padding = unit(c(5, 5), "pt"), fg_params = list(hjust = 1, x = 0.9)),
                        base_size = 9
                    )
                )
                out <- gtable_add_grob(
                    out,
                    grobs = segmentsGrob(x0 = unit(0, "npc"), y0 = unit(1, "npc"), x1 = unit(1, "npc"), y1 = unit(1, "npc"), gp = gpar(lwd = 2.0)),
                    t = 2, b = 2, l = 1, r = 2
                )
                out <- gtable_add_grob(
                    out,
                    grobs = segmentsGrob(x0 = unit(0, "npc"), y0 = unit(0, "npc"), x1 = unit(1, "npc"), y1 = unit(0, "npc"), gp = gpar(lwd = 2.0)),
                    t = ngene + 2, b = ngene + 2, l = 1, r = 2
                )
                out$widths <- unit(c(0, 1), 'npc')
                out$heights[length(out$heights)] <- unit(1, 'npc') - sum(out$heights[-length(out$heights)])
                out
            }
        )
        
        phase_bar <- ggplot(hdata[[1]]$data, aes(x = cell_name, y = 0, fill = phase)) +
            geom_raster() +
            scale_fill_manual(values = setNames(c('grey', brewer.pal(3, 'Dark2')), c('Not cycling', 'G1/S', 'Intermediate', 'G2/M'))) +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, 'pt'),
                axis.title = element_blank(),
                legend.justification = c(0, 1),
                plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
            )
        
        return(
            list(
                scatter = scatter,
                heatmap = plot_grid(
                    plot_grid(
                        blank_plot(),
                        blank_plot(),
                        phase_bar + theme(legend.position = 'none'),
                        ggplot() + theme(panel.background = element_blank()) + labs(y = 'G2/M genes'),
                        annots$g2m,
                        htmp_pair$g2m + theme(legend.position = 'none'),
                        ggplot() + theme(panel.background = element_blank()) + labs(y = 'G1/S genes'),
                        annots$g1s,
                        htmp_pair$g1s + theme(legend.position = 'none'),
                        blank_plot(),
                        blank_plot(),
                        ggplot() +
                            theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5, size = 11)) +
                            labs(title = 'Cells'),
                        nrow = 4,
                        ncol = 3,
                        rel_widths = c(0.5, 1, 8.5),
                        rel_heights = c(0.04, 0.45, 0.45, 0.06)
                    ),
                    plot_grid(get_legend(phase_bar), get_legend(htmp_pair$g1s), nrow = 2, ncol = 1),
                    nrow = 1,
                    ncol = 2,
                    rel_widths = c(0.84, 0.16)
                )
            )
        )
        
    } else {
        
        cts <- ccdata[, .(N = .N), by = cell_type][N >= min_pts, cell_type]
        if(length(cts) == 0) return(NULL)
        
        lim_x <- ccdata[cell_type %in% cts, c(min(g1s_score), max(g1s_score))]
        lim_y <- ccdata[cell_type %in% cts, c(min(g2m_score), max(g2m_score))]
        scatter <- slapply(
            cts,
            function(ct) ggplot(ccdata[cell_type == ct], aes(x = g1s_score, y = g2m_score, colour = phase)) +
                scale_x_continuous(limits = lim_x) +
                scale_y_continuous(limits = lim_y) +
                scale_colour_manual(values = setNames(c('grey', brewer.pal(3, 'Dark2')), c('Not cycling', 'G1/S', 'Intermediate', 'G2/M'))) +
                geom_point(alpha = 0.5) +
                annotate(
                    'text',
                    x = min(lim_x) + 0.85*diff(range(lim_x)),
                    y = min(lim_y) + 0.95*diff(range(lim_y)),
                    label = paste(percent(ccdata[cell_type == ct, nrow(.SD[phase != 'Not cycling'])/.N], accuracy = 0.1), 'cycling cells')
                ) +
                labs(x = 'G1/S score', y = 'G2/M score', colour = 'Phase') +
                theme_test()
        )
        
        htmps <- slapply(
            hdata,
            function(li) {
                
                htmp_pair <- slapply(
                    c('g1s', 'g2m'),
                    function(x) {
                        ggplot(li$data[gene %in% get(x)], aes(x = cell_name, y = gene, fill = exp_level)) +
                            geom_raster() +
                            scale_fill_gradientn(
                                colours = brewer.pal(9, 'Reds'),
                                limits = c(1, 5),
                                breaks = 1:5,
                                labels = c('1' = '\u2264 1', '2' = '2', '3' = '3', '4' = '4', '5' = '\u2265 5'),
                                oob = squish
                            ) +
                            scale_x_discrete(expand = c(0, 0)) +
                            scale_y_discrete(expand = c(0, 0)) +
                            geom_vline(xintercept = li$vline[-1] + 0.5) +
                            theme(
                                axis.text = element_blank(),
                                axis.ticks = element_blank(),
                                axis.ticks.length = unit(0, 'pt'),
                                axis.title = element_blank(),
                                panel.border = element_rect(fill = NA, size = 1),
                                legend.justification = c(0, 0.2),
                                plot.margin = unit(c(0, 5.5, 0, 0), 'pt')
                            ) +
                            labs(fill = 'Relative\nexpression\nlevel')
                    }
                )
                
                annots <- slapply(
                    c('g1s', 'g2m'),
                    function(sig) {
                        out <- tableGrob(
                            data.table(a = c(rev(tail(levels(li$data$gene)[levels(li$data$gene) %in% get(sig)], ngene)), '')),
                            theme = ttheme_minimal(
                                rowhead = list(padding = unit(c(0, 0), "pt"), fg_params = list(cex = 0)),
                                colhead = list(padding = unit(c(0, 0), "pt"), fg_params = list(cex = 0)),
                                core = list(padding = unit(c(5, 5), "pt"), fg_params = list(hjust = 1, x = 0.9)),
                                base_size = 9
                            )
                        )
                        out <- gtable_add_grob(
                            out,
                            grobs = segmentsGrob(
                                x0 = unit(0, "npc"),
                                y0 = unit(1, "npc"),
                                x1 = unit(1, "npc"),
                                y1 = unit(1, "npc"),
                                gp = gpar(lwd = 2.0)
                            ),
                            t = 2, b = 2, l = 1, r = 2
                        )
                        out <- gtable_add_grob(
                            out,
                            grobs = segmentsGrob(
                                x0 = unit(0, "npc"),
                                y0 = unit(0, "npc"),
                                x1 = unit(1, "npc"),
                                y1 = unit(0, "npc"),
                                gp = gpar(lwd = 2.0)
                            ),
                            t = ngene + 2, b = ngene + 2, l = 1, r = 2
                        )
                        out$widths <- unit(c(0, 1), 'npc')
                        out$heights[length(out$heights)] <- unit(1, 'npc') - sum(out$heights[-length(out$heights)])
                        out
                    }
                )
                
                phase_bar <- ggplot(li$data, aes(x = cell_name, y = 0, fill = phase)) +
                    geom_raster() +
                    scale_fill_manual(values = setNames(c('grey', brewer.pal(3, 'Dark2')), c('Not cycling', 'G1/S', 'Intermediate', 'G2/M'))) +
                    scale_x_discrete(expand = c(0, 0)) +
                    scale_y_continuous(expand = c(0, 0)) +
                    theme(
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.ticks.length = unit(0, 'pt'),
                        axis.title = element_blank(),
                        legend.justification = c(0, 1),
                        plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
                    )
                
                return(
                    plot_grid(
                        plot_grid(
                            blank_plot(),
                            blank_plot(),
                            phase_bar + theme(legend.position = 'none'),
                            ggplot() + theme(panel.background = element_blank()) + labs(y = 'G2/M genes'),
                            annots$g2m,
                            htmp_pair$g2m + theme(legend.position = 'none'),
                            ggplot() + theme(panel.background = element_blank()) + labs(y = 'G1/S genes'),
                            annots$g1s,
                            htmp_pair$g1s + theme(legend.position = 'none'),
                            blank_plot(),
                            blank_plot(),
                            ggplot() +
                                theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5, size = 11)) +
                                labs(title = 'Cells'),
                            nrow = 4,
                            ncol = 3,
                            rel_widths = c(0.5, 1, 8.5),
                            rel_heights = c(0.04, 0.45, 0.45, 0.06)
                        ),
                        plot_grid(get_legend(phase_bar), get_legend(htmp_pair$g1s), nrow = 2, ncol = 1),
                        nrow = 1,
                        ncol = 2,
                        rel_widths = c(0.84, 0.16)
                    )
                )
                
            }
        )
        
        return(list(scatter = scatter, heatmap = htmps))
        
    }

}
