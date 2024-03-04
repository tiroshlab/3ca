# The command line arguments should supply study name and cancer type (in that order).
r = commandArgs(trailingOnly = TRUE)





library(data.table)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(randomcoloR)
library(scales)
library(cowplot)
library(ggtext)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'))





paths <- apply(paths_table[as.list(r), .(cells, genes, expmat)], 1, as.list, simplify = FALSE)

if('data_cna.rds' %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1]))) {
    
    plot_data <- readRDS(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/data_cna.rds'))
    nullcond <- sapply(plot_data, function(x) ifelse(is.null(x), TRUE, all(sapply(x[names(x) != 'path'], is.null))))
    
    out <- lapply(which(!nullcond), function(i) {
        
        # To use in plot titles:
        if(sum(unique(paths_table[, .(study, cancer_type)])$study == r[1]) == 1) {
            if(paths_table[as.list(r), .N] > 1) {
                title_tail <- paste0(r[1], ' - ', paths_table[as.list(r)][i, if(group_name == '') paste('Group', group) else group_name])
            } else {
                title_tail <- r[1]
            }
        } else {
            if(paths_table[as.list(r), .N] > 1) {
                title_tail <- paste0(r[1], ', ', r[2], ' - ', paths_table[as.list(r)][i, if(group_name == '') paste('Group', group) else group_name])
            } else {
                title_tail <- paste0(r[1], ', ', r[2])
            }
        }
        
        cells <- suppressWarnings(fread(paths[[i]]$cells, na.strings = ''))
        cells$cell_name <- as.character(cells$cell_name)
        setkey(cells, cell_name)
        
        pdata <- plot_data[[i]]$cna_data[, .(gene, cell_name, cna, sample, chr)]
        pdata[, cell_type := do.call(`[`, list(cells, cell_name))$cell_type]
        pdata <- pdata[!is.na(cell_type)]
        chr_n <- pdata[, .(n = length(unique(gene))), keyby = chr][n > 20]
        pdata <- pdata[chr %in% chr_n$chr]
        pdata[, malignant := ifelse(cell_type == 'Malignant', 'm', 'nm')]
        
        # Sample cells so that we have no more than 2000, of which between 60 and 85% are malignant:
        grp_n <- pdata[, .(n = length(unique(cell_name))), keyby = .(malignant, sample)]
        grp_tot <- grp_n[, min(max(sum(malignant == 'm')/.N, 0.6), 0.85)]
        grp_tot <- 2000*c(m = grp_tot, nm = 1 - grp_tot)
        grp_tot_true <- grp_n[, .(N = sum(n)), keyby = malignant]$N
        if(any(grp_tot > grp_tot_true)) {wm <- which.min(grp_tot_true - grp_tot); grp_tot <- grp_tot*grp_tot_true[wm]/grp_tot[wm]}
        grp_n[, n_to_sample := { # Returns number of cells to take from each sample
            tot <- grp_tot[unique(malignant)] # Total allowed for this group
            baseline <- pmin(n, min(20, floor(tot/length(unique(sample))))) # Number of cells to keep (samples with too few cells are untouched)
            N <- n - baseline
            to_lose <- sum(n) - tot # Number of cells from each sample that we can lose
            N <- N - floor(to_lose*N/sum(N)) # Take number of cells from each sample proportional to the total number in that sample
            to_lose <- sum(N + baseline) - tot # Any left over due to the floor() function
            sub_places <- order(order(-N)) %in% 1:to_lose # Take these off evenly from the samples having the most cells
            N[sub_places] <- N[sub_places] - 1
            N + baseline
        }, by = malignant]
        set.seed(1012)
        pdata <- pdata[, .SD[cell_name %in% sample( # The actual sampling
            unique(cell_name),
            do.call(`[`, list(grp_n, .(unique(malignant), unique(sample))))$n_to_sample,
            replace = FALSE
        )], by = .(malignant, sample)]
        pdata[, gene := factor(gene, levels = unique(gene))]
        pdata[, cell_name := factor(cell_name, levels = unique(cell_name[order(cell_type != 'Malignant', sample)]))]
        pdata[, cell_num := as.numeric(cell_name)]
        
        # Manually choose chromosome labels:
        chr_lab <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '14', '16', '17', '18', '19', '20', '22', 'X')
        chr_lab <- unique(pdata[, .(gene, chr)])[,
            {inds <- chr_n[, floor(cumsum(n) - n/2)[chr %in% chr_lab]]; setNames(gsub('^0', '', chr[inds]), as.character(gene)[inds])}
        ]
        # To choose them automatically using some threshold for minimum number of genes (100 in this case):
        # chr_lab <- unique(pdata[, .(gene, chr)])[,
        #     {inds <- chr_n[, floor(cumsum(n) - n/2)[n >= 100]]; setNames(gsub('^0', '', chr[inds]), as.character(gene)[inds])}
        # ]
        
        cells_lab <- pdata[, .(n = length(unique(cell_num))), keyby = malignant][, floor(c(n[1]/2, n[1] + n[2]/2))]
        
        set.seed(1358)
        bar_smpl <- ggplot(pdata) +
            geom_raster(aes(x = 'a', y = cell_num, fill = sample)) +
            scale_fill_manual(name = 'Sample', values = pdata[, setNames(distinctColorPalette(length(unique(sample))), unique(sample))]) +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme(
                axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, 'pt'),
                legend.position = 'left',
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12),
                legend.justification = c(1, 0.5),
                panel.border = element_rect(fill = NA, colour = 'black', linewidth = 0.5),
                plot.margin = unit(c(5.5, 2, 5.5, 5.5), 'pt')
            ) +
            guides(fill = guide_legend(keywidth = unit(10, 'pt'), keyheight = unit(15, 'pt')))
        
        htmp <- ggplot(pdata) +
            geom_raster(aes(x = gene, y = cell_num, fill = cna)) +
            geom_vline(aes(xintercept = xi), data = data.table(xi = chr_n[1:(.N - 1), cumsum(n) + 0.5]), linewidth = 0.3) +
            geom_hline(aes(yintercept = yi), data = data.table(yi = pdata[cell_type == 'Malignant', length(unique(cell_name)) + 0.5]),
                linewidth = 0.5) +
            scale_fill_gradientn(colours = rev(brewer.pal(11, 'RdBu')), limits = c(-1, 1), breaks = c(-1, 0, 1), oob = squish) +
            scale_x_discrete(expand = c(0, 0), breaks = names(chr_lab), labels = chr_lab) +
            scale_y_continuous(expand = c(0, 0), sec.axis = dup_axis(breaks = cells_lab, labels = c('Malignant', 'TME'))) +
            theme(
                axis.text.x = element_text(size = 12),
                axis.text.y.left = element_blank(),
                axis.text.y.right = element_text(angle = -90, hjust = 0.5, size = 14),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, 'pt'),
                axis.title.x = element_text(size = 14),
                axis.title.y.left = element_blank(),
                axis.title.y.right = element_text(size = 14),
                legend.title = element_text(size = 12, margin = margin(b = 6)),
                legend.text = element_text(size = 10),
                legend.justification = c(0.5, 0),
                panel.border = element_rect(fill = NA, colour = 'black', linewidth = 0.5),
                plot.title = element_text(size = 14, face = 'bold'),
                plot.subtitle = element_text(size = 12),
                plot.margin = unit(c(5.5, 5.5, 5.5, 1), 'pt')
            ) +
            guides(fill = guide_colourbar(barwidth = unit(15, 'pt'), barheight = unit(70, 'pt'), ticks.colour = 'black', frame.colour = 'black')) +
            labs(
                x = 'Chromosome',
                y = 'Cells',
                fill = 'Inferred CNA\n(log2 ratio)',
                title = paste('Inferred CNAs in', title_tail),
                subtitle = paste(
                    'Reference cells:',
                    paste(plot_data[[i]]$ref_cells[, .N, by = cell_type][order(-N), gsub('_', ' ', cell_type)], collapse = ', ')
                )
            )
        
        return(list(heatmap = htmp, bar = bar_smpl))
        
    })
    
    if(sum(unique(paths_table[, .(study, cancer_type)])$study == r[1]) == 1) caption <- r[1] else caption <- paste0(r[1], ' (', r[2], ')')
    if(length(out) == 1) {
        caption <- paste0('Heatmap of inferred copy number alteration (CNA) values (quantified as log2 ratio) at each chromosomal position',
            ' for a representative subset of cells in ', caption, ', with colour bar (left) showing the corresponding samples.')
    } else {
        caption <- paste0('Heatmaps of inferred copy number alteration (CNA) values (quantified as log2 ratio) at each chromosomal position',
            ' for representative subsets of cells in ', caption, ', with colour bars (left) showing the corresponding samples.')
    }
    caption <- ggplot() +
        theme(panel.background = element_rect(fill = NA), plot.caption.position = 'plot',
            plot.caption = element_textbox_simple(hjust = 0, size = 12)) +
        labs(caption = caption)
    
    leg <- lapply(out, function(x) get_legend(x$bar))
    leg_width <- max(sapply(leg, function(x) as.numeric(x$widths[3])))
    w <- c(leg_width + 0.4, 0.7, 26.7)
    
    # Open the PNG before using plot_grid(), because this function calls ggdraw(), which opens a graphic device, leading to X11 error on the server.
    png(
        paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/CNAs.png'),
        width = sum(w), height = length(out)*18 + (length(out) - 1)*4 + 4, units = 'cm', res = 300
    )
    bar_grob <- lapply(out, function(x) ggplotGrob(x$bar + theme(legend.position = 'none')))
    htmp_grob <- lapply(out, function(x) ggplotGrob(x$heatmap))
    for(j in 1:length(bar_grob)) bar_grob[[j]]$heights[c(1, 7, 12)] <- unit(c(1.7, 14.9, 1.4), 'cm')
    for(j in 1:length(htmp_grob)) htmp_grob[[j]]$heights <- unit(c(0.2, 0, 0.8, 0.7, 0, 0, 14.9, 0.5, 0.7, 0, 0, 0.2), 'cm')
    plist <- lapply(1:length(out), function(j) {
        if(j == 1) {
            plot_grid(leg[[j]], bar_grob[[j]], htmp_grob[[j]], nrow = 1, ncol = 3, rel_widths = w)
        } else plot_grid(
            ggplot() + theme_void(),
            plot_grid(leg[[j]], bar_grob[[j]], htmp_grob[[j]], nrow = 1, ncol = 3, rel_widths = w),
            nrow = 2, ncol = 1, rel_heights = c(4, 18)
        )
    })
    if(length(plist) == 1) {
        print(plot_grid(plist[[1]], caption, nrow = 2, ncol = 1, rel_heights = c(18, 4)))
    } else {
        print(plot_grid(plotlist = c(plist, list(caption)), nrow = length(plist) + 1, ncol = 1, rel_heights = c(18, rep(22, length(plist) - 1), 4)))
    }
    dev.off()
    
}
