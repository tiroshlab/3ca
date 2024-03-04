library(data.table)
library(ggplot2)
library(magrittr)
library(plyr)
library(cowplot)
library(RColorBrewer)
library(scales)
library(randomcoloR)
library(gridExtra)
library(ggpubr)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'))
cts <- c('Malignant', 'B_cell', 'Endothelial', 'Epithelial', 'Fibroblast', 'Macrophage', 'T_cell')





pdf('temp.pdf') # To stop ggplotGrob() opening graphics device, which leads to X11 error on server

for(r in transpose(as.list(unique(paths_table[, .(study, cancer_type)])))) tryNULL({
    
    cat(r, '\n')
    
    if(!('data_dist.rds' %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])))) next
    
    data_dist <- readRDS(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/data_dist.rds'))
    nullcond <- sapply(data_dist, is.null)
    if(all(nullcond)) next
    
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
        
        # Pie chart of MPs for each cell type across all samples in the study:
        pies <- slapply(names(data_dist[[i]]$pie_data), function(ct) {
            pd <- data_dist[[i]]$pie_data[[ct]]
            pie <- ggplot(pd, aes(x = '', y = prop, fill = a_ct)) +
                geom_bar(width = 1, stat = 'identity') +
                coord_polar('y', start = 0, direction = -1) +
                scale_fill_manual(name = 'Meta-program', values = distinctColorPalette(pd[, length(unique(a_ct))])) +
                geom_text(data = pd[prop >= 0.05], aes(x = 1.3, y = 1 - cumprop - prop/2, label = percent(prop, accuracy = 0.1)), size = 3) +
                theme_minimal() +
                theme(
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.border = element_blank(),
                    panel.grid = element_blank(),
                    axis.ticks = element_blank(),
                    legend.text = element_text(size = 8),
                    legend.title = element_text(size = 10),
                    legend.key.size = unit(12, 'pt')
                )
            return(pie)
        })
        
        # Bar charts per cell type and sample showing proportion of cells of that type in that sample assigned to each meta-program:
        bars <- slapply(names(data_dist[[i]]$bar_data), function(ct) {
            slapply(names(data_dist[[i]]$bar_data[[ct]]), function(smpl) {
                bd <- data_dist[[i]]$bar_data[[ct]][[smpl]]
                ct_lab <- mapvalues(ct, cts, c('malignant cells', 'B cells', 'endothelial cells', 'epithelial cells', 'fibroblasts', 'macrophages',
                    'T cells'), warn_missing = FALSE)
                bar <- ggplot() +
                    geom_col(data = bd, aes(x = a_smpl, y = 100*prop), width = 0.7, colour = 'black', fill = 'lemonchiffon', linewidth = 0.3) +
                    coord_flip() +
                    scale_y_continuous(
                        expand = c(0, 0),
                        breaks = seq(0, 100, 5), # Next add 10% of max value to scale, plus extra to allow space for axis text:
                        limits = c(0, bd[, {ub <- 110*max(prop); if((ub %% 5)/ub < 1/25) 5*(ub %/% 5) + ub/25 else ub}])
                    ) +
                    scale_x_discrete(expand = c(0, 0.5)) + # 2nd expand arg counts from the centre of the bar, so needs to be 0.5
                    geom_vline(aes(xintercept = 2:nrow(bd) - 0.5), linetype = 'dotted', colour = 'grey80', linewidth = 0.3) +
                    theme_pubclean() +
                    theme(
                        axis.text = element_text(size = 8),
                        axis.title.x = element_text(size = 10, margin = margin(t = 7)),
                        axis.title.y = element_text(size = 10, margin = margin(r = 7)),
                        axis.ticks.x = element_line(linewidth = 0.3),
                        axis.ticks.y = element_blank(),
                        panel.grid.major.x = element_line(colour = 'grey92', linewidth = 0.3),
                        panel.grid.major.y = element_blank(),
                        panel.border = element_rect(fill = NA, colour = 'black', linewidth = 0.3),
                        plot.margin = unit(c(5.5, 0, 5.5, 5.5), 'pt'),
                        plot.title = element_text(size = 11, face = 'bold'),
                        plot.title.position = 'plot'
                    ) +
                    labs(x = 'Meta-program', y = paste0('% of ', ct_lab, '\nin tumour'), title = paste('Sample', smpl))
                return(bar)
            })
        })
        
        # Expression heatmaps for each cell type and sample:
        heatmaps <- slapply(names(data_dist[[i]]$heatmap_data), function(ct) {
            slapply(names(data_dist[[i]]$heatmap_data[[ct]]), function(smpl) {
                htmp <- ggplot(data = data_dist[[i]]$heatmap_data[[ct]][[smpl]], aes(x = cell_name, y = gene_num, fill = value)) +
                    geom_raster() +
                    facet_grid(rows = vars(a_smpl), scales = 'free') +
                    scale_fill_gradientn(limits = c(-4, 4), breaks = c(-4, 0, 4), colours = rev(brewer.pal(11, 'RdBu')), oob = squish) +
                    scale_y_continuous(expand = c(0, 0), sec.axis = dup_axis()) +
                    scale_x_discrete(expand = c(0, 0)) +
                    theme(
                        axis.text = element_blank(),
                        axis.title.x = element_text(size = 10),
                        axis.title.y.left = element_blank(),
                        axis.title.y.right = element_text(size = 10, margin = margin(l = 5)),
                        axis.ticks = element_blank(),
                        axis.ticks.length = unit(0, 'pt'),
                        strip.background = element_blank(),
                        strip.text = element_blank(),
                        panel.spacing = unit(0, 'pt'),
                        panel.border = element_rect(fill = NA, linewidth = 0.3),
                        legend.justification = c(0.5, 0),
                        legend.text = element_text(size = 8),
                        legend.title = element_text(size = 9, margin = margin(b = 5)),
                        plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
                    ) +
                    guides(fill = guide_colourbar(barwidth = unit(10, 'pt'), barheight = unit(50, 'pt'), ticks.colour = 'black',
                        frame.colour = 'black')) +
                    labs(x = 'Cells', y = 'Genes', fill = 'Expression\n(log2 ratio)')
                return(htmp)
            })
        })
        
        grobs <- slapply(names(bars), function(ct) {
            slapply(names(bars[[ct]]), function(smpl) {
                bar_grob <- ggplotGrob(bars[[ct]][[smpl]]); htmp_grob <- ggplotGrob(heatmaps[[ct]][[smpl]])
                n_facets <- data_dist[[i]]$heatmap_data[[ct]][[smpl]][, length(unique(a_smpl))] # Facets contribute 2*n_facets - 1 entries in heights
                htmp_grob$heights[c(3, 2*n_facets + c(7, 10))] <- unit(c(1.5, 0.6, 0.9), 'cm')
                bar_grob$heights[c(3, 8, 9, 12)] <- unit(c(1.5, 0.5, 0.8, 0.2), 'cm')
                return(list(bar = bar_grob, htmp = htmp_grob))
            })
        })
        
        return(list(pies = pies, bars = bars, heatmaps = heatmaps, grobs = grobs, title_tail = title_tail))
        
    })
    
    rmdlines <- c(
        '---\n',
        'title: ""\n',
        'header-includes:\n',
        ' \\renewcommand{\\familydefault}{\\sfdefault}\n', # Indent is important but must use spaces, NOT tabs!
        ' \\pagenumbering{gobble}\n',
        'geometry: margin=1cm\n',
        'output: pdf_document\n',
        'papersize: a4\n',
        '---\n',
        '\n',
        '```{r setup, include = FALSE}\n',
        'knitr::opts_chunk$set(echo = FALSE, warning = FALSE, error = FALSE, dev = "cairo_pdf")\n',
        '```\n'
    )
    for(i in 1:length(out)) {
        if(i > 1) rmdlines <- c(rmdlines, '\n\\newpage\n')
        rmdlines <- c(rmdlines, '\n# Meta-program distribution in ', out[[i]]$title_tail, '\n\n')
        for(j in 1:length(out[[i]]$pies)) {
            if(j > 1) rmdlines <- c(rmdlines, '\n\\newpage\n')
            ct <- names(out[[i]]$pies)[j]
            saveRDS(out[[i]]$pies[[ct]], paste0('dist_plots_', i, '_pies_', ct, '.rds'))
            ct_lab_title <- mapvalues(ct, cts, c('Malignant cells', 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts',
                'Macrophages', 'T cells'), warn_missing = FALSE)
            rmdlines <- c(
                rmdlines,
                '## ', ct_lab_title, '\n',
                '\n',
                '```{r}\n',
                'pie <- readRDS("dist_plots_', i, '_pies_', ct, '.rds")\n',
                '```\n',
                '\n',
                '### **A**\n',
                '\n',
                '```{r fig.align = "center", fig.height = ', 90/25.4, ', fig.width = ', 180/25.4, '}\n',
                'pie\n',
                '```\n',
                '### **B**\n',
                '\n'
            )
            for(smpl in names(out[[i]]$bars[[ct]])) {
                saveRDS(out[[i]]$grobs[[ct]][[smpl]], paste0('dist_grobs_', i, '_', ct, '_', smpl, '.rds'))
                rmdlines <- c(
                    rmdlines,
                    '```{r}\n',
                    'grobs <- readRDS("dist_grobs_', i, '_', ct, '_', smpl, '.rds")\n',
                    '```\n',
                    '\n',
                    '```{r fig.align = "center", fig.height = ', 90/25.4, ', fig.width = ', 200/25.4, '}\n',
                    'plot_grid(plotlist = grobs, nrow = 1, ncol = 2, rel_widths = c(90, 90))\n',
                    '```\n',
                    '\n'
                )
            }
            ct_lab <- mapvalues(ct, cts, c('malignant cells', 'B cells', 'endothelial cells', 'epithelial cells', 'fibroblasts', 'macrophages',
                'T cells'), warn_missing = FALSE)
            rmdlines <- c(
                rmdlines,
                '**A.** Pie chart showing the meta-program composition of ', ct_lab, ' in ', out[[i]]$title_tail, '. **B.** ',
                'Bar plots and heatmaps showing the meta-program composition of ', ct_lab, ' of each sample in ', out[[i]]$title_tail,
                '. For each sample, only those meta-programs to which at least 3% of ', ct_lab, ' were assigned are shown. ',
                'Heatmaps (right panels) show relative expression levels of meta-program genes in a representative sample of those ', ct_lab,
                ' which were assigned to at least one of these meta-programs.\n',
                '\n'
            )
        }
    }
    out_con <- file('temp.Rmd')
    writeLines(rmdlines, con = out_con, sep = '')
    close(out_con)
    rmarkdown::render('temp.Rmd', output_file = paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/Meta-programs.pdf'))
    file.remove('temp.Rmd')
    for(i in 1:length(out)) {
        for(ct in names(out[[i]]$pies)) {
            file.remove(paste0('dist_plots_', i, '_pies_', ct, '.rds'))
            for(smpl in names(out[[i]]$bars[[ct]])) {
                file.remove(paste0('dist_grobs_', i, '_', ct, '_', smpl, '.rds'))
            }
        }
    }
    
})

dev.off()
