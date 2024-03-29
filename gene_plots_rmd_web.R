n = commandArgs(trailingOnly = TRUE)





library(data.table)
library(magrittr)
library(ggplot2)
library(stringr)
library(plyr)
library(RColorBrewer)
library(scales)
library(cowplot)
library(matkot)

source('functions.R')

aliases <- fread('../data/alias_table.csv')





genes_n <- gsub('.csv$', '', gsub('_', '/', dir('../data/gene_plots/gene_plots_data_all_web')))
genes_n <- genes_n[cut(1:length(genes_n), 200, labels = FALSE) == n]
for(g in genes_n) {
    
    data_g <- fread(paste0('../data/gene_plots/gene_plots_data_all_web/', gsub('/', '_', g), '.csv'))
    data_g[, cell_type := gsub('_', ' ', cell_type)]
    
    data_ct <- data_g[study == 'all', -'study']
    data_study <- data_g[study != 'all']
    
    data_ct[, facet_y := 'c']
    data_ct <- rbind(
        data_ct,
        cbind(
            cancer_type = 'All cancers',
            data_ct[, .(ave = mean(ave[!is.na(ave)]), prop_pos = mean(prop_pos[!is.na(prop_pos)]), n_cell = sum(n_cell), n_sample = sum(n_sample),
                n_sample_thresh = sum(n_sample_thresh)), by = cell_type]
        )[, facet_y := 'a'],
        use.names = TRUE
    )
    data_ct[,
        c('facet_y', 'cell_type', 'cancer_type') := .(
            factor(facet_y, levels = c('c', 'a')),
            factor(cell_type, levels = c('Malignant', 'Macrophage', 'T cell', 'NK cell', 'B cell', 'Plasma', 'Dendritic', 'Mast',
                'Fibroblast', 'Pericyte', 'Endothelial', 'Epithelial', 'Oligodendrocyte')),
            factor(cancer_type, levels = sort(unique(cancer_type), decreasing = TRUE))
        )
    ]
    
    data_study[, cell_type := factor(cell_type, levels = c('Malignant', 'Macrophage', 'T cell', 'NK cell', 'B cell', 'Plasma',
        'Dendritic', 'Mast', 'Fibroblast', 'Pericyte', 'Endothelial', 'Epithelial', 'Oligodendrocyte'))]
    
    # Cut cancer types into groups for separate pages:
    ct_groups <- list(character(0))
    i <- 0
    j <- 1
    for(ct in data_study[, sort(unique(as.character(cancer_type)))]) {
        ns <- data_study[cancer_type == ct, length(unique(study))]
        if(30 - i >= ns) {
            ct_groups[[j]] <- c(ct_groups[[j]], ct)
            i <- i + ns
        } else {
            j <- j + 1
            ct_groups[[j]] <- ct
            i <- ns
        }
    }
    
    plots <- list(
        dotplot = ggplot(data_ct, aes(x = cell_type, y = cancer_type)) +
            geom_point(aes(fill = ave, size = 100*prop_pos), shape = 21, colour = 'black') +
            facet_grid(rows = vars(facet_y), scales = 'free', space = 'free') +
            scale_x_discrete(labels = c('Malignant' = expression(bold('Malignant')), parse = TRUE)) +
            scale_y_discrete(labels = c('All cancers' = expression(bold('All cancers')), parse = TRUE)) +
            scale_fill_gradientn(
                colours = brewer.pal(9, 'YlOrRd'),
                limits = c(0, 8),
                breaks = c(0, 2, 4, 6, 8),
                labels = c('0' = '0', '2' = '2', '4' = '4', '6' = '6', '8' = '\u2265 8'),
                oob = squish
            ) +
            scale_radius(limits = c(0, 100), range = c(1, 8), breaks = c(0, 20, 40, 60, 80, 100)) +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle = 55, hjust = 1),
                axis.title.x = element_blank(),
                axis.title.y = element_text(margin = margin(r = 10)),
                panel.grid = element_line(linewidth = 0.4),
                strip.background = element_blank(),
                strip.text = element_blank()
            ) +
            guides(fill = guide_colourbar(frame.colour = 'black')) +
            labs(y = 'Cancer type', fill = 'Mean', size = '% expressing\ncells'),
        dotplot_stdy = lapply(
            ct_groups,
            function(grp) {
                pdata <- data_study[cancer_type %in% grp]
                ggplot(pdata, aes(x = cell_type, y = study)) +
                    geom_point(aes(fill = ave, size = 100*prop_pos), shape = 21, colour = 'black') +
                    facet_grid(rows = vars(cancer_type), scales = 'free', space = 'free') +
                    scale_x_discrete(labels = c('Malignant' = expression(bold('Malignant')), parse = TRUE)) +
                    scale_fill_gradientn(
                        colours = brewer.pal(9, 'YlOrRd'),
                        limits = c(0, 8),
                        breaks = c(0, 2, 4, 6, 8),
                        labels = c('0' = '0', '2' = '2', '4' = '4', '6' = '6', '8' = '\u2265 8'),
                        oob = squish
                    ) +
                    scale_radius(limits = c(0, 100), range = c(1, 8), breaks = c(0, 20, 40, 60, 80, 100)) +
                    theme_bw() +
                    theme(
                        axis.text.x = element_text(angle = 55, hjust = 1),
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(margin = margin(r = 10)),
                        panel.grid = element_line(size = 0.4),
                        strip.background = element_rect(colour = NA, fill = NA),
                        strip.text = element_text(vjust = 0),
                        legend.justification = ifelse(data_study[cancer_type %in% grp, length(unique(study)) <= 14], 'top', 'center')
                    ) +
                    guides(fill = guide_colourbar(frame.colour = 'black')) +
                    labs(y = 'Study', fill = 'Mean', size = '% expressing\ncells')
            }
        )
    )
    
    saveRDS(plots, paste0('../data/gene_plots/rds_plots_web/', gsub('/', '_', g), '.rds'))
    
    rmdlines <- c('---\n', 'title: "**', g, '**"\n')
    
    if(g %in% aliases$symbol) rmdlines <- c(rmdlines, 'subtitle: "', paste(aliases[symbol == g, symbol_alt], collapse = ', '), '"\n')
    
    rmdlines <- c(
        rmdlines,
        'header-includes:\n',
        ' \\renewcommand{\\familydefault}{\\sfdefault}\n',
        ' \\pagenumbering{gobble}\n',
        'geometry: margin=1cm\n',
        'output: pdf_document\n',
        'papersize: a4\n',
        '---\n',
        '\n',
        '```{r setup, include = FALSE}\n',
        'knitr::opts_chunk$set(echo = FALSE, warning = FALSE, error = FALSE, dev = "cairo_pdf")\n',
        '```\n',
        '\n',
        '```{r}\n',
        'plots <- readRDS("../rds_plots_web/', gsub('/', '_', g), '.rds")\n', # This path has to be relative to where the Rmd file will be
        '```\n',
        '\n',
        '## **A**\n',
        '\n',
        '```{r fig.align="center", fig.height=', 165/25.4, ', fig.width=', 180/25.4, '}\n',
        'plots$dotplot\n',
        '```\n',
        '\n',
        '\\vspace*{\\fill}\n',
        '\n',
        '**A. Expression of ', g, ' per cell type and cancer type.** Plot showing the average expression level of ', g, ' and the percentage of ',
        'cells expressing ', g, ' in each cancer type and each of the most common cell types, namely those cell types constituting at least 1% of ',
        'at least 5 individual datasets. Expression levels are defined as log2(TPM/10). Average expression levels and percentages were measured ',
        'per cell type within each dataset, then averaged across studies within each cancer type.\n'
    )
    
    for(i in 1:length(plots$dotplot_stdy)) {
        plot_title <- ifelse(i == 1, '**B**', '**B (continued)**')
        plot_height <- 30 + 7.5*data_study[cancer_type %in% ct_groups[[i]], length(unique(study))] + 2*(length(ct_groups[[i]]) - 1)
        legend_space <- ifelse(i == length(plots$dotplot_stdy), 20, 0)
        rmdlines <- c(
            rmdlines,
            '\n',
            '\\newpage\n',
            '\n',
            '## ', plot_title, '\n',
            '\n',
            '```{r fig.align="center", fig.height=', (267 - legend_space)/25.4, ', fig.width=', 190/25.4, '}\n',
            'cowplot::plot_grid(plots$dotplot_stdy[[', i, ']], nrow = 2, ncol = 1, rel_heights = c(', plot_height, ', ',
            267 - legend_space -  plot_height, '))\n',
            '```\n'
        )
    }
    
    rmdlines <- c(
        rmdlines,
        '\n',
        '\\vspace*{\\fill}\n',
        '\n',
        '**B. Expression of ', g, ' per cell type and study.** Plot showing the average expression level of ', g, ' and the percentage of cells ',
        'expressing ', g, ' in each study and each of the cell types shown in A. Expression levels are defined as log2(TPM/10).\n'
    )
    
    out_con <- file(paste0('../data/gene_plots/rmds_web/', gsub('/', '_', g), '.Rmd'))
    writeLines(rmdlines, con = out_con, sep = '')
    close(out_con)
    
}
