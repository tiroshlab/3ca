library(data.table)
library(ggplot2)
library(magrittr)
library(Matrix)
library(stringr)
library(plyr)
library(cowplot)
library(RColorBrewer)
library(scales)
library(gtable)
library(grid)
library(gridExtra)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', encoding = 'UTF-8', key = c('study', 'cancer_type'))

ct_to_title <- function(ct) {
    if(grepl('ic$|al$|ar$|oid$|crine$|ant$', ct) | ct %in% c('Mast', 'Plasma', 'Immune', 'Stem', 'Schwann', 'Langerhans', 'Tuft', 'Stellate')) {
        return(paste0(ct, ' cells'))
    } else {
        return(paste0(gsub('_', ' ', ct), 's'))
    }
}





for(r in transpose(as.list(unique(paths_table[, .(study, cancer_type)])))) {
    
    cat(r, '\n')
    
    if(all(c('data_cc.RDS', 'Cell cycle.pdf') %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])) == c(TRUE, FALSE))) {
        
        plot_data <- readRDS(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/data_cc.RDS'))
        nullcond <- sapply(plot_data, function(x) ifelse(is.null(x), TRUE, all(sapply(x[names(x) != 'path'], is.null))))
        
        if(!all(nullcond)) {
            
            suff <- paste(gsub(' et al. ', '', r[1]), gsub('/| ', '-', r[2]), sep = '_')
            
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
                        title_tail <- paste0(r[1], ', ', r[2], ' - ', paths_table[as.list(r)][i, if(group_name == '') paste('Group', group) else
                            group_name])
                    } else {
                        title_tail <- paste0(r[1], ', ', r[2])
                    }
                }
                
                plot_out <- do.call(cc_plot, args = plot_data[[i]][c('ccdata', 'hdata', 'g1s', 'g2m')])
                plot_title <- paste0('Cell cycle activity in ', title_tail)
                return(list(plot_out = plot_out, plot_title = plot_title))
                
            })
            
            if(!all(sapply(out, function(x) is.null(x$plot_out)))) {
                
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
                    if(is.null(out[[i]]$plot_out)) next
                    saveRDS(out[[i]]$plot_out, paste0('comp_plots_', i, '_', suff, '.rds'))
                    if(i > 1) rmdlines <- c(rmdlines, '\n\\newpage\n')
                    rmdlines <- c(
                        rmdlines,
                        '\n',
                        '```{r}\n',
                        'plots <- readRDS("comp_plots_', i, '_', suff, '.rds")\n',
                        '```\n'
                    )
                    if(names(out[[i]]$plot_out$scatter)[1] == 'data') {
                        # In this case there's only one scatter plot ('data' is the 1st element of a ggplot object), and there will be a heatmap.
                        rmdlines <- c(
                            rmdlines,
                            '\n',
                            '### ', out[[i]]$plot_title, '\n',
                            '\n',
                            '### **A**\n',
                            '\n',
                            '```{r fig.align="center", fig.height=', 110/25.4, ', fig.width=', 170/25.4, ', out.width="70%", out.height="70%"}\n',
                            'plots$scatter\n',
                            '```\n',
                            '\n',
                            '### **B**\n',
                            '\n',
                            '```{r fig.align="center", fig.height=', 175/25.4, ', fig.width=', 210/25.4, ', out.width="80%", out.height="80%"}\n',
                            'plots$heatmap\n',
                            '```\n',
                            '\n',
                            '\\vspace*{\\fill}\n', # Aligns text with bottom of page
                            '\n',
                            '**A.** Scatter plot showing the scores for G1/S (x axis) and G2/M (y axis) cell cycle phases for individual cells ',
                            '(points). The color indicates the phase assigned to each cell, and the percentage refers to all cells which are ',
                            'assigned as either G1/S, G2/M or intermediate. **B.** Heatmap of relative expression levels in individual cells ',
                            '(columns) of the genes (rows) best distinguishing G1/S and G2/M phases. Cells are ordered by phase (indicated by ',
                            'the color bar, top) and, within each phase, by the maximum of the G1/S and G2/M scores.\n'
                        )
                    } else { # In this case every element is a separate scatter plot, corresponding to a cell type, possibly with no heatmap.
                        cts <- names(out[[i]]$plot_out$scatter)
                        if('Malignant' %in% cts) {cts <- c('Malignant', sort(cts[cts != 'Malignant']))} else {cts <- sort(cts)}
                        for(j in 1:length(out[[i]]$plot_out$scatter)) {
                            if(j > 1) rmdlines <- c(rmdlines, '\n\\newpage\n')
                            ct <- cts[j]
                            plot_title <- paste0(out[[i]]$plot_title, ' - ', ct_to_title(ct))
                            if(ct %in% names(out[[i]]$plot_out$heatmap)) {
                                rmdlines <- c(
                                    rmdlines,
                                    '\n',
                                    '### ', plot_title, '\n',
                                    '\n',
                                    '### **A**\n',
                                    '\n',
                                    '```{r fig.align="center", fig.height=',
                                    110/25.4,
                                    ', fig.width=',
                                    170/25.4,
                                    ', out.width="70%", out.height="70%"}\n',
                                    'plots$scatter[["', ct, '"]]\n',
                                    '```\n',
                                    '\n',
                                    '### **B**\n',
                                    '\n',
                                    '```{r fig.align="center", fig.height=',
                                    175/25.4,
                                    ', fig.width=',
                                    210/25.4,
                                    ', out.width="80%", out.height="80%"}\n',
                                    'plots$heatmap[["', ct, '"]]\n',
                                    '```\n',
                                    '\n',
                                    '\\vspace*{\\fill}\n', # Aligns text with bottom of page
                                    '\n',
                                    '**A.** Scatter plot showing the scores for G1/S (x axis) and G2/M (y axis) cell cycle phases for individual ',
                                    ct_to_title(ct),
                                    ' (points). The color indicates the phase assigned to each cell, and the percentage refers to all cells ',
                                    'which are assigned as either G1/S, G2/M or intermediate. **B.** Heatmap of relative expression levels in ',
                                    'individual ',
                                    ct_to_title(ct),
                                    ' (columns) of the genes (rows) best distinguishing G1/S and G2/M phases. Cells are ordered by phase ',
                                    '(indicated by the color bar, top) and, within each phase, by the maximum of the G1/S and G2/M scores.\n'
                                )
                            } else {
                                rmdlines <- c(
                                    rmdlines,
                                    '\n',
                                    '### ', plot_title, '\n',
                                    '\n',
                                    '```{r fig.align = "center", fig.height = ', 120/25.4, ', fig.width = ', 180/25.4,
                                    ', out.width = "70%", out.height = "70%"}\n',
                                    'plots$scatter[["', ct, '"]]\n',
                                    '```\n',
                                    '\n',
                                    'Scatter plot showing the scores for G1/S (x axis) and G2/M (y axis) cell cycle phases for individual ',
                                    ct_to_title(ct),
                                    ' (points). The color indicates the phase assigned to each cell, and the percentage refers to all cells which ',
                                    'are assigned as either G1/S, G2/M or intermediate.\n'
                                )
                            }
                        }
                    }
                }
                out_con <- file(paste0('temp_', suff, '.Rmd'))
                writeLines(rmdlines, con = out_con, sep = '')
                close(out_con)
                rmarkdown::render(
                    paste0('temp_', suff, '.Rmd'),
                    output_file = paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/Cell cycle.pdf')
                )
                file.remove(paste0('temp_', suff, '.Rmd'))
                for(i in 1:length(out)) file.remove(paste0('comp_plots_', i, '_', suff, '.rds'))
                
            }
            
        }
        
    }
    
}
