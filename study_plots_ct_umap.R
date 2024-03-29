# bsub -q new-medium -R "select[model=Intel_Skylake] rusage[mem=32000]" -oo log/study_plots_ct_umap.o -eo log/study_plots_ct_umap.e Rscript study_plots_ct_umap.R

library(data.table)
library(ggplot2)
library(magrittr)
library(Matrix)
library(stringr)
library(plyr)
library(cowplot)
library(RColorBrewer)
library(scales)
library(matkot)

try(library(randomcoloR), silent = TRUE)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', encoding = 'UTF-8', key = c('study', 'cancer_type'))





for(r in transpose(as.list(unique(paths_table[, .(study, cancer_type)])))) {
    
    cat(r, '\n')
    
    if(all(c('data_ct_umap.RDS', 'UMAP.pdf') %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])) == c(TRUE, FALSE))) {
        
        plot_data <- readRDS(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/data_ct_umap.RDS'))
        
        nullcond <- sapply(plot_data, function(x) ifelse(is.null(x), TRUE, is.null(x$data)))
        if(all(nullcond)) next
        
        if('randomcoloR' %in% .packages(TRUE)) {
            set.seed(9728)
            plot_colours <- lapply(
                which(!nullcond),
                function(i) {
                    pdata <- plot_data[[i]]$data
                    slapply(c('cell_type', 'sample')[c('cell_type', 'sample') %in% names(pdata)], function(vn) unique(pdata[[vn]]))
                }
            )
            plot_colours <- slapply(
                c('cell_type', 'sample'),
                function(vn) {
                    vn_unique <- unique(unlist(lapply(plot_colours, function(x) if(vn %in% names(x)) x[[vn]])))
                    if(!is.null(vn_unique)) setNames(distinctColorPalette(length(vn_unique)), vn_unique)
                }
            )
        }
        
        out <- lapply(which(!nullcond), function(i) {
            
            pdata <- plot_data[[i]]$data
            if(!any(c('cell_type', 'sample') %in% names(pdata))) return(NULL)
            
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
            
            if('cell_type' %in% names(pdata)) {
                if('Malignant' %in% pdata$cell_type) {
                    pdata[, cell_type := factor(cell_type, levels = c('Malignant', sort(unique(cell_type[cell_type != 'Malignant']))))]
                } else {
                    pdata[, cell_type := factor(cell_type, levels = sort(unique(cell_type)))]
                }
            }
            if('sample' %in% names(pdata)) {
                if(!any(is.na(as.numeric(unique(pdata$sample))))) {
                    pdata[, sample := factor(as.character(sample), levels = sort(as.numeric(unique(sample))))]
                } else pdata[, sample := factor(as.character(sample), levels = sort(as.character(unique(sample))))]
            }
            if('randomcoloR' %in% .packages(TRUE)) c(list(
                title_tail = title_tail,
                plots = slapply(
                    c('cell_type', 'sample')[c('cell_type', 'sample') %in% names(pdata)],
                    function(vn) ct_umap_plot(pdata, vn, colours = plot_colours[[vn]][levels(pdata[[vn]])],
                        legend_title = gsub('_', ' ', str_to_title(vn)))
                )
            )) else c(list(
                title_tail = title_tail,
                plots = slapply(
                    c('cell_type', 'sample')[c('cell_type', 'sample') %in% names(pdata)],
                    function(vn) ct_umap_plot(pdata, vn, legend_title = gsub('_', ' ', str_to_title(vn)))
                )
            ))
            
        })
        
        if(all(sapply(out, is.null))) next
        
        # A4 page is 210x297mm.  Here, I'm allowing 20mm for the title and 15mm each for the "A" and "B" labels.
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
            saveRDS(out[[i]]$plots, paste0('umap_plots_', i, '.rds'))
            if(i > 1) rmdlines <- c(rmdlines, '\n\\newpage\n')
            if(length(out[[i]]$plots) == 2) {
                rmdlines <- c(
                    rmdlines,
                    '\n',
                    '### UMAP in ', out[[i]]$title_tail, '\n',
                    '\n',
                    '```{r}\n',
                    'plots <- readRDS("umap_plots_', i, '.rds")\n',
                    '```\n',
                    '\n',
                    '### **A**\n',
                    '\n',
                    '```{r fig.align = "center", fig.height = ', 130/25.4, ', fig.width = ', 180/25.4, ', out.width = "80%", out.height = "80%"}\n',
                    'plots$cell_type\n',
                    '```\n',
                    '\n',
                    '### **B**\n',
                    '\n',
                    '```{r fig.align = "center", fig.height = ', 130/25.4, ', fig.width = ', 180/25.4, ', out.width = "80%", out.height = "80%"}\n',
                    'plots$sample\n',
                    '```\n',
                    '\n',
                    '\\vspace*{\\fill}\n', # Aligns text with bottom of page
                    '\n',
                    '**A.** UMAP plot of single cells in ', out[[i]]$title_tail, ', after quality control. Cells are colored by cell type. **B.** ',
                    'UMAP plot as in A, with cells colored by sample.\n'
                )
            } else {
                rmdlines <- c(
                    rmdlines,
                    '\n',
                    '### UMAP in ', out[[i]]$title_tail, '\n',
                    '\n',
                    '```{r}\n',
                    'plots <- readRDS("umap_plots_', i, '.rds")\n',
                    '```\n',
                    '\n',
                    '```{r fig.align = "center", fig.height = ', 130/25.4, ', fig.width = ', 180/25.4, ', out.width = "80%", out.height = "80%"}\n',
                    'plots[[1]]\n',
                    '```\n',
                    '\n',
                    'UMAP plot of single cells in ', out[[i]]$title_tail, ', after quality control. Cells are colored by ',
                    gsub('_', ' ', names(out[[i]]$plots)), '.\n'
                )
            }
        }
        out_con <- file('temp.Rmd')
        writeLines(rmdlines, con = out_con, sep = '')
        close(out_con)
        rmarkdown::render('temp.Rmd', output_file = paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/UMAP.pdf'))
        file.remove('temp.Rmd')
        for(i in 1:length(out)) file.remove(paste0('umap_plots_', i, '.rds'))
        
    }
    
}
