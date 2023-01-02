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

study_folder_map <- fread('../data/study_folder_map.csv')





for(r in transpose(as.list(unique(study_folder_map[study != 'Kinker et al. 2020', .(study, cancer_type)])))) {
    
    cat(r, '\n')
    
    if('data_ct_comp.RDS' %in% dir(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1]))) {
        
        plot_data <- readRDS(paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/data_ct_comp.RDS'))
        
        nullcond <- sapply(plot_data, function(x) ifelse(is.null(x), TRUE, all(sapply(x[names(x) != 'path'], is.null))))
        if(all(nullcond)) next
        
        out <- lapply(
            which(!nullcond),
            function(i) {
                
                # To use in plot titles:
                if(sum(!nullcond) > 1) {
                    title_tail <- str_split_fixed(
                        plot_data[[i]]$path$cells,
                        paste0(study_folder_map[study == r[1] & cancer_type == r[2], directory], '/'),
                        2
                    )[, 2]
                    if(grepl('/', title_tail)) {
                        title_tail <- str_split(title_tail, '/')[[1]][1]
                        if(sum(study_folder_map$study == r[1]) == 1) {
                            title_tail <- paste0(r[1], ', ', title_tail)
                        } else {
                            title_tail <- paste0(r[1], ' (', r[2], '), ', title_tail)
                        }
                    } else {
                        if(sum(study_folder_map$study == r[1]) == 1) {
                            title_tail <- paste0(r[1], ', Group ', i)
                        } else {
                            title_tail <- paste0(r[1], ' (', r[2], '), Group ', i)
                        }
                    }
                } else {
                    if(sum(study_folder_map$study == r[1]) == 1) {
                        title_tail <- r[1]
                    } else {
                        title_tail <- paste0(r[1], ' (', r[2], ')')
                    }
                }
                
                if(!is.null(plot_data[[i]]$data)) {
                    
                    pdata <- plot_data[[i]]$data
                    
                    marker_cond <- 'marker_plot_data' %in% names(pdata) &&
                        pdata$marker_plot_data[, length(unique(gene)) > 1 & length(unique(cell_type)) > 1]
                    if('marker_plot_data' %in% names(pdata) && pdata$marker_plot_data[, length(unique(gene)) < 10]) {
                        leg_arr <- 'horizontal'
                    } else {leg_arr <- 'vertical'}
                    
                    if('randomcoloR' %in% .packages(TRUE)) {
                        set.seed(9728)
                        plot_colours <- pdata$pie_data[, setNames(distinctColorPalette(.N), cell_type)]
                        plot_out <- do.call(
                            ct_comp_plot,
                            args = c(switch(marker_cond + 1, pdata['pie_data'], pdata), list(colours = plot_colours, legends_arrange = leg_arr))
                        )
                    } else {
                        plot_out <- do.call(
                            ct_comp_plot,
                            args = c(switch(marker_cond + 1, pdata['pie_data'], pdata), list(legends_arrange = leg_arr))
                        )
                    }
                    
                    if(marker_cond) {
                        plot_title <- paste0('Cell type composition and marker gene expression in ', title_tail)
                    } else {plot_title <- paste0('Cell type composition in ', title_tail)}
                    
                    return(list(plot_out = plot_out, plot_title = plot_title, leg_arr = leg_arr))
                    
                }
                
            }
        )
        
        # A4 page is 210x297mm. Allow 20mm for the title and 15mm each for the "A" and "B" labels.
        rmdlines <- c(
            '---\n',
            'title: ""\n',
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
            '```\n'
        )
        for(i in 1:length(out)) {
            if('marker_plot' %in% names(out[[i]]$plot_out)) {
                margin_x <- (ifelse(out[[i]]$leg_arr == 'vertical', 150, 125) - # Take 25mm off margin if legends are arranged horizontally
                        10*length(unique(plot_data[!nullcond][[i]]$data$marker_plot_data$cell_type)))/2
                n_genes <- length(unique(plot_data[!nullcond][[i]]$data$marker_plot_data$gene))
                h <- (log2(24/n_genes) + 6)*n_genes + 30
                # The following top-aligns legends in case of too few genes, else the tops of the legends get cut off.  It also depends on the length
                # of the x labels, as long labels squash the plot.  4 seems OK, but 5 is safer in case there's a longer label than "Oligodendrocyte".
                if(n_genes <= 5) out[[i]]$plot_out$marker_plot <- out[[i]]$plot_out$marker_plot + theme(legend.justification = 'top')
            }
            saveRDS(out[[i]]$plot_out, paste0('comp_plots_', i, '.rds'))
            if(i > 1) rmdlines <- c(rmdlines, '\n\\newpage\n')
            if('marker_plot' %in% names(out[[i]]$plot_out)) {
                rmdlines <- c(
                    rmdlines,
                    '\n',
                    '### ', out[[i]]$plot_title, '\n',
                    '\n',
                    '```{r}\n',
                    'plots <- readRDS("comp_plots_', i, '.rds")\n',
                    '```\n',
                    '\n',
                    '### **A**\n',
                    '\n',
                    '```{r fig.align = "center", fig.height = ', 117/25.4, ', fig.width = ', 210/25.4, ', out.width = "80%", out.height = "80%"}\n',
                    'plots$pie\n',
                    '```\n',
                    '\n',
                    '### **B**\n',
                    '\n',
                    '```{r fig.align = "center", fig.height = ', h/25.4, ', fig.width = ', 210/25.4, ', out.width = "80%", out.height = "80%"}\n',
                    'plots$marker_plot + theme(plot.margin = unit(c(0, ', margin_x, ', 0, ', margin_x, '), "mm"))\n',
                    '```\n',
                    '\n',
                    '\\vspace*{\\fill}\n', # Aligns text with bottom of page
                    '\n',
                    '**A.** Pie chart showing proportions of cell types after quality control. Cell types are ordered by abundance. **B.** Plot ',
                    'showing the average log2(TPM/10) expression levels of marker genes in each cell type and, for each gene, the percentage of ',
                    'cells of each type in which that gene was detected. Points are not shown for combinations where the percentage of expressing ',
                    'cells is below 50%. Only non-malignant cell types with at least 10 cells are shown, and these are ordered as in A.\n'
                )
            } else {
                rmdlines <- c(
                    rmdlines,
                    '\n',
                    '### ', out[[i]]$plot_title, '\n',
                    '\n',
                    '```{r}\n',
                    'plots <- readRDS("comp_plots_', i, '.rds")\n',
                    '```\n',
                    '\n',
                    '```{r fig.align = "center", fig.height = ', 117/25.4, ', fig.width = ', 210/25.4, ', out.width = "80%", out.height = "80%"}\n',
                    'plots$pie\n',
                    '```\n',
                    '\n',
                    'Pie chart showing proportions of cell types after quality control. Cell types are ordered by abundance.\n'
                )
            }
        }
        out_con <- file('temp.Rmd')
        writeLines(rmdlines, con = out_con, sep = '')
        close(out_con)
        rmarkdown::render('temp.Rmd', output_file = paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1], '/Cell types.pdf'))
        file.remove('temp.Rmd')
        for(i in 1:length(out)) file.remove(paste0('comp_plots_', i, '.rds'))
        
    }
    
}
