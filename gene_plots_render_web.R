n = commandArgs(trailingOnly = TRUE)

genes_n <- gsub('.csv$', '', dir('../data/gene_plots/gene_plots_data_all_web'))
genes_n <- genes_n[cut(1:length(genes_n), 500, labels = FALSE) == n]
for(g in genes_n) {
    if(paste0(g, '.Rmd') %in% dir('../data/gene_plots/rmds_web') & paste0(g, '.rds') %in% dir('../data/gene_plots/rds_plots_web')) {
        rmarkdown::render(
            paste0('/home/labs/tirosh/tyler/pan_cancer/data/gene_plots/rmds_web/', g, '.Rmd'),
            intermediates_dir = paste0('/home/labs/tirosh/tyler/pan_cancer/data/gene_plots/pdfs_web/tmp/', g),
            output_dir = '/home/labs/tirosh/tyler/pan_cancer/data/gene_plots/pdfs_web'
        )
    }
}
