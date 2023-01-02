library(data.table)
library(magrittr)
library(limma)
library(stringr)
library(plyr)
library(RCurl)
library(matkot)

cancer_types <- c('ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC',
    'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM')

hgnc_complete_set <- fread('../data/hgnc_complete_set.txt', key = 'ensembl_gene_id')[
    !(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1]),
    .(symbol = symbol, alias = alias_symbol, entrez_id = entrez_id)
][, entrez_id := as.character(entrez_id)]





for(ct in cancer_types) {
    
    cat(ct, '\n')
    
    # The CNV file URLs use one of three suffixes, so we need to choose whichever one exists:
    data_url <- sapply(
        c('-TP', '-TB', '-TM'),
        function(suff) paste0( 'http://gdac.broadinstitute.org/runs/analyses__2016_01_28/data/', ct, suff, '/20160128/gdac.broadinstitute.org_', ct,
            suff, '.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz'),
        USE.NAMES = FALSE
    )
    data_url <- data_url[sapply(data_url, url.exists)]
    
    download.file(data_url, destfile = 'tmp.tar.gz', quiet = TRUE)
    
    file_names <- untar('tmp.tar.gz', list = TRUE)
    untar('tmp.tar.gz', files = file_names[endsWith(file_names, 'all_thresholded.by_genes.txt')], exdir = 'tmp')
    file.remove('tmp.tar.gz')
    
    # Read in the data:
    cnv_data <- fread(paste0('tmp/', file_names[endsWith(file_names, 'all_thresholded.by_genes.txt')]))
    
    # Remove the created directory:
    unlink('tmp', recursive = TRUE)
    
    cnv_data[,
        c('Gene Symbol', 'Locus ID', 'Cytoband') := .(
            do.call(
                function(gl) {
                    gl[grep('C[0-9]+ORF[0-9]+', gl)] <- gsub('ORF', 'orf', gl[grep('C[0-9]+ORF[0-9]+', gl)])
                    gl_mapped <- alias2SymbolTable(gl)
                    gl_log <- gl %in% hgnc_complete_set$symbol
                    gl_mapped_log <- gl_mapped %in% hgnc_complete_set$symbol
                    out <- gl
                    out[!gl_log & gl_mapped_log] <- gl_mapped[!gl_log & gl_mapped_log]
                    out[!gl_log & !gl_mapped_log] <- sapply(
                        gl_mapped[!gl_log & !gl_mapped_log],
                        function(g) hgnc_complete_set[
                            grepl(paste0('^', g, '\\||\\|', g, '\\||\\|', g, '$'), alias),
                            switch((.N == 1) + 1, NA, symbol)
                        ]
                    )
                    return(out)
                },
                args = list(gl = cnv_data$`Gene Symbol`)
            ),
            NULL,
            NULL
        )
    ]
    
    cnv_data <- cnv_data[!is.na(`Gene Symbol`)]
    cnv_data <- cnv_data[, set_rownames(as.matrix(.SD), `Gene Symbol`), .SDcols = -'Gene Symbol']
    colnames(cnv_data) <- gsub('-', '\\.', colnames(cnv_data))
    
    fwrite(as.data.frame(cnv_data), paste0('~/TCGA_data/', ct, '/CNV_GISTIC_data.csv'), row.names = TRUE)
    
}





# Mutations annotations files were downloaded from GDC Data Portal, filtering by MAF, open access, Masked Somatic Mutation, MuTect2.

# These files have full sample IDs, but since they don't match the sample IDs in the expression data, I'll extract just the patient IDs.

setkey(hgnc_complete_set, entrez_id)

for(ct in cancer_types) {
    cat(ct, '\n')
    cont <- dir(paste0('~/TCGA_data/mut_data/', ct))
    mut_data <- fread(paste0('~/TCGA_data/mut_data/', ct, '/', cont[grep(ct, cont)]))[,
        c('Entrez_Gene_Id', 'patient_id', 'sample_type') := .(
            as.character(Entrez_Gene_Id),
            apply(str_split_fixed(Tumor_Sample_Barcode, '-', 4)[, 1:3], 1, paste, collapse = '.'),
            mapvalues(
                gsub('[A-Z]', '', str_split_fixed(Tumor_Sample_Barcode, '-', 5)[, 4]),
                c('01', '02', '03', '05', '06', '07', '11'),
                c('primary', 'recurrent', 'primary', 'primary_additional', 'metastatic', 'metastatic_additional', 'normal'),
                warn_missing = FALSE
            )
        )
    ][,
        gene := hgnc_complete_set[Entrez_Gene_Id, symbol]
    ][
        !is.na(gene),
        .(gene, patient_id, sample_type, Chromosome, Start_Position, End_Position, Strand, Variant_Classification, Variant_Type,
            Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, Mutation_Status)
    ]
    fwrite(mut_data, paste0('~/TCGA_data/', ct, '/Mutation_data.csv'))
}
