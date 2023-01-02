library(data.table)
library(magrittr)
library(stringr)
library(readxl)
library(plyr)
library(stringi)
library(matkot)

cancer_types <- c('ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC',
    'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM')

hgnc_complete_set <- fread('../data/hgnc_complete_set.txt')[!is.na(entrez_id)][, entrez_id := as.character(entrez_id)]
setkey(hgnc_complete_set, entrez_id)





for(ct in cancer_types) {
    
    cat(ct, '\b...')
    
    if(!(ct %in% dir('~/TCGA_data'))) dir.create(paste0('~/TCGA_data/', ct))
    
    if(paste0(ct, '_illuminahiseq_rnaseqv2_RSEM_genes.txt') %in% dir(paste0('~/TCGA_data/', ct))) {
        expmat <- fread(paste0(ct, '/', ct, '_illuminahiseq_rnaseqv2_RSEM_genes.txt'), showProgress = FALSE)
    } else {
        
        download.file(
            paste0(
                'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/',
                ct,
                '/20160128/gdac.broadinstitute.org_',
                ct,
                '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz'
            ),
            destfile = 'tmp.tar.gz',
            quiet = TRUE
        )
        
        file_names <- untar('tmp.tar.gz', list = TRUE)
        untar('tmp.tar.gz', files = file_names[endsWith(file_names, 'data.txt')], exdir = 'tmp')
        file.remove('tmp.tar.gz')
        
        # Rename file, else the file name is too long for R to handle:
        file.rename(paste0('tmp/', str_split_fixed(file_names[endsWith(file_names, 'data.txt')], '/', 2)[, 1]), 'tmp/tempdir')
        
        # Move the data file up one directory and delete tempdir:
        file.rename(
            paste0('tmp/tempdir/', str_split_fixed(file_names[endsWith(file_names, 'data.txt')], '/', 2)[, 2]),
            paste0('tmp/', str_split_fixed(file_names[endsWith(file_names, 'data.txt')], '/', 2)[, 2])
        )
        unlink('tmp/tempdir', recursive = TRUE)
        
        # Read in the data:
        expmat <- fread(paste0('tmp/', str_split_fixed(file_names[endsWith(file_names, 'data.txt')], '/', 2)[, 2]), showProgress = FALSE)
        
        # Remove the created directory:
        unlink('tmp', recursive = TRUE)
        
    }
    
    names(expmat)[1] <- 'id'
    
    # Extract scaled estimates:
    expmat <- expmat[, expmat[1] %in% c('gene_id', 'scaled_estimate'), with = FALSE][-1]
    
    # Make columns numeric:
    expmat[, names(expmat[, -'id']) := lapply(.SD, as.numeric), .SDcols = -'id']
    
    # Get up-to-date gene symbols from Entrez IDs:
    expmat[, id := hgnc_complete_set[str_split_fixed(id, '\\|', 2)[, 2], symbol]]
    expmat <- expmat[!is.na(id)]
    
    # Convert to matrix and split IDs for extraction of patient IDs and sample type codes:
    expmat <- expmat[, set_rownames(as.matrix(.SD), id), .SDcols = -'id']
    colnames(expmat) <- gsub('-', '\\.', colnames(expmat))
    
    sample_codes <- str_split_fixed(colnames(expmat), '\\.', 5)
    if(ct != 'LAML') { # For LAML, every sample is type '03', which is "primary blood derived"
        to_keep <- grep('^01|^02|^05|^06|^07|^11', sample_codes[, 4])
        expmat <- expmat[, to_keep]
        sample_codes <- sample_codes[to_keep, ]
    }
    
    # Convert to (something like) log TPM by multiplying by 1e+06 and taking log:
    expmat <- apply(expmat, 2, function(x) round(log2(x*1e+06 + 1), 4))
    
    # Construct metadata table:
    meta <- data.table(
        sample_id = colnames(expmat),
        patient_id = apply(sample_codes[, 1:3], 1, paste, collapse = '.'),
        cancer_type = ct,
        sample_type = mapvalues(
            str_split_fixed(sample_codes[, 4], '[A-Z]', 2)[, 1],
            c('01', '02', '03', '05', '06', '07', '11'),
            c('primary', 'recurrent', 'primary', 'primary_additional', 'metastatic', 'metastatic_additional', 'normal'),
            warn_missing = FALSE
        )
    )
    
    # Make sure it's in alphabetical order:
    setkey(meta, sample_id)
    expmat <- expmat[, meta$sample_id]
    
    # Put the genes in alphabetical order as well:
    expmat <- expmat[order(rownames(expmat)), ]
    
    # Write to file:
    fwrite(meta, paste0('~/TCGA_data/', ct, '/Cells.csv'))
    fwrite(as.data.frame(expmat), paste0('~/TCGA_data/', ct, '/Exp_data_TPM.csv'), row.names = TRUE)
    
    cat('Done!\n')
    
}





# Add purity data:

# First get all metadata and combine into one table:
meta <- rbindlist(lapply(cancer_types, function(ct) fread(paste0('~/TCGA_data/', ct, '/Cells.csv'))))

# Purity data from tcga_annotations file:

tcga_annotations <- fread('~/tcga_annotations.csv', showProgress = FALSE)[,
    .(id = individual_id, cancer_type = tumor_type, purity = purity)
][!is.na(purity)]

# Fix empty cancer_type entries in tcga_annotations:
setkey(meta, patient_id)
tcga_annotations[cancer_type == '', cancer_type := meta[sample_type == 'primary'][id, cancer_type]]

# The table has only patient IDs.  For patients with multiple samples, decide which sample to assign the purity value to by defining sample type:
tcga_annotations[,
    sample_type := meta[
        id,
        switch(('primary' %in% sample_type) + 1, switch((length(unique(sample_type)) > 1) + 1, unique(sample_type), NA), 'primary')
    ],
    by = id
]

tcga_annotations <- tcga_annotations[!is.na(sample_type)]
setnames(tcga_annotations, 'id', 'patient_id')
setcolorder(tcga_annotations, c('patient_id', 'cancer_type', 'sample_type', 'purity'))

# Purity data from the paper by Taylor et al. 2018 (https://doi.org/10.1016/j.ccell.2018.03.007):

purity_taylor <- as.data.table(read_xlsx('~/TCGA_data/taylor_1-s2.0-S1535610818301119-mmc2.xlsx', skip = 1))[
    Type %in% cancer_types,
    .(
        patient_id = apply(str_split_fixed(Sample, '-', 4)[, 1:3], 1, paste, collapse = '.'),
        cancer_type = Type,
        sample_type = mapvalues(
            str_split_fixed(Sample, '-', 4)[, 4],
            c('01', '02', '05', '06'), # There are no purity values for LAML
            c('primary', 'recurrent', 'primary_additional', 'metastatic'),
            warn_missing = FALSE
        ),
        purity = Purity
    )
][!is.na(purity)]

# Extra dataset for oesophageal and stomach (from TCGA oesophageal paper, https://doi.org/10.1038/nature20805):

purity_esca_stad <- as.data.table(read_xlsx('~/TCGA_data/ESCA/nature20805-s1-1.xlsx', skip = 1))[
    `Absolute extract purity` != 'NA',
    .(id = gsub('-', '\\.', as.character(barcode)), cancer_type = `Disease code`, purity = as.numeric(`Absolute extract purity`))
][!is.na(purity)]

# Decide which sample to assign the purity value to by defining sample type:
purity_esca_stad[,
    sample_type := meta[
        id,
        switch(('primary' %in% sample_type) + 1, switch((length(unique(sample_type)) > 1) + 1, unique(sample_type), NA), 'primary')
    ],
    by = id
] # All the sample types end up being 'primary'.

setnames(purity_esca_stad, 'id', 'patient_id')
setcolorder(purity_esca_stad, c('patient_id', 'cancer_type', 'sample_type', 'purity'))

# Merge these data tables:
purity_data <- Reduce(function(x, y) merge(x, y, all = TRUE), list(tcga_annotations, purity_taylor, purity_esca_stad))

rm(tcga_annotations)
rm(purity_taylor)
rm(purity_esca_stad)

# Average purity values for repeated patient_id/sample_type combinations:
purity_data <- purity_data[, .(purity = round(mean(purity), 2)), by = .(patient_id, cancer_type, sample_type)]

fwrite(purity_data, '~/TCGA_data/purity_data.csv')





# Add purity data to meta_data table:

setkey(purity_data, patient_id, sample_type)

for(ct in cancer_types) {
    meta <- fread(paste0('~/TCGA_data/', ct, '/Cells.csv'))
    meta[, purity := purity_data[.(meta$patient_id, meta$sample_type), purity]]
    fwrite(meta, paste0('~/TCGA_data/', ct, '/Cells.csv'))
}





# Subtypes data:

subtypes_data <- rbindlist(
    list(
        # Robertson et al. 2017, TCGA bladder cancer paper: https://doi.org/10.1016/j.cell.2017.09.007
        BLCA = as.data.table(read_xlsx('~/TCGA_data/BLCA/1-s2.0-S0092867417310565-mmc1.xlsx', sheet = 'Master table'))[,
            .(patient_id = gsub('-', '\\.', `Case ID`), cancer_type = 'BLCA', subtype = `mRNA cluster`)
        ],
        # TCGA 2012 breast cancer paper: https://doi.org/10.1038/nature11412
        BRCA = as.data.table(read_xls('~/TCGA_data/BRCA/Supplementary Tables 1-4.xls', skip = 1))[,
            .(patient_id = gsub('-', '\\.', `Complete TCGA ID`), cancer_type = 'BRCA', subtype = `PAM50 mRNA`)
        ][subtype != 'NA'],
        # TCGA 2017 oesophageal paper, https://doi.org/10.1038/nature20805
        ESCA = as.data.table(read_xlsx('~/TCGA_data/ESCA/nature20805-s1-1.xlsx', skip = 1))[
            `Disease code` == 'ESCA' & `Histological Type` != 'NA',
            .(patient_id = gsub('-', '\\.', barcode), cancer_type = 'ESCA', subtype = `Histological Type`)
        ],
        # Brennan et al. 2013, TCGA GBM paper: https://doi.org/10.1016/j.cell.2013.09.034
        GBM = as.data.table(read_xlsx('~/TCGA_data/GBM/1-s2.0-S0092867413012087-mmc7.xlsx', skip = 2))[
            !is.na(`Path Dx`),
            .(patient_id = gsub('-', '\\.', `Case ID`), cancer_type = 'GBM', subtype = `Expression\r\nSubclass`)
        ],
        # TCGA 2015 HNSCC paper: https://doi.org/10.1038/nature14129
        HNSC = as.data.table(read_xlsx('~/TCGA_data/HNSC/7.2.xlsx'))[,
            .(patient_id = as.character(Barcode), cancer_type = 'HNSC', subtype = as.character(RNA))
        ],
        # TCGA 2017 liver cancer paper: https://doi.org/10.1016/j.cell.2017.05.046
        LIHC = as.data.table(read_xlsx('~/TCGA_data/LIHC/1-s2.0-S0092867417306396-mmc3.xlsx', skip = 1))[
            `iCluster clusters` != 'NA',
            .(patient_id = gsub('-', '\\.', short.id), cancer_type = 'LIHC', subtype = gsub(':', ' ', `iCluster clusters`))
        ],
        # TCGA 2014 lung adenocarcinoma paper: https://doi.org/10.1038/nature13385
        LUAD = as.data.table(read_xlsx('~/TCGA_data/LUAD/nature13385-s2.xlsx', sheet = 7, skip = 4))[,
            .(
                patient_id = gsub('-', '\\.', `Tumor ID`),
                cancer_type = 'LUAD',
                subtype = mapvalues(expression_subtype, c('TRU', 'prox.-inflam', 'prox.-prolif.'), c('Bronchioid', 'Squamoid', 'Magnoid'))
            )
        ],
        # TCGA 2012 lung squamous paper: https://doi.org/10.1038/nature11404
        LUSC = as.data.table(read_xls('~/TCGA_data/LUSC/data.file.S7.5.clinical.and.genomic.data.table.xls', skip = 2))[,
            .(patient_id = gsub('LUSC', 'TCGA', gsub('-', '\\.', `Tumor ID`)), cancer_type = 'LUSC', subtype = get(names(.SD)[length(.SD)]))
        ],
        # Verhaak et al. 2012: https://doi.org/10.1172/JCI65833
        OV = as.data.table(read_xls('~/TCGA_data/OV/JCI65833sd1.xls', skip = 1, sheet = 1))[
            DATASET == 'TCGA-discovery',
            .(patient_id = gsub('-', '\\.', ID), cancer_type = 'OV', subtype = SUBTYPE)
        ],
        # TCGA 2017 PDAC paper: https://doi.org/10.1016/j.ccell.2017.07.007
        PAAD = as.data.table(read_xlsx('~/TCGA_data/PAAD/mmc2.xlsx', skip = 1))[,
            .(
                patient_id = sapply(`Tumor Sample ID`, function(tsid) paste(str_split_fixed(tsid, '-', 4)[, 1:3], collapse = '.')),
                cancer_type = 'PAAD',
                subtype = mapvalues(`mRNA Moffitt clusters (All 150 Samples) 1basal  2classical`, c(1, 2), c('Basal', 'Classical'))
            )
        ],
        # TCGA 2015 prostate cancer paper: https://doi.org/10.1016/j.cell.2015.10.025
        PRAD = as.data.table(read_xlsx('~/TCGA_data/PRAD/1-s2.0-S0092867415013392-mmc2.xlsx'))[,
            .(patient_id = gsub('-', '\\.', PATIENT_ID), cancer_type = 'PRAD', subtype = mRNA_cluster)
        ],
        # TCGA 2015 melanoma paper: https://doi.org/10.1016/j.cell.2015.05.044
        SKCM = as.data.table(
            read_xlsx('~/TCGA_data/SKCM/1-s2.0-S0092867415006340-mmc2.xlsx', sheet = 'Supplemental Table S1D', skip = 1)
        )[
            `RNASEQ-CLUSTER_CONSENHIER` != '-', # All *patient* IDs in this subset are unique.
            .(
                patient_id = apply(str_split_fixed(Name, '-', 4)[, 1:3], 1, paste, collapse = '.'),
                cancer_type = 'SKCM',
                subtype = `RNASEQ-CLUSTER_CONSENHIER`
            )
        ],
        # TCGA 2017 oesophageal paper, https://doi.org/10.1038/nature20805
        STAD = as.data.table(read_xlsx('~/TCGA_data/ESCA/nature20805-s1-1.xlsx', skip = 1))[
            `Disease code` == 'STAD',
            .(patient_id = gsub('-', '\\.', barcode), cancer_type = 'STAD', subtype = `Gastric classification`)
        ]
    )
)

setkey(subtypes_data, patient_id)

fwrite(subtypes_data, '~/TCGA_data/subtypes_data.csv')

for(ct in cancer_types) {
    meta <- fread(paste0('~/TCGA_data/', ct, '/Cells.csv'))
    meta[, subtype := subtypes_data[meta$patient_id, subtype]]
    fwrite(meta, paste0('~/TCGA_data/', ct, '/Cells.csv'))
}
