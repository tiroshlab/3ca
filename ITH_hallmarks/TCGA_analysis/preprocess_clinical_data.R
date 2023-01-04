library(data.table)
library(magrittr)
library(matkot)

# source('functions.R')

cancer_types <- c('ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC',
    'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM')

for(ct in cancer_types) {
    
    cat(ct, '\b...')
    
    if(ct %in% dir('~/TCGA_data') && 'All_CDEs.txt' %in% dir(paste0('~/TCGA_data/', ct))) {
        clinical_data <- fread(paste0('~/TCGA_data/', ct, '/All_CDEs.txt'), showProgress = FALSE)
    } else {
        
        download.file(
            paste0(
                'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/',
                ct,
                '/20160128/gdac.broadinstitute.org_',
                ct,
                '.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz'
            ),
            destfile = 'tmp.tar.gz',
            quiet = TRUE
        )
        
        file_names <- untar('tmp.tar.gz', list = TRUE)
        untar('tmp.tar.gz', files = file_names[endsWith(file_names, 'All_CDEs.txt')], exdir = 'tmp')
        file.remove('tmp.tar.gz')
        
        # Read in the data:
        clinical_data <- fread(paste0('tmp/', file_names[endsWith(file_names, 'All_CDEs.txt')]), showProgress = FALSE)
        
        # Remove the created directory:
        unlink('tmp', recursive = TRUE)
        
    }
    
    clinical_data <- transpose(clinical_data, make.names = 'bcr_patient_barcode', keep.names = 'bcr_patient_barcode')
    clinical_data <- clinical_data[, -'patient_id']
    setnames(clinical_data, 'bcr_patient_barcode', 'patient_id')
    clinical_data[, c('patient_id', 'cancer_type') := .(toupper(gsub('-', '\\.', patient_id)), ct)]
    clinical_data <- clinical_data[order(patient_id)]
    setcolorder(clinical_data, c('patient_id', 'cancer_type'))
    
    fwrite(clinical_data, paste0('~/TCGA_data/', ct, '/Clinical.csv'))
    
    cat('Done!\n')
    
}
