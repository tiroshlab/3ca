# This script extends the RNA subtype data already collected in preprocess_tcga_data.R.

library(data.table)
library(magrittr)
library(stringr)
library(plyr)
library(TCGAbiolinks)
library(matkot)

clin_test_data <- fread('../data/clin_test_data.csv')
subtypes_merged <- as.data.table(TabSubtypesCol_merged)[, -'color'][, samples := gsub('-', '\\.', samples)]
subtypes_merged[subtype != 'Not_Applicable', cancer_type := str_split_fixed(subtype, '\\.', 2)[, 1]]
subtypes_merged[subtype != 'Not_Applicable', subtype := str_split_fixed(subtype, '\\.', 2)[, 2]]
subtypes_merged[subtype %in% c('Nor_Applicable', 'Notassigned', '-'), subtype := NA]
setkey(subtypes_merged, samples)





# BLCA

meta <- fread('~/TCGA_data/BLCA/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'BLCA'] # Empty

as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes) # Empty

subtypes_query <- as.data.table(TCGAquery_subtype('BLCA'))[, patient := gsub('-', '\\.', patient)]
# This contains a "Histologic subtype" column, which classifies the tumours as "Papillary" or "Non-Papillary". I'm not sure this is useful. The
# "mRNA cluster" column seems to contain the same info as is already in my metadata file:
setkey(subtypes_query, patient)
meta[, all(subtypes_query[patient_id, `mRNA cluster`] == subtype)]

fwrite(
    meta[,
        .(
            patient_id = patient_id,
            expression_subtype = subtype,
            histologic_subtype = subtypes_query[patient_id, mapvalues(`Histologic subtype`, 'ND', NA)]
        )
    ] %>% unique,
    '~/TCGA_data/BLCA/subtypes.csv'
)





# BRCA

meta <- fread('~/TCGA_data/BRCA/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]
as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes) # All NAs
subtypes_query <- as.data.table(TCGAquery_subtype('BRCA'))[, patient := gsub('-', '\\.', patient)]
setkey(subtypes_query, patient)

# subtypes_merged gives no new information:
meta[!is.na(subtype), all(subtypes_merged[patient_id, subtype] == subtype)]
# I also checked the patient IDs where subtype is NA, and they don't have data in subtypes_merged.

meta[, subtype_pam50 := subtypes_query[patient_id, mapvalues(BRCA_Subtype_PAM50, c('LumA', 'LumB', 'Basal', 'Her2', 'Normal'),
    c('Luminal A', 'Luminal B', 'Basal-like', 'HER2-enriched', 'Normal-like'))]]

mmc2 <- as.data.table(read_xlsx('~/TCGA_data/BRCA/mmc2.xlsx', skip = 2))[, Case.ID := gsub('-', '\\.', Case.ID)]
setkey(mmc2, Case.ID)

meta[, c('mmc2_pam50', 'mmc2_pathology') := mmc2[patient_id, .(
    mapvalues(PAM50, c('LumA', 'LumB', 'Basal', 'Her2', 'Normal'), c('Luminal A', 'Luminal B', 'Basal-like', 'HER2-enriched', 'Normal-like')),
    `Final Pathology`
)]]
# PAM50 designations in the mmc2 table match those from subtypes_query:
meta[!is.na(subtype_pam50) & !is.na(mmc2_pam50), sum(subtype_pam50 == mmc2_pam50)/.N]
# This is not the case for the PAM50 designations I previously had:
meta[!is.na(subtype) & !is.na(mmc2_pam50), sum(subtype == mmc2_pam50)/.N]

# So mmc2 and subtypes_query seem more reliable, and therefore will take precedence over the "old" data.

fwrite(
    meta[,
        .(
            patient_id = patient_id,
            subtype_pam50 = ifelse(is.na(subtype_pam50), ifelse(is.na(mmc2_pam50), ifelse(is.na(subtype), NA, subtype), mmc2_pam50), subtype_pam50),
            subtype_pathologic = mmc2_pathology
        )
    ] %>% unique,
    '~/TCGA_data/BRCA/subtypes.csv'
)





# COAD

meta <- fread('~/TCGA_data/COAD/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'COAD'] # Empty

subtypes_mol <- as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes)[!is.na(samples), -c('color', 'patients')]
subtypes_mol[, samples := gsub('-', '\\.', samples)]
setkey(subtypes_mol, samples)

subtypes_query <- as.data.table(TCGAquery_subtype('COAD'))[, patient := gsub('-', '\\.', patient)]
setkey(subtypes_query, patient)

meta[, c('subtype_mol', 'subtype_query', 'subtype_msi', 'subtype_histologic') := c(
    .(subtypes_mol[patient_id, subtype]),
    subtypes_query[patient_id, .(expression_subtype, MSI_status, histological_type)]
)]

# All the columns from the subtypes_query table are mostly NA, so let's ignore them.

fwrite(unique(meta[, .(patient_id = patient_id, subtype = gsub('GI\\.', '', subtype_mol))]), '~/TCGA_data/COAD/subtypes.csv')





# ESCA

meta <- fread('~/TCGA_data/ESCA/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'ESCA'] # Empty

subtypes_mol <- as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes)[!is.na(samples), -c('color', 'patients')]
subtypes_mol[, samples := gsub('-', '\\.', samples)]
setkey(subtypes_mol, samples)

# All my samples are here, but I don't think I have the numbers to justify dividing any further:
all(meta$patient_id %in% subtypes_mol$samples)
subtypes_mol[, .N, by = subtype]

meta[, subtype_mol := subtypes_mol[patient_id, gsub('GI\\.', '', unique(subtype))], by = patient_id]

fwrite(
    unique(meta[, .(patient_id = patient_id, disease = subtype, subtype = ifelse(subtype == 'ESCC', NA, subtype_mol))]),
    '~/TCGA_data/ESCA/subtypes.csv'
)





# GBM

meta <- fread('~/TCGA_data/GBM/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'GBM', .N, by = subtype] # These are not what I expected

subtypes_mol <- as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes)[!is.na(samples), -c('color', 'patients')]
subtypes_mol[, unique(subtype)] # Also these

subtypes_query <- as.data.table(TCGAquery_subtype('GBM'))[, patient := gsub('-', '\\.', patient)]
setkey(subtypes_query, patient)

meta[, subtype_query := subtypes_query[patient_id, Original.Subtype]]

meta[!is.na(subtype_query) & !is.na(subtype), all(subtype == subtype_query)]

fwrite(
    unique(meta[, .(patient_id = patient_id, subtype = ifelse(is.na(subtype), subtype_query, subtype))]),
    '~/TCGA_data/GBM/subtypes.csv'
)





# HNSC

meta <- fread('~/TCGA_data/HNSC/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

meta[!is.na(subtype), all(subtype == subtypes_merged[cancer_type == 'HNSC'][patient_id, subtype])] # Nothing new here
meta[is.na(subtype), all(is.na(subtypes_merged[cancer_type == 'HNSC'][patient_id, subtype]))] # Nothing new here

subtypes_mol <- as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes)[,
    .(samples = gsub('-', '\\.', samples), subtype = gsub('HNSC\\.', '', subtype))
]
setkey(subtypes_mol, samples)

meta[!is.na(subtype), all(subtype == subtypes_mol[patient_id, subtype])] # Nothing new here
meta[is.na(subtype), all(is.na(subtypes_mol[patient_id, subtype]))] # Nor here

subtypes_query <- as.data.table(TCGAquery_subtype('HNSC'))[, patient := gsub('-', '\\.', patient)]
setkey(subtypes_query, patient)

meta[!is.na(subtype), all(subtype == subtypes_query[patient_id, RNA])] # Nothing new here
meta[is.na(subtype), all(is.na(subtypes_query[patient_id, RNA]))] # Nor here

fwrite(unique(meta[, .(patient_id, subtype)]), '~/TCGA_data/HNSC/subtypes.csv')

# Add HPV status:
subtypes <- fread('~/TCGA_data/HNSC/subtypes.csv', na.strings = '')
clin <- fread('~/TCGA_data/HNSC/Clinical.csv', na.strings = '', key = 'patient_id')
subtypes[, hpv_status := do.call(`[`, list(clin, patient_id))$hpv_status]
subtypes[hpv_status == 'indeterminate', hpv_status := NA]
fwrite(subtypes, '~/TCGA_data/HNSC/subtypes.csv')





# KICH

meta <- fread('~/TCGA_data/KICH/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

meta[, subtype_merged := subtypes_merged[cancer_type == 'KICH'][patient_id, subtype]]

as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes) # All NAs

subtypes_query <- as.data.table(TCGAquery_subtype('KICH'))[, patient := gsub('-', '\\.', patient)]
setkey(subtypes_query, patient)

meta[, subtype_query := subtypes_query[patient_id, Eosinophilic.vs.Classic]]

meta[, all(subtype_merged == subtype_query)] # TRUE

fwrite(unique(meta[, .(patient_id = patient_id, subtype = subtype_merged)]), '~/TCGA_data/KICH/subtypes.csv')





# KIRP

meta <- fread('~/TCGA_data/KIRP/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'KIRP'] # Empty
as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes) # Just cluster numbers, though one is also annotated as CIMP

subtypes_query <- as.data.table(TCGAquery_subtype('KIRP'))[, patient := gsub('-', '\\.', patient)]
setkey(subtypes_query, patient)

meta[, c('type1', 'type2') := subtypes_query[patient_id, .(tumor_type.KIRP.path., CDE_ID.3104937)]]
meta[grep('Unclassified', type1), type1 := NA]
meta[type2 %in% c('[Not Available]', '[Unknown]'), type2 := NA]
meta[, c('type1', 'type2') := .(str_extract(type1, 'Type [12]'), str_extract(type2, 'Type [12]'))]

meta[!is.na(type1) & !is.na(type2), sum(type1 == type2)/.N]
# The two type columns mosty agree, but not everywhere. I will choose tumor_type.KIRP.path. to take precedence, since it has more non-NA values.

fwrite(
    unique(meta[, .(patient_id = patient_id, subtype = ifelse(is.na(type1), ifelse(is.na(type2), NA, type2), type1))]),
    '~/TCGA_data/KIRP/subtypes.csv'
)





# LGG

meta <- fread('~/TCGA_data/LGG/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'LGG'] # Not sure about these
as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes) # Similar to above

subtypes_query <- as.data.table(TCGAquery_subtype('LGG'))[, patient := gsub('-', '\\.', patient)]
# IDH.codel.subtype and Original.Subtype columns give IDH and 1p/19q status, which I think is the most relevant. The Histology column defines as
# astrocytoma, oligodendroglioma or oligoastrocytoma, but I'm not sure how accurate this is - I think the molecular definition is final.
setkey(subtypes_query, patient)

meta[, c('molecular_subtype', 'histologic_subtype') := subtypes_query[patient_id, .(IDH.codel.subtype, Histology)]]

fwrite(unique(meta[, .(patient_id, molecular_subtype, histologic_subtype)]), '~/TCGA_data/LGG/subtypes.csv')





# LUAD

meta <- fread('~/TCGA_data/LUAD/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'LUAD'] # Just numbers
as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes) # Also numbers
as.data.table(TCGAquery_subtype('LUAD')) # expression_subtype column matches my subtype column

fwrite(unique(meta[, .(patient_id, subtype)]), '~/TCGA_data/LUAD/subtypes.csv')





# LUSC

meta <- fread('~/TCGA_data/LUSC/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'LUSC'] # Same as I already have
as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes) # Same
as.data.table(TCGAquery_subtype('LUSC')) # Same

fwrite(unique(meta[, .(patient_id, subtype)]), '~/TCGA_data/LUSC/subtypes.csv')





# OV

meta <- fread('~/TCGA_data/OV/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'OV'] # Empty
as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes) # Empty
as.data.table(TCGAquery_subtype('OV')) # No data

fwrite(unique(meta[, .(patient_id, subtype)]), '~/TCGA_data/OV/subtypes.csv')





# PAAD

meta <- fread('~/TCGA_data/PAAD/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'PAAD'] # Empty
as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes) # Empty
as.data.table(TCGAquery_subtype('PAAD')) # This has the same Moffitt subtypes that I already have.

fwrite(unique(meta[, .(patient_id, subtype)]), '~/TCGA_data/PAAD/subtypes.csv')





# READ

meta <- fread('~/TCGA_data/READ/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'READ'] # Empty

subtypes_mol <- as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes)[!is.na(samples), -c('color', 'patients')]
subtypes_mol[, samples := gsub('-', '\\.', samples)]
setkey(subtypes_mol, samples)

subtypes_query <- as.data.table(TCGAquery_subtype('READ'))[, patient := gsub('-', '\\.', patient)]
setkey(subtypes_query, patient)

meta[, c('subtype_mol', 'subtype_query', 'subtype_msi', 'subtype_histologic') := c(
    .(subtypes_mol[patient_id, subtype]),
    subtypes_query[patient_id, .(expression_subtype, MSI_status, histological_type)]
)]

# All the columns from the subtypes_query table are mostly NA, so let's ignore them.

fwrite(unique(meta[, .(patient_id = patient_id, subtype = gsub('GI\\.', '', subtype_mol))]), '~/TCGA_data/READ/subtypes.csv')





# SARC

meta <- fread('~/TCGA_data/SARC/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'SARC'] # Empty
as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes) # Empty

subtypes_query <- as.data.table(TCGAquery_subtype('SARC'))[, patient := gsub('-', '\\.', patient)] # Histologic subtypes
setkey(subtypes_query, patient)

meta[, histologic_subtype := subtypes_query[patient_id, histology]]

fwrite(unique(meta[, .(patient_id, histologic_subtype)]), '~/TCGA_data/SARC/subtypes.csv')





# SKCM

meta <- fread('~/TCGA_data/SKCM/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes) # Empty

meta[, subtype_merged := subtypes_merged[cancer_type == 'SKCM'][patient_id, subtype]]

subtypes_query <- as.data.table(TCGAquery_subtype('SKCM'))[, patient := gsub('-', '\\.', patient)]
setkey(subtypes_query, patient)

meta[, c('subtype_mutation', 'subtype_rna') := subtypes_query[patient_id, .(MUTATIONSUBTYPES, RNASEQ.CLUSTER_CONSENHIER)]]

# subtype_rna is the same as what I already have, and subtype_merged is the same as subtype_mutation.

fwrite(
    unique(meta[, .(patient_id = patient_id, expression_subtype = subtype, mutation_subtype = subtype_merged)]),
    '~/TCGA_data/SKCM/subtypes.csv'
)





# STAD

meta <- fread('~/TCGA_data/STAD/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'STAD'] # Empty

subtypes_mol <- as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes)[!is.na(samples), -c('color', 'patients')]
subtypes_mol[, samples := gsub('-', '\\.', samples)]
setkey(subtypes_mol, samples)

subtypes_query <- as.data.table(TCGAquery_subtype('STAD'))[, patient := gsub('-', '\\.', patient)]
setkey(subtypes_query, patient)

# WHO.Class subsumes Intestinal.Type.Subclass.

# This also has a 'Hypermutated' column.

meta[, c('subtype_mol', 'subtype_query', 'subtype_lauren', 'subtype_who') := c(
    .(subtypes_mol[patient_id, subtype]),
    subtypes_query[patient_id, .(Molecular.Subtype, Lauren.Class, WHO.Class)]
)]

meta[subtype_who == 'NA', subtype_who := NA]
meta[subtype_lauren == 'NA', subtype_lauren := NA]

# subtype and subtype_query have the same levels but don't exactly agree. I'll combine them, giving precedence to subtype, which has fewer NAs.
meta[is.na(subtype), subtype := subtype_query]

# subtype_mol doesn't agree so well with subtype and it has different levels, so I'll keep it as an alternative subtype definition.

fwrite(
    meta[, .(patient_id = patient_id, subtype = subtype, subtype_alt = gsub('GI\\.', '', subtype_mol), subtype_lauren = subtype_lauren,
        subtype_who = subtype_who)] %>% unique,
    '~/TCGA_data/STAD/subtypes.csv'
)





# THCA

meta <- fread('~/TCGA_data/THCA/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

subtypes_merged[cancer_type == 'THCA'] # Just numbers
as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes) # Empty

subtypes_query <- as.data.table(TCGAquery_subtype('THCA'))[, patient := gsub('-', '\\.', patient)]
setkey(subtypes_query, patient)

meta[, histologic_subtype := subtypes_query[patient_id, str_extract(histological_type, 'Classical|Follicular|Tall Cell')]]

# braf_genotype_lab doesn't have many non-missing values. BRAF_RAF_class agrees pretty well with BRAFV600E_RAS, has more values, and seems to
# better represent a subtype (BRAFV600E_RAS possibly reflects actual mutations), so I'll use BRAF_RAF_class.

meta[, molecular_subtype := subtypes_query[patient_id, mapvalues(BRAF_RAF_class, '', NA)]]

fwrite(unique(meta[, .(patient_id, histologic_subtype, molecular_subtype)]), '~/TCGA_data/THCA/subtypes.csv')





# UCEC

meta <- fread('~/TCGA_data/UCEC/Cells.csv', key = 'sample_type', na.strings = '')[sample_type != 'normal' & patient_id %in% clin_test_data$patient_id]

meta[, subtype_merged := subtypes_merged[cancer_type == 'UCEC'][patient_id, subtype]]

subtypes_mol <- as.data.table(TCGA_MolecularSubtype(gsub('\\.', '-', meta$patient_id))$subtypes)[!is.na(samples), -c('color', 'patients')]
subtypes_mol[, samples := gsub('-', '\\.', samples)]
setkey(subtypes_mol, samples)
meta[, subtype_mol := subtypes_mol[patient_id, subtype]]

subtypes_query <- as.data.table(TCGAquery_subtype('UCEC'))[, patient := gsub('-', '\\.', patient)]
setkey(subtypes_query, patient)

meta[, c('subtype_hist', 'subtype_msi', 'subtype_int') := subtypes_query[patient_id, .(histology, msi_status_7_marker_call, IntegrativeCluster)]]

# subtype_merged has only 2 non-NA values, which agree with subtype_mol, so ignore that. subtype_int has only 2 values which are not NA or
# "Notassigned", and they match both subtype_merged and subtype_mol, so ignore that as well. subtype_msi and subtype_hist have only 20 non-NA
# values, so I don't think they're that useful. Best stick to subtype_mol. The subtype_msi column also kind of agrees with subtype mol, where it's
# not NA. Note subtype_mol includes CN_HIGH and CN_LOW, where I think CN means copy number.

fwrite(unique(meta[, .(patient_id = patient_id, subtype = gsub('UCEC\\.', '', subtype_mol))]), '~/TCGA_data/UCEC/subtypes.csv')
