library(ggplot2)
library(ggrepel)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(magrittr)
library(matkot)

map.organs <- function(studies) {
    studies <- gsub('_cell_line(s)|_high_grade|_mouse_model', '', studies)
    cancer.names <- unlist(lapply(str_split(studies, '_'), function(x) if (length(x) == 2) x[2] else x[3]))
    organs.class <- list(
        'Brain' = c('Pituitary', 'Ependymoma', 'Glioma', 'gliomas', 'Glioblastoma', 'GBM', 'H3G34R', 'BrMs'),
        'Hematologic' = c('AML', 'MM', 'Leukemia', 'Myeloproliferative', 'CTCL', 'CLL'),
        'Skin' = c('Melanoma', 'SkinSCC'),
        'Sarcoma' = c('SyS', 'Ewing'),
        'Colorectal' = c('CRC', 'Colon'),
        'Pancreas' = c('PDAC', 'IPMN'),
        'Lung' = c('lungcancer', 'SCLC', 'LUAD'),
        'Head and neck' = c('HNSCC', 'LSCC', 'OPSCC', 'NasopharyngealCarcinoma', 'HTAN'),
        'NE' = c('NET', 'Neuroblastoma'),
        'Other' = c('SCHW', 'ESCC'),
        'Kidney' = c('RCC', 'WilmsTumor'),
        'Liver' = c('HCC', 'LiverMets'),
        'Lymphoma' = c('DLBCL', 'FL'),
        'Gyne' = c('HGSOC', 'Ovarian')
    )
    for (org.name in names(organs.class)) {cancer.names[cancer.names %in% organs.class[[org.name]]] <- org.name}
    return(cancer.names)
}

ct.list <- c('Fibroblast', 'B_cell', 'T_cell', 'Macrophage', 'Endothelial', 'Epithelial', 'Malignant')

ct.res <- lapply(ct.list, function(ct) {
    corr.df <- readRDS(paste0('../data/corr_data_nonCenteredGenes/', ct, '_pearson.rds'))
    study.names <- unlist(lapply(str_split(row.names(corr.df), '_(scRNAseq|snRNAseq).'), '[', 1))
    colnames(corr.df) <- study.names
    rownames(corr.df) <- study.names
    corr.pairs = reshape2::melt(corr.df)
    remove.studies <- c('Gao_2021_Breast_and_ATC', 'Gonzalez_2022_BrMs', 'Franses_2020_PDAC_CTC', 'Miyamoto_2015_Prostate_CTC',
        'Gao_2021_Breast_and_ATC', 'Jordan_2016_Breast_CTC', 'Bhan_2018_HCC_CTC', 'Massalha_2020_LiverMets', 'Griffiths_2021_Breast',
        'Jansky_2021_Neuroblastoma', 'Mercatelli_2021_Neuroblastoma')
    corr.pairs <- corr.pairs[(!corr.pairs$Var1 %in% remove.studies) & (!corr.pairs$Var2 %in% remove.studies), ]
    corr.pairs$same.study <- corr.pairs$Var1 == corr.pairs$Var2
    corr.pairs$Type1 <- map.organs(corr.pairs$Var1) # "Type" means cancer type
    corr.pairs$Type2 <- map.organs(corr.pairs$Var2)
    corr.pairs$same.type <- corr.pairs$Type1 ==  corr.pairs$Type2
    corr.pairs.ref <- corr.pairs
    type.list <- lapply(unique(corr.pairs$Type1), function(type) {
        corr.pairs = corr.pairs.ref[corr.pairs.ref$Type1==type | corr.pairs.ref$Type2==type, ]
        N.studies = length(unique(corr.pairs$Var1[corr.pairs$Type1==type]))
        N.same.tumor = 0 
        N.same.study = sum(corr.pairs$same.study==T)
        N.diff.study_same.type = sum(corr.pairs$same.study==F & corr.pairs$same.type==T)
        N.diff.study_diff.type = sum(corr.pairs$same.study==F & corr.pairs$same.type==F)
        same.tumor = 1
        same.study = mean(corr.pairs$value[corr.pairs$same.study==T])
        diff.study_same.type = mean(corr.pairs$value[corr.pairs$same.study==F & corr.pairs$same.type==T])
        diff.study_diff.type = mean(corr.pairs$value[corr.pairs$same.study==F & corr.pairs$same.type==F])
        all.pairs = mean(corr.pairs$value)
        means_vec = c(
            'N.studies'=N.studies,
            'N.same.tumor'=N.same.tumor, 'N.same.study'=N.same.study,
            'N.diff.study_same.type'=N.diff.study_same.type,
            'N.diff.study_diff.type'=N.diff.study_diff.type,
            'same.tumor'=same.tumor, 'same.study'=same.study,
            'diff.study_same.type'=diff.study_same.type,
            'diff.study_diff.type'=diff.study_diff.type,
            'all.pairs' = all.pairs
        )
        effect_vec = c('patient.specificity' = 1 - same.study, 'cancer.type.specificity' = diff.study_same.type - diff.study_diff.type)
        ct_vec = c(means_vec, effect_vec)
        ct_vec
    })
    names(type.list) <- unique(corr.pairs$Type1)
    type.df <- as.data.frame(t(data.frame(type.list)))
    type.df <- type.df[type.df$N.studies >= 3, ]
    type.df$cell.type <- ct
    type.df$type <- row.names(type.df)
    return(type.df)
}) %>% rbindlist

ct.res.long <- reshape2::melt(ct.res[, c("patient.specificity", "cancer.type.specificity", "cell.type", "type")])
ct.mean <- ct.res.long %>% group_by(variable, cell.type) %>% summarise_at(vars(value), list(name = mean))
ct.order <- ct.res.long %>% filter(variable=='cancer.type.specificity') %>% group_by(cell.type) %>% summarise_at(vars(value), list(name = mean))
ct.order <- ct.order$cell.type[order(ct.order$name)]
ct.res.long$cell.type <- factor(ct.res.long$cell.type, levels = ct.order)

fwrite(ct.res.long, '../3ca_manuscript/data/3f_ct_specificity.csv')
