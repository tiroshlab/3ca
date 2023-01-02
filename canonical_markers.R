library(data.table)

canonical_markers <- list(
    Adipocyte = c('ADIPOQ', 'FABP4', 'LEP'),
    Astrocyte = c('GFAP', 'GJA1'),
    B_cell = c('MS4A1', 'CD79A', 'CD79B'),
    Basophil = c('CCR3', 'IL3RA'),
    Dendritic = c('CCR7', 'CD86', 'CLEC10A'),
    Endothelial = c('VWF', 'CDH5', 'CLEC14A', 'CLDN5', 'ADGRL4'),
    Epithelial = c('EPCAM', 'KRT19', 'KRT7', 'MUC1'),
    Erythrocyte = c('HBA1', 'HBA2', 'HBB'),
    Fibroblast = c('COL1A1', 'COL1A2', 'COL3A1', 'LUM'),
    Immune = 'PTPRC',
    Keratinocyte = c('FLG', 'IVL', 'LORICRIN'),
    Lymphovascular = c('CCL21', 'TFF3'),
    Macrophage = c('AIF1', 'CD14', 'CD163'),
    Mast = c('CPA3', 'MS4A2', 'TPSB2'),
    Melanocyte = c('DCT', 'MLANA', 'PMEL', 'TYR', 'TYRP1'),
    Monocyte = c('CCR2', 'CSF1R'),
    Myocyte = c('ACTA1', 'DES', 'MYL1', 'TTN', 'NEB'),
    Neuron = c('ENO2', 'RBFOX3'),
    Neutrophil = c('AZU1', 'CTSG', 'ELANE', 'MPO'),
    NK_cell = c('NKG7', 'KLRD1', 'GZMB', 'KLRF1'),
    Oligodendrocyte = c('TF', 'PLP1', 'CLDN11', 'MAG', 'MBP', 'MOG'),
    OPC = c('OLIG1', 'OLIG2', 'NEU4'),
    Pericyte = c('RGS5', 'ESAM', 'MEF2C'),
    Plasma = c('JCHAIN', 'MZB1', 'TNFRSF17'),
    T_cell = c('CD2', 'CD3D', 'CD3E')
)

canonical_markers <- rbindlist(lapply(names(canonical_markers), function(ct) data.table(cell_type = ct, gene = canonical_markers[[ct]])))

fwrite(canonical_markers, '../data/canonical_markers.csv')
