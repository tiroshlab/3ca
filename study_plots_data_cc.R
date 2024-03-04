# The command line arguments should supply study name and cancer type (in that order).
r = commandArgs(trailingOnly = TRUE)





library(data.table)
library(ggplot2)
library(magrittr)
library(Matrix)
library(stringr)
library(plyr)
library(matkot)

source('functions.R')

paths_table <- fread('../data/paths_table.csv', key = c('study', 'cancer_type'))

hgnc_complete_set <- fread('../data/hgnc_complete_set.txt', key = 'symbol')
hgnc_complete_set <- hgnc_complete_set[!(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1])]
alias_table <- make_alias_table(hgnc_complete_set)

bins_nm <- fread('../data/bins_nm.csv')

# Study-specific bins for malignant cells and whole samples:
rdir <- paste0('../data/study_plots/', gsub('/', '-', r[2]), '/', r[1])
gene_ave <- fread(paste0(rdir, '/gene_ave.csv'))
gene_ave_list <- list(gene_ave[cell_type == 'Malignant' & n_cell >= 10, .(group, symbol, sample, cell_type, ave)])
if(unique(gene_ave[cell_type == 'Malignant', .(group, sample, n_cell)])[, sum(n_cell) >= 10]) gene_ave_list <- c(gene_ave_list, list(
    gene_ave[cell_type == 'Malignant', .(group = 'all', sample = 'all', cell_type = 'Malignant', ave = sum(ave*n_cell)/sum(n_cell)), by = symbol]
))
gene_ave_list <- c(gene_ave_list, list(gene_ave[,
    if(unique(.SD[, .(cell_type, n_cell)])[, sum(n_cell) >= 10]) .SD[, .(cell_type = 'all', ave = sum(ave*n_cell)/sum(n_cell)), by = symbol],
    by = .(group, sample)
]))
if(unique(gene_ave[, .(sample, cell_type, n_cell)])[, sum(n_cell) >= 10]) gene_ave_list <- c(gene_ave_list, list(
    gene_ave[, .(group = 'all', sample = 'all', cell_type = 'all', ave = sum(ave*n_cell)/sum(n_cell)), by = symbol]
))
gene_ave <- rbindlist(gene_ave_list, use.names = TRUE)
set.seed(7327) # Shuffle the genes before ranking so that the zero genes don't simply get ordered alphabetically:
gene_ave[, ave_rank := .SD[sample(1:.N), order(order(ave))[order(symbol)]]/.N, .SDcols = c('symbol', 'ave'), by = .(group, sample, cell_type)]
bins_m <- gene_ave[
    cell_type == 'Malignant',
    cbind(.SD[order(ave_rank)], bin = cut(1:.N, 15, labels = FALSE)),
    by = .(group, sample),
    .SDcols = -'cell_type'
]
bins_all <- gene_ave[
    cell_type == 'all',
    cbind(.SD[order(ave_rank)], bin = cut(1:.N, 15, labels = FALSE)),
    by = .(group, sample),
    .SDcols = -'cell_type'
]

genesets <- list(
    g1s = list(
        scandal = c("RRM2", "TYMS" , "UBE2T", "CDK1", "HMGB2", "MAD2L1", "PCNA", "UBE2C", "PBK", "TOP2A", "NUSAP1", "KIAA0101", "HIST1H4C",
            "MLF1IP", "GMNN", "BIRC5", "FAM64A", "RNASEH2A", "MELK", "CENPK", "PTTG1", "TK1", "TPX2", "TMEM106C", "CDCA5", "CKS1B", "CDC45", "MCM3",
            "CENPM", "AURKB", "PKMYT1", "KIF22","MCM4", "ASF1B", "GINS2", "MCM2", "NUF2", "CDKN3", "GGH", "NDC80", "FEN1", "RRM1", "PRC1" , "DUT",
            "RAD51AP1", "CKS2", "MCM7", "CCNE2", "ZWINT"),
        oligo = c('MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP',
            'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2',
            'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8'),
        neftel = c('RRM2', 'PCNA', 'KIAA0101', 'HIST1H4C', 'MLF1IP', 'GMNN', 'RNASEH2A', 'MELK', 'CENPK', 'TK1', 'TMEM106C', 'CDCA5', 'CKS1B',
            'CDC45', 'MCM3', 'CENPM', 'AURKB', 'PKMYT1', 'MCM4', 'ASF1B', 'GINS2', 'MCM2', 'FEN1', 'RRM1', 'DUT', 'RAD51AP1', 'MCM7', 'CCNE2',
            'ZWINT'),
        kinker = c('HIST1H4C', 'CLSPN', 'ATAD2', 'E2F1', 'HELLS', 'RRM2', 'HIST2H2AC', 'HIST1H1E', 'FAM111A', 'GINS2', 'CENPU', 'CDCA5', 'ASF1B',
            'CHAF1A', 'TCF19', 'HIST1H1D', 'KIAA0101', 'FEN1', 'HIST1H1B', 'PCNA', 'FBXO5', 'HIST1H1C', 'CDCA4', 'MYBL2', 'PKMYT1', 'FAM111B',
            'USP1', 'CDC6', 'EZH2', 'PSMC3IP', 'GMNN', 'HMGB2', 'ORC6', 'TYMS', 'DNAJC9', 'CDK1', 'UBE2T', 'SLBP', 'BRCA1', 'C21orf58', 'CDT1',
            'ESCO2', 'TEX30', 'ATAD5', 'CCNE1', 'RAD51AP1'),
        mps = c('PCNA', 'RRM2', 'FEN1', 'GINS2', 'TYMS', 'MCM3', 'GMNN', 'HIST1H4C', 'CLSPN', 'ATAD2', 'TK1', 'KIAA0101', 'DUT', 'HELLS', 'MCM7',
            'UBE2T', 'MCM4', 'CENPU', 'DHFR', 'ZWINT', 'ASF1B', 'MCM5', 'DNAJC9', 'RFC4', 'HMGB2', 'CDC6', 'RRM1', 'ORC6', 'CDK1', 'RAD51AP1',
            'RNASEH2A', 'CHAF1A', 'CENPK', 'CDCA5', 'SLBP', 'MCM6', 'TMEM106C', 'CENPM', 'MYBL2', 'E2F1', 'USP1', 'DNMT1', 'PKMYT1', 'MAD2L1',
            'PSMC3IP', 'CDCA4', 'RFC2', 'CDC45', 'UHRF1', 'MCM2')
    ),
    g2m = list(
        scandal = c("CCNB1", "UBE2C", "PTTG1", "CDC20", "CCNB2", "TOP2A", "FAM64A", "NUSAP1", "CDKN3", "PBK", "PLK1", "HMGB2", "TPX2", "BIRC5",
            "MAD2L1", "PRC1", "NUF2", "UBE2T", "CDK1", "CKS2", "CCNA2", "CKAP2", "KNSTRN", "RACGAP1", "CDCA3", "TROAP", "KIF2C", "AURKA", "CENPF",
            "KPNA2", "KIF20A", "ECT2", "BUB1", "CDCA8", "BUB1B", "TACC3", "NDC80", "TTK", "TUBA1C", "NCAPD2", "ARL6IP1", "KIF4A", "CKAP2L", "MZT1",
            "KIFC1", "KIF22", "TYMS", "SPAG5", "ANP32E", "KIF11", "PSRC1", "TUBB4B", "SMC4", "MXD3", "CDC25B", "OIP5", "GGH", "REEP4", "FOXM1",
            "TMPO", "GPSM2", "HMGB3", "ARHGAP11A", "RANGAP1", "H2AFZ"),
        oligo = c('HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF',
            'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP',
            'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA',
            'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA'),
        neftel = c('CCNB1', 'CDC20', 'CCNB2', 'PLK1', 'CCNA2', 'CKAP2', 'KNSTRN', 'RACGAP1', 'CDCA3', 'TROAP', 'KIF2C', 'AURKA', 'CENPF', 'KPNA2',
            'KIF20A', 'ECT2', 'BUB1', 'CDCA8', 'BUB1B', 'TACC3', 'TTK', 'TUBA1C', 'NCAPD2', 'ARL6IP1', 'KIF4A', 'CKAP2L', 'MZT1', 'KIFC1', 'SPAG5',
            'ANP32E', 'KIF11', 'PSRC1', 'TUBB4B', 'SMC4', 'MXD3', 'CDC25B', 'OIP5', 'REEP4', 'FOXM1', 'TMPO', 'GPSM2', 'HMGB3', 'ARHGAP11A',
            'RANGAP1', 'H2AFZ'),
        kinker = c('AURKA', 'CENPF', 'PLK1', 'TOP2A', 'UBE2C', 'ASPM', 'TPX2', 'CENPA', 'CKAP2', 'GTSE1', 'CCNB1', 'ARL6IP1', 'MKI67', 'CENPE',
            'CKS2', 'HMMR', 'DEPDC1', 'NUSAP1', 'PRC1', 'SGOL2', 'CCNA2', 'KPNA2', 'CDCA8', 'HMGB2', 'NUF2', 'KNSTRN', 'CDCA3', 'CEP55', 'KIF20B',
            'FAM83D', 'PIF1', 'CDC20', 'DLGAP5', 'KIF2C', 'PRR11', 'ARHGAP11A', 'KIF23', 'AURKB', 'CDK1', 'KIF14', 'FAM64A', 'CCNB2', 'PSRC1',
            'NEK2', 'CDCA2', 'BIRC5', 'TACC3', 'CKAP2L', 'HJURP', 'KIF4A', 'UBALD2', 'CDC25B', 'RACGAP1', 'SGOL1', 'ECT2', 'CCNF', 'UBE2S', 'ANLN',
            'CDKN3', 'DBF4', 'G2E3', 'PBK'),
        mps = c('TOP2A', 'UBE2C', 'HMGB2', 'NUSAP1', 'CENPF', 'CCNB1', 'TPX2', 'CKS2', 'BIRC5', 'PRC1', 'PTTG1', 'KPNA2', 'MKI67', 'CDC20', 'CDK1',
            'CCNB2', 'CDKN3', 'SMC4', 'NUF2', 'ARL6IP1', 'CKAP2', 'ASPM', 'PLK1', 'CKS1B', 'CCNA2', 'AURKA', 'MAD2L1', 'GTSE1', 'HMMR', 'UBE2T',
            'CENPE', 'CENPA', 'KIF20B', 'AURKB', 'CDCA3', 'CDCA8', 'UBE2S', 'KNSTRN', 'KIF2C', 'PBK', 'TUBA1B', 'DLGAP5', 'TACC3', 'STMN1',
            'DEPDC1', 'ECT2', 'CENPW', 'ZWINT', 'HIST1H4C', 'KIF23')
    )
)

genesets <- slapply(genesets, function(x) slapply(x, update_symbols_fast, alias_table))





paths <- apply(paths_table[as.list(r), .(cells, genes, expmat)], 1, as.list, simplify = FALSE)

out <- lapply(1:length(paths), function(paths_i) {
    
    p <- paths[[paths_i]]
    
    if(!endsWith(p$expmat, 'mtx')) return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
    
    cat('Preparing data\n')
    cells <- suppressWarnings(fread(p$cells, na.strings = ''))
    cells$cell_name <- as.character(cells$cell_name) # In case cell names are the same as row numbers
    if('sample' %in% names(cells)) cells$sample <- as.character(cells$sample)
    if(cells[, .N > length(unique(cell_name))]) cells[, cell_name := paste(cell_name, .I, sep = '_')]
    genes <- fread(p$genes, header = FALSE)$V1
    expmat <- readMM(p$expmat)
    expmat <- expmat[genes != '', ]
    genes <- genes[genes != '']
    expmat <- expmat[genes %in% names(table(genes))[table(genes) == 1], ]
    genes <- genes[genes %in% names(table(genes))[table(genes) == 1]]
    genes <- update_symbols_fast(genes, alias_table) # Update gene symbols
    rownames(expmat) <- genes
    colnames(expmat) <- cells$cell_name
    cells <- cells[col_nnz(expmat) >= 1000] # Remove low-complexity cells
    if(nrow(cells) < 30) return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
    expmat <- round(log_transform(1e+05*to_frac(expmat[, cells$cell_name])), 4) # Normalise to log TPM/10
    bins_nm <- bins_nm[symbol %in% rownames(expmat)]
    bins_m <- bins_m[symbol %in% rownames(expmat) & group %in% c(paths_i, 'all'), -'group']; setnames(bins_m, 'sample', 'group')
    bins_all <- bins_all[symbol %in% rownames(expmat) & group %in% c(paths_i, 'all'), -'group']; setnames(bins_all, 'sample', 'group')
    binslist <- list(bins_nm, bins_m, bins_all)[sapply(list(bins_nm, bins_m, bins_all), function(x) nrow(x) > 0)]
    if('cell_type' %in% names(cells)) {
        cells <- cells[!is.na(cell_type) & cell_type != 'Unassigned']
        cells <- cells[cell_type %in% cells[, .(N = .N), by = cell_type][N >= 30, cell_type]] # Only consider cell types with at least 30 cells
        expmat <- expmat[, cells$cell_name]
    }
    
    cat('Obtaining cell cycle signatures\n')
    g1s <- Reduce(union, genesets$g1s); g2m <- Reduce(union, genesets$g2m)
    int <- intersect(g1s, g2m)
    g1s <- g1s[!(g1s %in% int)]; g2m <- g2m[!(g2m %in% int)]
    # If too few genes are in expmat, it might be mouse data, so convert them (de-capitalising seems to work better than biomaRt):
    if(sum(g1s %in% rownames(expmat)) < 10 | sum(g2m %in% rownames(expmat)) < 10) {
        g1s <- str_to_title(g1s)
        if(sum(g1s %in% rownames(expmat)) >= 30) {
            g2m <- str_to_title(g2m)
            if(sum(g2m %in% rownames(expmat)) >= 30) {
                g1s <- g1s[g1s %in% rownames(expmat)]
                g2m <- g2m[g2m %in% rownames(expmat)]
            } else return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
        } else return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
    } else {
        g1s <- g1s[g1s %in% rownames(expmat)]
        g2m <- g2m[g2m %in% rownames(expmat)]
    }
    g1s <- g1s[g1s %in% Reduce(intersect, lapply(binslist, function(x) x$symbol))]
    g2m <- g2m[g2m %in% Reduce(intersect, lapply(binslist, function(x) x$symbol))]
    
    if(length(g1s) < 20 | length(g2m) < 20) return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
    
    # Consolidate bins tables (and add <group> column to <cells>, for easy matching with <bins>):
    cat('Consolidating bins tables\n')
    if('cell_type' %in% names(cells)) {
        if('sample' %in% names(cells)) {
            if(all(unique(cells$sample) %in% bins_m$group)) {
                bins <- rbind(bins_nm[group %in% cells$cell_type], bins_m[group != 'all'])
                cells[cell_type %in% c(unique(bins_nm$group), 'Malignant'), group := ifelse(cell_type == 'Malignant', sample, cell_type)]
            } else { # Use "all" group for samples with too few malignant cells:
                bins <- rbind(bins_nm[group %in% cells$cell_type], bins_m)
                cells[
                    cell_type %in% c(unique(bins_nm$group), 'Malignant'),
                    group := ifelse(cell_type != 'Malignant' | !(sample %in% bins_m$group), cell_type, sample)
                ]
            }
        } else {
            bins <- rbind(bins_nm[group %in% cells$cell_type], bins_m[group == 'all'])
            cells[cell_type %in% c(unique(bins_nm$group), 'Malignant'), group := cell_type]
        }
        if(!all(cells[cell_type != 'Malignant', unique(cell_type)] %in% bins_nm$group)) {
            bins <- rbind(bins, bins_all[group == 'all']) # Use "all" group for non-malignant cell types not listed in bins_nm
            cells[cell_type != 'Malignant' & !(cell_type %in% bins_nm$group), group := 'all']
        }
    } else if('sample' %in% names(cells)) {
        bins <- bins_all # Use "all" group for samples with too few cells
        cells[, group := ifelse(sample %in% bins_all$group, sample, 'all')]
    } else {
        bins <- bins_all[group == 'all']
        cells[, group := 'all']
    }
    setkey(bins, group, symbol)
    
    # Compute matrices of bin averages for the cell cycle signature genes, to be used for several computations:
    cat('Computing bin averages\n')
    # To compute scores, for each signature gene g we compute its expression relative to the average expression of a sample of n control genes.
    # Here we're just computing the control gene averages, to be used later for computing scores. We typically use n = 100, but in some cases the
    # control gene sets (bins) can have fewer than 100 genes. We set a minimum of 50, so that if there are between 50 and 100 genes there is no
    # sampling to do, while if there are fewer than 50 we'll return NA. Note this procedure is designed mainly for non-malignant cell types, for
    # which we use the same control gene bins for all datasets, and which tend to have low proliferation rates, so that some cell cycle genes are
    # in the lowest bins. Since the lowest bins contain less important genes, often only a fraction of them is contained in a given dataset, and it
    # can be fewer than 100. This is especially a problem after redefining gene bins (a later step, which we could argue should not be done for
    # non-malignant cell types), after which cell cycle genes are more likely to end up in the lowest bins.
    set.seed(9432)
    bin_aves <- slapply(c('g1s', 'g2m'), function(x) {
        xout <- slapply(get(x), function(g) {
            out <- cells[, {
                grp <- unique(group)
                bins_sub <- bins[grp][bin == bins[.(grp, g), bin]]
                if(bins_sub[, .N] < 50) {
                    .(cell_name = cell_name, ave = as.numeric(NA))
                } else {
                    bin_sample <- bins_sub[, if(.N <= 100) symbol else sample(symbol, 100, replace = FALSE)]
                    .(cell_name = cell_name, ave = colMeans(expmat[bin_sample, cell_name, drop = FALSE]))
                }
            }, by = group]
            setkey(out, cell_name)
            return(out[colnames(expmat), setNames(ave, cell_name)])
        })
        xout <- set_colnames(Reduce(cbind, xout), names(xout))
        navec <- apply(xout, 1, function(r) sum(!is.na(r))) < 20
        xout[navec, ] <- NA
        return(list(xout = t(xout), na_ids = names(navec)[navec]))
    }) # Output is a matrix. Some columns (cells/groups) may contain NAs or be all NA, depending on how many sig genes had enough control genes.
    if(any(sapply(bin_aves, is.null))) return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
    
    # Remove groups with all NAs in either matrix in bin_aves (meaning they don't have enough control genes for either g1s or g2m signatures):
    na_ids_all <- unlist(lapply(bin_aves, `[[`, 'na_ids'), use.names = FALSE)
    if(length(na_ids_all) > 0) {cells <- cells[!(cell_name %in% na_ids_all)]; expmat <- expmat[, cells$cell_name]}
    bin_aves <- slapply(bin_aves, function(x) x$xout[, !(colnames(x$xout) %in% na_ids_all)])
    
    # Filter cell cycle genes:
    cat('Filtering cell cycle genes\n')
    scores <- slapply(c('g1s', 'g2m'), function(x) {
        xout <- expmat[get(x), ] - bin_aves[[x]] # Some columns (cells/groups) in xout may contain NAs or be all NA.
        xout[xout > 3] <- 3; xout[xout < -3] <- 3 # Cap relative expression levels to reduce influence of extreme values.
        colMeans(xout, na.rm = TRUE) # Vector of scores - values for some cells/groups may be NA/NaN.
    })
    cor_ids <- colnames(expmat)[!is.na(scores$g1s) & !is.na(scores$g2m)] # Compute correlations over cells with non-NA values for both g1s and g2m.
    g1s_cor <- cor(as.matrix(t(expmat[g1s, cor_ids])), data.table(g1s = scores$g1s[cor_ids], g2m = scores$g2m[cor_ids]))
    g1s_cor <- as.data.table(g1s_cor, keep.rownames = 'gene')
    g2m_cor <- cor(as.matrix(t(expmat[g2m, cor_ids])), data.table(g1s = scores$g1s[cor_ids], g2m = scores$g2m[cor_ids]))
    g2m_cor <- as.data.table(g2m_cor, keep.rownames = 'gene')
    g1s <- g1s_cor[g1s > g2m + 0.1, .(gene = gene, r = rowMeans(data.table(diff = order(order(g2m - g1s)), abs = order(order(-g1s)))))]
    g1s <- g1s[order(r), head(gene, max(20, min(50, .N)))]
    g2m <- g2m_cor[g2m > g1s + 0.1, .(gene = gene, r = rowMeans(data.table(diff = order(order(g1s - g2m)), abs = order(order(-g2m)))))]
    g2m <- g2m[order(r), head(gene, max(20, min(50, .N)))]
    if(length(g1s) < 20 | length(g2m) < 20) return(list(ccdata = NULL, hdata = NULL, g1s = NULL, g2m = NULL, path = p))
    for(x in c('g1s', 'g2m')) bin_aves[[x]] <- bin_aves[[x]][get(x), ]
    
    # Initial classification of cycling cells:
    cat('Initial classification of cycling cells\n')
    
    # The following computes mean- and count-based G1/S and G2/M scores and stores them in a data table 'ccdata':
    ccdata_vars <- c('cell_name', 'cell_type', 'sample', 'group')
    ccdata <- copy(cells[, ccdata_vars[ccdata_vars %in% names(cells)], with = FALSE])
    ccdata <- ccdata[,
        c('g1s_score', 'g1s_count', 'g2m_score', 'g2m_count') := lapply(c('g1s', 'g2m'), function(x) {
            xout <- expmat[get(x), ] - bin_aves[[x]] # Some columns (cells/groups) in xout may contain NAs or be all NA.
            xout[xout > 3] <- 3; xout[xout < -3] <- 3 # Cap relative expression levels to reduce influence of extreme values.
            # Return vectors of mean- and count-based scores - values for some cells/groups may be NA/NaN:
            list(colMeans(xout, na.rm = TRUE), colMeans(expmat[get(x), cell_name] > bin_aves[[x]], na.rm = TRUE))
        }) %>% unlist(recursive = FALSE)
    ]
    setkey(ccdata, cell_name)
    
    # Define cycling cells relative to null distributions generated by bootstrapping:
    varnames <- c('g1s_count', 'g1s_score', 'g2m_count', 'g2m_score')
    if('cell_type' %in% names(cells)) {
        null_data <- cells[,
            do.call(merge, lapply(c('g1s', 'g2m'), function(x) {
                cell_ids <- .SD[ # If cells have NA in score or count columns, it means their columns in bin_aves are all NA, so exclude these.
                    apply(do.call(`[`, list(ccdata, cell_name))[, paste0(x, c('_score', '_count')), with = FALSE], 1, function(y) all(!is.na(y))),
                    cell_name
                ]
                if(length(cell_ids) < 30) return(NULL)
                sig <- get(x)
                expmat_sub <- expmat[sig, cell_ids]
                bin_aves_sub <- bin_aves[[x]][, cell_ids]
                lapply(1:1000, function(i) return(data.table(trial = i, ind = 1:100)[,
                    (paste0(x, c('_count', '_score'))) := {
                        smpls <- slapply(sig, function(g) sample(cell_ids, 100, replace = TRUE))
                        expmat_smpl <- sapply(sig, function(g) expmat_sub[g, smpls[[g]]])
                        bin_ave_smpl <- sapply(sig, function(g) bin_aves_sub[g, smpls[[g]]])
                        scmat <- expmat_smpl - bin_ave_smpl
                        scmat[scmat > 3] <- 3; scmat[scmat < -3] <- -3
                        list(rowMeans(expmat_smpl > bin_ave_smpl, na.rm = TRUE), rowMeans(scmat, na.rm = TRUE))
                    }
                ])) %>% rbindlist
            })),
            by = cell_type
        ]
        null_max <- null_data[, slapply(.SD, max), .SDcols = varnames, keyby = .(cell_type, trial)]
        ccdata[,
            (paste0(varnames, '_n')) := lapply(varnames, function(vn) sum(get(vn) > do.call(function(x) null_max[x, get(vn)], list(cell_type)))),
            by = cell_name
        ]
        ccdata[,
            (paste0(varnames, '_p')) := lapply(
                varnames,
                function(vn) if(all(is.na(get(paste0(vn, '_n'))))) NA else {
                    ctprob <- do.call(function(x) ccdata[cell_type == x, sum(get(paste0(vn, '_n')))/(.N*1000)], list(unique(cell_type)))
                    .SD[, .(y = binom.test(get(paste0(vn, '_n')), 1000, ctprob, alternative = 'greater')$p.value), by = cell_name]$y
                }
            ),
            by = cell_type
        ]
    } else {
        null_data <- cells[,
            do.call(merge, lapply(c('g1s', 'g2m'), function(x) {
                cell_ids <- .SD[
                    apply(do.call(`[`, list(ccdata, cell_name))[, paste0(x, c('_score', '_count')), with = FALSE], 1, function(y) all(!is.na(y))),
                    cell_name
                ]
                if(length(cell_ids) < 30) return(NULL)
                sig <- get(x)
                expmat_sub <- expmat[sig, cell_ids]
                bin_aves_sub <- bin_aves[[x]][, cell_ids]
                lapply(1:1000, function(i) return(data.table(trial = i, ind = 1:100)[,
                    (paste0(x, c('_count', '_score'))) := {
                        smpls <- slapply(sig, function(g) sample(cell_ids, 100, replace = TRUE))
                        expmat_smpl <- sapply(sig, function(g) expmat_sub[g, smpls[[g]]])
                        bin_ave_smpl <- sapply(sig, function(g) bin_aves_sub[g, smpls[[g]]])
                        scmat <- expmat_smpl - bin_ave_smpl
                        scmat[scmat > 3] <- 3; scmat[scmat < -3] <- -3
                        list(rowMeans(expmat_smpl > bin_ave_smpl, na.rm = TRUE), rowMeans(scmat, na.rm = TRUE))
                    }
                ])) %>% rbindlist
            }))
        ]
        null_max <- null_data[, slapply(.SD, max), .SDcols = varnames, by = trial]
        ccdata[, (paste0(varnames, '_n')) := lapply(varnames, function(vn) sum(get(vn) > null_max[[vn]])), by = cell_name]
        ccdata[,
            (paste0(varnames, '_p')) := lapply(
                varnames,
                function(vn) if(is.na(get(paste0(vn, '_n')))) NA else {
                    vnprob <- ccdata[, sum(get(paste0(vn, '_n')))/(.N*1000)]
                    binom.test(get(paste0(vn, '_n')), 1000, vnprob, alternative = 'greater')$p.value
                }
            ),
            by = cell_name
        ]
    }
    ccdata[, (paste0(varnames, '_padj')) := lapply(varnames, function(vn) p.adjust(get(paste0(vn, '_p')), 'BH'))]
    
    ccdata[,
        c('phase_count', 'phase_score') := lapply(c('count', 'score'), function(vt) ifelse(
            (get(paste0('g1s_', vt, '_padj')) < 0.05 | get(paste0('g2m_', vt, '_padj')) < 0.05),
            'Cycling',
            'Not cycling'
        ))
    ]
    ccdata[, phase_cons := ifelse(phase_count == 'Cycling' & phase_score == 'Cycling', 'Cycling', 'Not cycling')]
    
    # ggplot(ccdata[cell_type == 'Malignant'], aes(x = g1s_score, y = g2m_score)) + geom_point(aes(colour = phase_cons)) + theme_test()
    
    cat('Redefining cell cycle gene bins\n')
    
    if('cell_type' %in% names(cells)) {
        ct_change <- ccdata[, .(prop_cyc = sum(phase_cons == 'Cycling')/.N), by = cell_type][prop_cyc > 0.1, cell_type]
        grp_change <- cells[cell_type %in% ct_change, unique(group)]
    } else if(ccdata[, sum(phase_cons == 'Cycling')/.N > 0.1]) {
        grp_change <- cells[, unique(group)]
    } else {
        grp_change <- character()
    }
    ccbins <- ccdata[phase_cons == 'Not cycling' & group %in% grp_change, .(
        symb = do.call(`[`, list(bins, unique(group)))$symbol,
        r = order(order(rowMeans(expmat[do.call(`[`, list(bins, unique(group)))$symbol, cell_name, drop = FALSE])))/nrow(expmat)
    ), by = .(grp = group)][symb %in% c(g1s, g2m)]
    
    if(nrow(ccbins) > 0) {
        
        setkey(ccbins, grp, symb)
        bin_range <- bins[group %in% ccbins$grp, .(rr = range(ave_rank)), keyby = .(group, bin)]
        ccbins <- ccbins[, bin_range[grp, .(rm = r %between% rr), by = bin], by = .(grp, symb)][,
            .(bin = ifelse(any(rm), bin[rm], 1L)),
            keyby = .(grp, symb)
        ]
        
        # Redefine bins in a copy of the bins table first, so we can check for cases of not enough control genes:
        bins_ccbins <- copy(bins)
        bins_ccbins[symbol %in% c(g1s, g2m) & group %in% ccbins$grp, bin := ccbins[.(group, symbol), bin], by = group]
        
        # Recompute matrices of bin averages for cell cycle signature genes using the updated bins, only for those groups whose bins were updated:
        cat('Recomputing bin averages\n')
        set.seed(3215)
        bin_aves_ccbins <- slapply(c('g1s', 'g2m'), function(x) {
            xout <- slapply(get(x), function(g) {
                out <- cells[group %in% ccbins$grp, {
                    grp <- unique(group)
                    bins_sub <- bins_ccbins[grp][bin == bins_ccbins[.(grp, g), bin]]
                    if(bins_sub[, .N] < 50) {
                        .(cell_name = cell_name, ave = as.numeric(NA))
                    } else {
                        bin_sample <- bins_sub[, if(.N <= 100) symbol else sample(symbol, 100, replace = FALSE)]
                        .(cell_name = cell_name, ave = colMeans(expmat[bin_sample, cell_name, drop = FALSE]))
                    }
                }, by = group]
                setkey(out, cell_name)
                return(out[colnames(expmat)[colnames(expmat) %in% cell_name], setNames(ave, cell_name)])
            })
            xout <- set_colnames(Reduce(cbind, xout), names(xout))
            navec <- apply(xout, 1, function(r) sum(!is.na(r))) < 20
            xout[navec, ] <- NA
            return(list(xout = t(xout), na_ids = names(navec)[navec]))
        }) # Output 2 matrices. Some columns (cells/groups) may contain NAs or be all NA - those that are all NA are stored in na_ids.
        
        # Don't redefine bins for groups with all NAs in bin_aves_ccbins (not enough control genes) - remove these from ccbins and bin_aves_ccbins:
        na_ids_all <- unlist(lapply(bin_aves_ccbins, `[[`, 'na_ids'), use.names = FALSE)
        ccbins <- ccbins[!(grp %in% cells[cell_name %in% na_ids_all, unique(group)])]
        bin_aves_ccbins <- slapply(bin_aves_ccbins, function(x) x$xout[, !(colnames(x$xout) %in% na_ids_all)])
        
    }
    
    if(nrow(ccbins) > 0) {
        
        # Redefine bins for the groups remaining in ccbins:
        bins[symbol %in% c(g1s, g2m) & group %in% ccbins$grp, bin := ccbins[.(group, symbol), bin], by = group]
        for(x in c('g1s', 'g2m')) {bin_aves[[x]][, colnames(bin_aves[[x]]) %in% colnames(bin_aves_ccbins[[x]])] <- bin_aves_ccbins[[x]]}
        
        # Final cycling classifications for groups with updated bins, using cells initially classified as non-cycling to define null distribution:
        cat('Final classification of cycling cells\n')
        
        ccdata[
            group %in% ccbins$grp,
            c('g1s_score', 'g1s_count', 'g2m_score', 'g2m_count') := lapply(c('g1s', 'g2m'), function(x) {
                xout <- expmat[get(x), cell_name] - bin_aves[[x]][, cell_name] # Some columns (cells/groups) in xout may contain NAs or be all NA.
                xout[xout > 3] <- 3; xout[xout < -3] <- 3 # Cap relative expression levels to reduce influence of extreme values.
                # Return vectors of mean- and count-based scores - values for some cells/groups may be NA/NaN:
                list(colMeans(xout, na.rm = TRUE), colMeans(expmat[get(x), cell_name] > bin_aves[[x]][, cell_name], na.rm = TRUE))
            }) %>% unlist(recursive = FALSE)
        ]
        
        # Define cycling cells relative to null distributions generated by bootstrapping:
        if('cell_type' %in% names(cells)) {
            ct_change <- ccdata[group %in% ccbins$grp, unique(cell_type)]
            null_ids <- ccdata[cell_type %in% ct_change & !(group %in% ccbins$grp & phase_cons == 'Cycling'), cell_name]
            null_data <- cells[
                cell_name %in% null_ids,
                do.call(merge, lapply(c('g1s', 'g2m'), function(x) {
                    cell_ids <- .SD[
                        apply(do.call(`[`, list(ccdata, cell_name))[, paste0(x, c('_score', '_count')), with = FALSE], 1, function(y) all(!is.na(y))),
                        cell_name
                    ]
                    if(length(cell_ids) < 30) return(NULL)
                    sig <- get(x)
                    expmat_sub <- expmat[sig, cell_ids]
                    bin_aves_sub <- bin_aves[[x]][, cell_ids]
                    lapply(1:1000, function(i) return(data.table(trial = i, ind = 1:100)[,
                        (paste0(x, c('_count', '_score'))) := {
                            smpls <- slapply(sig, function(g) sample(cell_ids, 100, replace = TRUE))
                            expmat_smpl <- sapply(sig, function(g) expmat_sub[g, smpls[[g]]])
                            bin_ave_smpl <- sapply(sig, function(g) bin_aves_sub[g, smpls[[g]]])
                            scmat <- expmat_smpl - bin_ave_smpl
                            scmat[scmat > 3] <- 3; scmat[scmat < -3] <- -3
                            list(rowMeans(expmat_smpl > bin_ave_smpl, na.rm = TRUE), rowMeans(scmat, na.rm = TRUE))
                        }
                    ])) %>% rbindlist
                })),
                by = cell_type
            ]
            if(nrow(null_data) > 0) {
                null_max <- null_data[, slapply(.SD, max), .SDcols = varnames, keyby = .(cell_type, trial)]
                ccdata[
                    cell_type %in% null_data$cell_type, # Sometimes cell types in ct_change don't end up in null_data (return NULL in lapply)
                    (paste0(varnames, '_n')) := lapply(
                        varnames,
                        function(vn) sum(get(vn) > do.call(function(x) null_max[x, get(vn)], list(cell_type)))
                    ),
                    by = cell_name
                ]
                ccdata[
                    cell_type %in% null_data$cell_type,
                    (paste0(varnames, '_p')) := lapply(
                        varnames,
                        function(vn) if(all(is.na(get(paste0(vn, '_n'))))) NA else {
                            ctprob <- do.call(function(x) ccdata[cell_type == x, sum(get(paste0(vn, '_n')))/(.N*1000)], list(unique(cell_type)))
                            .SD[, .(y = binom.test(get(paste0(vn, '_n')), 1000, ctprob, alternative = 'greater')$p.value), by = cell_name]$y
                        }
                    ),
                    by = cell_type
                ]
            }
        } else {
            # We're under the condition that nrow(ccbins) > 0, so we know that at least some groups have updated bins. Since there is no cell_type
            # column, the null data and p values are defined across all cells, hence we need to recompute null data using all cells.
            null_data <- cells[,
                do.call(merge, lapply(c('g1s', 'g2m'), function(x) {
                    cell_ids <- .SD[
                        apply(do.call(`[`, list(ccdata, cell_name))[, paste0(x, c('_score', '_count')), with = FALSE], 1, function(y) all(!is.na(y))),
                        cell_name
                    ]
                    if(length(cell_ids) < 30) return(NULL)
                    sig <- get(x)
                    expmat_sub <- expmat[sig, cell_ids]
                    bin_aves_sub <- bin_aves[[x]][, cell_ids]
                    lapply(1:1000, function(i) return(data.table(trial = i, ind = 1:100)[,
                        (paste0(x, c('_count', '_score'))) := {
                            smpls <- slapply(sig, function(g) sample(cell_ids, 100, replace = TRUE))
                            expmat_smpl <- sapply(sig, function(g) expmat_sub[g, smpls[[g]]])
                            bin_ave_smpl <- sapply(sig, function(g) bin_aves_sub[g, smpls[[g]]])
                            scmat <- expmat_smpl - bin_ave_smpl
                            scmat[scmat > 3] <- 3; scmat[scmat < -3] <- -3
                            list(rowMeans(expmat_smpl > bin_ave_smpl, na.rm = TRUE), rowMeans(scmat, na.rm = TRUE))
                        }
                    ])) %>% rbindlist
                }))
            ]
            if(nrow(null_data) > 0) {
                null_max <- null_data[, slapply(.SD, max), .SDcols = varnames, by = trial]
                ccdata[, (paste0(varnames, '_n')) := lapply(varnames, function(vn) sum(get(vn) > null_max[[vn]])), by = cell_name]
                ccdata[,
                    (paste0(varnames, '_p')) := lapply(
                        varnames,
                        function(vn) if(is.na(get(paste0(vn, '_n')))) NA else {
                            vnprob <- ccdata[, sum(get(paste0(vn, '_n')))/(.N*1000)]
                            binom.test(get(paste0(vn, '_n')), 1000, vnprob, alternative = 'greater')$p.value
                        }
                    ),
                    by = cell_name
                ]
            }
        }
        ccdata[, (paste0(varnames, '_padj')) := lapply(varnames, function(vn) p.adjust(get(paste0(vn, '_p')), 'BH'))]
        
        ccdata[,
            c('phase_count', 'phase_score') := lapply(c('count', 'score'), function(vt) ifelse(
                (get(paste0('g1s_', vt, '_padj')) < 0.05 | get(paste0('g2m_', vt, '_padj')) < 0.05),
                'Cycling',
                'Not cycling'
            ))
        ]
        ccdata[, phase_cons := ifelse(phase_count == 'Cycling' & phase_score == 'Cycling', 'Cycling', 'Not cycling')]
        
    }
    
    ccdata[,
        c('phase_g1s', 'phase_g2m') := .(all(c(g1s_score_padj, g1s_count_padj) < 0.05), all(c(g2m_score_padj, g2m_count_padj) < 0.05)),
        by = cell_name
    ]
    if('cell_type' %in% names(ccdata)) {
        ccdata[,
            c('lr_g1s', 'lr_g2m') := .(
                glm(phase_g1s ~ g1s_score, family = binomial(), .SD)$fitted.values > 0.5,
                glm(phase_g2m ~ g2m_score, family = binomial(), .SD)$fitted.values > 0.5
            ),
            by = cell_type
        ]
    } else {
        ccdata[,
            c('lr_g1s', 'lr_g2m') := .(
                glm(phase_g1s ~ g1s_score, family = binomial(), .SD)$fitted.values > 0.5,
                glm(phase_g2m ~ g2m_score, family = binomial(), .SD)$fitted.values > 0.5
            )
        ]
    }
    ccdata[, cycling := ifelse(lr_g1s | lr_g2m, 'Cycling', 'Not cycling')]
    ccdata[, c('phase_g1s', 'phase_g2m') := NULL]
    
    # Re-centre the scores relative to means of the non-cycling cells:
    if('cell_type' %in% names(ccdata)) {
        ccdata[,
            c('g1s_score', 'g2m_score') := lapply(.(g1s_score, g2m_score), function(x) x - mean(x[cycling == 'Not cycling'])),
            by = cell_type
        ]
    } else {
        ccdata[, c('g1s_score', 'g2m_score') := lapply(.(g1s_score, g2m_score), function(x) x - mean(x[cycling == 'Not cycling']))]
    }
    
    # ggplot(ccdata[cell_type == 'Malignant'], aes(x = g1s_score, y = g2m_score)) + geom_point(aes(colour = cycling)) + theme_test()
    
    # Define separate data for heatmap:
    cat('Constructing heatmap data\n')
    annot_cols <- c('cell_name', 'g1s_score', 'g2m_score', 'cycling')
    if('cell_type' %in% names(ccdata)) {
        hdata <- ccdata[,
            merge(.SD, as.data.table(apply(expmat[c(g1s, g2m), cell_name], 1, function(x) x - mean(x)), keep.rownames = 'cell_name')),
            .SDcols = annot_cols,
            by = cell_type
        ][, -'cell_type']
    } else {
        hdata <- ccdata[,
            merge(.SD, as.data.table(apply(expmat[c(g1s, g2m), ], 1, function(x) x - mean(x)), keep.rownames = 'cell_name')),
            .SDcols = annot_cols
        ]
    }
    
    cat('Done!\n\n')
    return(list(ccdata = ccdata, hdata = hdata, g1s = g1s, g2m = g2m, path = p))
    
})

if(!all(sapply(out, function(x) is.null(x$ccdata)))) {
    ct <- gsub('/', '-', r[2])
    if(!(ct %in% dir('../data/study_plots'))) {dir.create(paste0('../data/study_plots/', ct))}
    if(!(r[1] %in% dir(paste0('../data/study_plots/', ct)))) {dir.create(paste0('../data/study_plots/', ct, '/', r[1]))}
    saveRDS(out, paste0('../data/study_plots/', ct, '/', r[1], '/data_cc.RDS'))
}
