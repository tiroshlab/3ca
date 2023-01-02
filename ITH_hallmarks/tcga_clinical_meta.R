library(data.table)
library(magrittr)
library(ggplot2)
library(stringr)
library(randomcoloR)
library(plyr)
library(scales)
library(latex2exp)
library(RColorBrewer)
library(cowplot)
library(ggrepel)
library(matkot)





# Combine results with and without deconvolution into the same heatmaps:

clin_test_res <- rbind(fread('../data/clin_test_res.csv'), fread('../data/clin_test_res_deconv.csv')[, cancer_type := paste(cancer_type, '(d)')])
clin_test_res[, pval_adj := p.adjust(pval, method = 'BH')]

setkey(clin_test_res, cancer_type, meta_program)
cancer_types <- sort(unique(clin_test_res$cancer_type))
mps <- unique(clin_test_res$meta_program)

pdata <- slapply(
    c('OS', 'PFI', 'stage', 'grade', 't_stage', 'm_stage', 'ln', 'ther'),
    function(v) {
        htmp_data <- clin_test_res[var_name == v][CJ(cancer_types, mps), -c('var_name', 'pval')]
        htmp_data[is.na(eff), eff := 0]
        htmp_data[is.na(pval_adj), pval_adj := 1]
        clust_data <- dcast(htmp_data[, -'pval_adj'], meta_program ~ cancer_type)[,
            set_rownames(as.matrix(.SD), meta_program),
            .SDcols = -'meta_program'
        ]
        clust_data <- clust_data[apply(clust_data, 1, function(x) !all(x == 0)), apply(clust_data, 2, function(x) !all(x == 0))]
        clust_x <- colSums(abs(clust_data))
        clust_x <- data.table(ctd = names(clust_x), ct = gsub(' \\(d\\)', '', names(clust_x)), vald = clust_x)[, val := mean(vald), by = ct]
        clust_x_temp <- clust_x[, unique(.SD)[, .(cancer_type = ct, r = order(order(-val))/.N)], .SDcols = c('ct', 'val')]
        setkey(clust_x_temp, cancer_type)
        clust_x[, r := clust_x_temp[ct, r]]
        clust_y <- list(labels = rownames(clust_data), order = order(-rowSums(abs(clust_data))))
        htmp_data <- htmp_data[cancer_type %in% colnames(clust_data) & meta_program %in% rownames(clust_data)]
        return(list(data = htmp_data, clust_x = clust_x[, .(ctd, ct, r)], clust_y = clust_y))
    }
)

mps_data <- rbindlist(lapply(names(pdata), function(v) with(pdata[[v]]$clust_y, data.table(v = v, mp = labels, r = order(order)/length(order)))))
setkey(mps_data, v, mp)
mps_data <- mps_data[CJ(c('OS', 'PFI', 'stage', 'grade', 't_stage', 'm_stage', 'ln', 'ther'), mps)]
mps_data[, rmean := mean(r[!is.na(r)]), by = mp]
mps_ord <- unique(mps_data[, .(mp, rmean)])[is.finite(rmean), mp[order(rmean)]]

cts_data <- rbindlist(lapply(pdata, `[[`, 'clust_x'))
cts_data_temp <- cts_data[, .(r = mean(r)), keyby = ct]
cts_ord <- unique(cts_data[, .(ctd = ctd, ct = ct, r = do.call(`[`, list(cts_data_temp, ct))$r)])[order(ctd), ctd[order(r)]]

# Overall survival:
os_data <- copy(pdata$OS)$data
os_data[, c('cancer_type', 'meta_program') := .(factor(cancer_type, levels = cts_ord), factor(meta_program, levels = mps_ord))]
os_data[, signif := ifelse(pval_adj < 0.05, ifelse(pval_adj < 0.01, ifelse(pval_adj < 0.001, '***', '**'), '*'), '')]
htmp_os <- ggplot(os_data, aes(x = cancer_type, y = meta_program, fill = exp(eff))) +
    geom_tile(colour = 'grey90') +
    geom_text(aes(label = signif), hjust = 0.5, vjust = 0.5, size = 3) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(colours = rev(brewer.pal(9, "PiYG")), limits = c(0.5, 2.5),
        values = rescale(c(0, 1, 2, 3, 4, 7, 10, 13, 16), to = c(0, 1)), oob = squish, na.value = '#F7F7F7', name = 'Hazard\nratio',
        guide = guide_colourbar(frame.colour = 'black', ticks.colour = 'black')) +
    theme(panel.border = element_rect(size = 0.5, fill = NA), axis.text.x = element_text(angle = 55, hjust = 1),
        legend.title = element_text(size = 13)) +
    labs(title = 'Overall survival', x = 'Cancer type', y = NULL)

# LN mets and therapy resistance:
ln_data <- copy(pdata$ln)$data
setkey(ln_data, cancer_type, meta_program)
ther_data <- copy(pdata$ther)$data
setkey(ther_data, cancer_type, meta_program)
mps_ln_ther <- union(ln_data$meta_program, ther_data$meta_program)
ln_data <- ln_data[CJ(unique(cancer_type), mps_ln_ther)]
ln_data[, c('cancer_type', 'meta_program') := .(factor(cancer_type, levels = cts_ord), factor(meta_program, levels = mps_ord))]
ln_data[, signif := ifelse(is.na(pval_adj), '', ifelse(pval_adj < 0.05, ifelse(pval_adj < 0.01, ifelse(pval_adj < 0.001, '***', '**'), '*'), ''))]
ther_data <- ther_data[CJ(unique(cancer_type), mps_ln_ther)]
ther_data[, c('cancer_type', 'meta_program') := .(factor(cancer_type, levels = cts_ord), factor(meta_program, levels = mps_ord))]
ther_data[, signif := ifelse(is.na(pval_adj), '', ifelse(pval_adj < 0.05, ifelse(pval_adj < 0.01, ifelse(pval_adj < 0.001, '***', '**'), '*'), ''))]
htmp_ln <- ggplot(ln_data, aes(x = cancer_type, y = meta_program, fill = eff)) +
    geom_tile(colour = 'grey90') +
    geom_text(aes(label = signif), hjust = 0.5, vjust = 0.5, size = 3) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(colours = rev(brewer.pal(9, "PiYG")), limits = c(-0.5, 0.5), breaks = c('-0.5' = -0.5, '0' = 0, '0.5' = 0.5), oob = squish,
        na.value = '#F7F7F7', name = expression(Delta*'(score)'), guide = guide_colourbar(frame.colour = 'black', ticks.colour = 'black')) +
    theme(panel.border = element_rect(size = 0.5, fill = NA), axis.text.x = element_text(angle = 55, hjust = 1),
        legend.title = element_text(size = 13)) +
    labs(title = 'Lymph node metastasis', x = 'Cancer type', y = NULL)
htmp_ther <- ggplot(ther_data, aes(x = cancer_type, y = meta_program, fill = eff)) +
    geom_tile(colour = 'grey90') +
    geom_text(aes(label = signif), hjust = 0.5, vjust = 0.5, size = 3) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(colours = rev(brewer.pal(9, "PiYG")), limits = c(-0.5, 0.5), breaks = c('-0.5' = -0.5, '0' = 0, '0.5' = 0.5), oob = squish,
        na.value = '#F7F7F7', name = expression(Delta*'(score)'), guide = guide_colourbar(frame.colour = 'black', ticks.colour = 'black')) +
    theme(panel.border = element_rect(size = 0.5, fill = NA), axis.text.x = element_text(angle = 55, hjust = 1),
        legend.title = element_text(size = 13)) +
    labs(title = 'Therapy resistance', x = 'Cancer type', y = NULL)

# Make grobs and fix widths and heights:
htmp_os_grob <- ggplotGrob(htmp_os)
htmp_os_grob$widths <- unit(c(0.3, 0, 0, 4.5, 0.5*os_data[, length(unique(cancer_type))], 0, 0, 0.5, 2, 0, 0.3), 'cm')
htmp_os_grob$heights <- unit(c(0.3, 0, 0.8, 0, 0, 0, 0.4*os_data[, length(unique(meta_program))], 2.5, 0.6, 0, 0, 0.3), 'cm')
htmp_ln_grob <- ggplotGrob(htmp_ln + theme(legend.position = 'none'))
htmp_ln_grob$widths <- unit(c(0.3, 0, 0, 4.5, 0.5*ln_data[, length(unique(cancer_type))], 0, 0, 0, 0.3), 'cm')
htmp_ln_grob$heights <- unit(c(0.3, 0, 0.8, 0, 0, 0, 0.4*ln_data[, length(unique(meta_program))], 2.5, 0.6, 0, 0, 0.3), 'cm')
htmp_ther_grob <- ggplotGrob(htmp_ther + theme(axis.text.y = element_blank()))
htmp_ther_grob$widths <- unit(c(0.3, 0, 0, 0, 0.5*ther_data[, length(unique(cancer_type))], 0, 0, 0.5, 2, 0, 0.3), 'cm')
htmp_ther_grob$heights <- unit(c(0.3, 0, 0.8, 0, 0, 0, 0.4*ther_data[, length(unique(meta_program))], 2.5, 0.6, 0, 0, 0.3), 'cm')

pdf('../data/clin_htmps_sub.pdf', width = sum(c(htmp_ln_grob$widths, htmp_ther_grob$widths))/2.54,
    height = sum(c(htmp_os_grob$heights, htmp_ln_grob$heights))/2.54)
plot_grid(
    plot_grid(
        htmp_os_grob,
        ggplot() + theme_void(),
        nrow = 1,
        ncol = 2,
        rel_widths = c(
            as.numeric(sum(htmp_os_grob$widths)),
            sum(c(htmp_ln_grob$widths, htmp_ther_grob$widths)) - as.numeric(sum(htmp_os_grob$widths))
        )
    ),
    plot_grid(
        htmp_ln_grob,
        htmp_ther_grob,
        nrow = 1,
        ncol = 2,
        rel_widths = c(sum(htmp_ln_grob$widths), sum(htmp_ther_grob$widths))
    ),
    nrow = 2,
    ncol = 1,
    rel_heights = c(sum(htmp_os_grob$heights), sum(htmp_ln_grob$heights))
)
dev.off()





# Summary showing pan-cancer consistency of most common MPs:

clin_test_res <- rbind(fread('../data/clin_test_res.csv'), fread('../data/clin_test_res_deconv.csv')[, cancer_type := paste(cancer_type, '(d)')])
clin_test_res[, pval_adj := p.adjust(pval, method = 'BH')]

study_contrib <- readRDS('../data/study_contribution_per_MP.RDS')
names(study_contrib) <- gsub('  ', ' ', names(study_contrib))
names(study_contrib)[41] <- 'MP41 Unassigned'

study_tcga_map <- fread('../data/study_tcga_map.csv', na.strings = '', key = 'study')

mps_tab <- lapply(
    names(study_contrib),
    function(mp_name) study_tcga_map[study_contrib[[mp_name]], .(meta_program = mp_name, study = study, cancer_type = cancer_type)]
) %>% rbindlist

common_mps <- mps_tab[, .(n_study = .N, n_ct = length(unique(cancer_type[!is.na(cancer_type)]))), by = meta_program][n_ct >= 5 & n_study >= 10]
common_mps <- common_mps$meta_program

pdata_scatter <- clin_test_res[meta_program %in% common_mps, .(s = sum(sign(eff)), s_sig = sum(sign(eff[pval_adj < 0.05]))), by = .(meta_program)]

set.seed(8331)
pdf('../data/clin_summary_scatter.pdf', width = 14/2.54, height = 11.5/2.54)
ggplot(pdata_scatter, aes(x = s_sig, y = s)) +
    theme_test() +
    geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5, colour = 'lightgrey') +
    geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5, colour = 'lightgrey') +
    geom_point(shape = 21, fill = 'tomato', stroke = 0.3, size = 3) +
    geom_text_repel(data = pdata_scatter[grep('Cell Cycle|Stress \\(in vitro\\)|Proteasomal|Hypoxia', meta_program)], aes(label = meta_program),
        colour = 'steelblue') +
    geom_text_repel(data = pdata_scatter[grep('Stress$|PDAC', meta_program)], aes(label = meta_program),
        colour = 'steelblue', nudge_x = 14, nudge_y = 8) +
    geom_text_repel(data = pdata_scatter[grep('MHC-II \\(I\\)', meta_program)], aes(label = meta_program),
        colour = 'steelblue', nudge_x = 1, nudge_y = -1) +
    labs(x = expression(sum()*italic(sign)*'(significant effects)'), y = expression(sum()*italic(sign)*'(all effects)'))
dev.off()
