#!/usr/bin/R5

############
# Figure 5 #
############

# This script will create plots for figure 5 from single-cell data

.libPaths('~/R/newer_packs_library/3.5/')
library(Seurat)
source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')

## Destiny Folder ##
root <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma/results'
setwd(paste0(root, '/a1_final_figures/Figure_5'))
setwd(paste0(root, '/a1_final_figures/Figure_6_Th2diseasecomp_IL9'))

mycells <- theObjectSavedIn(paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/clustering/zetInfo/clustCells13PCs_30Ks_0.06667JD.RData'))
gr.cols <- read.csv("/mnt/BioHome/ciro/asthma/info/AS3_final_groupColours.csv",row.names=1,stringsAsFactors=F,header=1)

newlabs <- c("ACT1", "ACT2", "ACT3", "TH1", "TH2", "THIFNR", "TH17")
names(newlabs) <- 0:6
gorder <- c("TH2", "TH1", "TH17", "THIFNR", "ACT3", "ACT2", "ACT1")
mycells@meta.data$Cell_type <- factor(c(newlabs[mycells@meta.data$RNA_snn_res.0.4]), gorder)

### 5.B Volcano TH2 AS_AL vs AR ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- mycells@meta.data
expdata <- cts2cpm(mycells@assays$RNA@counts)
compid <- "AS_ALvsAR"
wgrcol <- "orig.diseasegroup"
gr <- "4"
res <- readfile(paste0(root, '/mast/AS3EsmTeff_daf_30p/comprs/dg', gr, "/", compid, '/results_', compid, '_mastlog2cpm.csv'), stringsAsFactors = FALSE, check.names = FALSE)
wgr <- unlist(strsplit(compid, "vs"))
thesecells <- getsubset(list(c("Cell_type", newlabs[gr]), c(wgrcol, wgr)), annot, v = TRUE)
suffy <- ""

cat(newlabs[gr], "\n")
res$gene <- gsub("'", "", res$gene_name)
rownames(res) <- res$gene
res <- res[, !duplicated(colnames(res))]
grps <- make_list(x = annot[thesecells, ], colname = wgrcol, grouping = TRUE)
grps <- sort(grps)
table(grps)
# ## Checking stats are correct
# expdata <- mycells@assays$RNA@data[, thesecells]
# stattab <- get_stat_report(mat = expdata, groups = grps, rname = res$gene, moments = c("mn"), v = TRUE)
# head(stattab)
# head(res[, grepl("_mean", colnames(res))])
# head(log2(stattab[res$gene, ] + 1))
## Correcting
# expdata <- cts2cpm(mycells@assays$RNA@counts[, thesecells])
# stattab <- get_stat_report(mat = log2(expdata[, thesecells] + 1), groups = grps, rname = res$gene, moments = c("mn", "p"), v = TRUE)
stattab <- get_stat_report(mat = expdata[, thesecells], groups = grps, rname = res$gene, moments = c("mn", "p"), v = TRUE)
stattab[, grepl("_mean", colnames(stattab))] <- log2(stattab[, grepl("_mean", colnames(stattab))] + 1)
head(stattab)
summary(stattab)
head(res[, grepl("_mean", colnames(res))])
summary(res[, grepl("_mean", colnames(res))])
res <- cbind(stattab, res)
head(res)
res <- res[, !duplicated(colnames(res))]

fname <- paste0(root, '/mast/AS3EsmTeff_daf_30p/summary/clusters/padj0.05_FC0/RNA_snn_res.0.4_SummaryDEGsTable_suas.csv')
res_spec <- readfile(fname, stringsAsFactors = FALSE, check.name = FALSE)
res_spec <- res_spec[!is.na(res_spec$group), ]
rownames(res_spec) <- res_spec$gene_name
genes <- getsubset(c('group', gr), res_spec, v = TRUE)
res$TH2_specific <- rownames(res) %in% sub("'", "", genes)
head(res)

tres <- res[, c("gene_name", "group", "log2FoldChange", "padj", colnames(stattab))]
head(tres)
write.csv(tres, file = paste0("results_", newlabs[gr], '_', compid, suffy, ".csv"), row.names = FALSE)

res$Fold <- res$log2FoldChange
res$FDR <- -log10(res$padj)
res$pct <- res[, paste0(wgr[1], "_exprFrac")]
res$pct[res$group == wgr[2]] <- res[res$group == wgr[2], paste0(wgr[2], "_exprFrac")]
tvar <- res$padj <= 0.05 & abs(res$log2FoldChange) >= 0.25
table(tvar)
res$pct[!tvar] <- 0 # size for only significant genes
res$mean <- res[, paste0(wgr[1], "_mean")]
res$mean[res$group == wgr[2]] <- res[res$group == wgr[2], paste0(wgr[2], "_mean")]
degs <- getDEGenes(x =res, pv = 0.05, fc = 0.25, v = TRUE)
res$mean[!res$gene %in% degs] <- NA

# p <- ggplot(res, aes(x = Fold, y = FDR)) +
#   xlim(c(-max(abs(x$Fold)), max(abs(x$Fold)))) +
#   geom_point(data = x, aes(color = Bmean, size = pct))
# pdf(paste0('volcano_manual_', newlabs[gr], '_', compid, '.pdf'), width = 10, height = 10)
# print(p)
# dev.off()

# all genes or a subset
if(TRUE){
  tvar <- res
  fsufix <- paste0(newlabs[gr], '_', compid, suffy, "_all_genes")
}else{
  tvar <- res[res$gene_name %in% genes, ]
  fsufix <- paste0(newlabs[gr], '_', compid, suffy, "")
}
source("~/scripts/functions/volcano_variant_color.R")
p <- try(volplot(
  x = tvar,
  pvalth = 0.05,
  lfcth = 0.25,
  pvaltype = 'padj',
  lfctype = 'log2FoldChange',
  col_feature = "mean",
  gene_name = "gene_name",
  size_feature = "pct",
  ngenes = 10,
  return_plot = TRUE,
  v = TRUE
))
pdf(paste0('nvolcano_colours_', fsufix, '_gray.pdf'), width = 10, height = 10)
print(p)
dev.off()
pdf(paste0('nvolcano_colours_', fsufix, '_grayblank.pdf'), width = 10, height = 10)
print(shut_it(p))
dev.off()

ttvar <- tvar
ttvar$FDR <- trans(ttvar$FDR, 10)
ttvar$log2FoldChange <- trans(ttvar$log2FoldChange, 4, 0.05)
tvar[bordering(tvar, cnames="FDR", 6), c("FDR", "log2FoldChange")]
ttvar[bordering(ttvar, cnames="FDR", 6), c("FDR", "log2FoldChange")]
tvar[bordering(tvar, cnames="log2FoldChange", 6), c("FDR", "log2FoldChange")]
ttvar[bordering(ttvar, cnames="log2FoldChange", 6), c("FDR", "log2FoldChange")]
yticks <- c(0, 5, 10, 15)
xticks <- c(-4, -2, -1, 0, 1, 2, 4)
source("~/scripts/functions/volcano_variant_color.R")
p <- try(volplot(
  x = ttvar,
  pvalth = trans(-log10(0.05)),
  lfcth = trans(0.25),
  pvaltype = 'FDR',
  lfctype = 'log2FoldChange',
  col_feature = "mean",
  size_feature = "pct",
  gene_name = "gene_name",
  do_fdr_log = FALSE,
  ngenes = 10,
  return_plot = TRUE,
  v = TRUE
)) + scale_x_continuous(limits = c(-4.2, 4.2), breaks = trans(xticks, 4, 0.05), labels = xticks) +
  scale_y_continuous(limits = c(0, NA), breaks = trans(yticks, 10), labels = yticks)
pdf(paste0('nvolcano_colours_', fsufix, '_gray_breaks.pdf'), width = 10, height = 10)
print(p)
dev.off()
pdf(paste0('nvolcano_colours_', fsufix, '_gray_breaksblank.pdf'), width = 10, height = 10)
print(shut_it(p))
dev.off()

### 5.X Heatmap TH2 AS_AL/AR ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FNAMES=(`ls`)
# for FNAME in ${FNAMES[@]}; do
#   mv ${FNAME} ${FNAME/heatmap_/}
# done
dir.create("heatmaps")
library(pheatmap)
coulrange = c('blue', 'black', 'yellow')
mypalette <- colorRampPalette(coulrange, space = 'Lab')
bounds <- 1
palettebreaks <- seq(-bounds, bounds, 0.1) # Heatmap setting

gr <- "4"
annot <- mycells@meta.data
annot$sCell_type <- "NCT"
annot[rownames(mycells_th2@meta.data), ]$sCell_type <- mycells_th2@meta.data[, "Cell_type"]
cdg <- "sCell_type"
dg <- paste0("TH2_", c("0", "1"))
clord <- rlord <- TRUE

# annot <- mycells@meta.data
# cdg <- "orig.diseasegroup"
# dg <- unlist(strsplit(compid, "vs"))[1]
# genes <- getDEGenes(x = res, pv = 0.05, fc = 0.5, upreg = ifelse(length(dg) == 2, NA, dg == "AS_AL"), v = TRUE)

genes <- c("EIF3J", "EIF5B", "CALM3", "FKBP1A", "ATP13A3", "PSMD13", "UBE2S", "DUSP4", "BTLA", "TGFBR3",
  "IL17RB", "IL1RL1", "CSF2", "IL13", "IL5", "IL9", "GZMB", "PPARG", "ZEB2", "CD28", "RAB27A", "CTLA4",
  "MAP3K8", "DUSP6", "ICOS", "IL4", "IL21", "IL3", "IL31", "CFLAR", "BCL2A1", "TNFSF14", "NFKBID",
  "FOSL2", "IRF4", "STAT4", "GATA3", "RUNX3", "SATB1", "PDCD1", "NIPA1", "NEDD9", "BIRC3")
genes <- getfound(genes, rownames(mycells), v = TRUE)
cutoff <- c(0.1, 10)
# suffix <- paste0("_", newlabs[gr], "_", paste0(dg, collapse = "n"))
filters <- list(c("Cell_type", newlabs[gr]), c(cdg, dg), c("orig.diseasegroup", "AS_AL", "AR"), c("orig.donor", "-1838"))
suffix <- paste0(unlist(sapply(filters, function(x) x[-1] )), collapse = "_")

# ### 5.X TH2 polyfunction ###
# genes <- c("IL2", "IL5", "IL4", "IL13", "IL21", "IL9", "IL3", "IL31")
# # genes <- rownames(mycells)[grepl("^IL", rownames(mycells))]
# cutoff <- c(0.1, 10)
# suffix <- paste0("_", newlabs[gr], "_", paste0(dg, collapse = "n"), "_polyfunction")

annot <- remove.factors(annot[order(annot[, cdg]), ])
thesecells <- getsubset(filters, annot, v = TRUE)
tvar <- sapply(annot, function(x) length(table(x)) )
annot <- annot[thesecells, (tvar < 25 & tvar > 1) & !grepl("clusters|res|lib|ST|ND|gem", colnames(annot))]
head(annot)

matex <- round(cts2cpm(mycells@assays$RNA@counts[, thesecells]), 2)
# matex <- log2(cts2cpm(saveymat[, thesecells]) + 1)

matex <- matex[genes, ]
stattab <- get_stat_report(mat = matex, moments = c("mn", "p"), v = TRUE)
dim(matex)
tvar <- stattab[, 1] > cutoff[1] & stattab[, 2] > cutoff[2]
cat(sum(tvar), "of", nrow(matex), "\n")
genes[!tvar]
matex <- matex[, ]
dim(matex)
matext <- matex

headmat(matex)
# # write.csv(matex, file = paste0('heatmaps/', newlabs[gr], suffix, '_matrix.csv'))
# qmatf <- qlucore_format(mat = matex, metadata = annot, v = TRUE)
# write.csv(qmatf, file = paste0('heatmaps/', newlabs[gr], suffix, '_qlucore.csv'))
# matex <- log2(matex + 1)

# # no more modification on columns
# matcolms <- annot
# colnames(matcolms) <- colnames(annot)

# average
# matex <- t(scale(apply(matex, 1, function(vec) tapply(vec, annot[, cdg], mean, na.rm = T) )))
# mean per disease
clor <- c(Donor = "orig.donor", Disease = "orig.diseasegroup", Cluster = "sCell_type")[1]
annot$cl_donor <- paste0(annot[, cdg], "_", annot[, clor]); table(annot$cl_donor)
tvar <- annot; tvar <- tvar[order(tvar$cl_donor), ]
matex <- t(apply(matex[, rownames(tvar)], 1, function(vec) tapply(vec, annot[, "cl_donor"], mean, na.rm = T) ))
matcolms <- remove.factors(tvar[!duplicated(tvar$cl_donor), ])
rownames(matcolms) <- matcolms$cl_donor
matcolms <- matcolms[rev(with(matcolms, order(sCell_type, orig.diseasegroup))), -ncol(matcolms)]
if(names(clor) == "Disease") matcolms <- matcolms[, c("orig.diseasegroup", "sCell_type")]
if(names(clor) == "Cluster") matcolms <- matcolms[, c("sCell_type"), drop = FALSE]
suffix <- paste0(suffix, "_meanPer", names(clor))

# add disease avarage!!!!!!!!!!!!!! Christ!
clor <- c(Donor = "orig.donor", Disease = "orig.diseasegroup", Cluster = "sCell_type")[2]
annot$cl_donor <- paste0(annot[, cdg], "_", annot[, clor]); table(annot$cl_donor)
tvar <- annot; tvar <- tvar[order(tvar$cl_donor), ]
tmp <- t(apply(matext[, rownames(tvar)], 1, function(vec) tapply(vec, annot[, "cl_donor"], mean, na.rm = T) ))
matex <- cbind(tmp, matex)
nmatcolms <- rbind(matcolms, mat_names(rev(colnames(tmp)), colnames(matcolms)))[, c("sCell_type", "orig.diseasegroup")]
nmatcolms[rev(colnames(tmp)), ] <- data.frame(sCell_type = rep(c("TH2_1", "TH2_0"), each = 2), orig.diseasegroup = rev(colnames(tmp)), stringsAsFactors = F)
matcolms <- nmatcolms
gr.cols[names(table(annot$cl_donor)), ] <- gr.cols[sub("TH2_._", "", names(table(annot$cl_donor))), ]
suffix <- paste0(suffix, "_diseaseadded")

# add cluster avarage
clor <- c(Donor = "orig.donor", Disease = "orig.diseasegroup", Cluster = "sCell_type")[3]
annot$cl_donor <- paste0(annot[, cdg], "_", annot[, clor]); table(annot$cl_donor)
tvar <- annot; tvar <- tvar[order(tvar$cl_donor), ]
tmp <- t(apply(matext[, rownames(tvar)], 1, function(vec) tapply(vec, annot[, "cl_donor"], mean, na.rm = T) ))
matex <- cbind(tmp, matex)
nmatcolms <- rbind(matcolms, mat_names(rev(colnames(tmp)), colnames(matcolms)))[, c("sCell_type", "orig.diseasegroup")]
nmatcolms[rev(colnames(tmp)), ] <- data.frame(sCell_type = c("TH2_1", "TH2_0"), orig.diseasegroup = c("TH2_1", "TH2_0"), stringsAsFactors = F)
matcolms <- nmatcolms
matex <- matex[, rev(rownames(matcolms))]
suffix <- paste0(suffix, "_clusteradded")

fname <- paste0(root, '/a1_final_figures/Figure_5_TH2co-exp_and_sub-clustering/correlations/saver_correlation_clustTH2_specTH2_quantile_pct5_genes214.csv')
list.files(fname)
clusts <- read.csv(fname, row.names= 2)
clusts[genes, 2] <- as.character(clusts[genes, 2])
# # adding cluster row colour and ordering hclust
# matrows <- clusts[genes, 2, drop = FALSE]
# rlord <- TRUE
# # order genes on module
# matrows <- clusts[rownames(matrows) %in% rownames(matex), 2, drop = FALSE]
# matex <- matex[rownames(matrows), ]
# rlord <- FALSE
# suffix <- paste0(suffix, "_modor")
# selected module order
thisorder <- c("3" = "1", "5" = "2", "1" = "3", "2" = "4", "4" = "5")
clusts$gregs <- thisorder[clusts$cluster]
matrows <- clusts[unlist(orderthis(clusts, order_by = "gregs", grupos = names(thisorder), v = TRUE)), ]
head(clusts)
head(matrows)
matrows <- matrows[rownames(matrows) %in% rownames(matex), 3, drop = FALSE]
matex <- matex[rownames(matrows), ]
rlord <- FALSE
suffix <- paste0(suffix, "_gregor_cols")
tcols <- v2cols(1:5)[names(thisorder)]
names(tcols) <- paste0("C", 1:5)
matrows[, 1] <- paste0("C", matrows[, 1])
gr.cols[names(tcols), "colour"] <- tcols

head(matcolms)
dopca <- FALSE
if(dopca){
  pc1 <- orderthis(annot = matcolms, order_by = 'pca', cname = cdg, mat = matex, v = TRUE)
  pc1 <- unlist(pc1)
  matex <- matex[, pc1]
  matcolms <- matcolms[pc1, ]

  # pc1 <- prcomp_irlba(matex, n = min(c(3, dim(matex) - 1)) )$x[, 1]
  # names(pc1) <- genes
  # matex <- matex[rev(names(sort(pc1, decreasing = T))), ]
  headmat(matex)
  mysufix <- paste0("", suffix, "_pca")
}else if(all(c(cdg, "orig.diseasegroup") %in% colnames(matcolms)) && !grepl("clusteradded", suffix)){
  setorderv(matcolms, c(cdg, "orig.diseasegroup"))
  clord <- FALSE
  matex <- matex[, rownames(matcolms)]
  mysufix <- paste0("", suffix, "_diseaseor")
}else{
  clord <- FALSE
  mysufix <- paste0("", suffix, "")
}

# ## this is the same order pheatmap gives!!
# d <- dist(as.matrix(matexz))
# hc <- hclust(d)
# head(hc$labels[hc$order])

matcouls <- lapply(matcolms, v2cols, gr.cols)
matcouls <- c(matcouls, lapply(matrows, v2cols, gr.cols))
# manually scale
# matexz <- t(scale(t(matex), scale = ncol(matex) > 2))
# matexz <- t(scale(t(matex), center = ncol(matex) > 2))
matexz <- t(scale(t(matex)))
matexz[matexz < (-bounds)] <- -bounds
matexz[matexz > (bounds)] <- bounds

x <- pheatmap(
  mat               = matexz,
  cluster_rows      = rlord,
  cluster_cols      = clord,
  color             = mypalette(length(palettebreaks)-1),
  breaks            = palettebreaks,
  scale             = 'none',
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  main              = "Cluster-specific genes",
  fontsize          = 10,
  fontsize_col      = 5,
  fontsize_number   = 5,
  annotation_col    = matcolms,
  annotation_row    = matrows,
  annotation_colors = matcouls,
  annotation_legend = T,
  annotation_names_col = F,
  annotation_names_row = F,
  drop_levels       = TRUE,
  # gaps_col          = 1:(length(gorder)-1), # gaps per type
  # gaps_row          = grgenessep,
  filename          = paste0('heatmaps/', mysufix, '.pdf'),
  width = 12, height = 12
)
try(dev.off(), silent = T)
ogenes <- x$tree_row$labels[x$tree_row$order]
write.csv(ogenes, file = paste0('heatmaps/', mysufix, '.csv'))

### 5.X Coexpression TH2 AS_AL/AR ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thisdata <- matex
datatype <- "log2CPM"

# savey <- theObjectSavedIn(paste0(root, '/../raw/10X/AS3Esm/AGGR/AS3Esm/outs/imputed/saver_matrix.RData'))
# savey <- as.matrix(savey)
# savey <- sweep(savey, 2, colSums(savey), '/') * 1e6
# saveymat <- log2(savey + 1)
# # saver_obj <- theObjectSavedIn(paste0(root, '../../raw/10X/AS3Esm/AGGR/AS3Esm/outs/imputed/combined.RData'))
# # saver_obj_ss <- list(estimate = saver_obj[[1]][tmp, mysamples], se = saver_obj[[2]][mygenes, mysamples])
# # saver_cor_gene <- cor.genes(saver_obj)
thisdata <- saveymat
datatype <- 'saver'

mysufix <- paste0(dg, collapse = "n")
fname <- paste0('heatmaps/correlation_', newlabs[gr], '_', mysufix, "_", datatype)
pdf(paste0(fname, '.pdf'), width = 10, height = 10)
void <- geneCorr(
  expdata = thisdata,
  cortype = 'spearman',
  brks = 'quantile',
  mmain = paste('Correlation:\nCluster', commas(dg), 'cells'),
  sublab = paste0('Genes: ', commas(dg), '-specific'),
  keepgenes = genes,
  samples = thesecells,
  # cols = showgenes,
  verbose = TRUE,
  rank_genes = T
)
dev.off()
write.csv(void, file = paste0(fname, '.csv'))

# cluname <- paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/th2_0.1/clustering/zetInfo/clustCells6PCs_30Ks_0.06667JD.RData')
# cluname <- paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/th2p_0.1_pctfilt_first/clustering/zetInfo/clustCells15PCs_30Ks_0.06667JD.RData')
# cluname <- paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/th2pv_0.1/clustering/zetInfo/clustCells5PCs_30Ks_0.06667JD.RData')
cluname <- paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/th2pv_0.1/clustering/zetInfo/clustCells14PCs_30Ks_0.06667JD.RData')
mycells_th2 <- theObjectSavedIn(cluname)
mycells_th2@reductions$tsne <- mycells_th2@reductions$tsne50
mycells_th2@meta.data$RNA_snn_res.0.3 <- paste0("TH2_", mycells_th2@meta.data$RNA_snn_res.0.2)

### 5.X TH2 t-SNE ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xcols <- v2cols(levels(factor(mycells_th2@meta.data$Cell_type)), gr.cols)
pdf(paste0(basename(dirnames(cluname, 3)), "_R0.2_tSNE.pdf"))
DimPlot(mycells_th2, cols = xcols, group.by = "Cell_type")
dev.off()

### 5.X TH2 subcluster disease composition ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resolut <- 'Cell_type'
annot <- mycells_th2@meta.data[getsubset(c("orig.diseasegroup", "AS_AL", "AR"), mycells_th2@meta.data, v = TRUE), ]
table(annot$orig.diseasegroup)
thesecells <- sample_grp(annot, cname = "orig.diseasegroup")
p <- get_props(metadata = annot[thesecells, ], group_by = 'orig.diseasegroup', resolution = resolut,
  couls = gr.cols, reverseit = TRUE, ssize = c(3, 1, 5), v = TRUE)
dirnames <- function(x, n = 1){ for(i in 1:n) x <- dirname(x); return(x) }
fname <- paste0(basename(dirnames(cluname, 3)), "_TH2_subclustering_sampled")
pdf(paste0("clusters_pies_", fname, ".pdf"))
print(p$pies);
dev.off()
pdf(paste0("clusters_pies_", fname, "blank.pdf"), height = 7, width = 7)
print(shut_it(p$pies));
dev.off()
grps <- rev(make_list(mycells_th2@meta.data[thesecells, ], colname = "orig.diseasegroup"))
xcols <- v2cols(levels(factor(mycells_th2@meta.data$Cell_type)), gr.cols)
p <- lapply(grps, function(x) DimPlot(mycells_th2, pt.size= 2, cells = x, group.by = resolut, cols = xcols) )
pdf(paste0("clusters_tsne_", fname, ".pdf"), width = 12)
plot_grid(plotlist = p, ncol = 2)
dev.off()
p <- lapply(grps, function(x) shut_it(DimPlot(mycells_th2, pt.size= 2,cells = x, group.by = resolut, cols = xcols)) )
pdf(paste0("clusters_tsne_", fname, "blank.pdf"), width = 12)
plot_grid(plotlist = p, ncol = 2)
dev.off()

p <- get_props(metadata = annot[thesecells, ], group_by = resolut, resolution = 'orig.diseasegroup',
  couls = gr.cols, reverseit = TRUE, ssize = c(3, 1, 5), v = TRUE)
dirnames <- function(x, n = 1){ for(i in 1:n) x <- dirname(x); return(x) }
fname <- paste0(basename(dirnames(cluname, 3)), "_TH2_subclustering_sampled")
pdf(paste0("clusters_pies_", sub("TH2", "disease", fname), ".pdf"))
print(p$pies);
dev.off()
pdf(paste0("clusters_pies_", sub("TH2", "disease", fname), "blank.pdf"), height = 7, width = 7)
print(shut_it(p$pies));
dev.off()

### 5.X Proportion of cluster per patient per disease in TH2 subclustering ###%%
annot$diseaseclust <- paste0(annot$orig.diseasegroup, annot$orig.donor)
for(tmp in c("orig.donor", "orig.diseasegroup", "diseaseclust")[1]){
  tvar <- c('RNA_snn_res.0.3', tmp)
  if(1) tvar <- rev(tvar)
  fname <- paste0(basename(dirnames(cluname, 3)), "_",
    sub("orig.|RNA_snn_", "", tvar[1]), "_vs_", sub("orig.|RNA_snn_", "", tvar[2]))
  p <- get_props(metadata = annot[thesecells, ], group_by = tvar[1], resolution = tvar[2],
    couls = gr.cols, reverseit = TRUE, ssize = c(3, 1, 5), v = TRUE)
  pdf(paste0(fname, ".pdf"))
  print(p$pies);
  dev.off()
}

### 5.X IL9 subcluster genes ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Scatter plots of IL9 subcluster genes with Th2 AS-AL cells
fname <- paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/th2_0.1/markers/6PCs_RNA_snn_res.0.4_MAST/result_MAST_DEGsTAT_LFC0.25_QVAL0.05.csv')
res <- readfile(fname, stringsAsFactors = FALSE, check.names = FALSE)
head(res)
genes_tab <- res[res$cluster == "1", ]
rownames(genes_tab) <- sub("'", "", genes_tab$gene_name)
genes <- getDEGenes(genes_tab, pv = 0.05, fc = 0.4, pvtype = "p_val_adj", lfc.type = "avg_logFC", v = TRUE)

dir.create("TH2_subclust_scatters")
tsne <- Embeddings(mycells_th2, reduction = "tsne")
xmax <- tsne[bordering(tsne), ]
for(gg in genes){
  cat(gg, "\n")
  p <- FeaturePlot(mycells_th2, features = gg, cells = getsubset(c("orig.diseasegroup", "AS_AL"), mycells_th2@meta.data, v = TRUE))
  p <- p + xlim(c(min(tsne[, 1]), max(tsne[, 1]))) + ylim(c(min(tsne[, 2]), max(tsne[, 2])))
  pdf(paste0("TH2_subclust_scatters/AS_AL_", gg,".pdf")); print(p); dev.off()
  p <- FeaturePlot(mycells_th2, features = gg, cells = getsubset(c("orig.diseasegroup", "AR"), mycells_th2@meta.data, v = TRUE))
  p <- p + xlim(c(min(tsne[, 1]), max(tsne[, 1]))) + ylim(c(min(tsne[, 2]), max(tsne[, 2])))
  pdf(paste0("TH2_subclust_scatters/AR_", gg,".pdf")); print(p); dev.off()
}

### 5.X ILs comparison within AS_AL ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- mycells@meta.data
expdata <- cts2cpm(mycells@assays$RNA@counts)
void <- add_gene_tag(lgenes = c("IL9", "IL3", "IL5"), annot = annot, mat = expdata, tag = c('tag', 'P', 'N'), v = TRUE)
annot <- cbind(annot, void[rownames(annot), ])
thisgene <- paste0("IL", 9)
compid <- paste0(thisgene, "Pvs", thisgene, "N")
wgrcol <- "tag_IL9"
gr <- "4"
res <- readfile(paste0(root, '/mast/AS3EsmTeff_daf_30p/comprs/TH2_IL9_AS_AL/', compid, '/results_', compid, '_mastlog2cpm.csv'), stringsAsFactors = FALSE, check.names = FALSE)
headmat(res)
res <- res[res$gene != thisgene, ]
headmat(res)
wgr <- unlist(strsplit(compid, "vs"))
thesecells <- getsubset(list(c("Cell_type", newlabs[gr]), c("orig.diseasegroup", "AS_AL"), c(wgrcol, wgr)), annot, v = TRUE)
suffy <- "_AS_AL"

### 5.X ILs comparison within 0 ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- mycells_th2@meta.data
expdata <- cts2cpm(mycells_th2@assays$RNA@counts)
void <- add_gene_tag(lgenes = c("IL9", "IL3", "IL5"), annot = annot, mat = expdata, tag = c('tag', 'p', 'n'), v = TRUE)
annot <- cbind(annot, void[rownames(annot), ])
thisgene <- paste0("IL", 9)
compid <- paste0(thisgene, "pvs", thisgene, "n")
wgrcol <- "tag_IL9"
gr <- "4"
res <- readfile(paste0(root, '/mast/th2pv_0.1/comprs/within0/', compid, '/results_', compid, '_mastlog2cpm.csv'), stringsAsFactors = FALSE, check.names = FALSE)
headmat(res)
res <- res[res$gene != thisgene, ]
headmat(res)
wgr <- unlist(strsplit(compid, "vs"))
thesecells <- getsubset(list(c("Cell_type", "TH2_0"), c("orig.diseasegroup", "AS_AL", "AR")), annot, v = TRUE)
suffy <- "_TH0nly"

### 5.X AS_AL vs AR comparison within 0 ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- mycells_th2@meta.data
expdata <- cts2cpm(mycells_th2@assays$RNA@counts)
compid <- "AS_ALvsAR"
wgrcol <- "orig.diseasegroup"
gr <- "4"
res <- readfile(paste0(root, '/mast/th2pv_0.1/comprs/within0/', compid, '/results_', compid, '_mastlog2cpm.csv'), stringsAsFactors = FALSE, check.names = FALSE)
headmat(res)
wgr <- unlist(strsplit(compid, "vs"))
thesecells <- getsubset(list(c("Cell_type", "TH2_0"), c("orig.diseasegroup", "AS_AL", "AR")), annot, v = TRUE)
suffy <- "_within0"

### 5.X ILs comparison ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- mycells@meta.data
expdata <- cts2cpm(mycells@assays$RNA@counts)
void <- add_gene_tag(lgenes = c("IL9", "IL3", "IL5"), annot = annot, mat = expdata, v = TRUE)
annot <- cbind(annot, void[rownames(annot), ])
thisgene <- paste0("IL", 9)
compid <- paste0(thisgene, "+vs", thisgene, "-")
fsufix <- paste0("", thisgene)
wgrcol <- "tag_IL9"
gr <- "4"
res <- readfile(paste0(root, '/mast/AS3EsmTeff_daf_30p/comprs/per_gene/', compid, '/results_', compid, '_mastlog2cpm.csv'), stringsAsFactors = FALSE, check.names = FALSE)
headmat(res)
res <- res[res$gene != thisgene, ]
headmat(res)
wgr <- unlist(strsplit(compid, "vs"))
thesecells <- getsubset(list(c("Cell_type", newlabs[gr]), c(wgrcol, wgr)), annot, v = TRUE)
suffy <- ""
pdf(paste0('volcano_', fsufix, '.pdf'), width = 10, height = 10)
void <- try(volplot(res,
  pvalth = 0.05,
  lfcth = 0.25,
  pvaltype = 'padj',
  lfctype = 'log2FoldChange',
  gene_name = 'gene',
  interact = FALSE,
  ngenes = 10,
  legends = TRUE,
  titl = paste("Volcano Plot:", sub("vs", " vs ", compid)),
  subt = NULL,
  # check_genes = markers,
  gname = "Markers",
  only_degs = FALSE,
  return_plot = FALSE,
  clipp = ifelse(thisgene == "IL3", 2, FALSE),
  v = TRUE
))
dev.off()

### donors tSNEs ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sobject <- mycells
dset <- "donors_tsnes"
sobject <- mycells_th2[, getsubset(c("orig.diseasegroup", "AS_AL", "AR"), mycells_th2@meta.data, v = TRUE)]
dset <- "donors_tsnes_th2"

diseaseorder <- c('AS_AL', 'AR', 'AS_NA', 'HC')
axmax <- max(abs(Embeddings(sobject, reduction = "tsne")))
nsmall <- min(table(sobject@meta.data$orig.diseasegroup))
sobject@meta.data$ddisease <- paste(sobject@meta.data$orig.donor, "-", sobject@meta.data$orig.diseasegroup)
idorder <- paste0("D", 1:24)
names(idorder) <- c("1014", "1020", "1196", "1197", "1437", "1441", "1016", "1178", "1362",
  "1442", "1447", "1663", "1175", "1176", "1179", "1198", "1448", "1838", "1562",
  "1565", "2103", "2133", "2184", "2194")
sobject@meta.data$ddisease <- idorder[sobject@meta.data$orig.donor]
sampleit <- TRUE
compileplots <- list()
dir.create(dset)
for(gr in diseaseorder[diseaseorder %in% sobject@meta.data$orig.diseasegroup]){
  thesecells <- getsubset(c("orig.diseasegroup", gr), sobject@meta.data, v = T)
  fname <- gr
  if(isTRUE(sampleit)){
    fname <- paste0("sampled_", gr)
    set.seed(27); thesecells <- sample(thesecells, nsmall)
  }
  ddf <- FetchData(sobject, vars = c("tSNE_1", "tSNE_2", "Cell_type", "ddisease"))
  ddf <- ddf[thesecells, ]
  p <- ggplot(ddf, aes(x = tSNE_1, y = tSNE_2, colour = Cell_type)) + geom_point() +
    facet_wrap(~ ddisease, ncol = 2) + scale_colour_manual(values = v2cols(ddf$Cell_type, gr.cols)) +
    xlim(c(-axmax, axmax)) + ylim(c(-axmax, axmax)) + drawsq + theme_classic() +
    theme(legend.position = "none", strip.text.x = element_text(face = 'bold', size = 23, hjust = 0),
      strip.background = element_rect(fill = "#FFFFFF", linetype = 0),
      axis.text = element_blank(), axis.title = element_blank())
  compileplots[[fname]] <- p
  pdf(paste0(dset, "/", fname, ".pdf"), height = 12, width = 9)
  print(p);
  dev.off()
  pdf(paste0(dset, "/", fname, "blank.pdf"), height = 12, width = 9)
  print(p + shut_up)
  dev.off()
  compileplots[[paste0(fname, "blank")]] <- p + shut_up
}

### TH2 vs IL9 numbers ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thesecells <- getsubset(c("orig.diseasegroup", "AS_AL", "AR"), mycells_th2@meta.data, v = TRUE)
annot <- mycells_th2@meta.data[thesecells, ]
annot <- mycells_th2@meta.data[mycells_th2@meta.data$orig.diseasegroup %in% c("AS_AL", "AR"), ]
thesecells <- rownames(annot)
expdata <- cts2cpm(mycells_th2@assays$RNA@counts[, thesecells])
void <- add_gene_tag(lgenes = "IL9", annot = annot, mat = expdata, v = TRUE)
annot <- cbind(annot, void[rownames(annot), , drop = FALSE])
prop.table(table(annot[, c("tag_IL9", "Cell_type")]), 2)
prop.table(table(annot[, c("tag_IL9", "orig.diseasegroup")]), 2)
annot$tmp <- paste0(annot$Cell_type, "_", annot$orig.diseasegroup)
mytab <- melt(table(annot$tag_IL9, annot$orig.donor, annot$orig.diseasegroup, annot$Cell_type))
mytab_t <- dcast.data.table(setDT(mytab), Var2 + Var3 + Var4 ~ Var1)
write.csv(mytab_t, file = "IL9_groups_proportions.csv")

### TH2 vs IL9 numbers ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expdata <- cts2cpm(mycells_th2@assays$RNA@counts)
ddf <- mycells_th2@reductions$tsne50@cell.embeddings[, c("tSNE_1", "tSNE_2")]
ddf <- cbind(ddf, mycells_th2@meta.data[rownames(ddf), "Cell_type", drop = FALSE])
ddf <- cbind(ddf, log2(t(expdata["IL9", rownames(ddf), drop = FALSE]) + 1))
ddf[ddf$IL9 < 1, "IL9"] <- NA
summary(ddf[ddf$IL9 < 1, "IL9"])
headmat(ddf)
couls <- c("#ffc100","#ff7400", "#ff0000","#a10000", "#670000")
brks <- round(c(min(ddf$IL9, na.rm = TRUE), 10, 12.5, 15, max(ddf$IL9, na.rm = TRUE)), 1)
brks <- c(8, 10, 12.5, 15, 17)
pdf("IL9coloured_fixed.pdf")
ggplot(ddf, aes(x = tSNE_1, y = tSNE_2, color = IL9)) + geom_point() +
  scale_colour_gradientn(colors = couls, na.value = "#bebebe", breaks = brks)
dev.off()
pdf("IL9coloured_ellipse_fixed.pdf")
ggplot(ddf, aes(x = tSNE_1, y = tSNE_2, color = IL9)) + geom_point() +
  scale_colour_gradientn(colors = couls, na.value = "#bebebe", breaks = brks) +
  stat_ellipse(aes(group = Cell_type))
dev.off()

### Chilli plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lgenes <- list(
  TH2_0_IL9P = c("IL9", "IL1RL1", "IL5", "DUSP6", "CD109", "GZMB", "ZEB2", "EFHD2", "SEC61G",
  "SEC61B", "PHLDA1", "FKBP1A", "FKBP2", "MAF")
)
lgenes <- list(
  TH2_0_IL9P = c("IL9", "IL1RL1", "IL5", "DUSP6", "CD109", "GZMB", "ZEB2", "EFHD2",
  "SEC61B", "PHLDA1", "FKBP2", "MAF")
)
mergegroups <- NULL
annot <- remove.factors(mycells_th2@meta.data)
thesecells <- getsubset(c("Cell_type", "TH2_0"), annot, v = TRUE)
expdata <- cts2cpm(mycells_th2@assays$RNA@counts)[, thesecells]
annot <- annot[thesecells, ]
void <- add_gene_tag(lgenes = c("IL9", "IL3", "IL5"), annot = annot, mat = expdata, tag = c('tag', 'p', 'n'), v = TRUE)
annot <- cbind(annot, void[rownames(annot), ])
table(annot$Cell_type, annot$tag_IL9)
annot$Cell_type <- annot$tag_IL9
get_stat_report(expdata, groups = make_list(annot, colname= "Cell_type", grouping = T), rnames = lgenes[[1]], v = T)
sum(expdata["IL9", getsubset(c("Cell_type", "IL9p"), annot, v = TRUE)] > 0)
brkies <- list(filname = '%+cells', brks = c(1, 10, 20, 30, 40, 50))
# colies <- c("#fffeee", "#ffc100","#ff7400", "#ff0000","#a10000", "#670000")
colies <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
brkies <- NULL
nc <- 2
extramods <- coord_cartesian(ylim = c(5, 20)) #+ theme(panel.spacing = unit(2, "lines"))
# Use fig2.R at 2.D Chilli plot
