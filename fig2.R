#!/usr/bin/R5

############
# Figure 2 #
############

# This script will create plots for figure 2 from single-cell data

library(Seurat)
source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')

## Destiny Folder ##
root <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma/results'
setwd(paste0(root, '/a1_final_figures/Figure_2'))

mycells <- theObjectSavedIn(paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/clustering/zetInfo/clustCells13PCs_30Ks_0.06667JD.RData'))
gr.cols <- read.csv("/mnt/BioHome/ciro/asthma/info/AS3_final_groupColours.csv",row.names=1,stringsAsFactors=F,header=1)

newlabs <- c("ACT1", "ACT2", "ACT3", "TH1", "TH2", "THIFNR", "TH17")
names(newlabs) <- 0:6
gorder <- c("TH2", "TH1", "TH17", "THIFNR", "ACT3", "ACT2", "ACT1")
mycells@meta.data$Cell_type <- factor(c(newlabs[as.character(mycells@meta.data$RNA_snn_res.0.4)], gorder)

fname <- paste0(root, '/mast/AS3EsmTeff_daf_30p/summary/clusters/padj0.05_FC0/RNA_snn_res.0.4_SummaryDEGsTable_suas.csv')
res <- readfile(fname, stringsAsFactors = FALSE, check.name = FALSE)
res <- res[!is.na(res$group), ]
# modifying the label names for the clusters
void <- sapply(1:length(newlabs), function(x){
  y <- newlabs[x]
  colnames(res) <<- gsub(paste0("^", names(y)), y, colnames(res))
  colnames(res) <<- gsub(paste0("\\(", names(y), "\\)"), paste0("(", y, ")"), colnames(res))
  res$group <<- gsub(paste0("^", names(y), "$"), y, res$group)
  res$group <<- gsub(paste0("^", names(y), "n"), paste0(y, "n"), res$group)
  res$group <<- gsub(paste0("n", names(y), "$"), paste0("n", y), res$group)
  NULL
})

### 1.A Volcano plot ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### 2.A Density plot ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mycellsa <- theObjectSavedIn(paste0(root, '/clustering_seurat/AS3Esm_df/clustering/zetInfo/clustCells20PCs_30Ks_0.06667JD.RData'))
fname <- paste0(root, '/mast/AS3EsmTeff_sng_amb/comprs/celltype/TeffvsTNeg/activatedTeff_110g/activatedTeff_110g_signature.csv')
sig_tab <- readfile(fname, row.names = 1)
samples <- c(rownames(mycellsa@meta.data[mycellsa@meta.data$orig.celltype == "TNeg", ]), rownames(mycellsa@meta.data))
df <- sig_tab[rownames(sig_tab) %in% samples, , drop = FALSE]
df$celltype <- mycellsa@meta.data[rownames(df), "orig.celltype"]
p <- ggplot(df, aes_string(x = "activatedTeff_110g1", fill = "celltype")) +
  geom_density(aes(y = ..scaled..), alpha = 0.5) +
  scale_fill_manual(values = v2cols(unique(df[, "celltype"]), gr.cols)) +
  labs(title = "Activation score - Teff vs TNeg", x = "Level of activation", y = "Density")
pdf("density.pdf", width = 12, height = 8)
p
dev.off()
pdf("densityblank.pdf", width = 12, height = 8)
p + shut_up
dev.off())
pdf("densitylineblank.pdf", width = 12, height = 8)
p + shut_up + geom_vline(xintercept = -0.273)
dev.off()

### 2.B Heatmap plot Teff vs TNeg signature ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname <- paste0(root, '/mast/AS3EsmTeff_sng_amb/comprs/celltype/TeffvsTNeg/activatedTeff_110g/used_genes.csv')
system(paste("cp", fname, "./activation_genes.csv"))
genes <- readfile(fname, row.names = 1, stringsAsFactor = FALSE)[, 1]
expdata <- cts2cpm(mycellsa@assays$RNA@counts[, rownames(df)])
df <- df[order(df[, 1]), ]
head(df)
samples <- unlist(sapply(unique(df$celltype), function(x){
  tdf <- df[df$celltype == x, ]
  rownames(tdf[order(tdf[, 1]), ])
}))

matcolms <- data.frame(Group = df[samples, ]$celltype)
rownames(matcolms) <- colnames(matex) # column names and type
matcouls <- list(Group = v2cols(unique(df$celltype), gr.cols, v = TRUE)) # colours for type
matex <- expdata[genes, samples]

coulrange = c('blue', 'black', 'yellow')
palettebreaks <- seq(-2, 2, 0.1) # Heatmap setting
mypalette <- colorRampPalette(coulrange, space = 'Lab')

# manually scale
matexz <- t(scale(t(matex)))
matexz[matexz < (-2)] <- -2
matexz[matexz > (2)] <- 2
library(pheatmap)
x <- pheatmap(
  mat               = matexz,
  cluster_rows      = 0,
  cluster_cols      = F,
  color             = mypalette(length(palettebreaks)-1),
  breaks            = palettebreaks,
  scale             = 'none',
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  main              = "Cluster-specific genes",
  fontsize          = 10,
  fontsize_col      = 5,
  fontsize_number   = 5,
  annotation_col    = matcolms,
  # annotation_row    = matrows,
  annotation_colors = matcouls,
  annotation_legend = T,
  annotation_names_col = F,
  annotation_names_row = F,
  drop_levels       = TRUE,
  # gaps_col          = 1:(length(gorder)-1), # gaps per type
  # gaps_row          = grgenessep,
  filename          = paste0('heatmap_teff_tneg_activation.pdf'),
  width = 12, height = 12
)
try(dev.off(), silent = T)

## Copying gene list
fname <- paste0(root, '/mast/AS3EsmTeff_sng_amb/comprs/celltype/TeffvsTNeg/activatedTeff_110g/used_genes.csv')

### 2.X Heatmap plot cluster-specific ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# data from single-cell
coldata_ <- mycells@meta.data
coldata_$Cell_type <- as.character(coldata_$Cell_type)
coldata_$Disease <- as.character(coldata_$orig.diseasegroup)
# samples <- rownames(coldata_[with(coldata_, order(Cell_type, Disease)), ])
samples <- orderthis(annot = coldata_, order_by = 'Cell_type', grupos = gorder, lgenes = 'lol', v = TRUE)
unique(res$group)
rownames(res) <- gsub("'", "", res$gene_name)
genes <- gsub("'", "", res[!is.na(res$group), ]$gene_name)

# Get necessary variables
groups_col = v2cols(gorder, gr.cols)
gr.genes.col = unique(res$colour_key)
names(gr.genes.col) <- unique(res$group)
gr.genes.col[gorder] <- groups_col
gr.genes.size = table(res$group)
gr.genes.size <- c(gr.genes.size[gorder], rev(sort(gr.genes.size[!names(gr.genes.size) %in% gorder])))
geneorder <- orderthis(annot = res, order_by = 'group', grupos = names(gr.genes.size), lgenes = 'lol', v = TRUE)
res <- res[geneorder, ]
coulrange = c('blue', 'black', 'yellow')

# plotting settings
matcouls <- list(Specific = gr.genes.col, Group = groups_col) # colours for type
grgenessep <- cumsum(gr.genes.size) # to separate gene groups
palettebreaks <- seq(-2, 2, 0.1) # Heatmap setting
mypalette <- colorRampPalette(coulrange, space = 'Lab')
matrows <- data.frame(Specific = rep(names(gr.genes.size), gr.genes.size))
rownames(matrows) <- gsub("'", "", res$gene_name)

# Plot means per group
matex <- res[, paste0(gorder, "_mean")]
colnames(matex) <- sub("_mean", "", colnames(matex))
matcolms <- data.frame(Group = gorder)
rownames(matcolms) <- colnames(matex) # column names and type

# manually scale
matexz <- t(scale(t(matex)))
matexz[matexz < (-2)] <- -2
matexz[matexz > (2)] <- 2
library(pheatmap)
x <- pheatmap(
  mat               = matexz,
  cluster_rows      = 0,
  cluster_cols      = F,
  color             = mypalette(length(palettebreaks)-1),
  breaks            = palettebreaks,
  scale             = 'none',
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
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
  gaps_col          = 1:(length(gorder)-1), # gaps per type
  gaps_row          = grgenessep,
  filename          = paste0('heatmap_byGroups_means_new.pdf'),
  width = 12, height = 12
)
try(dev.off(), silent = T)

### 1.X Heatmap plot cluster-specific bottom ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library("UpSetR")
digr <- list()
void <- apply(res, 1, function(x){
  cat('.')
  y <- unlist(strsplit(x[2], "n"))
  for(gr in y){
    digr[[gr]] <<- c(digr[[gr]], x[3])
  }
}); cat("\n")
digr <- lapply(digr, unique)

# bin_tab <- data.table::melt(digr)
# bin_tab <- table(bin_tab[, 1], bin_tab[, 2])
# bin_tab <- as.data.frame.matrix(bin_tab, stringsAsFactors = F, check.names = F)
# head(bin_tab)
# # bin_tab <- bin_tab[rowSums(bin_tab) == 1, gorder]
# bin_tab$n <- sample(1:nrow(bin_tab))
# head(bin_tab)
# myqueries <- lapply(head(gorder, 7), function(x){
#   list(query = intersects, params = list(x), active = T, color = groups_col[x])
# })
# pdf("degs_bars_comb.pdf", onefile = F)
# upset(bin_tab,
#   sets = gorder,
#   sets.bar.color = gr.cols[gorder, ],
#   keep.order = TRUE,
#   point.size = 3,
#   queries = myqueries
# )
# dev.off()
digr2 <- overlap_calc(groups_list = digr, sharedmax = 1)
library("ComplexHeatmap")
m <- make_comb_mat(digr)
thiscols <- c(groups_col, rep("black", 12))
# UpSet(m, set_order = gorder, comb_col = thiscols, pt_size = unit(5, "mm"))
ss = set_size(m)
myorder <- order(comb_degree(m), -comb_size(m))
myorder <- c(1:7, myorder[-c(1:7)][1:3])
m <- m[, myorder]
myorder <- c(1:10)
pdf("degs_bars_reordered_comb.pdf", onefile = F, height = 7, width = 10)
UpSet(m, set_order = gorder, comb_col = thiscols, pt_size = unit(5, "mm"),
  comb_order = myorder,
  top_annotation = HeatmapAnnotation(
    "Intersection Size" = anno_barplot(
      comb_size(m),
      border = FALSE,
      gp = gpar(fill = thiscols),
      height = unit(10, "cm")
    ),
    annotation_name_side = "left",
    annotation_name_rot = 90
  ),
  left_annotation = rowAnnotation(
    "Cluster Size" = anno_barplot(-ss,
      baseline = 0,
      axis_param = list(
        at = c(0, -100, -200, -300),
        labels = c(0, 100, 200, 300),
        labels_rot = 0),
      border = FALSE,
      gp = gpar(fill = thiscols[1:7]),
      width = unit(4, "cm")
    ),
    set_name = anno_text(
      set_name(m),
      location = 0.5,
      just = "center",
      width = max_text_width(set_name(m)) + unit(4, "mm"))
    ),
    right_annotation = NULL,
    show_row_names = FALSE
)
dev.off()

### 1.X Heatmap plot cluster-specific sheets ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Adding the rest of the genes
expdata <- cts2cpm(mycells@assays$RNA@counts)
dim(res)
dim(expdata)
genes2add <- rownames(expdata)[!rownames(expdata) %in% gsub("'", "", res$gene_name) ]
resrest <- data.frame(mat_names(genes2add, colnames(res)), stringsAsFactors = FALSE, check.names = FALSE)
resrest$gene_name <- paste0("'", genes2add)
# samples <- colnames(resrest)#[-c(1:33)]
# all(samples %in% rownames(coldata_)) # are all samples in the meta data?
group_list <- make_list(x = coldata_, colname = "Cell_type", grouping = TRUE)
stattab <- get_stat_report(mat = expdata, groups = group_list, rnames = genes2add, v = TRUE)
# Proving if calculations are consistent
stattabdone <- get_stat_report(mat = expdata, groups = group_list, rnames = gsub("'", "", res$gene_name), v = TRUE)
head(stattabdone[gsub("'", "", res$gene_name), ])
head(res[, colnames(stattabdone)])
takencols <- colnames(resrest)[colnames(resrest) %in% colnames(stattab)]
resrest[, takencols] <- stattab[rownames(resrest), takencols]
# resrest[, samples] <- expdata[rownames(resrest), samples]
res_all_genes <- rbind(res, resrest)[, -1]
headmat(res_all_genes)
write.csv(res_all_genes, file = "heatmapp_celltype_genes_all_info_sc.csv")
write.csv(res_all_genes[, c("group", takencols)], file = "heatmapp_celltype_genes_only_stats_and_group_sc.csv")
sgenes <- rownames(expdata[gsub("'", "", res_all_genes[res_all_genes$Bmean > 1, ]$gene_name), ])
qlucoref <- qlucore_format(mat = expdata, metadata = coldata_[samples, ], rnames = sgenes, v = TRUE)
write.csv(qlucoref, file = "qlucore_file_celltype_genes_mean_gt1.csv")

dir.create("clusters_tables")
coldata_ <- mycells@meta.data
## expression tables per cell type cluster
for(gr in as.character(unique(coldata_$Cell_type))[-7]){
  cat(gr, "\n")
  genes <- gsub("'", "", res[which(res$group == gr), ]$gene_name)
  samples <- getsubset(c("Cell_type", gr), coldata_, v = TRUE)
  samples <- names(head(sort(colSums(expdata[genes, samples] > 0), decreasing = TRUE), 1000))
  write.csv(expdata[genes, samples], file = paste0("clusters_tables/", gr, "_specific_genes.csv"))
  newmat <- expdata[rowMeans(expdata[, samples]) > 1, samples]
  write.csv(newmat, file = paste0("clusters_tables/", gr, "_genes_mean_gt1.csv"))
  qlucoref <- qlucore_format(mat = newmat, metadata = coldata_, cnames = samples, v = TRUE)
  write.csv(qlucoref, file = paste0("clusters_tables/", gr, "_genes_mean_gt1_qlucore.csv"))
}
# Now for TNegs
mycellsa <- theObjectSavedIn(paste0(root, '/clustering_seurat/AS3Esm_df/clustering/zetInfo/clustCells20PCs_30Ks_0.06667JD.RData'))
gr <- "TNeg"
cat(gr, "\n")
fname <- paste0(root, '/mast/AS3EsmTeff_sng_amb/comprs/celltype/TeffvsTNeg/results_TeffvsTNeg_mastlog2cpm.csv')
resy <- read.csv(fname, stringsAsFactors = FALSE)
genes <- gsub("'", "", resy[which(resy$group == gr), ]$gene_name)
genes <- genes[genes %in% rownames(expdatatneg)]
fname <- paste0(root, '/mast/AS3EsmTeff_sng_amb/comprs/celltype/TeffvsTNeg/activatedTeff_110g/non_activated_tnegs.csv')
samples <- read.csv(fname, stringsAsFactors = FALSE)[, 2]
expdatatneg <- cts2cpm(mycellsa@assays$RNA@counts[, samples])
samples <- names(head(sort(colSums(expdatatneg[genes, ] > 0), decreasing = TRUE), 1000))
write.csv(expdatatneg[genes, samples], file = paste0("clusters_tables/", gr, "_specific_genes.csv"))
newmat <- expdatatneg[rowMeans(expdatatneg[, samples]) > 1, samples]
write.csv(newmat, file = paste0("clusters_tables/", gr, "_genes_mean_gt1.csv"))
qlucoref <- qlucore_format(mat = newmat, metadata = mycellsa@meta.data, cnames = samples, v = TRUE)
write.csv(qlucoref, file = paste0("clusters_tables/", gr, "_genes_mean_gt1_qlucore.csv"))
annot <- mycellsa@meta.data[samples, ]
table(annot$seurat_clusters)

## mega table
fnames <- list.files("clusters_tables", pattern = "_specific_genes", full.name = TRUE)
fnames <- fnames[!grepl("all_clusters_top_cells", fnames)]
gr <- "all_clusters_top_cells"
cellnames <- sapply(fnames, function(x){
  colnames(readfile(x, check.names = FALSE, row.names = 1))
})
sapply(cellnames, length)
samples <- unlist(cellnames)
head(samples)
names(samples) <- NULL
head(samples)
expdataall <- as.matrix(mycellsa@assays$RNA@counts[, samples])
expdatatneg <- cts2cpm(expdataall)
genes <- gsub("'", "", res$gene_names)
write.csv(expdatatneg[genes, samples], file = paste0("clusters_tables/", gr, "_specific_genes_fixed.csv"))
newmat <- expdatatneg[rowMeans(expdatatneg[, samples]) > 1, samples]
write.csv(newmat, file = paste0("clusters_tables/", gr, "_genes_mean_gt1_fixed.csv"))
# conames <- colnames(mycells@meta.data)
# conames[conames %in% colnames(mycellsa@meta.data)]
annot <- remove.factors(data.frame(data.table::rbindlist(list(
  mycells@meta.data[unlist(cellnames[-8]), ],
  mycellsa@meta.data[cellnames[[8]], ]), fill = TRUE
), check.names = FALSE))
rownames(annot) <- annot$cell_id
annot$Cell_type[is.na(annot$Cell_type)] <- "TNeg"
sapply(annot, function(x) any(is.na(x)) )
annot[is.na(annot)] <- "TNeg"
table(annot$Cell_type)
table(annot$orig.r88gnb_clusters)
annot <- annot[, -c(10:30)]
annot <- annot[, !colnames(annot) %in% c("RNA_snn_res.0.6", "seurat_clusters")]
qlucoref <- qlucore_format(mat = newmat, metadata = annot, cnames = samples, v = TRUE)
write.csv(qlucoref, file = paste0("clusters_tables/", gr, "_genes_mean_gt1_qlucore_fixed.csv"))
write.csv(samples, file = paste0("clusters_tables/", gr, "_genes_mean_gt1_qlucore_fixed_samples.csv"))

### 2.C t-SNE plot ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p <- DimPlot(mycells, reduction = "tsne75", pt.size = 1.5,
  group.by = "Cell_type", cols = v2cols(levels(mycells@meta.data$Cell_type), gr.cols)) +
  drawsq
# p <- squareplot(p)
pdf("tSNE.pdf", width = 10, height = 9)
p
dev.off()
pdf("tSNEblank.pdf", width = 10, height = 9)
p + shut_up
dev.off()
metdata <- FetchData(mycells, var = c("tSNE_1", "tSNE_2", "Cell_type"))
p <- ggplot(metdata, aes(x = tSNE_1, y = tSNE_2, colour = Cell_type)) +
  geom_density_2d(h = c(2, 2)) +
  scale_colour_manual(values = v2cols(as.character(metdata$Cell_type), gr.cols))
pdf("tSNE_contour.pdf", width = 10, height = 9)
p
graphics.off()
pdf("tSNE_contourblank.pdf", width = 10, height = 9)
p + shut_up
dev.off()

### 2.X histogram and pie chart ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- mycells@meta.data
annot$Cell_type <- factor(annot$Cell_type, gorder)
p <- get_props(metadata = annot, group_by = 'Cell_type', resolution = 'orig.ident',
  couls = gr.cols, reverseit = TRUE, v = TRUE)
pdf("pie.pdf")
p$pies
dev.off()
pdf("pieblank.pdf")
shut_it(p$pies)
dev.off()

df <- data.frame(table(annot[, c("Cell_type", "orig.ident")]), stringsAsFactors = FALSE)[, -2]
pp <- ggplot(df, aes(x = Cell_type, y = Freq, fill = Cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = gr.cols[levels(df[, 1]), ]) +
  theme(axis.text.x = element_text(angle = 45, face = "bold", hjust = 1)) +
  labs(title = "Size per cell type", y = "Number of cells", x = "Cell type")
# pdf("histogram.pdf", height = 7, width = 10)
# pp + drawsq
# dev.off()
# pdf("histogramblank.pdf", height = 7, width = 10)
# pp + drawsq + shut_up
# dev.off()

mygrob <- ggplotGrobf(p$pies + shut_up)
pdf("histogram_pie.pdf", height = 7, width = 10)
pp +
annotation_custom(
  grob = mygrob,
  xmin = -1.3,
  xmax = 6,
  ymin = 4000,
  ymax = 11000
)
dev.off()

mygrob <- ggplotGrobf(shut_it(p$pies))
pdf("histogram_pieblank.pdf", height = 7, width = 10)
pp + shut_up +
annotation_custom(
  grob = mygrob,
  xmin = -1.2,
  xmax = 6,
  ymin = 4000,
  ymax = 10800
)
dev.off()

newgorder <- rev(names(sort(table(annot$Cell_type))))
annot$Cell_type <- factor(annot$Cell_type, newgorder)
p <- get_props(metadata = annot, group_by = 'Cell_type', resolution = 'orig.ident',
  couls = gr.cols, reverseit = TRUE, v = TRUE)
pdf("pie_ordered.pdf")
p$pies
dev.off()
p$pies$layers <- p$pies$layers[1]
pdf("pie_orderedblank.pdf")
p$pies + shut_up
dev.off()

### 2.D Chilli plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- remove.factors(mycells@meta.data)
expdata <- cts2cpm(mycells@assays$RNA@counts)
# grps <- c("TH2", "TH1", "TH17")
# lgenes <- make_list(res[getsubset(c("group", grps), df = res, v= TRUE), ], colname = "group")
lgenes <- list(
  TH2 = c("IL13", "IL4", "IL1RL1", "IL5", "PLA2G16", "ICOS", "GATA3", "IL17RB", "GADD45G", "IL21"),
  TH1 = c("ZBED2", "SEMA7A", "CXCR3", "FASLG", "IFNG", "PRF1", "KLRG1", "CD27", "CD226", "NFATC1"),
  TH17 = c("IL17F", "IL17A", "CCL20", "CTSH", "IL22", "CCR6"),
  THIFNR = c("MX1", "IFI6", "ISG15", "IFI44L", "ISG20", "IFIT3", "IFITM1", "OAS1", "OAS3", "IFIT1")
)
lgenes <- list(
  TH2 = c("IL13", "IL4", "IL1RL1", "IL5", "PLA2G16", "ICOS", "GATA3", "IL17RB", "GADD45G", "IL21"),
  TH1 = c("XCL2", "SEMA7A", "CXCR3", "FASLG", "IFNG", "PRF1", "KLRG1", "XCL1", "CD226", "NFATC1"),
  TH17 = c("IL17F", "IL17A", "CCL20", "CTSH", "IL22", "CCR6"),
  THIFNR = c("MX1", "IFI6", "ISG15", "IFI44L", "ISG20", "IFIT3", "IFITM1", "OAS1", "OAS3", "IFIT1")
)
lgenes <- list(
  TH2 = c("IL5", "IL4", "IL13", "IL1RL1", "IL17RB", "GATA3"),
  TH1 = c("IFNG", "CXCR3", "NFATC1", "FASLG", "XCL1", "KLRG1"),
  TH17 = c("IL17A", "IL17F", "CCR6", "IL22", "CTSH", "CCL20"),
  THIFNR = c("MX1", "IFI6", "ISG15", "IFI44L", "IFIT3", "OAS1")
)
brkies <- list(filname = '%+cells', brks = c(1, 10, 20, 30, 40, 50))
colies <- c("#fffeee", "#ffc100","#ff7400", "#ff0000","#a10000", "#670000")
extramods <- NULL
mergegroups <- list(ACT = c("ACT1", "ACT2", "ACT3"))
nc <- 1

for(gr in names(lgenes)){
  cat(gr, "\n")
  dname <- paste0(gr, "_chilli", ifelse(length(mergegroups), "_2g/", "/"))
  dir.create(dname)
  genes <- if(FALSE) lgenes[[gr]] else lgenes[gr]
  dname <- paste0(dname, ifelse(is.list(genes), "grid_", ""))
  # annot$tmp <- factor(ifelse(annot$Cell_type == gr, gr, "REST"), c(gr, "REST"))
  if(!is.null(mergegroups)){
    annot$tmp <- ifelse(annot$Cell_type %in% mergegroups[[1]], names(mergegroups), annot$Cell_type)
    annot$tmp <- factor(annot$tmp, c(gorder[!gorder %in% mergegroups[[1]]], names(mergegroups)))
  }else{
    annot$tmp <- annot$Cell_type
  }
  for(g in genes){
    cat(" ", g, "\n")
    fname <- paste0(dname, ifelse(length(g) < 2, g, paste0(length(g), "genes")))
    fname <- paste0(fname, ifelse(is.null(brkies), "_pct100", tail(brkies[[2]], 1)))
    g <- getfound(g, rownames(expdata), v = TRUE)
    # source('/mnt/BioHome/ciro/scripts/functions/myplots_functions.R')
    p <- vlnplot(
      cpmdata = expdata,
      metadata = annot,
      gg = g,
      orderby = ifelse(length(mergegroups), "tmp", "Cell_type"),
      noncero = TRUE,
      # couls = c("#fffeee", "#ffe080", "#ffc100", "#ff9a00", "#ff7400", "#ff4d00", "#ff0000", "#EE0000", "#a10000", "#670000"),
      couls = colies,
      brks = brkies,
      datatype = expression(bold('Log'[2]*'(CPM + 1)')),
      ncolp = nc,
      log2t = TRUE,
      legendt = 'posiv',
      cuof = 0,
      tags = FALSE,
      rotx = 45,
      return_plot = TRUE,
      v = TRUE
    )
    pdf(paste0(fname, ".pdf"), 5, 15)
    print(p) #draw_vline(p)
    dev.off()
    pdf(paste0(fname, "blank.pdf"), 5, 15)
    print(p + shut_up)
    dev.off()
    if(!is.null(extramods)){
      pdf(paste0(fname, "_extramods.pdf"), 5, 15)
      print(p + extramods)
      dev.off()
      pdf(paste0(fname, "_extramodsblank.pdf"), 5, 15)
      print(p + extramods + shut_up)
      dev.off()
    }
  }
}

### 2.X Scatter/coexpression plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- mycells@meta.data
expdata <- cts2cpm(mycells@assays$RNA@counts)
dtype <- "cpm"
gslist <- list(
  list(genes = c("IFNG", "XCL1", "XCL2", "KLRG1"),
  ssamples = c('Cell_type', 'TH1')),
  list(genes = c("IL17A", "CTSH", "IL22", "CCL20"),
  ssamples = c('Cell_type', 'TH17'))
)
dir.create('scatters_coexp')
for(gslisty in gslist){
  cat(gslisty$ssamples[2], "\n")
  samples <- getsubset(gslisty$ssamples, annot, v = T)
  genes <- getfound(gslisty$genes[-1], rownames(expdata), v = TRUE)
  # ddf <- remove.factors(melt(t(expdata[genes, samples, drop = FALSE])))
  # colnames(ddf)[3] <- "CPM"
  # ddf <- cbind(ddf, unlist(expdata[gslisty$genes[1], ddf[, 1]]))
  # colnames(ddf)[4] <- gslisty$genes[1]
  # ddf <- cbind(ddf, annot[ddf[, 1], , drop = FALSE])
  # headmat(ddf)
  # ddf[, gslisty$genes[1]] <- log2(ddf[, gslisty$genes[1]] + 1)
  # ddf$CPM <- log2(ddf$CPM + 1)
  # p <- ggplot(ddf, aes_string(x = gslisty$genes[1], y = "CPM", color = "orig.diseasegroup")) +
  #   geom_jitter(size = 3) + facet_wrap(~ Var2, scales = 'fixed') +
  #   scale_colour_manual(values = v2cols(ddf$orig.diseasegroup, gr.cols)) +
  #   labs(x = gslisty$genes[1], y = "log2(CPM + 1)") +
  #   drawsq
  # fname <- paste0('scatters_coexp/grid_', paste0(c(dtype, gslisty$ssamples[2]), collapse = "_"), '.pdf')
  # pdf(fname)
  # print(p)
  # dev.off()
  # pdf(sub("\\.pdf", "blank.pdf", fname))
  # print(p + shut_up)
  # dev.off()

# Individual
  for(gg in lapply(genes, function(x) c(gslisty$genes[1], x) )){
    cat(commas(gg), "\n")
    p <- get_densities(mat = expdata[, samples], cuof = 1, genes = gg, subname = "10_", return_plot = TRUE)$scatter
    p <- p + geom_point(color = "black") + theme(legend.position = "none")
    fname <- paste0('scatters_coexp/', paste0(c(dtype, gslisty$ssamples[2], gg), collapse = "_"),'.pdf')
    pdf(fname)
    print(p)
    dev.off()
    pdf(sub("\\.pdf", "blank.pdf", fname))
    print(shut_it(p))
    dev.off()
  }
}
