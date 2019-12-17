#!/usr/bin/R5

############
# Figure 4 #
############

# This script will create plots for figure 4 from single-cell data

library(Seurat)
source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')

## Destiny Folder ##
root <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma/results'
setwd(paste0(root, '/a1_final_figures/Figure_4'))
setwd(paste0(root, '/a1_final_figures/Figure_3_Diseasecomps+THIFNR'))

mycells <- theObjectSavedIn(paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/clustering/zetInfo/clustCells13PCs_30Ks_0.06667JD.RData'))
gr.cols <- read.csv("/mnt/BioHome/ciro/asthma/info/AS3_final_groupColours.csv",row.names=1,stringsAsFactors=F,header=1)

diseaseorder <- c('AS_AL', 'AR', 'AS_NA', 'HC')
newlabs <- c("ACT1", "ACT2", "ACT3", "TH1", "TH2", "THIFNR", "TH17")
names(newlabs) <- 0:6
gorder <- c("TH2", "TH1", "TH17", "THIFNR", "ACT3", "ACT2", "ACT1")
mycells@meta.data$Cell_type <- factor(c(newlabs[mycells@meta.data$RNA_snn_res.0.4]), gorder)

### 4.A Pie/t-SNE charts ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sobject <- mycells
sobject@reductions$tsne <- sobject@reductions$tsne75
mycelltypes <- levels(sobject@meta.data$Cell_type)
ct <- "Teff"

dnames <- paste0(c('disease_pies_', 'disease_tsnes_'), ct, "/")
sapply(dnames, dir.create)
axmax <- max(abs(Embeddings(sobject, reduction = "tsne")))
nsmall <- min(table(sobject@meta.data$orig.diseasegroup))
sampleit <- TRUE
compileplots <- list()
for(gr in diseaseorder){
  thesecells <- getsubset(c("orig.diseasegroup", gr), sobject@meta.data, v = T)
  fname <- gr
  if(isTRUE(sampleit)){
    fname <- paste0("sampled_", gr)
    set.seed(27); thesecells <- sample(thesecells, nsmall)
  }

  p <- DimPlot(sobject, cells = thesecells, cols = v2cols(mycelltypes, gr.cols), group.by = "Cell_type", pt.size = 1.2, , reduction = "tsne")
  p <- p + xlim(c(-axmax, axmax)) + ylim(c(-axmax, axmax)) + drawsq
  compileplots[[fname]] <- p
  # pdf(paste0(dnames[1], fname, ".pdf"), height = 7, width = 8)
  # print(p);
  # dev.off()
  # pdf(paste0(dnames[1], fname, "blank.pdf"), height = 7, width = 7)
  # print(p + shut_up)
  # dev.off()
  compileplots[[paste0(fname, "blank")]] <- p + shut_up

  annot <- sobject@meta.data[thesecells, ]
  p <- get_props(metadata = annot, group_by = 'Cell_type', resolution = 'orig.diseasegroup', couls = gr.cols, reverseit = TRUE, v = TRUE)
  compileplots[[paste0(fname, "pie")]] <- p$pies
  # pdf(paste0(dnames[2], fname, ".pdf"))
  # print(p$pies);
  # dev.off()
  p$pies$layers <- p$pies$layers[1]
  # pdf(paste0(dnames[2], fname, "blank.pdf"), height = 7, width = 7)
  # print(p$pies + shut_up);
  # dev.off()
  compileplots[[paste0(fname, "pieblank")]] <- p$pies + shut_up
}
prefix <- gsub("pies_.*/$", ifelse(isTRUE(sampleit), "sampled_", "_"), dnames[1])
str(compileplots, max.level = 1)
names(compileplots)
tmpp <- compileplots[seq(1, length(compileplots), 4)]
tmpp <- lapply(tmpp, function(x) x + NoLegend())
pdf(paste0(prefix, "compiled_tsnes.pdf"), height = 7, width = 28)
plot_grid(plotlist = tmpp, ncol = length(tmpp))
dev.off()
tmpp <- compileplots[seq(2, length(compileplots), 4)]
tmpp <- lapply(tmpp, function(x) x + NoLegend() + NoAxes() )
pdf(paste0(prefix, "compiled_tsnesblank.pdf"), height = 7, width = 28)
plot_grid(plotlist = tmpp, ncol = length(tmpp))
dev.off()
tmpp <- compileplots[seq(3, length(compileplots), 4)]
tmpp <- lapply(tmpp, function(x) x + NoLegend())
pdf(paste0(prefix, "compiled_pies.pdf"), height = 7, width = 28)
plot_grid(plotlist = tmpp, ncol = length(tmpp))
dev.off()
tmpp <- compileplots[seq(4, length(compileplots), 4)]
tmpp <- lapply(tmpp, function(x) x + NoLegend() + NoAxes() )
pdf(paste0(prefix, "compiled_piesblank.pdf"), height = 7, width = 28)
plot_grid(plotlist = tmpp, ncol = length(tmpp))
dev.off()

### 4.B Pie charts per cluster disease proportion ###%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnames <- paste0("clusters_pies_", ct, "/")
dir.create(dnames)
annot <- sobject@meta.data
annot$Cell_type <- as.character(annot$Cell_type)
sampleit <- TRUE
nsmall <- min(table(annot$orig.diseasegroup))
thesecells <- lapply(c("AR", "AS_AL", "AS_NA", "HC"), function(gr){
  sample(getsubset(c("orig.diseasegroup", gr), annot, v = T), nsmall)
})
# annot <- annot[sample_grp(annot, "orig.diseasegroup"), ]
annot <- annot[unlist(thesecells), ]
compileplots <- list()
for(gr in mycelltypes){
  thesecells <- getsubset(c("Cell_type", gr), annot, v = T)
  fname <- ifelse(isTRUE(sampleit), paste0("sampled_", gr), gr)
  p <- get_props(metadata = annot[thesecells, ], group_by = 'orig.diseasegroup', resolution = 'Cell_type', couls = gr.cols, reverseit = TRUE, v = TRUE)
  compileplots[[paste0(fname, "pie")]] <- p$pies
  pdf(paste0(dnames, fname, ".pdf"))
  print(p$pies);
  dev.off()
  p$pies$layers <- p$pies$layers[1]
  pdf(paste0(dnames, fname, "blank.pdf"), height = 7, width = 7)
  print(p$pies + shut_up);
  dev.off()
  compileplots[[paste0(fname, "pieblank")]] <- p$pies + shut_up
}
prefix <- gsub("pies_.*/", ifelse(isTRUE(sampleit), "sampled_", "_"), dnames[1])
tmpp <- compileplots[seq(1, length(compileplots), 2)]
tmpp <- lapply(tmpp, function(x) x + NoLegend())
pdf(paste0(prefix, "compiled_pies.pdf"), height = 6, width = 24)
plot_grid(plotlist = tmpp, ncol = length(tmpp))
dev.off()
tmpp <- compileplots[seq(2, length(compileplots), 2)]
tmpp <- lapply(tmpp, function(x) x + NoLegend() + NoAxes() )
pdf(paste0(prefix, "compiled_piesblank.pdf"), height = 7, width = 28)
plot_grid(plotlist = tmpp, ncol = length(tmpp))
dev.off()

tvary <- names(compileplots)
tvary[seq(2, length(compileplots), 2)] <- NA
tmpp <- lapply(tvary, function(x) if(!is.na(x)) compileplots[[x]] + NoLegend() else NULL )
pdf(paste0(prefix, "compiled_pies_check.pdf"), height = 7, width = 32)
plot_grid(plotlist = tmpp, ncol = length(tmpp), rel_widths = rep(c(1, 0.2), length(tmpp)))
dev.off()
tvary <- names(compileplots)
tvary[seq(1, length(compileplots), 2)] <- NA
tmpp <- lapply(tvary, function(x) if(!is.na(x)) compileplots[[x]] + NoLegend() else NULL )
pdf(paste0(prefix, "compiled_piesblankcheck.pdf"), height = 7, width = 32)
plot_grid(plotlist = tmpp, ncol = length(tmpp), rel_widths = rep(c(0.2, 1), length(tmpp)))
dev.off()

### 4.X Disease differences ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mycells@meta.data$orig.allergy <- ifelse(mycells@meta.data$orig.diseasegroup %in% c("AS_AL", "AR"), "AL", "NAL")
mycells@meta.data$orig.asthma <- ifelse(mycells@meta.data$orig.diseasegroup %in% c("AS_AL", "AS_NA"), "AS", "NAS")
table(mycells@meta.data[, c("Cell_type", "orig.diseasegroup")])
table(mycells@meta.data[, c("Cell_type", "orig.allergy")])
table(mycells@meta.data[, c("Cell_type", "orig.asthma")])

cnames <- c("orig.allergy", "orig.asthma")
cnames <- c("orig.diseasegroup")
annot <- remove.factors(mycells@meta.data[getsubset(c("Cell_type", "-TH2"), mycells@meta.data, v = TRUE), ])
minc <- sapply(cnames, function(x){
  y <- table(annot[, c("Cell_type", x)])
  print(y)
  min(y)
})
void <- lapply(unique(annot$Cell_type), function(x){
  cellstab <<- annot[getsubset(c("Cell_type", x), annot, v = TRUE), ]
  thesecells <- lapply(names(minc), function(cname){
    lapply(unique(cellstab[, cname]), function(gr){
      set.seed(27)
      scells <- sample(getsubset(c(cname, gr), cellstab, v = TRUE), minc[cname])
      if(cname == names(minc)[1]) cellstab <<- cellstab[!rownames(cellstab) %in% scells, ]
      scells
    })
  })
})
str(void)
sampled_cells <- unlist(void)
head(sampled_cells)
sum(duplicated(sampled_cells))
table(annot[sampled_cells[duplicated(sampled_cells)], c("Cell_type", "orig.asthma")])
table(annot[sampled_cells, c("Cell_type", "orig.asthma")])
table(annot[sampled_cells, c("Cell_type", "orig.allergy")])
table(annot[sampled_cells, c("orig.allergy", "orig.asthma")])
write.csv(sampled_cells, file = "sampled_clusters_per_disease.csv")

minc <- min(table(annot$Cell_type))
void <- lapply(unique(annot$Cell_type), function(gr){
    sample(getsubset(c("Cell_type", gr), annot, v = TRUE), minc)
})
str(void)
sampled_cells <- unlist(void)
head(sampled_cells)
write.csv(sampled_cells, file = "sampled_clusters_no_th2.csv")

minc <- min(table(annot$orig.diseasegroup))
void <- lapply(unique(annot$orig.diseasegroup), function(gr){
    sample(getsubset(c("orig.diseasegroup", gr), annot, v = TRUE), minc)
})
str(void)
sampled_cells <- unlist(void)
head(sampled_cells)
write.csv(sampled_cells, file = "sampled_disease_no_th2.csv")

### 4.X Combined volcanos ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname <- paste0(root, '/mast/AS3EsmTeff_daf_30p/summary/clusters/padj0.05_FC0/RNA_snn_res.0.4_SummaryDEGsTable_suas.csv')
res_spec <- readfile(fname, stringsAsFactors = FALSE, check.name = FALSE)
res_spec <- res_spec[!is.na(res_spec$group), ]
res_spec <- res_spec[, !duplicated(colnames(res_spec))]
rownames(res_spec) <- gsub("'", "", res_spec$gene_name)
# couls <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','red2','#b30000', '#670000')
# brks <- round(seq(min(res[, col_feature], na.rm = TRUE), max(res[, col_feature], na.rm = TRUE), length.out = length(couls)), 1)
# colbar <- scale_color_gradientn(name = "group", colours = couls, na.value = 'gray',
#   breaks = brks, labels = brks, limits = c(0, max(brks)), guide = "colorbar")
dir.create('crater')
axmax <- 2.5
thres <- 0.5 # log2FoldChange threshold
lgenes <- list()
for(gr in names(newlabs)[4:7]){
  # gr = names(newlabs)[6]
  cat(newlabs[gr], "\n")
  thesecells <- getsubset(c("Cell_type", newlabs[gr]), mycells@meta.data, v = TRUE)
  fstat <- list(
    Asthma = readfile(paste0(root, '/mast/AS3EsmTeff_daf_30p/comprs/dg', gr,'/ASvsNAS/results_ASvsNAS_mastlog2cpm.csv'), stringsAsFactors = FALSE, check.names = FALSE),
    Allergy = readfile(paste0(root, '/mast/AS3EsmTeff_daf_30p/comprs/dg', gr,'/ALvsNAL/results_ALvsNAL_mastlog2cpm.csv'), stringsAsFactors = FALSE, check.names = FALSE)
  )
  # intergenes <- sub("'", "", unique(unlist(sapply(fstat, "[", 7))))
  intergenes <- unlist(sapply(fstat, "[", 7))
  intergenes <- intergenes[duplicated(intergenes)]
  intergenes <- fstat[[1]]$gene_name[fstat[[1]]$gene_name %in% intergenes]
  fstat <- lapply(fstat, function(x){
    y <- x
    rownames(y) <- x$gene_name
    y[intergenes, ]
  })
  cts <- mycells@assays$RNA@counts[, thesecells]
  expdata <- cts2cpm(cts)
  myfcs <- lapply(fstat, "[", 3)
  # fstat[[1]][, 3, drop = FALSE]
  mydatafc <- data.frame(myfcs, stringsAsFactors = F)
  names(mydatafc) <- names(myfcs)
  rownames(mydatafc) <- gsub("'", "", intergenes)
  head(mydatafc)
  mydatafc$gene_name <- rownames(mydatafc)
  mydatafc$padj <- apply(data.frame(lapply(fstat, "[", 6)), 1, min, na.rm = TRUE)
  mydatafc$significance <- -log10(mydatafc$padj)
  mydatafc <- mydatafc[complete.cases(mydatafc), ]
  mydatafc$pct <- get_stat_report(mat = expdata, moments = "p", rnames = mydatafc$gene_name, v = TRUE)
  tvar <- mydatafc$padj > 0.05 | rowSums(abs(mydatafc[, 1:2]) <= thres) == 2
  max(abs(mydatafc[tvar, 1:2]))
  table(tvar)
  lgenes[[newlabs[gr]]] <- rownames(mydatafc[tvar, ])
  mydatafc$pct[tvar] <- 0
  mydatafc$Mean <- rowMeans(expdata[rownames(mydatafc), ])
  mydatafc$Mean[tvar] <- NA
  tvar <- mydatafc$Mean > 10
  table(tvar)
  # mydatafc <- mydatafc[tvar, ]## Genes over 10 CPM of mean expression
  mydatafc$Log2_Mean <- log2(mydatafc$Mean)
  mydatafc$Log2_Mean[mydatafc$Log2_Mean > 12] <- 12
  tvar <- unique(c(bordering(mydatafc, 1:2), bordering(mydatafc[!is.na(mydatafc$Mean), ], 1:2), "TNFSF10", "CXCL10"))
  mydatafc_repel <- mydatafc[tvar, ]
  p <- ggplot(mydatafc, aes(x = Asthma, y = Allergy, color = Log2_Mean, size = pct)) +
    geom_point() + scale_color_gradientn(colours = c('#fffff0', '#670000'), na.value = 'gray') + ##CB4154
    geom_hline(yintercept = c(-thres, thres), linetype = "dashed", alpha = 0.4) +
    geom_vline(xintercept = c(-thres, thres), linetype = "dashed", alpha = 0.4) +
    labs(title = 'Allergy - Asthma differences',
      subtitle = paste('Threshold: FC >', thres), x = 'Asthma', y = 'Allergy',
      color = 'Expression', size = '%+cells') +
    drawsq
  p <- squareplot(p, axmax)
  fname <- paste0('crater/', newlabs[gr], 'gene_FC', thres, '_update')
  pdf(paste0(fname, '.pdf'), 10, 10)
  print(p + geom_text_repel(data = mydatafc_repel, aes(label = gene_name), color = 'black'))
  dev.off()
  pdf(paste0(fname, 'blank.pdf'), 10, 10)
  print(p + shut_up)
  dev.off()
  p$data <- p$data[p$data$gene_name %in% getsubset(c("group", gr), res_spec, v = TRUE), ]
  tvar <- unique(c(bordering(p$data, 1:2, 20), bordering(p$data[!is.na(mydatafc$Mean), ], 1:2, 20), "TNFSF10", "CXCL10"))
  mydatafc_repel <- p$data[tvar, ]
  mydatafc_repel <- mydatafc_repel[mydatafc_repel$gene_name %in% getsubset(c("group", gr), res_spec, v = TRUE), ]
  fname <- paste0('crater/', newlabs[gr], 'gene_FC', thres, '_specific_update')
  pdf(paste0(fname, '.pdf'), 10, 10)
  print(p + geom_text_repel(data = mydatafc_repel, aes(label = gene_name), color = 'black'), inherit.aes = FALSE)
  dev.off()
  pdf(paste0(fname, 'blank.pdf'), 10, 10)
  print(p + shut_up)
  dev.off()
  write.csv(mydatafc, file = paste0('crater/', newlabs[gr], 'gene_FC', thres, '_specific.csv'))
}
## Checking DEGs
str(lgenes)
lgenes_no_spec <- list()
dir.create('disease_specific_degs')
for(gr in names(newlabs)[4:7]){
  lgenes_no_spec[[newlabs[gr]]] <- set_ops(lgenes[[newlabs[gr]]],
    rownames(res_spec[res_spec$group == gr, ]), addlen = TRUE, rename = c("disease", "specific"))
  head(lgenes_no_spec[[newlabs[gr]]])
  write.csv(lgenes_no_spec[[newlabs[gr]]], file = paste0("disease_specific_degs/", newlabs[gr], ".csv"))
}
str(lgenes_no_spec)

### 4.X Chilli/scatter plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- mycells@meta.data
if(!exists("expdata")) expdata <- cts2cpm(mycells@assays$RNA@counts)
genes <- c("TNFSF10", "CXCL10", paste0("TNFRSF10", LETTERS[1:5]), "IFNAR1", "IFNAR2", "IFNGR1")
dname <- c("chilli/", "scatters/")

genes <- getfound(genes, rownames(expdata), v = TRUE)
sapply(dname, dir.create)
mergegroups <- list(ACT = c("ACT1", "ACT2", "ACT3"))
annot$tmp <- ifelse(as.character(annot$Cell_type) %in% mergegroups[[1]], names(mergegroups), as.character(annot$Cell_type))
annot$tmp <- factor(annot$tmp, c(gorder[!gorder %in% mergegroups[[1]]], names(mergegroups)))
for(g in genes){
  cat(" ", g, "\n")
  p <- vlnplot(
    cpmdata = expdata,
    metadata = annot,
    gg = g,
    orderby = "tmp",
    noncero = TRUE,
    # couls = c("#fffeee", "#ffe080", "#ffc100", "#ff9a00", "#ff7400", "#ff4d00", "#ff0000", "#EE0000", "#a10000", "#670000"),
    couls = c("#fffeee", "#ffc100","#ff7400", "#ff0000","#a10000", "#670000"),
    brks = list(filname = '%+cells', brks = c(1, 10, 20, 30, 40, 50)),
    datatype = 'log2(CPM + 1)',
    log2t = TRUE,
    legendt = 'posiv',
    cuof = 0,
    tags = FALSE,
    rotx = FALSE,
    return_plot = TRUE,
    v = FALSE
  ) + drawsq
  pdf(paste0(dname[1], g, "_update.pdf"))
  print(p)
  dev.off()
  pdf(paste0(dname[1], g, "blank.pdf"))
  print(p + shut_up)
  dev.off()

  p <- FeaturePlot(mycells, features = g, reduction = 'tsne', #pt.size = 1.2,
    cols = c("#f7f7f7", "#ffc100","#ff7400", "#ff0000","#a10000", "#670000")
    # cols = c("#ffe080", "#ffc100", "#ff9a00", "#ff7400", "#ff4d00", "#ff0000", "#EE0000", "#a10000", "#670000")
  )
  pdf(paste0(dname[2], g, "_update.pdf"))
  print(p)
  dev.off()
}

####
sobject <- mycells
sobject <- mycells_tneg
source('/mnt/BioHome/ciro/scripts/functions/signature.R')
lol <- get_signature(
  mat = sobject,
  genes = c("TNFRSF10A", "TNFRSF10B", "TNFRSF10D"),
  catg = "Cell_type",
  sname = "ifnrs_tneg",
  smethods = 1,
  # feat2see = c('orig.celltype', 'TNeg'), # reference group
  # track_genes = c("IL5", "IL4", "IL13", "IL3", "IL9", "IL17A", "IL17F", "IFNG", "IL22", "CXCL10", "TNFSF10", "IL17RB"),
  grpct = 5, # percentage of reference keeping
  filtercells = FALSE,
  # gcouls = gr.cols,
  path = './',
  reversing = FALSE,
  v = TRUE
)
p <- FeaturePlot(mycells_tneg, features = c("TNFRSF10A", "TNFRSF10B", "TNFRSF10D"), reduction = 'tsne',
  cols = c("#f7f7f7", "#ffc100","#ff7400", "#ff0000","#a10000", "#670000")
)
pdf(paste0("ifnrs_tneg_update.pdf"))
print(p)
dev.off()

### 4.X Scatter plots of TNegs with TNFSF10 ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tanals <- c("AS3EsmTNeg_daf_30p_0.01m", "AS3EsmTNeg_daf_30p_0.1m")[2]
for(tanal in tanals){
  rpath <- paste0(root, '/clustering_seurat/', tanal, '/clustering/zetInfo/')
  sfile <- list.files(rpath, pattern = "clustCells")
  sfiles <- sfile[grepl("RData", sfile)]
  sfile <- "clustCells16PCs_30Ks_0.06667JD.RData"
  for(sfile in sfiles){
    mycells_tneg <- theObjectSavedIn(paste0(rpath, sfile))
    mycells@reductions$tsne <- mycells@reductions$tsne75
    g <- "TNFSF10"
    p <- FeaturePlot(mycells_tneg, features = g, reduction = 'tsne',
      cols = c("#ffe080", "#ffc100", "#ff9a00", "#ff7400", "#ff4d00", "#ff0000", "#EE0000", "#a10000", "#670000")
    )
    pdf(paste0("TNFSF10_", tanal, "_", sub(".*(..PCs).*", "\\1", sfile), ".pdf"))
    print(p)
    dev.off()
  }
}

### 4.X Signature ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mycells_tneg <- theObjectSavedIn(paste0(root, '/clustering_seurat/AS3EsmTNeg_daf_30p_0.1m/clustering/zetInfo/clustCells16PCs_30Ks_0.06667JD.RData'))
mycells_tneg@reductions$tsne <- mycells_tneg@reductions$tsne75
mycells_tneg@meta.data$Cell_type <- paste0("TNeg_", mycells_tneg@meta.data$RNA_snn_res.0.6)

fname <- paste0(root, '/mast/AS3EsmTeff_daf_30p/summary/clusters/padj0.05_FC0/RNA_snn_res.0.4_SummaryDEGsTable_suas.csv')
res <- readfile(fname, stringsAsFactors = FALSE, check.name = FALSE)
res <- res[!is.na(res$group), ]
newlabs <- c("ACT1", "ACT2", "ACT3", "TH1", "TH2", "THIFNR", "TH17")
names(newlabs) <- 0:6

fname <- paste0(root, '/clustering_seurat/AS3EsmTreg_daf_30p_0.1m/markers/20PCs_RNA_snn_res.0.2_MAST/result_MAST_DEGsTAT_LFC0.25_QVAL0.05.csv')
res <- readfile(fname, stringsAsFactors = FALSE, check.name = FALSE)
res$group <- res$cluster
newlabs <- paste0("treg", 0:2)
names(newlabs) <- 0:2

for(gr in names(newlabs)[-6]){
  genes2 <- sub("'", "", res[res$group == gr, "gene_name"])
  gactiv <- paste0(newlabs[gr], "_tneg")
  mycatg <- 'RNA_snn_res.0.6'

  if(get_version(mycells_tneg) < 3) sobject <- UpdateSeuratObject(mycells_tneg) else sobject <- mycells_tneg
  cat(gactiv, '\n')
  source('/mnt/BioHome/ciro/scripts/functions/signature.R')
  lol <- get_signature(
    mat = sobject,
    genes = genes2,
    catg = mycatg,
    sname = gactiv,
    smethods = 1,
    # feat2see = c('orig.celltype', 'TNeg'), # reference group
    # track_genes = c("IL5", "IL4", "IL13", "IL3", "IL9", "IL17A", "IL17F", "IFNG", "IL22", "CXCL10", "TNFSF10", "IL17RB"),
    grpct = 5, # percentage of reference keeping
    filtercells = FALSE,
    # gcouls = gr.cols,
    path = './',
    reversing = FALSE,
    v = TRUE
  )
}

#### 4.B Disease bar plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Treg too - with and without the dots; stats
dir.create("barplots")
# for(gr in levels(annot$Cell_type)){
annot <- mycells@meta.data
annot$orig.diseasegroup <- factor(annot$orig.diseasegroup, diseaseorder)
str(annot)
# thesecells <- getsubset(c("Cell_type", gr), annot, v = TRUE)
thesecells <- sample_grp(annot, cname = "orig.diseasegroup")
myannot <- annot[thesecells, ]
gr <- c("Cell_type", "TH2", "TH1", "TH17", "THIFNR")
# annot$orig.donor <- factor(annot$orig.donor)
source('/mnt/BioHome/ciro/scripts/functions/myplots_functions.R')
p1 <- bar_dis(annot = myannot, cnames = c("orig.donor", "orig.diseasegroup", "Cell_type"),
  subsamp = gr, cols = gr.cols)
pdf(paste0('barplots/', paste0(gr, collapse = "_"), '.pdf'), width = 10, height = 10)
lapply(p1, plot)
dev.off()
p1 <- bar_dis(annot = myannot, cnames = c("orig.donor", "orig.diseasegroup", "Cell_type"),
  subsamp = gr, cols = gr.cols, test = NULL, ncolp = 1)
pdf(paste0('barplots/', paste0(gr, collapse = "_"), '_notest.pdf'), width = 4, height = 16)
print(p1)
dev.off()
p1 <- bar_dis(annot = myannot, cnames = c("orig.donor", "orig.diseasegroup", "Cell_type"), cols = gr.cols,
  subsamp = gr, plotdots = FALSE, test = NULL, ncolp = 1)
pdf(paste0('barplots/', paste0(gr, collapse = "_"), '_notest_nodots.pdf'), width = 5, height = 20)
print(p1)
dev.off()

#### Correlation of positive cells ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sobject <- mycells_treg
subsetc <- c("Cell_type", "TREG_1")
# DimPlot(mycells_treg, group.by = "Cell_type", label = TRUE, reduction = "tsne"); dev.off()

sobject <- mycells
subsetc <- c("Cell_type", "THIFNR")

sobject <- mycells_tneg
sobject@meta.data$Celltype <- "TNEG"
subsetc <- c("Celltype", "TNEG")

thesecells <- getsubset(subsetc, sobject@meta.data, v = TRUE)
expdata <- log2(cts2cpm(sobject@assays$RNA@counts[, thesecells]) + 1); datatype <- ""
# expdata <- cts2cpm(sobject@assays$RNA@counts[, thesecells]); datatype <- ""
# expdata <- thisdata[, thesecells]; datatype <- "saver"

tgenes <- c("TNFSF10", "IFI6", "ISG15", "MX1")
ddf <- FetchData(sobject, vars = c("tSNE_1", "tSNE_2", "Cell_type"), cells = thesecells)
ddf <- cbind(ddf, t(expdata[tgenes, thesecells]))

dgenes <- remove.factors(create_pairs(data.frame(g = tgenes, stringsAsFactors = FALSE))[, 2:1])
dgenes <- dgenes[dgenes[, 1] %in% "TNFSF10", ]
pp <- apply(dgenes, 1, function(x){
  # x <- unlist(dgenes[1, ])
  sddf <- ddf[rowSums(ddf[, x] > 0) > 1, ]
  # source('/mnt/BioHome/ciro/scripts/functions/myplots_functions.R')
  suffy <- paste0(c(subsetc, datatype), collapse = "")
  # plot_corr(sddf, x[1], x[2], return_plot = FALSE, suffix = suffy, v = TRUE)
  p <- get_densities(mat = expdata, genes = x, log2t = FALSE, cuof = 0, return_plot = TRUE)$scatter
  pdf(paste0("corrplot_", paste0(x, collapse = "_"), suffy, ".pdf"))
  print(p)
  dev.off()
  pdf(paste0("corrplot_", paste0(x, collapse = "_"), suffy, "blank.pdf"))
  print(shut_it(p, lays = "nsity|oin|ext") + geom_point(shape = 1, size = 2))
  dev.off()
  ggplot(sddf, aes_string(x = x[1], y = x[2])) + geom_point(shape = 1)
})
