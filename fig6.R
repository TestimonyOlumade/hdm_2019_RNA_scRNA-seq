#!/usr/bin/R5

############
# Figure 5 #
############

# This script will create plots for figure 5 from single-cell data

library(Seurat)
source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')

## Destiny Folder ##
root <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma/results'
setwd(paste0(root, '/a1_final_figures/Figure_6'))
setwd(paste0(root, '/a1_final_figures/Figure_4_Tregs'))

mycells <- theObjectSavedIn(paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/clustering/zetInfo/clustCells13PCs_30Ks_0.06667JD.RData'))
gr.cols <- read.csv("/mnt/BioHome/ciro/asthma/info/AS3_final_groupColours.csv",row.names=1,stringsAsFactors=F,header=1)

newlabs <- c("ACT1", "ACT2", "ACT3", "TH1", "TH2", "THIFNR", "TH17")
names(newlabs) <- 0:6
gorder <- c("TH2", "TH1", "TH17", "THIFNR", "ACT3", "ACT2", "ACT1")
mycells@meta.data$Cell_type <- factor(c(newlabs[mycells@meta.data$RNA_snn_res.0.4]), gorder)

### 6.X Volcano THIFNR AS_NA vs HC ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### 6.X Heatmap THIFNR AS_NA/HC ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## THIFNR is not balanaced so we need to make AS_NA vs HC
compid <- "AS_NAvsHC"
gr <- "5"
# Go to fig5.R and, run 5.B Volcano, 5.X Heatmap and 5.X Coexpression
fsufix <- paste0("THIFNR_", compid)
res <- readfile(paste0(root, '/mast/AS3EsmTeff_daf_30p/comprs/dg5/AS_NAvsHC/results_AS_NAvsHC_mastlog2cpm.csv'), stringsAsFactors = FALSE, check.names = FALSE)

# ### Treg business
# mycells_treg <- theObjectSavedIn(paste0(root, '/clustering_seurat/AS3EsmTreg_daf_30p_0.1m/clustering/zetInfo/clustCells20PCs_30Ks_0.06667JD.RData'))
# mycells_treg@reductions$tsne <- mycells_treg@reductions$tsne200
# mycells_treg@meta.data$Cell_type <- paste0("TREG_", mycells_treg@meta.data$RNA_snn_res.0.2)

mycelltypes <- unique(mycells_treg@meta.data$Cell_type)
sobject <- mycells_treg
ct <- "Treg"
# Go to fig4.R "Pie/t-SNE charts"

### 6.X t-SNE ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names(mycells_treg@reductions)
p1 <- DimPlot(mycells_treg, reduction = "tsne", group.by = "Cell_type", cols = v2cols(mycelltypes, gr.cols), pt.size = 1) +
  theme(legend.position = "none")
pdf("tSNE_tregs.pdf")
print(p1)
dev.off()
pdf("tSNE_tregsblank.pdf")
print(shut_it(p1))
dev.off()

### 6.X Chilli/scatter plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- mycells@meta.data
expdata <- cts2cpm(mycells@assays$RNA@counts)
genes <- c("TNFSF10", "CXCL10", paste0("TNFRSF10", LETTERS[1:5]))
# Go to fig4.R "Chilli/scatter plots"

annot <- mycells_treg@meta.data
if(!exists("expdata")) expdata <- cts2cpm(mycells_treg@assays$RNA@counts)
lgenes <- list(
  TREGINFR = c("ISG15", "IFI6", "MX1", "IFIT3", "SAMD9L", "OAS1", "TNFSF10")
)
brkies <- list(filname = '%+cells', brks = c(1, 10, 20, 30, 40, 50))
colies <- c("#fffeee", "#ffc100","#ff7400", "#ff0000","#a10000", "#670000")
brkies <- NULL
colies <- NULL
extramods <- ylim(c(5, 15))
mergegroups <- NULL

### 3.X qlucore file ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- remove.factors(mycells_treg@meta.data)
# annot <- annot[getsubset(c("orig.diseasegroup", "AR", "AS_AL"), annot, v = TRUE), ]
thesecells <- sample_grp(annot, cname = "Cell_type", maxln = -1000)

annot$sampled <- rownames(annot) %in% thesecells
tvar <- sapply(annot, function(x) length(table(x)) )
myannot <- annot[, tvar > 1 & tvar < 200 & !names(tvar) %in% c("tmp", "cluster_number")]
myannot <- annot[thesecells, c("Cell_type", "orig.experiment", "orig.diseasegroup", "sampled", "orig.class", "orig.donor")]
head(myannot)

tpm <- cts2cpm(mycells_treg@assays$RNA@counts)
void <- add_gene_tag(lgenes = c("CXCL10", "TNFSF10"), annot = myannot, mat = tpm, v = TRUE)
myannot <- cbind(myannot, void[rownames(myannot), ])

thresh <- 1
genes <- rownames(tpm[rowMeans(tpm[, rownames(myannot)]) > thresh, ])
dname <- "qlucore_treg_subclustering"
dir.create(dname)
qf <- qlucore_format(mat = tpm, metadata = myannot, rnames = genes, v = TRUE)
write.csv(qf, file = paste0(dname, '/qlucore_mean_cpm_greater_than', thresh, '.csv'))

qfa <- qlucore_format(mat = tpm, metadata = myannot, v = TRUE)
write.csv(qfa, file = paste0(dname, '/qlucore_all_genes.csv'))

### 6.X Heatmap cluster-specific Treg ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system("cp ~/large/asthma/results/mast/AS3EsmTreg_daf_30p_0.1m/summary/clusters_padj0.05_FC0/RNA_snn_res.0.2_byGroups_meanHeatmap_suas.pdf treg_byGroups_meanHeatmap_suas.pdf")
system("cp ~/large/asthma/results/mast/AS3EsmTreg_daf_30p_0.1m/summary/clusters_padj0.05_FC0/RNA_snn_res.0.2_SummaryDEGsTable_suas.csv treg_SummaryDEGsTable_suas.csv")

### Donors plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sobject <- mycells_treg
dset <- "donors_tsnes_treg"
# got ot fig5.R "donors tSNEs"

### Gene densities ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expdata <- cts2cpm(mycells@assays$RNA@counts)

gr <- "TH1"
genes <- c("IFNG", "XCL1", "XCL2", "KLRG1")
gr <- "TH17"
genes <- c("IL17A", "IL17F", "IL22", "CCL20", "CTSH")

thesecells <- getsubset(c("Cell_type", gr), mycells@meta.data, v = TRUE)
genesdf <- create_pairs(data.frame(genes, stringsAsFactors = FALSE))[, 1:2]
dname <- paste0("correlation_", gr, "/")
dir.create(dname)
void <- apply(genesdf, 1, function(x){
  p <- get_densities(mat = expdata[, thesecells], genes = x, return_plot = TRUE)$scatter
  pdf(paste0(dname, paste0(x, collapse = "_"), ".pdf"))#, 5, 5)
  print(p)
  dev.off()
  pdf(paste0(dname, paste0(x, collapse = "_"), "blank.pdf"))#, 5, 5)
  print(shut_it(p, "ext|ensity"))
  dev.off()
})
for(x in genes){
  pdf(paste0(dname, x, "_violin.pdf"), width = 3.5, height = 7)
  print(vlnplot(cpmdata = expdata[, thesecells], gg = x, plotdots = FALSE, v = FALSE, log2t = TRUE,
    datatype = 'log2(TPM + 1)', return_plot = TRUE, noncero = TRUE))
  dev.off()
}

### 4.X Combined volcanos ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sobject <- mycells_treg
fname <- paste0(root, '/mast/AS3EsmTreg_daf_30p_0.1m/summary/clusters_padj0.05_FC0/RNA_snn_res.0.2_SummaryDEGsTable_suas.csv')
res_spec <- readfile(fname, stringsAsFactors = FALSE, check.name = FALSE)
res_spec <- res_spec[!is.na(res_spec$group), ]
res_spec <- res_spec[, !duplicated(colnames(res_spec))]
rownames(res_spec) <- gsub("'", "", res_spec$gene_name)
dname <- "treg_craters"
dir.create(dname)
axmax <- 4
thres <- 0.5 # log2FoldChange threshold
lgenes <- list()
mylabs <- unique(sobject@meta.data$Cell_type)
names(mylabs) <- 0:(length(mylabs) - 1)
for(gr in names(mylabs)){
  # gr = names(mylabs)[1]
  cat(mylabs[gr], "\n")
  thesecells <- getsubset(c("Cell_type", mylabs[gr]), sobject@meta.data, v = TRUE)
  fstat <- list(
    Asthma = readfile(paste0(root, '/mast/AS3EsmTreg_daf_30p_0.1m/comprs/disease_cl', gr,'/ASvsNAS/results_ASvsNAS_mastlog2cpm.csv'), stringsAsFactors = FALSE, check.names = FALSE),
    Allergy = readfile(paste0(root, '/mast/AS3EsmTreg_daf_30p_0.1m/comprs/disease_cl', gr,'/ALvsNAL/results_ALvsNAL_mastlog2cpm.csv'), stringsAsFactors = FALSE, check.names = FALSE)
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
  cts <- sobject@assays$RNA@counts[, thesecells]
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
  lgenes[[mylabs[gr]]] <- rownames(mydatafc[tvar, ])
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
  fname <- paste0(dname, '/', mylabs[gr], '_FC', thres, '')
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
  fname <- paste0(fname, '_specific')
  pdf(paste0(fname, '.pdf'), 10, 10)
  print(p + geom_text_repel(data = mydatafc_repel, aes(label = gene_name), color = 'black'), inherit.aes = FALSE)
  dev.off()
  pdf(paste0(fname, 'blank.pdf'), 10, 10)
  print(p + shut_up)
  dev.off()
  write.csv(mydatafc, file = paste0(fname, '.csv'))
}
## Checking DEGs
str(lgenes)
lgenes_no_spec <- list()
dir.create('disease_specific_degs')
for(gr in names(mylabs)[4:7]){
  lgenes_no_spec[[mylabs[gr]]] <- set_ops(lgenes[[mylabs[gr]]],
    rownames(res_spec[res_spec$group == gr, ]), addlen = TRUE, rename = c("disease", "specific"))
  head(lgenes_no_spec[[mylabs[gr]]])
  write.csv(lgenes_no_spec[[mylabs[gr]]], file = paste0("disease_specific_degs/", mylabs[gr], ".csv"))
}
str(lgenes_no_spec)
