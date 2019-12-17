#!/usr/bin/R5

############
# Figure 3 #
############

# This script will create plots for figure 3 from single-cell data
# This script is to create gene correlations from AS3EsmTeff_daf_30p

.libPaths('~/R/newer_packs_library/3.5')
library(SAVER)
library(Seurat)
source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')

## Destiny Folder ##
root <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma/results'
setwd(paste0(root, '/a1_final_figures/Figure_3'))
setwd(paste0(root, '/a1_final_figures/Figure_5_TH2co-exp_and_sub-clustering'))

mycells <- theObjectSavedIn(paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/clustering/zetInfo/clustCells13PCs_30Ks_0.06667JD.RData'))
gr.cols <- read.csv("/mnt/BioHome/ciro/asthma/info/AS3_final_groupColours.csv",row.names=1,stringsAsFactors=F,header=1)

newlabs <- c("ACT1", "ACT2", "ACT3", "TH1", "TH2", "THIFNR", "TH17")
names(newlabs) <- 0:6
gorder <- c("TH2", "TH1", "TH17", "THIFNR", "ACT3", "ACT2", "ACT1")
mycells@meta.data$Cell_type <- factor(c(newlabs[mycells@meta.data$RNA_snn_res.0.4]), gorder)

### 3.A Correlations ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dname <- "mycorrelations"
dir.create(dname)
### choose data ###
# raw
thisdata <- mycells@assays$RNA@counts
datatype <- 'raw'

# cts <- cts2cpm(mycells@assays$RNA@counts)
# mat <- log2(cts + 1)
thisdata <- mat
datatype <- 'cpm'

# mag <- theObjectSavedIn(paste0(root, '../magic/magic_AS3Esm.RData'))
thisdata <- mag
datatype <- 'magic'

# savey <- theObjectSavedIn(paste0(root, '/../raw/10X/AS3Esm/AGGR/AS3Esm/outs/imputed/saver_matrix.RData'))
# savey <- cts2cpm(savey)
# saveymat <- log2(savey + 1)
# # saver_obj <- theObjectSavedIn(paste0(root, '../../raw/10X/AS3Esm/AGGR/AS3Esm/outs/imputed/combined.RData'))
# # saver_obj_ss <- list(estimate = saver_obj[[1]][tmp, mysamples], se = saver_obj[[2]][mygenes, mysamples])
# # saver_cor_gene <- cor.genes(saver_obj)
thisdata <- saveymat
datatype <- 'saver'

# knn <- theObjectSavedIn(paste0(root, '../magic/knn_smoothed_AS3Esm.RData'))
thisdata <- knn
datatype <- 'knn_smooth'

### Get genes file ###
resolut <- "RNA_snn_res.0.4"
tag <- c('clusters')[1]
co <- 'padj0.05_FC0'
tvar <- paste0(root, '/mast/AS3EsmTeff_daf_30p/summary/', tag,'/', co,'/')
list.files(tvar)
clust <- read.csv(paste0(tvar, resolut, '_SummaryDEGsTable_suas.csv'), stringsAsFactors = F, check.names = F)
headmat(clust)
clust <- clust[!is.na(clust$group), ]
rownames(clust) <- sub("'", "", clust$gene_name)

### Groups and cutoffs ###
plus_hvg <- FALSE # adding highly variable genes
pct <- 5 # filter by pct expressing
whichy <- 1
cl <- as.character(c(4, 6, 3, 1, 5)[whichy])
gr <- as.character(c(4, 6, 3, 1, 5)[whichy])
mybrks <- c('seq', 'quantile', 0.3)[2]
if(!is.na(as.numeric(mybrks))) mybrks <- as.numeric(mybrks)
mysufix <- paste0('clust', paste0(newlabs[cl], collapse = '_'),
                  '_spec', paste0(newlabs[gr], collapse = '_'),
                  '_', ifelse(is.numeric(mybrks), paste0('cut', mybrks), mybrks))
mysamples <- getsubset(c(resolut, cl), mycells@meta.data, v = T)
mysufix

### Get genes to use ###
## group
mygenes <- getsubset(c('group', gr), clust, v = T)

# adding marker genes
# loading markers from subclustering
base <- "/clustering_seurat/AS3EsmTeff_daf_30p/th2pv_0.1"
markname <- paste0(root, base, "/markers/14PCs_RNA_snn_res.0.2_MAST/result_MAST_DEGsTAT_LFC0.25_QVAL0.05.csv")
res_seu <- readfile(markname, row.names=1,stringsAsFactors=F,header=1)
head(res_seu)
res_seu$gene <- gsub("'", "", res_seu$gene_name)
if(FALSE){
  mygenes <- unique(c(mygenes, res_seu$gene)); mysufix <- paste0(mysufix, "_markers")
}

## highly variable
if(plus_hvg){
  mysufix <- paste0(mysufix, "_hvg")
  mycells_sub <- mycells[, mysamples]
  mycells_sub <- FindVariableFeatures(mycells_sub)
  hvftab <- HVFInfo(mycells_sub)
  hvftab <- hvftab[order(-hvftab[, 3]), ]
  hvftab <- hvftab[hvftab[, 1] > 0.1, ]
  head(HVFInfo(mycells_sub)[VariableFeatures(mycells_sub), ])
  sum(VariableFeatures(mycells_sub) %in% mygenes)
  sum(rownames(head(hvftab, Inf)) %in% mygenes)
  sum(VariableFeatures(mycells_sub) %in% rownames(head(hvftab, 2000)))
  maxgene <- max(which(rownames(hvftab) %in% mygenes))
  if(maxgene < length(mygenes) && nrow(hvftab) > length(mygenes)) maxgene <- length(mygenes)
  mygenes <- rownames(hvftab)[1:maxgene]
}


## Filter by pct of cells expressing
if(isTRUE(pct > 0)){
  # mygenes <- mygenes[paste0("'", mygenes) %in% clust[clust[, paste0(cl, "_exprFrac")] >= pct, "gene_name"]]
  pcts <- rowMeans(mycells@assays$RNA@counts[mygenes, mysamples] > 0) * 100
  summary(pcts)
  mygenes <- mygenes[pcts >= pct]
  mysufix <- paste0(mysufix, "_pct", pct)
}
mysufix <- paste0(mysufix, "_genes", length(mygenes))

### highlight some genes ####
# specified ones
lgenes <- list(
  th_treg = c("TNFRSF4", "WARS", "LTA", "BCL2L1", "PSMD11", "CD70", "EIF5A2", "FASLG", "EIF2C3", "PSMB2", "E2F6", "FARS2"),
  th = c("MIR155HG", "TNF", "CSF2", "IL5", "IFNG", "IL17A", "BCL2A1", "IL2", "IL4", "IL13", "CD40LG"),
  treg = c("IL2RA", "FOXP3", "CTLA4", "TNFRSF1B", "IL1R2", "HLA-DRA", "CORO2A", "PRDM1", "IKZF2", "TNFRSF8", "LRRC32"),
  th_fig_g_h = c("IL5", "IL13", "IL4", "IL9", "IL31", "IL1RL1", "IL17RB", "ISG15", "IFI6", "MX1", "IFI44", "IFNG")
)
tmp <- unique(c(unlist(lgenes), "IL3", "GMCSF", "IL21", "CD109", "STAT4", "TIGIT", "TNFSF10"))
showgenes <- rep('red', length(tmp))
addleg <- FALSE

# # to add a legend group colours
# grps <- head(unique(clust$group), Inf)
# tmp <- as.vector(unlist(sapply(grps, function(x) head(clust[clust$group == x, 'gene_name']) )))
# tmp <- sub("'", "", tmp)
# showgenes <- as.vector(unlist(sapply(grps, function(x) head(clust[clust$group == x, 'colour_key']) )))
# addleg <- TRUE
####

names(showgenes) <- tmp
showgenes <- showgenes[getfound(names(showgenes), mygenes, v = TRUE)]

# means <- clust[mygenes, grepl(paste0("^", cl, ""), colnames(clust)), drop = FALSE]
# head(means)
# means[c("NIPA1", "CMC2", "UFM1"), ]
# head(means[means[, 1] < 0.1, , drop = FALSE])
# summary(means)

fname <- paste0(dname, "/", datatype, '_correlation_', mysufix)
# # mydata <- thisdata[mygenes, mysamples]
# mydata <- t(prcomp(thisdata[mygenes, mysamples])$x[1:3, ])
# wss <- (nrow(mydata)-1) * sum(apply(mydata, 2, var))
# for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers = i)$withinss)
# require(vegan)
# # fit <- cascadeKM(scale(mydata, center = TRUE,  scale = TRUE), 1, 10, iter = 1000)
# fit <- cascadeKM(mydata, 1, 10, iter = 1000)
# calinski.best <- as.numeric(which.max(fit$results[2,]))
# cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")
# library(apcluster)
# d.apclus <- apcluster(negDistMat(r=2), mydata)
# cat("affinity propogation optimal number of clusters:", length(d.apclus@clusters), "\n")
# library(cluster)
# cgap <- clusGap(mydata, kmeans, 10, B = 100)
#
# pdf(paste0(fname, '_brk_chooseCLS.pdf'), width = 10, height = 10)
# plot(1:15, wss, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares")
# plot(fit, sortg = TRUE, grpmts.plot = TRUE)
# heatmap(d.apclus)
# plot(d.apclus, mydata)
# plot(1:10, cgap$Tab[, 3])
# dev.off()
source('/mnt/BioHome/ciro/scripts/functions/myplots_functions.R')
pdf(paste0(fname, '_brk.pdf'), width = 10, height = 10)
corrgenes <- geneCorr(
  expdata = thisdata,
  cortype = 'spearman',
  brks = 0.5,
  mmain = paste('Correlation:\nCluster', commas(newlabs[cl]), 'cells'),
  sublab = paste0('Genes: ', commas(newlabs[gr]), '-specific'),
  keepgenes = mygenes,
  samples = mysamples,
  cols = showgenes,
  verbose = TRUE,
  clust_cols = 5,
  rank_genes = T
)
if(addleg){
  legend('bottomleft', legend = grps,
  fill = unique(clust$colour_key)[1:length(grps)], cex=.6, bty='n',
  title = 'Genes groups', title.col = 'red')
}
dev.off()
write.csv(corrgenes, file = paste0(fname, '.csv'))

# overlap with marker genes
restype <- c("MAST", "S")[1]
if(restype == "MAST"){
  fnamey <- paste0(root, '/mast/th2pv_0.1/comprs/clusters2/0vs1/results_0vs1_mastlog2cpm.csv')
  res_mast <- readfile(fnamey)
  res_mast <- res_mast[getDEGenes(res_mast, pv = 0.05, v = TRUE), ]
  res_mast$gene <- gsub("'", "", res_mast$gene_name)
  res_mast$clustert <- paste0(restype, res_mast$group)
  resy <- res_mast
}else{
  resy <- res_seu
  resy$clustert <- paste0(restype, resy$cluster)
}

corrgenes$clustert <- paste0("M", corrgenes$cluster)
lgenes2ov <- c(make_list(resy, colname = "clustert", col_objects = "gene"),
  make_list(corrgenes, colname = "clustert", col_objects = "gene_name"))
str(lgenes2ov)
lgenesov <- overlap_calc(lgenes2ov, v = TRUE)
lgenesov <- lgenesov[sapply(lgenesov, length) > 0]
str(lgenesov)
lgenesovdf <- rbind(sapply(lgenesov, length), vlist2df(lgenesov))
head(lgenesovdf)
write.csv(lgenesovdf, file = paste0("modules_markersby", restype, "_overlaps.csv"), row.names = FALSE)

### Gephi files ###
grgmods <- c("1" = "3", "2" = "4")
modulo <- "4"
newsufix <- mysufix
if(!is.null(modulo)){
  newsufix <- paste0(mysufix, "_mod", paste0(grgmods[modulo], collapse = "n"))
}
selgenes <- corrgenes[getsubset(c("cluster", modulo), corrgenes, v = TRUE), "gene_name"]
dgephi <- "mygephi_input"
dir.create(dgephi)
thisdata[selgenes, mysamples] -> datExpr
library(WGCNA)

powers = c(c(1:10), seq(from = 12, to=60, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf(paste0(dgephi, "/", datatype, "_", newsufix, "_powerd_selection.pdf"), 9, 5)
par(mfrow = c(1,2));
cex1 = 0.5;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
cat("Predicted", min(which(-sign(sft$fitIndices[,3])*sft$fitIndices[,2] > 0.9)), "\n")

chosenpw <- switch(newsufix,
  cpm = 10,
  saver = 30,
  knn_smooth = 5,
  clustTH2_specTH2_quantile_genes240 = 35,
  clust4_spec4_quantile_pct5_genes495 = 36,
  clustTH2_specTH2_quantile_pct5_genes214 = 35,
  clustTH17_specTH17_quantile_pct5_genes14 = 20,
  30
)
chosenpw <- ifelse(chosenpw > 30, 30, 5)
chosenpw <- 3
TOM = TOMsimilarityFromExpr(t(datExpr), power = chosenpw);
TOM[1:5, 1:5]
summary(c(as.matrix(TOM)))

modProbes = rownames(datExpr)
modTOM = TOM
dimnames(modTOM) = list(modProbes, modProbes)
tedge <- exportNetworkToCytoscape(modTOM, weighted = TRUE, threshold = quantile(TOM,  prob = 0.5))$edgeData
head(tedge)
round(range(tedge$weight), 2)
dim(tedge)
modules <- paste0(datatype, "_Cluster", newlabs[cl], "_Genes", newlabs[gr])
cyt = exportNetworkToCytoscape(modTOM,
  edgeFile = NULL,
  nodeFile = NULL,
  # edgeFile = paste0("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt"),
  # nodeFile = paste0("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt"),
  weighted = TRUE, threshold = quantile(modTOM,  prob = 0.5))
str(cyt)

edgdata <- cyt$edgeData[, c(1:3, 5, 6)]
colnames(edgdata) <- c('fromNode','toNode','weight','Source','Target')
if(sum(levels(edgdata[,4]) == 'NA') > 0 || sum(levels(edgdata[,5]) == 'NA') > 0){
  edgdata[,4] <- edgdata[,1]
  edgdata[,5] <- edgdata[,2]
}
noddata <- cyt$nodeData
if(sum(levels(noddata[,2]) == 'NA') > 0){
  noddata[,2] <- noddata[,1]
}
if(sum(is.na(noddata[,3])) > 0){
  noddata[,3] <- paste0("genes", gr)
}
noddata$Label <- noddata[,2]
noddata <- noddata[,c(1,2,4,3)]
colnames(noddata) <- c('nodeName','ID','Label','Module')
dim(edgdata)
dim(noddata)
lapply(list(edgdata, noddata), function(x) summary(x) )
cat("Correlations:", min(edgdata$weight), "to", max(edgdata$weight), "\n")
write.csv(edgdata, file = paste0(dgephi, '/', newsufix, '_p', chosenpw, '_edges.csv'), row.names = F)
write.csv(noddata, file = paste0(dgephi, '/', newsufix, '_p', chosenpw, '_nodes.csv'), row.names = F)

## ggraph (g-giraffe)
library(ggraph)
library(tidygraph)
library(plotly)
library(htmlwidgets)
# head(highschool)
# #   from to year
# # 1    1 14 1957
# # 2    1 15 1957
# gdata <- highschool
# graph <- as_tbl_graph(gdata) %>%
#   mutate(Popularity = centrality_degree(mode = 'in'))
# graph
thesegenes <- switch(newlabs[gr],
  TH2 = c("IL5", "IL13", "IL4", "IL9", "IL31", "IL1RL1", "IL17RB", "GZMB", "ZEB2", "TGFBR3", "PPARG", "IL21", "CTLA4", "GATA3", "IL3"),
  names(showgenes)
)

## Checking correlation trends
thesegenes <- getfound(thesegenes, rownames(thisdata), v = TRUE)
melted1 <- data.frame(t(thisdata[head(thesegenes, 5), mysamples]))
mycols <- v2cols(mycells_th2@meta.data$Cell_type)
pdf(paste0(fname, '_scatters.pdf'), width = 12, height = 12)
plot(melted1, col = mycols[mycells_th2@meta.data$Cell_type])
dev.off()
nmycols <- mycols; nmycols[2] <- NA
pdf(paste0(fname, '_scatters_0.pdf'), width = 12, height = 12)
plot(melted1, col = nmycols[mycells_th2@meta.data$Cell_type])
dev.off()
nmycols <- mycols; nmycols[1] <- NA
pdf(paste0(fname, '_scatters_1.pdf'), width = 12, height = 12)
plot(melted1, col = nmycols[mycells_th2@meta.data$Cell_type])
dev.off()

df2tbl <- function(eddata, sgenes = NULL, addata = NULL, v = FALSE){
  eddata <- remove.factors(eddata)
  ggenes <- unique(unlist(eddata[, 1:2]))
  mynodes <- data.frame(N = 1:length(ggenes), name = ggenes, row.names = ggenes)
  if(!is.null(sgenes)) mynodes$highlight <- mynodes$name %in% sgenes
  if(!is.null(addata))
    mynodes <- cbind(mynodes, addata[rownames(mynodes), , drop = FALSE])
  myedges <- data.frame(from = mynodes[eddata[, 1], "N"],
                              to = mynodes[eddata[, 2], "N"])#,
                          # weight = eddata[, 3])
  tbl_graph(nodes = mynodes, edges = myedges)
}

# weight: strength between two nodes (genes) in terms of the correlation value obtained
# from TOMsimilarity function. This value should be always between 0 and 1.
# The higher value refers to a strong connection or co-expression of genes.

tcorrgenes <- corrgenes
rownames(tcorrgenes) <- tcorrgenes[, 1]
tcorrgenes$gene_name = ifelse(tcorrgenes$gene_name %in% thesegenes, tcorrgenes$gene_name, NA)
tcorrgenes$markers <- sapply(unique(res_seu$gene), function(x) commas(res_seu[res_seu$gene == x, "cluster"]) )[tcorrgenes$gene_name]
lvls <- c("NONE", names(table(tcorrgenes$markers)))
tcorrgenes$markers <- ifelse(is.na(tcorrgenes$markers), "NONE", tcorrgenes$markers)
tcorrgenes$markers <- factor(tcorrgenes$markers, lvls)
tcorrgenes$clusterz <- ifelse(tcorrgenes$cluster == 1, "3", NA)
# extrafont::font_import(); extrafont::loadfonts() # before usint theme_graph
dir.create(paste0(fname, '_graphs'))
for(qa in c(10, 50, seq(80, 95, 5), 99)){
  # qa <- 90
  adj <- quantile(modTOM,  prob = qa / 100)
  cat(names(adj), "\n")
  gdata <- exportNetworkToCytoscape(modTOM, weighted = TRUE, threshold = adj)$edgeData
  # tvar <- table(gdata[, 1])
  # gdata <- gdata[as.character(gdata[, 1]) %in% names(tvar[tvar > 1]), ]
  gdata <- gdata[as.character(gdata$fromNode) %in% rownames(tcorrgenes[tcorrgenes$cluster == 1, ]) & as.character(gdata$toNode) %in% rownames(tcorrgenes[tcorrgenes$cluster == 1, ]), ]
  if(nrow(gdata) == 0) next
  gdata[bordering(gdata, cnames = 3), ]

  keptgenes <- getfound(selgenes, unique(unlist(gdata[, 1:2])), v = TRUE)
  nname <- paste0(fname, '_graphs/adj', qa, '_', length(keptgenes), 'genes_annot', length(thesegenes), '')
  ggdata <- df2tbl(gdata, addata = tcorrgenes, sgenes = thesegenes, v = TRUE)
  graph <- ggdata %>%
      mutate(Popularity = centrality_degree(mode = 'in'))
  graph

  set.seed(27)# plot using ggraph
  g1 <- ggraph(graph, layout = 'kk') +
    geom_edge_fan(aes(alpha = stat(index)), show.legend = FALSE)
    # guide_geom(guide = geom_node_point, = guide_colorbar(barwidth = 0.8, barheight = 6))
  # if('markers' %in% colnames(g1$data)){
  #   g1 <- g1 + geom_node_point(aes(size = Popularity, colour = factor(cluster), shape = markers), alpha = 0.9)
  # }else{
    g1 <- g1 + geom_node_point(aes(size = Popularity, colour = factor(cluster)), alpha = 0.9)
  # }
  pdf(paste0(nname, '_nameless.pdf'))
  print(g1)
  dev.off()
  g1 <- g1 + geom_node_text(aes(label = name, size = 0.1))
  pdf(paste0(nname, '_names.pdf'))
  print(g1)
  dev.off()
  # p <- ggplotly(g1)
  # setwd(paste0(fname, '_graphs/'))
  # saveWidget(as_widget(p), paste0(basename(nname), '.html'))
  # setwd(paste0(root, '/a1_final_figures/Figure_3'))

  # Zoom
  yls <- c(min(g1$data$y), g1$data[g1$data$name == "ACSL4", "y"])
  xls <- c(g1$data[g1$data$name == "HCST", "x"], g1$data[g1$data$name == "MT-CO1", "x"])
  g1 <- ggraph(graph, layout = 'kk') +
    geom_edge_fan(aes(alpha = stat(index)), show.legend = FALSE) +
    geom_node_point(aes(size = Popularity, colour = factor(clusterz)), alpha = 0.9, na.rm = TRUE)
  pdf(paste0(nname, '_nameless_zoom.pdf'))
  print(g1)
  dev.off()
  pdf(paste0(nname, '_nameless_zoomcut.pdf'), width = 5, height = 3)
  print(g1 + coord_cartesian(ylim = yls))
  dev.off()
}

### 3.X Chilli/scatter plots teff ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- mycells@meta.data
if(!exists("expdata")) expdata <- cts2cpm(mycells@assays$RNA@counts)
genes <- c("IL5", "IL13", "IL4", "IL9", "IL31", "IL17RB", "IL3",
  "IL1RL1", "IL10", "IL21", "CSF2", "GZMB", "PTGER4", "RANKL" = "TNFSF11",
  "LIGHT" = "TNFSF14", "TGFBB3", "NRF2" = "NFE2L2", "SDE4", "ICOS", "CD200", "CD109")
dname <- c("chilli_teff/", "scatters_teff/")
# go to 4.X Chilli/scatter plots in fig.4.R

### 3.X qlucore file ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluname <- paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/th2pv_0.1/clustering/zetInfo/clustCells14PCs_30Ks_0.06667JD.RData')
# cluname <- paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/th2_0.1/clustering/zetInfo/clustCells6PCs_30Ks_0.06667JD.RData')
mycells_th2 <- theObjectSavedIn(cluname)
mycells_th2@reductions$tsne <- mycells_th2@reductions$tsne50
mycells_th2@meta.data$Cell_type <- paste0("TH2_", mycells_th2@meta.data$RNA_snn_res.0.2)
annot <- remove.factors(mycells_th2@meta.data)
annot <- annot[getsubset(c("orig.diseasegroup", "AR", "AS_AL"), annot, v = TRUE), ]
minc <- min(table(annot$orig.diseasegroup))
thesecells <- sample_grp(annot, "orig.diseasegroup")
annot$sampled <- rownames(annot) %in% thesecells
tvar <- sapply(annot, function(x) length(table(x)) )
myannot <- annot[, tvar > 1 & tvar < 200 & !names(tvar) %in% c("tmp", "cluster_number")]
myannot <- annot[, c("Cell_type", "orig.experiment", "orig.diseasegroup", "sampled", "orig.class", "orig.donor")]
head(myannot)

tpm <- cts2cpm(mycells_th2@assays$RNA@counts)
void <- add_gene_tag(lgenes = c("IL9", "IL5", "IL4", "IL31", "IL3"), annot = myannot, mat = tpm, v = TRUE)
myannot <- cbind(myannot, void[rownames(myannot), ])

thresh <- 1
genes <- rownames(tpm[rowMeans(tpm[, rownames(myannot)]) > thresh, ])
dir.create("qlucore_th2_subclustering")
qf <- qlucore_format(mat = tpm, metadata = myannot, rnames = genes, v = TRUE)
write.csv(qf, file = paste0('qlucore_th2_subclustering/qlucore_mean_cpm_greater_than', thresh, '.csv'))

qfa <- qlucore_format(mat = tpm, metadata = myannot, v = TRUE)
write.csv(qfa, file = 'qlucore_th2_subclustering/qlucore_all_genes.csv')

### LIGER ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GSE89809 - extract and check comaprison AS_AL-AR
# TH9_PPARG_Th2_ScImm_2019 - in subclusters
# TFH = in subcluster
source('/mnt/BioHome/ciro/scripts/gsea/gsea_liger.R')
ctag <- "2vsglobal"
markname <- paste0(root, "/mast/th2pv_0.1/comprs/global/", ctag, "/results_", ctag, "_mastlog2cpm.csv")
resdea <- readfile(markname, row.names=1,stringsAsFactors=F,header=1)
ctag <- "TH2vsTNeg"
markname <- paste0(root, "/mast/AS3EsmTeff_daf_30p/comprs/activation/", ctag, "/results_", ctag, "_mastlog2cpm.csv")
resdea <- readfile(markname, row.names=1,stringsAsFactors=F,header=1)
ctag <- "0vs1"
markname <- paste0(root, "/mast/th2pv_0.1/comprs/clusters2/", ctag, "/results_", ctag, "_mastlog2cpm.csv")
resdea <- readfile(markname, row.names=1,stringsAsFactors=F,header=1)
void <- gsea_liger(
  res = resdea,
  gene_name = "gene_name",
  pvtype = 'padj',
  lfc.type = 'log2FoldChange',
  gsea_file = "~/asthma/info/gsea_lists_extended.csv",
  myseed = 27,
  path = paste0("gsea_liger_", ctag),
  v = TRUE
)

### Volcano cluster 0 vs 1 TH2 PC14 resolution 0.2 ###%%%%%%%%%%%%%%%%%%%%%%%%%%
# FC -5,5 with a threshold of 0.5 colour like Fig5 TH2s (slide 14)
### 5.B Volcano TH2 AS_AL vs AR ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compid <- "0vs1"
gr <- "4"
cat(newlabs[gr], "\n")
wgr <- unlist(strsplit(compid, "vs"))
res <- readfile(paste0(root, '/mast/th2pv_0.1/comprs/clusters2/', compid, '/results_', compid, '_mastlog2cpm.csv'), stringsAsFactors = FALSE, check.names = FALSE)
res$gene <- gsub("'", "", res$gene_name)
rownames(res) <- res$gene
res <- res[, !duplicated(colnames(res))]

thesecells <- getsubset(list(c("RNA_snn_res.0.2", wgr)), mycells_th2@meta.data, v = TRUE)
grps <- make_list(x = mycells_th2@meta.data[thesecells, ], colname = "RNA_snn_res.0.2", grouping = TRUE)
grps <- sort(grps)
# ## Checking stats are correct
# expdata <- mycells_th2@assays$RNA@data[, thesecells]
# stattab <- get_stat_report(mat = expdata, groups = grps, rname = res$gene, moments = c("mn"), v = TRUE)
# head(stattab)
# head(res[, grepl("_mean", colnames(res))])
# head(log2(stattab[res$gene, ] + 1))
## Correcting
expdata <- log2(cts2cpm(mycells_th2@assays$RNA@counts[, thesecells]) + 1)
# expdata <- cts2cpm(mycells_th2@assays$RNA@counts[, thesecells])
stattab <- get_stat_report(mat = expdata, groups = grps, rname = res$gene, moments = c("mn", "p"), v = TRUE)
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

res$Fold <- res$log2FoldChange
res$FDR <- -log10(res$padj)
res$pct <- res[, paste0(wgr[1], "_exprFrac")]
res$pct[res$group == wgr[2]] <- res[res$group == wgr[2], paste0(wgr[2], "_exprFrac")]
tvar <- res$padj <= 0.05 & abs(res$log2FoldChange) >= 0.5
table(tvar)
res$pct[!tvar] <- 0 # size for only significant genes
res$mean <- res[, paste0(wgr[1], "_mean")]
res$mean[res$group == wgr[2]] <- res[res$group == wgr[2], paste0(wgr[2], "_mean")]
degs <- getDEGenes(x =res, pv = 0.05, fc = 0.5, v = TRUE)
res$mean[!res$gene %in% degs] <- NA

# all genes or a subset
if(TRUE){
  tvar <- res
  fsufix <- paste0(newlabs[gr], '_', compid, "_all_genes")
}else{
  tvar <- res[res$gene_name %in% genes, ]
  fsufix <- paste0(newlabs[gr], '_', compid, "")
}
source("~/scripts/functions/volcano_variant_color.R")
# tvar$log2FoldChange[tvar$log2FoldChange >= 4] <- 4
# tvar$log2FoldChange[tvar$log2FoldChange <= -4] <- -4
# tvar$padj[-log10(tvar$padj) >= 25] <- 10**(-(25))
# tvar <- tvar[tvar$FDR > 25, ]
p <- try(volplot(
  x = tvar,
  pvalth = 0.05,
  lfcth = 0.5,
  pvaltype = 'padj',
  lfctype = 'log2FoldChange',
  col_feature = "mean",
  size_feature = "pct",
  gene_name = "gene_name",
  ngenes = 25,
  return_plot = TRUE,
  v = TRUE
))
pdf(paste0('nvolcano_colours_', fsufix, '_graynolimits.pdf'), width = 10, height = 10)
print(p)
dev.off()
pdf(paste0('nvolcano_colours_', fsufix, '_gray.pdf'), width = 10, height = 10)
print(p + coord_cartesian(ylim = c(0, 25), xlim = c(-4, 4)))
dev.off()
pdf(paste0('nvolcano_colours_', fsufix, '_grayblank.pdf'), width = 10, height = 10)
print(shut_it(p + coord_cartesian(ylim = c(0, 25), xlim = c(-4, 4))))
dev.off()
pdf(paste0('nvolcano_colours_', fsufix, '_gray_freeY.pdf'), width = 10, height = 10)
print(p + coord_cartesian(xlim = c(-4, 4)))
dev.off()
pdf(paste0('nvolcano_colours_', fsufix, '_gray_freeYblank.pdf'), width = 10, height = 10)
print(shut_it(p + coord_cartesian(xlim = c(-4, 4))))
dev.off()
# print(p + scale_x_continuous(breaks = c(-10, -4, -2, -1, 0, 1, 2, 4, 10), limits = c(-10.1, 10.1)) +
#   scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 60, 70)))

# modified
ttvar <- tvar
ttvar$FDR <- trans(ttvar$FDR, 25)
ttvar$log2FoldChange <- trans(ttvar$log2FoldChange, 4, 0.05)
tvar[bordering(tvar, cnames="FDR", 6), c("FDR", "log2FoldChange")]
ttvar[bordering(ttvar, cnames="FDR", 6), c("FDR", "log2FoldChange")]
tvar[bordering(tvar, cnames="log2FoldChange", 6), c("FDR", "log2FoldChange")]
ttvar[bordering(ttvar, cnames="log2FoldChange", 6), c("FDR", "log2FoldChange")]
yticks <- c(0, 5, 10, 15, 20, 25, 60, 70)
xticks <- c(-10, -4, -2, -1, 0, 1, 2, 4, 10)
p <- try(volplot(
  x = ttvar,
  pvalth = -log10(0.05),
  lfcth = 0.5,
  pvaltype = 'FDR',
  lfctype = 'log2FoldChange',
  col_feature = "mean",
  size_feature = "pct",
  gene_name = "gene_name",
  do_fdr_log = FALSE,
  ngenes = 25,
  return_plot = TRUE,
  v = TRUE
)) + scale_x_continuous(breaks = trans(xticks, 4, 0.05), labels = xticks) +
  scale_y_continuous(limits = c(0, NA), breaks = trans(yticks, 25), labels = yticks)
pdf(paste0('nvolcano_colours_', fsufix, '_gray_breaks.pdf'), width = 10, height = 10)
print(p)
dev.off()
pdf(paste0('nvolcano_colours_', fsufix, '_gray_breaksblank.pdf'), width = 10, height = 10)
print(shut_it(p))
dev.off()

### TH2 TSNE PC14 resolution 0.2 ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TSNE, 0 is dark red
cluname <- paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/th2pv_0.1/clustering/zetInfo/clustCells14PCs_30Ks_0.06667JD.RData')
mycells_th2 <- theObjectSavedIn(cluname)
mycells_th2@reductions$tsne <- mycells_th2@reductions$tsne50
mycells_th2@meta.data$Cell_type <- paste0("TH2_", mycells_th2@meta.data$RNA_snn_res.0.2)
p <- DimPlot(mycells_th2, group.by = "Cell_type", cols = v2cols(mycells_th2@meta.data$Cell_type, gr.cols), reduction = "tsne") +
  drawsq + theme(legend.position = "none", axis.text = element_text(face = 'bold', size = 15),
  axis.title = element_text(face = 'bold', size = 15))
pdf('th2_clusters2.pdf')
p
dev.off()
pdf('th2_clusters2blank.pdf')
shut_it(p)
dev.off()

## Disease bars
thesecells <- getsubset(c("orig.diseasegroup", "AS_AL", "AR"), mycells_th2@meta.data, v = TRUE)
annot <- mycells_th2@meta.data[thesecells, ]
source('/mnt/BioHome/ciro/scripts/functions/myplots_functions.R')
p1 <- bar_dis(annot, cnames = c("orig.donor", "orig.diseasegroup", "Cell_type"), cols = gr.cols, filt = 0)
pdf('th2_cluster_disease_bars_norm.pdf', width = 5, height = 10)
lapply(p1, plot)
dev.off()
pdf('th2_cluster_disease_bars_normblank.pdf', width = 5, height = 10)
shut_it(bar_dis(annot, cnames = c("orig.donor", "orig.diseasegroup", "Cell_type"), cols = gr.cols, test = NULL))
dev.off()

p1 <- bar_dis(annot, cnames = c("orig.donor", "Cell_type", "orig.diseasegroup"), cols = gr.cols, filt = 0)
pdf('th2_disease_cluster_bars_norm.pdf', width = 5, height = 10)
lapply(p1, plot)
dev.off()

### 3.X TH2 subcluster disease composition ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resolut <- 'RNA_snn_res.0.2'
annot <- mycells_th2@meta.data[getsubset(c("orig.diseasegroup", "AS_AL", "AR"), mycells_th2@meta.data, v = TRUE), ]
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
grps <- make_list(mycells_th2@meta.data[thesecells, ], colname = "orig.diseasegroup")
p <- lapply(grps, function(x) DimPlot(mycells_th2, cells = x, group.by = resolut, cols = v2cols(mycells@meta.data[, resolut], gr.cols)) )
pdf(paste0("clusters_tsne_", fname, ".pdf"), width = 12)
plot_grid(plotlist = p, ncol = 2)
dev.off()

### 3.X chilli plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lgenes <- list(
  m1n2 = c("EIF3J", "EIF5B", "CALM3", "FKBP1A", "PTPN11", "ATP13A3", "PSMD13", "UBE2S", "DUSP4"),
  m3 = c("TGFBR3", "IL17RB", "IL1RL1", "CSF2", "IL13", "IL5", "IL9", "GZMB", "PPARG", "ZEB2", "GFI1", "RAB27A", "PLA2G16"),
  m4 = c("ICOS", "IFNAR2", "TNFSF8", "TNFRSF9", "PTGER4", "SDC4", "IL4", "IL10", "IL21", "IL3","IL31", "CFLAR", "BCL2A1"),
  m5 = c("IRF4", "STAT4", "GATA3", "RUNX3", "SATB1", "CREM", "CD109", "SLAMF1", "PDCD1")
)
lgenes <- list(
  mil9 = c("IL5", "IL9", "CSF2", "IL3", "IL17RB", "IL1RL1", "GZMB"),
  mil4 = c("IL4", "IL2", "IL31", "IL21", "CD40LG", "CD200", "TNF"),
  mil9_2 = c("IL5", "IL9", "GZMB", "TGFBR3", "PDCD1", "BTLA")
)[3]
mergegroups <- NULL
annot <- remove.factors(mycells_th2@meta.data)
expdata <- cts2cpm(mycells_th2@assays$RNA@counts)[, rownames(annot)]
get_stat_report(expdata, groups = make_list(annot,colname= "Cell_type", grouping = T), rnames = lgenes[[1]], v = T)
brkies <- list(filname = '%+cells', brks = c(1, 10, 20, 30, 40, 50))
colies <- c("#fffeee", "#ffc100","#ff7400", "#ff0000","#a10000", "#670000")
# colies <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
# brkies <- NULL
extramods <- coord_cartesian(ylim = c(5, 20)) #+ theme(panel.spacing = unit(2, "lines"))
# Use fig2.R at 2.D Chilli plot
