#!/usr/bin/R

############
# Figure 1 #
############

# This script will create plots for figure 1 from bulk data

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')

#### Parameters ####
root <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma/results'
setwd(paste0(root, '/a1_final_figures/Figure_1'))
counts.vis = '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma/raw/expMatrices/bulk/tpm_AS3E4d_4013.csv' # Path to counts for visualisation
counts.raw = '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma/raw/expMatrices/bulk/raw_AS3E4d_4013.tsv'
annot.f = '/mnt/BioHome/ciro/asthma/info/deseq2/AS3E4d_coldata.csv' # Annotation file: samples, batches, conditions, etc
gcolour = '/mnt/BioHome/ciro/asthma/info/AS3_final_groupColours.csv' # path to colors for groups in conditions (default='', random colors)
gorder = '/mnt/BioHome/ciro/asthma/info/AS3_final_groupOrder.csv' # Groups heatmap, x_order and y_order columns (default='', group size)

expdata <- read.csv(counts.vis, row.names=1, header=1, stringsAsFactors=F, check.names=F)
counts_ <- read.csv(counts.raw,sep='\t', row.names=1, header=1, stringsAsFactors=F, check.names=F)
coldata_ <- read.csv(annot.f, row.names=2,stringsAsFactors=F, check.names=F)
coldata_ <- data.frame(apply(coldata_,c(1,2),sub,pattern = '^$',replacement = 'no.class'),stringsAsFactors=F,check.names=F)
gr.cols <- read.csv(gcolour,row.names=1,stringsAsFactors=F,header=1)
gr.order <- read.csv(gorder,stringsAsFactors=F,header=1)

### Supplementary table 2 cell type DEGs ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname <- paste0(root, '/Bulk_RNASeq_Analysis/Asthma3_HDM_specific/Hg19/AS3E4_individual_tech_dups/Cell_type/summary/padj0.01_FC2/Cell_type_SummaryDEGsTable_BetasMLE_SSet_suas.csv')
res <- readfile(fname, stringsAsFactors = FALSE, check.name = FALSE)
rownames(res) <- sub("'", "", res$gene_name)
str(res[, 1:22])
mystats <- get_stat_report(
  mat = expdata,
  groups = make_list(coldata_, "Cell_type", grouping = TRUE),
  moments = c("mn", "p"),
  expr_cutoff = 10,
  v = TRUE
)
str(mystats)
head(mystats[rownames(res), ])
head(res[, grepl("_mean|expr", colnames(res))])
mytab <- summ_tables(ltabs = list(res[, 1:21], mystats), v = TRUE)
head(mytab)
write.csv(mytab, file = "bulk_celltypes_stats.csv")

### 1.D Heatmap plot ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname <- paste0(root, '/Bulk_RNASeq_Analysis/Asthma3_HDM_specific/Hg19/AS3E4_individual_tech_dups/Cell_type/summary/padj0.01_FC2/Cell_type_SummaryDEGsTable_BetasMLE_SSet_suas.csv')
res <- readfile(fname, stringsAsFactors = FALSE, check.name = FALSE)
samples <- rownames(coldata_[with(coldata_, order(Cell_type, Disease)), ])
unique(res$group)
genes <- gsub("'", "", res[!is.na(res$group), ]$gene_name)

## NMF
library(NMF)
titl <- "Group-specific genes: Cell type"
matex <- expdata[genes, samples]
# matcolms <- data.frame(Group = grps)
annoc <- coldata_[samples, c("Disease", "Cell_type")]
annog <- res[!is.na(res$group), c("group"), drop = FALSE]; rownames(annog) <- genes
anncolist <- lapply(coldata_[samples, c("Disease", "Cell_type"), drop = F], function(x){
  v2cols(select = x, sour = gr.cols, fw = "gg", v = F)
})
anncolist$group <- v2cols(select = unique(res[!is.na(res$group), ]$group), sour = gr.cols, fw = "gg", v = F)

# manually scale
maxz <- 1.5
matexz <- t(scale(t(matex)))
matexz[matexz < (-maxz)] <- -maxz
matexz[matexz > (maxz)] <- maxz

# pdf(paste0('heatmap_celltype_genes.pdf'), width = 10, height = 10, onefile = FALSE)
# aheatmap(matexz, annCol = annoc, annRow = annog, main = titl, annColors = anncolist,
#   scale = 'none', Rowv = NA, Colv = NA, distfun = 'spearman', subsetCol = NULL,
#   col = rev(colorRampPalette(colors = c('yellow', 'black', 'blue'))(20)))
# graphics.off()

## pheatmap
library(pheatmap)
palettebreaks <- seq(-maxz, maxz, 0.1)
mypalette <- colorRampPalette(c('blue', 'black', 'yellow'), space = 'rgb')
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
  main              = titl,
  fontsize          = 10,
  fontsize_col      = 5,
  fontsize_number   = 5,
  annotation_col    = annoc,
  annotation_row    = annog,
  annotation_colors = anncolist,
  annotation_legend = T,
  annotation_names_col = F,
  annotation_names_row = F,
  drop_levels       = TRUE,
  gaps_col          = cumsum(table(annoc$Cell_type)), # gaps per type
  gaps_row          = cumsum(table(annog$group)[unique(annog$group)]),
  filename          = paste0('heatmapp_celltype_genes_', maxz, '.pdf'),
  width = 10, height = 10
)
file.remove('Rplots.pdf')

### 1.D Heatmap plot sheet ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Adding the rest of the genes
dim(res)
dim(expdata)
genes2add <- rownames(expdata)[!rownames(expdata) %in% gsub("'", "", res$gene_name) ]
resrest <- data.frame(mat_names(genes2add, colnames(res)), stringsAsFactors = FALSE, check.names = FALSE)
resrest$gene_name <- paste0("'", genes2add)
samples <- colnames(resrest)[-c(1:33)]
all(samples %in% rownames(coldata_)) # are all samples in the meta data?
group_list <- make_list(x = coldata_, colname = "Cell_type", grouping = TRUE)
stattab <- get_stat_report(mat = expdata, groups = group_list, rnames = genes2add, v = TRUE)
# Proving if calculations are consistent
stattabdone <- get_stat_report(mat = expdata, groups = group_list, rnames = gsub("'", "", res$gene_name), v = TRUE)
head(stattabdone[gsub("'", "", res$gene_name), ])
head(res[, colnames(stattabdone)])
takencols <- colnames(resrest)[colnames(resrest) %in% colnames(stattab)]
resrest[, takencols] <- stattab[rownames(resrest), takencols]
resrest[, samples] <- expdata[rownames(resrest), samples]
res_all_genes <- rbind(res, resrest)[, -1]
headmat(res_all_genes)
tail(res_all_genes[, 1:10])
write.csv(res_all_genes, file = "heatmapp_celltype_genes_all_info.csv")
write.csv(res_all_genes[, c("group", "gene_name", takencols)], file = "heatmapp_celltype_genes_only_stats_and_group.csv")
sgenes <- rownames(expdata[gsub("'", "", res_all_genes[res_all_genes$Bmean > 1, ]$gene_name), ])
qlucoref <- qlucore_format(mat = expdata, metadata = coldata_[samples, ], rnames = sgenes, v = TRUE)
write.csv(qlucoref, file = "qlucore_file_celltype_genes_mean_gt1.csv")

### 1.X SEM plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 32 genes
coldata_$Treg_dg_tneg <- coldata_$Treg_dg
coldata_$Treg_dg_tneg[coldata_$Cell_type == 'TNeg'] <- 'TNeg'
lgenes <- list(sems = c('IL2', 'IL13', 'IL5', 'IL4', 'IL9', 'IL31', 'IL17F', 'IL22', 'TNF', 'IFNG', 'CSF2', 'CCL20', 'CXCL10'))

lgenes <- list(
  th_treg = c("TNFRSF4", "WARS", "LTA", "BCL2L1", "PSMD11", "CD70", "EIF5A2", "FASLG", "EIF2C3", "PSMB2", "E2F6", "FARS2"),
  th = c("MIR155HG", "TNF", "CSF2", "IL5", "IFNG", "IL17A", "IL17F", "BCL2A1", "IL2", "IL4", "IL13", "CD40LG"),
  treg = c("IL2RA", "FOXP3", "CTLA4", "TNFRSF1B", "IL1R2", "HLA-DRA", "CORO2A", "PRDM1", "IKZF2", "TNFRSF8", "LRRC32", "TIGIT"),
  th_fig_g_h = c("IL5", "IL13", "IL4", "IL9", "IL31", "IL1RL1", "IL17RB", "IL17F", "IFI6", "MX1", "IFI44", "IFNG")
)[4]
# lgenes <- list(th_treg = unlist(lgenes[1:3])) # height = 18
lgenes <- list(treg_diseaseDEGS = head(bordering(mydatafc, 1:2), Inf))
lgenes <- list(treg_diseaseDEGS = head(ifis, Inf))
lgenes <- list(treg_diseaseDEGS = head(bordering(mydatafc, 1:2, 3), Inf))
dir.create('sem_plots')
for(gname in names(lgenes)){
  cat(gname, "\n")
  fname <- paste0('sem_plots/', gname)
  if(1){ # For disease differences
    fname <- paste0(fname, '_disease')
    ssamples <- list(c('Cell_type', 'TNeg', switch(gname, th = 'Teff', treg = 'Treg',
      th_treg = c('Teff', 'Treg'), th_fig_g_h = c('Teff'), treg_diseaseDEGS = 'Treg')))
    if(gname == 'th_fig_g_h') ssamples <- list(ssamples[[1]][-2])
    cname <- switch(gname, th = 'Teff_dg_tneg', treg = 'Treg_dg_tneg', th_treg = 'Disease',
      th_fig_g_h = 'Disease', treg_diseaseDEGS = 'Disease')
    grps <- c('TNeg', 'AS_AL', 'AR', 'AS_NA', 'HC')
  }else{ # For cell types differences
    fname <- paste0(fname, '_celltype')
    ssamples <- list(c('Cell_type', 'Teff', 'TNeg', 'Treg'))
    cname <- 'Cell_type'
    grps <- c('TNeg', 'Teff', 'Treg')
  }

  genes <- getfound(lgenes[[gname]], rownames(expdata), v = TRUE)
  # source('/mnt/BioHome/ciro/scripts/functions/myplots_functions.R')
  p1 <- pSEM(
    mygene = genes,
    edata = expdata,
    metadata = coldata_[getsubset(ssamples, coldata_, v = T), ],
    group_by = cname,
    colour_by = cname,
    mylevels = grps,
    sourcols = gr.cols,
    ptsize = 3,
    group_name = NULL,
    ctstype = "TPM",
    legendpos = "right",
    return_plot = TRUE,
    v = TRUE
  )
  pdf(paste0(fname, length(genes), '.pdf'), height = 13, width = 14)
  print(p1)
  dev.off()
  pdf(paste0(fname, length(genes), 'blank.pdf'), height = 13, width = 14)
  print(shut_it(p1))
  dev.off()
}
# pdf(paste0(fname, '_individual.pdf'), height = 7, width = 7)
# for(gg in genes){
#   p1 <- pSEM(
#     mygene = gg,
#     edata = expdata,
#     metadata = coldata_[getsubset(ssamples, coldata_, v = T), ],
#     group_by = cname,
#     colour_by = cname,
#     mylevels = grps,
#     sourcols = gr.cols,
#     group_name = NULL,
#     ctstype = "TPM",
#     legendpos = "right",
#     return_plot = TRUE,
#     v = TRUE
#   )
#   print(p1)
#   print(p1 + shut_up)
# }
# dev.off()

### Table 2 Teff disease differences ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fstat <- list(
  Asthma = data.frame(theObjectSavedIn(paste0(root, '/Bulk_RNASeq_Analysis/Asthma3_HDM_specific/Hg19/AS3E4_individual_tech_dups/Teff_asthma/ASvsNAS/ASvsNAS_DEGs_BetasMLE_SSet.RData'))),
  Allergy = data.frame(theObjectSavedIn(paste0(root, '/Bulk_RNASeq_Analysis/Asthma3_HDM_specific/Hg19/AS3E4_individual_tech_dups/Teff_allergy/ALvsNAL/ALvsNAL_DEGs_BetasMLE_SSet.RData')))
)
sum(!is.na(fstat[['Asthma']]$padj))
fstat <- lapply(names(fstat), function(x){
  y <- fstat[[x]][, c("log2FoldChange", "padj")]
  colnames(y) <- paste0(x, "_", colnames(y))
  y$gene_name <- paste0("'", rownames(y))
  y
})
str(fstat)
str(coldata_)
mystats <- lapply(c("Allergy", "Asthma", "Teff_dg"), function(x){
  get_stat_report(
    mat = expdata,
    groups = make_list(coldata_[getsubset(c("Cell_type", "Teff"), coldata_, v = TRUE), ], x, grouping = TRUE),
    moments = c("mn", "p"),
    expr_cutoff = 10,
    v = TRUE
  )
})
str(mystats)
mytab <- summ_tables(ltabs = c(fstat, mystats), v = TRUE)
mytab <- mytab[order(-mytab$Allergy_log2FoldChange), ]
head(mytab)
sum(!is.na(mytab$Asthma_padj))
sum(!is.na(fstat[[1]]$Asthma_padj))

### 1.G Combined volcano ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thesecells <- getsubset(c("Cell_type", "Teff"), coldata_, v = TRUE)
thesecellst <- rownames(coldata_[order(coldata_[, "Disease"]), ])
thesecells <- thesecellst[thesecellst %in% thesecells]
library(DESeq2) # may need R's DESeq2 and not R5's
fstat <- list(
  Asthma = data.frame(theObjectSavedIn(paste0(root, '/Bulk_RNASeq_Analysis/Asthma3_HDM_specific/Hg19/AS3E4_individual_tech_dups/Teff_asthma/ASvsNAS/ASvsNAS_DEGs_BetasMLE_SSet.RData'))),
  Allergy = data.frame(theObjectSavedIn(paste0(root, '/Bulk_RNASeq_Analysis/Asthma3_HDM_specific/Hg19/AS3E4_individual_tech_dups/Teff_allergy/ALvsNAL/ALvsNAL_DEGs_BetasMLE_SSet.RData')))
)
myfcs <- sapply(fstat, "[", 2)
mydatafc <- data.frame(myfcs, stringsAsFactors = F, row.names = rownames(fstat[[1]]))
names(mydatafc) <- names(myfcs)
names(mydatafc) <- gsub(".log2FoldChange", "", names(myfcs))
mydatafc$padj <- apply(data.frame(sapply(fstat, "[", 7)), 1, min, na.rm = TRUE)
head(mydatafc)
mydatafc$padj[is.infinite(mydatafc$padj)] <- 1
mydatafc$significance <- -log10(mydatafc$padj)
mydatafc <- mydatafc[complete.cases(mydatafc), ]
mydatafc$Mean <- rowMeans(expdata[rownames(mydatafc), thesecells])
thesegenes <- rownames(mytab[apply(mytab[, grep("exprFrac", colnames(mytab))[1:4]], 1, max) > 50, ])
tvar <- mydatafc$Mean > 10 & rownames(mydatafc) %in% thesegenes
sum(tvar)
mydatafc <- mydatafc[tvar, ]## Genes over 10 TPM of mean expression
mydatafc$gene_name <- rownames(mydatafc)
mydatafc$Log2_Mean <- log2(mydatafc$Mean)
mydatafc$Log2_Mean[mydatafc$Log2_Mean > 5] <- 5
# mydatafc_repel <- mydatafc[rowSums(abs(mydatafc[, 1:2]) >= 1.2, na.rm = T) > 0, ]
# mydatafc <- mydatafc[rowSums(expdata[genes, ] > 1, na.rm=T) > (ncol(expdata) / 4), ]
borgenes <- unique(c(rownames(mydatafc)[grepl("^IFI", rownames(mydatafc))], bordering(mydatafc, 1:2)))
mydatafc_repel <- mydatafc[borgenes, ]
thres <- 1 # log2FoldChange threshold
p <- ggplot(mydatafc, aes(x = Asthma, y = Allergy, color = Log2_Mean, size = significance)) +
  geom_point() + scale_color_gradientn(colours = c('#ffffff', '#670000')) + ##CB4154
  geom_hline(yintercept = c(-thres, thres), linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = c(-thres, thres), linetype = "dashed", alpha = 0.4) +
  labs(title = 'Allergy - Asthma differences',
    subtitle = paste('Threshold: FC >', thres), x = 'Asthma', y = 'Allergy',
    color = 'Expression', size = 'Significance') +
  drawsq + scale_size_area() + scale_size(range = c(0, 10))
p <- squareplot(p)
fname <- paste0('fc_volcano_gene_FC_filtered50perCompGroup_', thres)
if(1){
  p <- p + guides(size = guide_legend(keywidth = 0.5, keyheight = 0.5, default.unit = "inch"))
  fname <- paste0(fname, '_update')
}; set.seed(27)
pdf(paste0(fname, '.pdf'), 10, 10)
print(p + geom_text_repel(data = mydatafc_repel, aes(label = gene_name), color = 'black'))
dev.off()
pdf(paste0(fname, 'blank.pdf'), 10, 10)
print(p + shut_up)
dev.off()

smytab <- mytab[rownames(mytab) %in% rownames(mydatafc), ]
smytab <- smytab[order(-smytab$Allergy_log2FoldChange), ]
dim(smytab)
head(smytab)
write.csv(smytab, file = "bulk_teff_disease_stats.csv")

### 1.F Teff t-SNE plot ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# for Teff, TNeg
# selectss <-  c('Cell_type',c('Teff','Treg','TNeg')[c(1,3)])
# pcs <- 3
# pxs <- 10
# More combinations without TNegs
selectss <-  c('Cell_type',c('Teff','Treg','TNeg')[c(1)])
pcs <- seq(3, 15, 2)
pxs <- seq(2, 15, 2)

samples <- getsubset(selectss, coldata_, v = T)
keep_genes <- mostVars(counts_[, samples], 200, "AR")
keep_genes <- mostVars(counts_[, samples], 300, "AR")
cv2 <- keep_genesl$vars/keep_genesl$means^2
keep_genes <- keep_genesl$genes
# ordgenes <- sort(keep_genesl$vars, decreasing = TRUE)
# which(names(ordgenes) %in% keep_genesl$genes)
file.remove('Rplots.pdf')
cts <- log2(expdata[keep_genes, samples] + 1)
# summary(rowMeans(expdata[keep_genes, samples])) # min mean = 2.582, min pct = 53.70

pca.m <- prcomp(cts, center = TRUE)
pca.mat <- data.frame(pca.m$rotation, stringsAsFactors = FALSE, check.names = FALSE)

dir.create('teff_tsnes')
for(pc in pcs){
  for(px in pxs){
    fname <- paste0('teff_tsnes_', paste0(pc, "PCs_", px), 'PX_', length(keep_genes), 'genes.pdf')
    cat(fname, '\n')
    set.seed(27)
    tsne_data <- Rtsne::Rtsne(dist(pca.mat[, 1:pc]), perplexity = px, pca = FALSE, theta = 0, is_distance = TRUE, max_iter = 1000)
    ddf <- data.frame(x.var = tsne_data$Y[, 1], y.var = tsne_data$Y[,2])
    rownames(ddf) <- rownames(pca.mat)
    ddf$col<-  factormix(coldata_[samples, 'Teff_dg_tneg'])

    gg.d <- ggplot(ddf, aes(x.var, y.var) ) +
        geom_point(size=5, shape=19, aes( colour=col)) +
        scale_colour_manual(values = v2cols(levels(ddf$col), gr.cols)) + theme_classic() +
        theme(panel.background = element_blank(), panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), plot.background = element_blank(),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          legend.position = "right", legend.text = element_text(size = 10)) +
        labs( x= "tSNE 1", y= "tSNE 2", colour = "Group") +
        drawsq

    gg.d <- squareplot(gg.d)
    pdf(fname)
    print(gg.d)
    dev.off()
    pdf(sub("PX.pdf", "PXblank.pdf", fname))
    print(gg.d + shut_up)
    dev.off()
  }
}

### 1.C PCA/t-SNE plot ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selectss <-  c('Cell_type',c('Teff','Treg','TNeg')[c(1:3)])
samples <- getsubset(selectss, coldata_, v = T)
keep_genes <- mostVars(counts_[, samples], 200, "AR")
file.remove('Rplots.pdf')
length(keep_genes)
cts <- log2(expdata[keep_genes, samples] + 1)
pca.m <- prcomp(cts, center = TRUE)
pca.mat <- data.frame(pca.m$rotation, stringsAsFactors = FALSE, check.names = FALSE)

set.seed(27)
tsne_data <- Rtsne::Rtsne(dist(pca.mat[, 1:5]), perplexity = 30, pca = FALSE, theta = 0, is_distance = TRUE, max_iter = 1000)

ddf <- data.frame(x.var = tsne_data$Y[, 1], y.var = tsne_data$Y[,2])
rownames(ddf) <- rownames(pca.mat)
ddf$col<-  factormix(coldata_[samples, 'Cell_type'])

gg.d <- ggplot(ddf, aes(x.var, y.var) ) +
    geom_point(size=5, shape=19, aes( colour=col)) +
    scale_colour_manual(values = v2cols(levels(ddf$col), gr.cols)) + theme_classic() +
    theme(panel.background = element_blank(), panel.border = element_blank(),panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), plot.background = element_blank(),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.position = "right", legend.text = element_text(size = 10)) +
    labs( x= "tSNE 1", y= "tSNE 2", colour = "Group") +
    drawsq

gg.d <- squareplot(gg.d)
pdf('tsne_celltypes.pdf')
gg.d
dev.off()
pdf('tsne_celltypesblank.pdf')
gg.d + shut_up
dev.off()

pca.mat$col<-  factormix(coldata_[samples, 'Cell_type'])
# pca.mat$col<-  factormix(coldata_[samples, 'Disease'])
p <- ggplot(pca.mat, aes(x = PC1, PC2)) +
  geom_point(size=5, shape=19, aes( colour=col)) +
  scale_colour_manual(values = v2cols(levels(pca.mat$col), gr.cols)) + theme_classic() +
  theme(panel.background = element_blank(), panel.border = element_blank(),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), plot.background = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right", legend.text = element_text(size = 10)) +
  labs( x= "PC 1", y= "PC 2", colour = "Group") +
  drawsq
# p <- squareplot(p)
pdf('pca.pdf')
p
dev.off()
pdf('pcablank.pdf')
p + shut_up
dev.off()

# pdf('pca_func.pdf')
# plotPCA(sca_obj = cts[keep_genes, ], coldata = coldata_[samples, ], condition = "Cell_type", colours = gr.cols)
# # pcaplot(x = cts[keep_genes, ], main.f = "Cell types", coldata = coldata_[samples, ],
# #   batchCol = "donor_id", cd = "Cell_type", grcols = gr.cols, v = TRUE)
# dev.off()

# pdf('pca_3d.pdf')
# pca3d(pca = pca.m, group = gr.cols[pca.mat$col, ], radius = 3)
# dev.off()

### 1.H Scatter plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir.create('scatters')
genes <- c("IL5", "IL4", "IL13", "IL17A")
ssamples <- c('Cell_type', 'Teff', 'TNeg', 'Treg')[1:2]
samples <- getsubset(ssamples, coldata_, v = T)
ddf <- remove.factors(melt(t(expdata[genes, samples])))
colnames(ddf)[3] <- "TPM"
ddf <- cbind(ddf, t(expdata["IFNG", ddf[, 1]]))
ddf <- cbind(ddf, coldata_[ddf[, 1], , drop = FALSE])
head(ddf)
ddf$IFNG <- log2(ddf$IFNG + 1)
ddf$TPM <- log2(ddf$TPM + 1)
p <- ggplot(ddf, aes(x = IFNG, y = TPM, color = Disease)) +
  geom_jitter(size = 3) + facet_wrap(~ Var2, scales = 'fixed') +
  scale_colour_manual(values = v2cols(ddf$Disease, gr.cols)) +
  labs(x = "IFNG", y = "log2(TPM + 1)") +
  drawsq
pdf(paste0('scatters/', paste0(ssamples[-1], collapse = "_"),'.pdf'))
print(p)
dev.off()
pdf(paste0('scatters/', paste0(ssamples[-1], collapse = "_"),'blank.pdf'))
print(p + shut_up)
dev.off()

# Individual
ddf <- data.frame(t(expdata[c("IFNG", genes), samples]))
ddf <- log2(ddf + 1)
axmax <- max(abs(ddf))
ddf <- cbind(samples = samples, ddf)
ddf <- cbind(ddf, coldata_[ddf[, 1], , drop = FALSE])
for(gg in genes[2]){
  cat(gg, "\n")
  p <- ggplot(ddf, aes_string(x = "IFNG", y = gg, color = "Disease")) +
    geom_jitter(size = 7) +
    scale_colour_manual(values = v2cols(ddf$Disease, gr.cols)) +
    labs(x = "IFNG", y = "log2(TPM + 1)") +
    drawsq
  p <- p + xlim(c(0, axmax)) + ylim(c(0, axmax))
  pdf(paste0('scatters/', gg,'.pdf'))
  print(p)
  dev.off()
  pdf(paste0('scatters/', gg,'blank.pdf'))
  print(p + shut_up)
  dev.off()
}

### Bulk samples ### -----------------------------------------------------------
raws <- readfile("/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma/raw/expMatrices/bulk/raw_AS3E4d_4013.tsv", header = 1, stringsAsFactors =FALSE, check.names = F)
eannot <- readfile("/mnt/BioHome/ciro/asthma/info/deseq2/AS3E4d_coldata.csv",  stringsAsFactors =FALSE, check.names = F)

colnames(raws)[!colnames(raws) %in% eannot[, "sample_name"]]
