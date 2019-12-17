#!/usr/bin/R5

####################
# Extended figures #
####################

# This script will create plots for Extended Figures

.libPaths('~/R/newer_packs_library/3.5/')
library(Seurat)
source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')

## Destiny Folder ##
root <- '/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma/results'
setwd(paste0(root, '/a1_final_figures/Extended_figures'))

gr.cols <- read.csv("/mnt/BioHome/ciro/asthma/info/AS3_final_groupColours.csv",row.names=1,stringsAsFactors=F,header=1)
mycells <- theObjectSavedIn(paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/clustering/zetInfo/clustCells13PCs_30Ks_0.06667JD.RData'))

newlabs <- c("ACT1", "ACT2", "ACT3", "TH1", "TH2", "THIFNR", "TH17")
names(newlabs) <- 0:6
gorder <- c("TH2", "TH1", "TH17", "THIFNR", "ACT3", "ACT2", "ACT1")
mycells@meta.data$Cell_type <- factor(c(newlabs[as.character(mycells@meta.data$RNA_snn_res.0.4)]), gorder)

mycells_treg <- theObjectSavedIn(paste0(root, '/clustering_seurat/AS3EsmTreg_daf_30p_0.1m/clustering/zetInfo/clustCells20PCs_30Ks_0.06667JD.RData'))
mycells_treg@reductions$tsne <- mycells_treg@reductions$tsne200
mycells_treg@meta.data$Cell_type <- paste0("TREG_", mycells_treg@meta.data$RNA_snn_res.0.2)

cluname <- paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/th2pv_0.1/clustering/zetInfo/clustCells14PCs_30Ks_0.06667JD.RData')
mycells_th2 <- theObjectSavedIn(cluname)
mycells_th2@reductions$tsne <- mycells_th2@reductions$tsne50
mycells_th2@meta.data$Cell_type <- paste0("TH2_", mycells_th2@meta.data$RNA_snn_res.0.2)

## Doublets
annot <- theObjectSavedIn('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma/results/demuxlet/asthma_libraries/reports/aggregated_cells_deconvoluted.RData')
mycellsa <- theObjectSavedIn('/mnt/BioAdHoc/Groups/vd-vijay/cramirez/asthma/results/clustering_seurat/AS3Esm/clustering/zetInfo/clustCells50PCs_30Ks_0.06667JD.RData')

thesecols <- v2cols(names(table(mycellsa@meta.data$orig.celltype)), gr.cols)
# DimPlot(mycellsa, cells = thesecells, reduction = "tsne", group.by = "celltype")
p <- DimPlot(mycellsa, reduction = "tsne", group.by = "orig.celltype")+
  scale_color_manual(values = thesecols,
    labels = c(expression("HDM"^"-"), expression("HDM"^"+ "*T[H]), expression("HDM"^"+ "*T[REG]))) +
  guides(colour = guide_legend(override.aes = list(size = 6)))
pdf("all_cells.pdf", width = 8, height = 7)
set_big_text(p, sizey = 24)
dev.off()

dim(mycellsa@meta.data)
dim(annot)
# mycellsa@meta.data <- cbind_repcol(mycellsa@meta.data, annot) # it already have DECONVOLUTION
mycellsa@meta.data <- annot[rownames(mycellsa@meta.data), ]
# mycellsa@meta.data$orig.class <- factor(mycellsa@meta.data$orig.class)
mycellsa@meta.data$orig.class[mycellsa@meta.data$orig.class == "void"] <- "Filtered"
table(mycellsa@meta.data$orig.class, useNA = 'always')

couls <- c("blue", "red", "orange", "gray")
couls <- c("gray", "red", "gray", "gray")
names(couls) <- names(table(mycellsa@meta.data$orig.class))
p <- DimPlot(mycellsa, reduction = "tsne", group.by = "orig.class", cols = couls, pt.size = 0.3)
# ddf <- FetchData(mycellsa, vars = c("tSNE_1", "tSNE_2", "orig.class"))
# p <- ggplot(ddf, aes(tSNE_1, tSNE_2, color = orig.class))
# p <- p + geom_point(data = ddf[ddf$orig.class == "DBL", ], aes(tSNE_1, tSNE_2), color = "red", size = 0.3)
pdf("demuxlet_classification.pdf", width = 8, height = 7)
set_big_text(p, sizey = 24)
dev.off()
pdf("demuxlet_classificationblank.pdf")
print(shut_it(p))
dev.off()

mycellsa@meta.data$ct_disease <- do.call('paste', c(mycellsa@meta.data[, c("orig.celltype", "orig.diseasegroup", "orig.experiment", "gem")], sep = "_"))
ddf <- melt(table(mycellsa@meta.data$ct_disease, mycellsa@meta.data$origlib))
ddf[ddf$value > 0, ]
cnames <- rev(c("orig.diseasegroup", "orig.celltype", "orig.experiment", "origlib", "ct_disease"))
mytabs <- rbindlist(lapply(cnames, function(x){
  cat(x, "\n")
  mytb <- table_pct(mycellsa@meta.data, cnames = c(x, "orig.class"))
  mytb <- data.frame(cbind(Group = rownames(mytb), mytb)[, 1:(ncol(mytb) + 1)])
  if(tail(cnames, 1) == x) mytb else mytb[-nrow(mytb), ]
}))
pdf("demuxlet_proportions.pdf", width = 8, height = 15)
plot(tableGrob(mytabs))
dev.off()
write.csv(mytabs, file = "demuxlet_proportions.csv", row.names = FALSE)

### Activation signature ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annot <- mycellsa@meta.data
expdata <- as.matrix(expm1(mycellsa@assays$RNA@data[, rownames(annot)]) * 100)
allstats <- get_stat_report(
  mat = expdata,
  groups = make_list(x = annot, colname = "orig.celltype", grouping = TRUE),
  # rnames = genes,
  moments = c("mn", "p"),
  v = TRUE
)
# thesecells <- getsubset(c("orig.celltype", "TNeg", celltype), annot, v = TRUE)

celltype <- "Teff"
fname <- list.files(paste0(root, "/mast/AS3EsmTeff_sng_amb/comprs"), pattern = paste0("results_", celltype, "vsTNeg_mastlog2cpm.csv"), recursive = TRUE, full.name = TRUE)
sgenes <- as.character(read.csv(paste0(dirname(fname), "/Teff_oActivation_110g/used_genes.csv"))[, 2])
titl <- expression("HDM"^"+ "*T[H]*" vs. HDM"^"-"*" cells")

celltype <- "Treg"
fname <- list.files(paste0(root, "/mast/AS3EsmTeff_sng_amb/comprs"), pattern = paste0("results_", celltype, "vsTNeg_mastlog2cpm.csv"), recursive = TRUE, full.name = TRUE)
sgenes <- as.character(read.csv(paste0(dirname(fname), "/signature_genes/activatedTreg_64g/a1_used_genes.csv"))[, 2])
titl <- expression("HDM"^"+ "*T[REG]*" vs. HDM"^"-"*" cells")

res <- readfile(fname, stringsAsFactors = FALSE)[, -1]#[, 2:8]
rownames(res) <- sub("'", "", res$gene_name)
colnames(res)[7:ncol(res)] <- paste0(colnames(res)[7:ncol(res)], "_dea")
headmat(res)
genes <- getfound(rownames(res), rownames(expdata), v = TRUE)

modres <- cbind(res[genes, ], allstats[genes, ])
modres$lmean <- log2(modres[, paste0(celltype, "_mean")] + 1)
modres$logpadj <- -log10(modres[, "padj"])
modres$logpadj[is.infinite(modres$logpadj)] <- max(modres[is.finite(modres$logpadj), ]$logpadj)
headmat(modres); tailmat(modres)
tvar <- paste0(c(celltype, "TNeg"), "_mean")
modres$means <- round(log2(ifelse(modres$group == celltype, modres[, tvar[1]], modres[, tvar[2]]) + 1), 1)
tvar <- paste0(c(celltype, "TNeg"), "_exprFrac")
modres$pcts <- round(ifelse(modres$group == celltype, modres[, tvar[1]], modres[, tvar[2]]), 1)
headmat(modres); tailmat(modres)
mysignames <- getDEGenes( # not same result, Nov 12 2019
  modres, pv = 1e-50, fc = 2,
  further = paste0(celltype, "_meanlog2.CPM...1._dea >= 0.75 & ", celltype, "_exprFrac_dea >= 37.5"),
  # further = paste0(celltype, "_mean >= 0.75 & ", celltype, "_exprFrac >= 37.5"),
  # further = paste0("lmean > 7.878 & ", celltype, "_exprFrac >= 37.5"),
  v = TRUE
)
# summary(modres[sgenes, "lmean"])
# summary(modres[sgenes, ])

source('/mnt/BioHome/ciro/scripts/functions/volcano_variant_color.R')
modres$Mean <- modres[, paste0(celltype, '_meanlog2.CPM...1._dea')]
modres$Percentage <- modres[, paste0(celltype, "_exprFrac_dea")]
modres[which.max(modres$Percentage), ]$Percentage <- 100
modres$Significance <- ifelse(modres$logpadj > 300, 300, modres$logpadj)
modres[!rownames(modres) %in% sgenes, ]$Percentage <- NA
p <- volplot(
  x = modres,
  pvalth = 0.5,
  lfcth = 2,
  pvaltype = 'padj',
  lfctype = 'log2FoldChange',
  col_feature = 'Percentage',
  # check_genes = sgenes,
  ngenes = 0,
  grepel = FALSE,
  return_plot = TRUE,
  v = TRUE
)
p$labels$subtitle <- sub(" \\n.*", "", p$labels$subtitle)
coulsy <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','red2','#b30000', '#670000')
p <- p + geom_point(data = p$data[!is.na(p$data$Percentage), ], aes(colour = Percentage)) +
  scale_color_gradientn(name = "Percentage", colours = coulsy, limits = c(0, 100))
p <- p + labs(title = titl, subtitle = NULL, y = "Significance", x = expression(log[2]*" fold change"))
pdf(paste0("signature_padj_lfc_", celltype, ".pdf"), 8, 8)
set_big_text(p, sizey = 24) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
dev.off()

modres$Percentage <- modres[, paste0(celltype, "_exprFrac_dea")]
modres[!rownames(modres) %in% sgenes, ]$Significance <- NA
p <- volplot(
  x = modres,
  pvalth = 37.5,
  lfcth = 0.75,
  pvaltype = 'Percentage',
  do_fdr_log = FALSE,
  lfctype = 'Mean',
  col_feature = "Significance",
  # check_genes = sgenes,
  ngenes = 0,
  grepel = FALSE,
  return_plot = TRUE,
  v = TRUE
)
p <- p + geom_point(data = p$data[!is.na(p$data$Significance), ], aes(colour = Significance))
p <- p + labs(title = titl, subtitle = NULL)
p$labels$subtitle <- sub(" \\n.*", "", p$labels$subtitle)
pdf(paste0("signature_pct_mean_", celltype, ".pdf"), 8, 8)
set_big_text(p, sizey = 24) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
dev.off()

### Activation ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thesecells <- getsubset(c("orig.celltype", "-TNeg"), mycellsa@meta.data, v = TRUE)
# pdf("filtered_before_cells.pdf", 10, 10)
# DimPlot(mycellsa, cells = thesecells, reduction = "tsne", group.by = "orig.celltype")
# dev.off()
thesecells <- getfound(c(colnames(mycells), colnames(mycells_treg)), colnames(mycellsa), v = TRUE)
mycellsa@meta.data$celltype <- mycellsa@meta.data$orig.celltype
mycellsa@meta.data$celltype[!rownames(mycellsa@meta.data) %in% thesecells] <- "Filtered"
thesecols <- v2cols(names(table(mycellsa@meta.data$celltype)), gr.cols)
thesecols['Filtered'] <- "#BEBEBE"
# DimPlot(mycellsa, cells = thesecells, reduction = "tsne", group.by = "celltype")
p <- DimPlot(mycellsa, reduction = "tsne", group.by = "celltype") +
  scale_color_manual(values = thesecols, labels = c("Filtered", expression(T[H]), expression(T[REG]))) +
  guides(colour = guide_legend(override.aes = list(size = 6)))
pdf("filtered_cells.pdf", width = 8, height = 7)
set_big_text(p, sizey = 24)
dev.off()

## Teff activation
celltype <- "Teff"
fname <- list.files(paste0(root, "/mast/AS3EsmTeff_sng_amb/comprs"), pattern = paste0("results_TeffvsTNeg_mastlog2cpm.csv"), recursive = TRUE, full.name = TRUE)
fname <- paste0(dirname(fname), "/Teff_oActivation_110g/Teff_oActivation_110g_signature.csv")
activ <- read.csv(fname, row.names = 1)
all(colnames(mycellsa) %in% rownames(activ))
# ddf <- FetchData(mycellsa, vars = c("tSNE_1", "tSNE_2", "orig.celltype", "orig.diseasegroup"))[rownames(activ), ]
# ddf$Activation <- activ[, 1]
sobject <- mycellsa[, getfound(colnames(mycellsa), rownames(activ), v = TRUE)]
activ$celltype <- sobject@meta.data[rownames(activ), "orig.celltype"]
activ_teff <- activ
sobject@meta.data$Activation <- activ[rownames(sobject@meta.data), 1]
brks <- c("+0.4", "0", "-0.4")
p <- FeaturePlot(sobject, features = "Activation", reduction = "tsne") + labs(title = expression(T[H]*" activation")) +
  scale_color_gradientn(colours = c('#ffdf32', '#ffe311', 'red2', 'red3', '#a30404', '#710117'),
    breaks = as.numeric(brks), labels = brks)
pdf("teff_activation.pdf", 8, 8)
set_big_text(p, sizey = 24)
dev.off()
pdf("teff_activation2.pdf", 8, 8)
set_big_text(p + scale_color_gradientn(colours = c('#ffdf32', 'red')), sizey = 24)
dev.off()
pdf("teff_activation3.pdf", 8, 8)
set_big_text(p + scale_color_gradientn(colours = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')), sizey = 24)
dev.off()

## Treg activation
celltype <- "Treg"
fname <- list.files(paste0(root, "/mast/AS3EsmTeff_sng_amb/comprs"), pattern = paste0("results_TregvsTNeg_mastlog2cpm.csv"), recursive = TRUE, full.name = TRUE)
list.files(paste0(dirname(fname), "/signature_genes/activatedTreg_64g"))
fname <- paste0(dirname(fname), "/signature_genes/activatedTreg_64g/activatedTreg_64g_signature.csv")
activ <- read.csv(fname, row.names = 1)
all(colnames(mycellsa) %in% rownames(activ))
# ddf <- FetchData(mycellsa, vars = c("tSNE_1", "tSNE_2", "orig.celltype", "orig.diseasegroup"))[rownames(activ), ]
# ddf$Activation <- activ[, 1]
sobject <- mycellsa[, getfound(colnames(mycellsa), rownames(activ), v = TRUE)]
activ$celltype <- sobject@meta.data[rownames(activ), "orig.celltype"]
sobject@meta.data$Activation <- activ[rownames(sobject@meta.data), 1]
brks <- c("+0.8", "+0.4", "0", "-0.4", "-0.8")
p <- FeaturePlot(sobject, features = "Activation", reduction = "tsne") + labs(title = expression(T[REG]*" activation")) +
  scale_color_gradientn(colours = c('#ffdf32', '#ffe311', 'red2', 'red3', '#a30404', '#710117'),
  breaks = as.numeric(brks), labels = brks)
pdf("treg_activation.pdf", 8, 8)
set_big_text(p, sizey = 24)
dev.off()
pdf("treg_activation2.pdf", 8, 8)
set_big_text(p + scale_color_gradientn(colours = c('#ffdf32', 'red')), sizey = 24)
dev.off()
pdf("treg_activation3.pdf", 8, 8)
set_big_text(p + scale_color_gradientn(colours = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')), sizey = 24)
dev.off()

thesecells <- intersect(rownames(activ_teff), rownames(activ))
colnames(activ) <- c("Activation", "celltype")
colnames(activ_teff) <- c("Activation", "celltype")
activ_comb <- rbind(activ_teff[!rownames(activ_teff) %in% thesecells, ], activ[!rownames(activ) %in% thesecells, ])
mymins_df <- activ_teff[thesecells, ]
mymins_df$Activation <- apply(cbind(activ_teff[thesecells, "Activation"], activ[thesecells, "Activation"]), 1, min)
activ_comb <- rbind(activ_comb, mymins_df)
table(activ_comb$celltype)

sobject <- mycellsa[, getfound(colnames(mycellsa), rownames(activ_comb), v = TRUE)]
sobject@meta.data$Activation <- activ_comb[rownames(sobject@meta.data), 1]
brks <- c("+0.8", "+0.4", "0", "-0.4", "-0.8")
p <- FeaturePlot(sobject, features = "Activation", reduction = "tsne") + labs(title = expression("HDM"^"+"*" activation")) +
  scale_color_gradientn(colours = c('#ffdf32', '#ffe311', 'red2', 'red3', '#a30404', '#710117'),
  breaks = as.numeric(brks), labels = brks)
pdf("combined_activation.pdf", 8, 8)
set_big_text(p, sizey = 24)
dev.off()
pdf("combined_activation2.pdf", 8, 8)
set_big_text(p + scale_color_gradientn(colours = c('#ffdf32', 'red')), sizey = 24)
dev.off()
pdf("combined_activation3.pdf", 8, 8)
set_big_text(p + scale_color_gradientn(colours = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')), sizey = 24)
dev.off()

### Clustering parameters ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cluname <- paste0(root, '/clustering_seurat/AS3EsmTeff_daf_30p/th2pv_0.1/clustering/zetInfo/clustCells14PCs_30Ks_0.06667JD.RData')
mycells_th2 <- theObjectSavedIn(cluname)
mycells_th2@reductions$tsne <- mycells_th2@reductions$tsne50
mycells_th2@meta.data$Cell_type <- paste0("TH2_", mycells_th2@meta.data$RNA_snn_res.0.2)

summobj <- list(sobject = mycells, title = expression(T[H] * " elbow"), elbow = 13, name = "teff")
summobj <- list(sobject = mycells_treg, title = expression(T[REG] * " elbow"), elbow = 20, name = "treg")
summobj <- list(sobject = mycells_th2, title = expression(T[H] * "2 elbow"), elbow = 4, name = "th2")

pcs.comp <- length(summobj$sobject@reductions$pca@stdev)
p <- ElbowPlot(object = summobj$sobject, ndims = pcs.comp) +
  ggtitle(summobj$title) +
  geom_vline(xintercept = summobj$elbow, linetype = "dotted", color = "gray")
fname <- paste0('clust_params_', summobj$name, "_set_", summobj$elbow, 'PCs.pdf')
pdf(fname, width = 7, height = 5); set_big_text(p, sizey = 24); dev.off()

# top variable genes`
top10 <- head(x = VariableFeatures(object = summobj$sobject), 20)
plot1 <- VariableFeaturePlot(summobj$sobject)
p <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size = 5) +
  theme(legend.position = c(0, 75), legend.direction = "vertical")
fname <- paste0('hvg_plot_', summobj$name, '.pdf')
pdf(fname, width = 7, height = 5.5); set_big_text(p, sizey = 24); dev.off()

# cumulative variance
p <- pct_variance(
  object = summobj$sobject,
  smethod = 'vst',
  cname = NULL,
  cname_ord = NULL,
  cutoff = 0.1,
  v = TRUE
)
p$plot$labels$subtitle <- gsub(".*mean.*: (.*%).*", "Percentage of variance explained: \\1", p$plot$labels$subtitle)
p$plot$labels$title <- "Highly Variable Genes"
p$plot$labels$caption <- NULL
p$plot$labels$x <- gsub(",.*", "", p$plot$labels$x)
p$plot$labels$y <- "Percentage"
# p$plot$layers <- p$plot$layers[-3]
fname <- paste0('hvg_cumulative_', summobj$name, '.pdf')
pdf(fname, width = 8, height = 6.5); set_big_text(p[[1]], sizey = 24); dev.off()

### Grid of scatter plots ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expdata <- expm1(mycells@assays$RNA@data) * 100; dtype <- ""
expdata <- saveymat; dtype <- "_saver"
expdata <- mycells@assays$RNA@data; dtype <- "_seurat"
gr <- "TH1"
genes <- c("IFNG", "XCL1", "XCL2", "KLRG1")
gr <- "TH17"
genes <- c("IL17F", "IL22", "CCL20", "CTSH")

plothem <- TRUE

thesecells <- getsubset(c("Cell_type", gr), mycells@meta.data, v = TRUE)
genesdf <- remove.factors(create_pairs(data.frame(genes, stringsAsFactors = FALSE))[, 1:2])
dname <- paste0("correlations_", gr, "", dtype, "/")
dir.create(dname)
source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
void <- apply(genesdf, 1, function(x){ # x = unname(unlist(genesdf[1, ]))
  gg <- rev(genes[genes %in% x])
  p <- get_densities(mat = expdata[, thesecells],
    log2t = FALSE, genes = gg, return_plot = TRUE, cuof = 1, v = T
  )$scatter
  # cat("Layers:", sapply(p$layers, function(x) class(x$geom)[1] ), "\n")
  tvar <- sapply(p$layers, function(x) !grepl("ensi|oint", class(x$geom)[1]) )
  p$layers <- p$layers[tvar]
  p <- p + stat_smooth(method = "lm", col = "red", se = TRUE)
  p <- p + geom_point(size = 3)
  if(plothem){
    pdf(paste0(dname, paste0(x, collapse = "_"), ".pdf"))#, 5, 5)
    print(p)
    dev.off()
    # pdf(paste0(dname, paste0(x, collapse = "_"), "blank.pdf"))#, 5, 5)
    # print(shut_it(p, "ext|ensity|mooth"))
    # dev.off()
  }
  p
})
names(void) <- do.call('paste', c(genesdf, sep = "_"))
voidv <- lapply(genes, function(x){ # x = unname(unlist(genesdf[1, ]))
  p <- vlnplot(cpmdata = expdata[, thesecells], gg = x, plotdots = FALSE, v = FALSE, log2t = TRUE,
    datatype = 'log2(CPM + 1)', return_plot = TRUE, noncero = TRUE) + ylim(c(0, 16))
  if(plothem){
    pdf(paste0(dname, x, "_violin.pdf"), width = 3.5, height = 7)
    print(p)
    dev.off()
    pdf(paste0(dname, x, "_violinblank.pdf"), width = 3.5, height = 7)
    print(shut_it(p))
    dev.off()
  }
  p
})
names(voidv) <- genes
pp <- list(); pmat <- mat_names(genes, genes)
layouty <- lapply(genes, function(x){
  lapply(genes, function(y){
    if(x == y){
      pmat[x, y] <<- "Violin"
      pp[[paste0(x, "_", y, "_Violin")]] <<- voidv[[x]] + labs(subtitle = NULL)
    }else if(which(genes == x) < which(genes == y)){
      pmat[y, x] <<- "Scatter"
      tvar <- strsplit(names(void), "_")
      tvar <- tvar[[which(sapply(tvar, function(gg) all(gg %in% c(x, y)) ))]]
      pp[[paste0(x, "_", y, "_Scatter")]] <<- void[[paste0(tvar, collapse = "_")]]
    }else{ pp[[paste0(x, "_", y, "_NULL")]] <<- "NULL" }
    return(NULL)
  })
})
pp <- lapply(pp, function(x) if(class(x)[1] != "character") x )
names(pp); pmat
# rel_wds <- ifelse(grepl("iolin", names(pp)), 0.5, 1)
pdf(paste0(dname, "_grid.pdf"), width = 25, height = 25)
plot_grid(plotlist = pp)#, rel_widths = rel_wds)
dev.off()
pps <- lapply(pp, function(x) if(!is.null(x)) shut_it(x) )
pdf(paste0(dname, "_gridblank.pdf"), width = 20, height = 20)
plot_grid(plotlist = pps)#, rel_widths = rel_wds)
dev.off()

### Table 2 Treg disease differences ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# now disease differences
fstat <- list(
  Asthma = data.frame(theObjectSavedIn(paste0(root, '/Bulk_RNASeq_Analysis/Asthma3_HDM_specific/Hg19/AS3E4_individual_tech_dups/Treg_asthma/ASvsNAS/ASvsNAS_DEGs_BetasMLE_SSet.RData'))),
  Allergy = data.frame(theObjectSavedIn(paste0(root, '/Bulk_RNASeq_Analysis/Asthma3_HDM_specific/Hg19/AS3E4_individual_tech_dups/Treg_allergy/ALvsNAL/ALvsNAL_DEGs_BetasMLE_SSet.RData')))
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
mystats <- lapply(c("Asthma", "Allergy", "Treg_dg"), function(x){
  get_stat_report(
    mat = expdata,
    groups = make_list(coldata_[getsubset(c("Cell_type", "Treg"), coldata_, v = TRUE), ], x, grouping = TRUE),
    moments = c("mn", "p"),
    expr_cutoff = 10,
    v = TRUE
  )
})
str(mystats)
mytab <- summ_tables(ltabs = c(fstat, mystats), v = TRUE)
mytab <- mytab[order(-mytab$Asthma_log2FoldChange), ]
head(mytab)
tail(mytab)
bordering(mytab[, 2:3])
sum(!is.na(mytab$Asthma_padj))
sum(!is.na(fstat[[1]]$Asthma_padj))

### Treg like figure 1 F and G ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thesecells <- getsubset(c("Cell_type", "Treg"), coldata_, v = TRUE)
thesecellst <- rownames(coldata_[order(coldata_[, "Disease"]), ])
# thesecellst <- rownames(coldata_[order(coldata_[, "Asthma"]), ])
# t(expdata["CXCR3", thesecells])
# tapply(unlist(expdata["CXCR3", thesecells]), INDEX = coldata_[thesecells, "Disease"], mean)
# tapply(unlist(expdata["CXCR3", thesecells]), INDEX = coldata_[thesecells, "Allergy"], mean)
# tapply(unlist(expdata["CXCR3", thesecells]), INDEX = coldata_[thesecells, "Asthma"], mean)
thesecells <- thesecellst[thesecellst %in% thesecells]
library(DESeq2) # may need R's DESeq2 and not R5's
fstat <- list(
  Asthma = data.frame(theObjectSavedIn(paste0(root, '/Bulk_RNASeq_Analysis/Asthma3_HDM_specific/Hg19/AS3E4_individual_tech_dups/Treg_asthma/ASvsNAS/ASvsNAS_DEGs_BetasMLE_SSet.RData'))),
  Allergy = data.frame(theObjectSavedIn(paste0(root, '/Bulk_RNASeq_Analysis/Asthma3_HDM_specific/Hg19/AS3E4_individual_tech_dups/Treg_allergy/ALvsNAL/ALvsNAL_DEGs_BetasMLE_SSet.RData')))
)
myfcs <- sapply(fstat, "[", 2)
mydatafc <- data.frame(myfcs, stringsAsFactors = F, row.names = rownames(fstat[[1]]))
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
ifis <- rownames(mydatafc)[grepl("^IFI", rownames(mydatafc))]
borgenes <- unique(c(ifis, bordering(mydatafc, 1:2)))
mydatafc_repel <- mydatafc[borgenes, ]
thres <- 1 # log2FoldChange threshold
# mydatafc[ifis, ]
# rowSums(mydatafc[ifis, 1:2] > thres)
# mydatafc[ifis, ][mydatafc[ifis, 3] < 0.5, ]
p <- ggplot(mydatafc, aes(x = Asthma, y = Allergy, color = Log2_Mean, size = significance)) +
  geom_point() + scale_color_gradientn(colours = c('#ffffff', '#670000')) + ##CB4154
  geom_hline(yintercept = c(-thres, thres), linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = c(-thres, thres), linetype = "dashed", alpha = 0.4) +
  labs(title = 'Allergy - Asthma differences',
    subtitle = paste('Threshold: FC >', thres), x = 'Asthma', y = 'Allergy',
    color = 'Expression', size = 'Significance') +
  drawsq + scale_size_area() + scale_size(range = c(0, 10), limits = c(0, 20))
p <- squareplot(p)
fname <- paste0('treg_fc_volcano_gene_FC_filtered50perCompGroup_', thres)
if(1){
  p <- p + guides(size = guide_legend(keywidth = 0.5, keyheight = 0.5, default.unit = "inch"))
  # fname <- paste0(fname, '_update')
}; set.seed(27)
library(ggrepel)
pdf(paste0(fname, '.pdf'), 10, 10)
print(p + geom_text_repel(data = mydatafc_repel, aes(label = gene_name), color = 'black'))
dev.off()
pdf(paste0(fname, 'blank.pdf'), 10, 10)
print(p + shut_up)
dev.off()

smytab <- mytab[rownames(mytab) %in% rownames(mydatafc), ]
smytab <- smytab[order(-smytab$Asthma_log2FoldChange), ]
dim(smytab)
head(smytab)
write.csv(smytab, file = "bulk_treg_disease_stats.csv")

### 7.X Treg t-SNE plot ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# for Treg, TNeg
selectss <-  c('Cell_type',c('Teff','Treg','TNeg')[c(2,3)])
pcs <- 3
pxs <- 10
coldata_$Treg_dg_tneg <- coldata_$Treg_dg
coldata_$Treg_dg_tneg[coldata_$Cell_type == 'TNeg'] <- 'TNeg'

samples <- getsubset(selectss, coldata_, v = T)
keep_genesl <- mostVars(counts_[, samples], 200, "AR")
keep_genes <- keep_genesl$genes
file.remove('Rplots.pdf')
cts <- log2(expdata[keep_genes, samples] + 1)
# summary(rowMeans(expdata[keep_genes, samples])) # min mean = 2.582, min pct = 53.70

pca.m <- prcomp(cts, center = TRUE)
pca.mat <- data.frame(pca.m$rotation, stringsAsFactors = FALSE, check.names = FALSE)

# dir.create('treg_tsnes')
for(pc in pcs){
  for(px in pxs){
    fname <- paste0('teff_tsnes_', paste0(pc, "PCs_", px), 'PX_', length(keep_genes), 'genes.pdf')
    cat(fname, '\n')
    set.seed(27)
    tsne_data <- Rtsne::Rtsne(dist(pca.mat[, 1:pc]), perplexity = px, pca = FALSE, theta = 0, is_distance = TRUE, max_iter = 1000)
    ddf <- data.frame(x.var = tsne_data$Y[, 1], y.var = tsne_data$Y[,2])
    rownames(ddf) <- rownames(pca.mat)
    ddf$col<-  factormix(coldata_[samples, 'Treg_dg_tneg'])

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
