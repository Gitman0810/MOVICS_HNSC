# https://xlucpu.github.io/MOVICS/MOVICS-VIGNETTE.html
setwd("~/bulk_work/HNSC_MOVICS")
library(MOVICS)
#load(system.file("extdata", "hnsc.tcga.RData", package = "MOVICS", mustWork = TRUE))
load("hnsc.tcga.8.19.2.Rdata")  
#---------------
#------------------------------4.2.1 GET Module------------------------------
# extract raw count data for downstream analyses
count     <- hnsc.tcga$count

# extract fpkm data for downstream analyses
fpkm      <- hnsc.tcga$fpkm

# extract maf for downstream analysis
maf       <- hnsc.tcga$maf

# extract segmented copy number for downstream analyses
segment   <- hnsc.tcga$segment

# extract survival information
surv.info <- hnsc.tcga$clin.info

mRNA.expr=hnsc.tcga$mRNA.expr
lncRNA.expr=hnsc.tcga$lncRNA.expr
miRNA.expr=hnsc.tcga$miRNA.expr
meth.beta=hnsc.tcga$meth.beta

load("./mo.data.8.19.Rdata")
optk.hnsc <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,F,T), # note: 第四个数据是体细胞突变，是一个二元矩阵
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.name    = "CLUSTER NUMBER OF TCGA-HNSC")
save(optk.hnsc,file="./optk.hnsc.8.27.Rdata")
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 2,
                         type        = c("gaussian", "gaussian", "gaussian","gaussian", "binomial"))
moic.res.list <- append(moic.res.list, 
                        list("iClusterBayes" = iClusterBayes.res))
save(mo.data,moic.res.list, file = "./moic.res.list.rda")
load("./moic.res.list.rda")
cmoic.hnsc <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.path = getwd(),
                               fig.name      = "CONSENSUS HEATMAP", 
                               distance      = "euclidean", #euclidean
                               linkage       = "average",showID = F)
getSilhouette(sil      = cmoic.hnsc$sil,
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)
indata <- mo.data
indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))

# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,T,F)) # no scale for mutation

feat   <- iClusterBayes.res$feat.res
feat1  <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA.expr"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "miRNA.expr"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
feat5  <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3,
               feat4, feat5)
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#fdbb2d")
miRNA.col=c("#40e0d0","#ff8c00","#ff0080")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col,miRNA.col, 
                   meth.col, mut.col)

getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","miRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","miRNA.FPKM","M value","Mutated"),
             clust.res     = iClusterBayes.res$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")

getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = moic.res.list$COCA$clust.res, # cluster results
             clust.dend    = moic.res.list$COCA$clust.dend, # show dendrogram for samples
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF COCA")

annCol    <- surv.info[,c("grade", "stage", "age"), drop = FALSE]

annColors <- list(age    = circlize::colorRamp2(breaks = c(min(annCol$age),
                                                           median(annCol$age),
                                                           max(annCol$age)), 
                                                colors = c("#AAAA00", "#555555", "#0000AA")),
                  grade  = c("G1" = "#b2d8d8",
                             "G2"   = "#66b2b2",
                             "G3"   = "#008080",
                             "G4"   = "#004c4c",
                             "GX" = "gray"),
                  stage = c("I"    = "#ffdbac",
                             "II"    = "#e0ac69",
                             "III"    = "#c68642",
                             "IV"    = "#8d5524", 
                             "TX"    = "black"))

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","miRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","miRNA.FPKM","M value","Mutated"),
             clust.res     = cmoic.hnsc$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F,F), # show no dendrogram for features
             annRow        = annRow, # no selected features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF CONSENSUSMOIC (annRow)",
             fig.path      = "./heatmap")
#12*15

#pca/tsne
tinyarray::draw_pca(mo.data$mRNA.expr,group,color = c("#2EC4B6", "#E71D36"),addEllipses = F)
tinyarray::draw_pca(mo.data$lncRNA.expr,group,color = c("#2EC4B6", "#E71D36"),addEllipses = F)
tinyarray::draw_pca(mo.data$miRNA.expr,group,color = c("#2EC4B6", "#E71D36"),addEllipses = F)
tinyarray::draw_pca(mo.data$meth.beta,group,color = c("#2EC4B6", "#E71D36"),addEllipses = F)

tinyarray::draw_tsne(mo.data$mRNA.expr,group,color = c("#2EC4B6", "#E71D36"))
tinyarray::draw_tsne(mo.data$lncRNA.expr,group,color = c("#2EC4B6", "#E71D36"))
tinyarray::draw_tsne(mo.data$miRNA.expr,group,color = c("#2EC4B6", "#E71D36"))
tinyarray::draw_tsne(mo.data$meth.beta,group,color = c("#2EC4B6", "#E71D36"))

tinyarray::draw_pca(fpkm,group,color = c("#2EC4B6", "#E71D36"),addEllipses = F)
