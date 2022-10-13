library(MOVICS)
runDEA(dea.method = "edger",
       expr       = count, # raw count data
       moic.res   = cmoic.hnsc,
       prefix     = "TCGA-HNSC",# prefix of figure name
       res.path = "./DEA/") 
# run DEA with DESeq2
runDEA(dea.method = "deseq2",
       expr       = count,
       moic.res   = cmoic.hnsc,
       prefix     = "TCGA-HNSC",# prefix of figure name
       res.path = "./DEA/") 
runDEA(dea.method = "limma",
       expr       = fpkm, # normalized expression data
       moic.res   = cmoic.hnsc,
       prefix     = "TCGA-HNSC",# prefix of figure name
       res.path = "./DEA/") 
marker.up.edgeR <- runMarker(moic.res      = cmoic.hnsc,
                             dea.method    = "edger", # name of DEA method
                             prefix        = "TCGA-HNSC", # MUST be the same of argument in runDEA()
                             dat.path      = "./DEA/", # path of DEA files
                             res.path      = "./markers/", # path to save marker files
                             p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                             p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                             dirct         = "up", # direction of dysregulation in expression
                             n.marker      = 500, # number of biomarkers for each subtype
                             doplot        = TRUE, # generate diagonal heatmap
                             norm.expr     = count, # use normalized expression as heatmap input
                             annCol        = annCol, # sample annotation in heatmap
                             annColors     = annColors, # colors for sample annotation
                             show_rownames = FALSE, # show no rownames (biomarker name)
                             fig.name      = "UPREGULATED BIOMARKER HEATMAP - edgeR")
                             fig.name      = "UPREGULATED BIOMARKER HEATMAP - deseq2")
plot.marker=fread("./DEA/consensusMOIC_TCGA-HNSC_edger_test_result.CS2_vs_Others.txt")
plot.marker=as.data.frame(plot.marker)
rownames(plot.marker)=plot.marker$id
plot.marker1 <- add_regulate(plot.marker, log2FC_name = "log2fc",
                     fdr_name = "padj",log2FC = 1, fdr = 0.01)
colnames(plot.marker1)=c("row","fc","log2FoldChange","pvalue","padj","regulate")
ggvolcano(plot.marker1, x = "log2FoldChange", y = "padj",
          label = "row", label_number = 10, output = FALSE,colors = c("#2EC4B6", "gray","#E71D36"))

load("~/yau_41613.Rdata")
hnsc.yau=yau_41613
yau.ntp.pred <- runNTP(expr       = hnsc.yau$mRNA.expr,
                       templates  = marker.up.edgeR$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR GSE41613",
                       fig.path =   "./ntpheatmap") 
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = hnsc.yau$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR GSE41613 (edgeR)",
                     fig.path =   "./ntpheatmap") 
agree.yau <- compAgree(moic.res  = yau.ntp.pred,
                       subt2comp = hnsc.yau$clin.info[, "PAM50", drop = FALSE],
                       doPlot    = TRUE,
                       fig.name  = "YAU PREDICTEDMOIC WITH PAM50")
print(agree.yau)
yau.pam.pred <- runPAM(train.expr  = fpkm,
                       moic.res    = cmoic.hnsc,
                       test.expr   = hnsc.yau$mRNA.expr)
print(yau.pam.pred$IGP)

surv.yau <- compSurv(moic.res         = yau.pam.pred,
                     surv.info        = hnsc.yau$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR GSE41613") 
tcga.ntp.pred <- runNTP(expr      = fpkm,
                        templates = marker.up$templates,
                        doPlot    = FALSE) 
tcga.pam.pred <- runPAM(train.expr  = fpkm,
                        moic.res    = cmoic.hnsc,
                        test.expr   = fpkm)
runKappa(subt1     = cmoic.hnsc$clust.res$clust,
         subt2     = tcga.ntp.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "NTP",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and NTP")
runKappa(subt1     = cmoic.hnsc$clust.res$clust,
         subt2     = tcga.pam.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "PAM",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and PAM")
runKappa(subt1     = yau.ntp.pred$clust.res$clust,
         subt2     = yau.pam.pred$clust.res$clust,
         subt1.lab = "NTP",
         subt2.lab = "PAM",
         fig.name  = "CONSISTENCY HEATMAP FOR GSE41613")
