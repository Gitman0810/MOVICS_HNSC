maf <- read.maf(maf = "/home/data/refdir/database/tcga_data/TCGA_SNV/mutec2/TCGA.HNSC.mutect.1aa33f25-3893-4f37-a6a4-361c9785d07e.DR-10.0.somatic.maf.gz", isTCGA = T)
maf_clin=cmoic.hnsc$clust.res
maf_clin$samID=str_sub(cmoic.hnsc$clust.res$samID,1,12)
maf_clin$clust=ifelse(maf_clin$clust==1,"CS1","CS2")
maf@clinical.data$Subtype=maf_clin$clust[match(maf@clinical.data$Tumor_Sample_Barcode,maf_clin$samID)]
# 从临床数据中提取性别对应的"Tumor_Sample_Barcode"
cs1_maf <- subset(maf_clin, clust=="CS1")$samID
cs2_maf <- subset(maf_clin, clust=="CS2")$samID
maf_all=c(cs1_maf,cs2_maf)
# 使用subsetMaf构建男性和女性的MAF对象
cs1_maf <- subsetMaf(maf=maf, tsb=cs1_maf, isTCGA=TRUE)
cs2_maf <- subsetMaf(maf=maf, tsb=cs2_maf, isTCGA=TRUE)
maf_all=subsetMaf(maf=maf, tsb=maf_all, isTCGA=TRUE)
# 使用mafCompare比较差异突变基因
fvsm <- mafCompare(m1=cs1_maf, m2=cs2_maf, m1Name="CS1", m2Name="CS2", minMut=5)
# 结果保存到文件"female_vs_male.tsv"
write.table(fvsm$results, file="female_vs_male.tsv", quote=FALSE, row.names=FALSE, sep="\t")

#------------------------------------------------------------------------------------------------------------
library(maftools)
geneCloud(input=maf_all, top=10)


vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

#1. 绘制森林图
forestPlot(mafCompareRes=fvsm, pVal=0.01, color=c("#E71D36","#2EC4B6"), geneFontSize=0.8)
#2. 比较突变率
coBarplot(m1 = cs1_maf, m2 = cs2_maf, m1Name = "CS1", m2Name = "CS2",colors = vc_cols)
#3. 瀑布图
subtypecolors = c("#2EC4B6","#E71D36")
names(subtypecolors) = c("CS1","CS2")
ancolors = list(Subtype = subtypecolors)
oncoplot(maf = maf_all, colors = vc_cols, top = 20,clinicalFeatures = c("Subtype"),sortByAnnotation = TRUE,
         annotationColor = ancolors
         )


output_cs1 <- somaticInteractions(maf=cs1_maf, top=50, pvalue=c(0.05, 0.01))
output_cs2 <- somaticInteractions(maf=cs2_maf, top=50, pvalue=c(0.05, 0.01))

#临床富集分析
clin_enrich <- clinicalEnrichment(maf=maf_all, clinicalFeature="Subtype")
#save(clin_enrich ,file = "./mutations/clin_enrich.Rdata")
load("./mutations/clin_enrich.Rdata")
plotEnrichmentResults(enrich_res=clin_enrich, pVal=0.01)


#药物基因互作
druggability_cs1 <- drugInteractions(maf=cs1_maf)
druggability_cs2 <- drugInteractions(maf=cs2_maf)
#save(druggability_cs1,druggability_cs2,file="./mutations/druggability_cs1&2.Rdata")
load("./mutations/druggability_cs1&2.Rdata")

#2. 药物和基因间互作
kras <- drugInteractions(genes="TP53", drugs=TRUE)
kras[,.(Gene, interaction_types, drug_name, drug_claim_name)]

#致癌信号通路
#OncogenicPathways函数可以分析TCGA中已知的10个致癌信号通路中基因突变的数量和比例。包括cell cycle, Hippo, Myc, Notch, Nrf2, PI-3-Kinase/Akt, RTK-RAS, TGFβ signaling, p53 and β-catenin/Wnt：
OncogenicPathways1(maf=cs1_maf)
OncogenicPathways(maf=cs2_maf)


PlotOncogenicPathways(maf = cs1_maf, pathways = "TP53")
PlotOncogenicPathways(maf = cs2_maf, pathways = "TP53")


#生存分析
maf.hnsc <- read.maf(maf = "/home/data/refdir/database/tcga_data/TCGA_SNV/mutec2/TCGA.HNSC.mutect.1aa33f25-3893-4f37-a6a4-361c9785d07e.DR-10.0.somatic.maf.gz", isTCGA = T)
load("../../m5c/TCGA-HNSC_sur_model.Rdata")
meta$samID=str_sub(meta$ID,1,12)
maf.hnsc@clinical.data$time=meta$time[match(maf.hnsc@clinical.data$Tumor_Sample_Barcode,meta$samID)]*30
maf.hnsc@clinical.data$event=meta$event[match(maf.hnsc@clinical.data$Tumor_Sample_Barcode,meta$samID)]
mafSurvival(maf = maf.hnsc, genes = c('TP53','TTN','FAT1','CDKN2A','CSMD3'), time = 'time', Status = 'event', isTCGA = TRUE)

