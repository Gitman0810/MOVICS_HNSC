library(MOVICS)
load("./Rdata/moic.res.list.rda")
cmoic.hnsc <- getConsensusMOIC(moic.res.list = moic.res.list, fig.name      = "CONSENSUS HEATMAP", distance      = "euclidean", linkage       = "average")
surv.hnsc <- compSurv(moic.res         = cmoic.hnsc,
                      surv.info        = surv.info,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")
print(surv.hnsc)
clin.hnsc <- compClinvar(moic.res      = cmoic.hnsc,
                         var2comp      = surv.info, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("gender","grade","stage","race","fustat"), #被认为是分类变量的features
                         nonnormalVars = "futime", # 被认为使用非参数检验的 feature(s)
                         exactVars     = c("stage"), # 考虑使用精确测试的feature(s)
                         doWord        = TRUE, # generate .docx file in local path
                         includeNA=F,
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")

print(clin.hnsc$compTab)
mut.hnsc <- compMut(moic.res     = cmoic.hnsc,
                    mut.matrix   = hnsc.tcga$mut.status, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
                    p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                    innerclust   = TRUE, # perform clustering within each subtype
                    annCol       = annCol, # same annotation for heatmap
                    annColors    = annColors, # same annotation color for heatmap
                    width        = 6, 
                    height       = 2,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")

print(mut.hnsc)
head(maf)
tmb.hnsc <- compTMB(moic.res     = cmoic.hnsc,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")

head(tmb.hnsc$TMB.dat)

colnames(segment) <- c("sample","chrom","start","end","value")

alldrugs=c("Erlotinib", "Rapamycin", "Sunitinib", "PHA-665752", "MG-132", "Paclitaxel", "Cyclopamine", "AZ628", "Sorafenib", "VX-680", "Imatinib", "TAE684", "Crizotinib", "Saracatinib", "S-Trityl-L-cysteine", "Z-LLNle-CHO", "Dasatinib", "GNF-2", "CGP-60474", "CGP-082996", "A-770041", "WH-4-023", "WZ-1-84", "BI-2536", "BMS-536924", "BMS-509744", "CMK", "Pyrimethamine", "JW-7-52-1", "A-443654", "GW843682X", "MS-275", "Parthenolide", "KIN001-135", "TGX221", "Bortezomib", "XMD8-85", "Roscovitine", "Salubrinal", "Lapatinib", "GSK269962A", "Doxorubicin", "Etoposide", "Gemcitabine", "Mitomycin C", "Vinorelbine", "NSC-87877", "Bicalutamide", "QS11", "CP466722", "Midostaurin", "CHIR-99021", "AP-24534", "AZD6482", "JNK-9L", "PF-562271", "HG-6-64-1", "JQ1", "JQ12", "DMOG", "FTI-277", "OSU-03012", "Shikonin", "AKT inhibitor VIII", "Embelin", "FH535", "PAC-1", "IPA-3", "GSK-650394", "BAY 61-3606", "5-Fluorouracil", "Thapsigargin", "Obatoclax Mesylate", "BMS-754807", "Lisitinib", "Bexarotene", "Bleomycin", "LFM-A13", "GW-2580", "AUY922", "Phenformin", "Bryostatin 1", "Pazopanib", "LAQ824", "Epothilone B", "GSK1904529A", "BMS345541", "Tipifarnib", "BMS-708163", "Ruxolitinib", "AS601245", "Ispinesib Mesylate", "TL-2-105", "AT-7519", "TAK-715", "BX-912", "ZSTK474", "AS605240", "Genentech Cpd 10", "GSK1070916", "KIN001-102", "LY317615", "GSK429286A", "FMK", "QL-XII-47", "CAL-101", "UNC0638", "XL-184", "WZ3105", "XMD14-99", "AC220", "CP724714", "JW-7-24-1", "NPK76-II-72-1", "STF-62247", "NG-25", "TL-1-85", "VX-11e", "FR-180204", "Tubastatin A", "Zibotentan", "YM155", "NSC-207895", "VNLG/124", "AR-42", "CUDC-101", "Belinostat", "I-BET-762", "CAY10603", "Linifanib ", "BIX02189", "CH5424802", "EKB-569", "GSK2126458", "KIN001-236", "KIN001-244", "KIN001-055", "KIN001-260", "KIN001-266", "Masitinib", "MP470", "MPS-1-IN-1", "BHG712", "OSI-930", "OSI-027", "CX-5461", "PHA-793887", "PI-103", "PIK-93", "SB52334", "TPCA-1", "TG101348", "Foretinib", "Y-39983", "YM201636", "Tivozanib", "GSK690693", "SNX-2112", "QL-XI-92", "XMD13-2", "QL-X-138", "XMD15-27")
drug.hnsc <- compDrugsen(moic.res    = cmoic.hnsc,
                         norm.expr   = fpkm[,cmoic.hnsc$clust.res$samID], # double guarantee sample order
                         drugs       = alldrugs, # a vector of names of drug in GDSC
                         tissueType  = "all", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED IC50",
                         fig.path    ="./drugs") 

head(drug.hnsc$Cisplatin)

surv.info$stage <- factor(surv.info$stage, levels = c("I","II","III","IV"))
surv.info$grade <- factor(surv.info$grade, levels = c("GX","G1","G2","G3","G4"))  #原文

agree.hnsc <- compAgree(moic.res  = cmoic.hnsc,
                        subt2comp = surv.info[,c("stage","grade")],
                        doPlot    = TRUE,
                        box.width = 0.2,
                        fig.name  = "AGREEMENT OF CONSENSUSMOIC WITH GRADE AND STAGE")

print(agree.hnsc)
