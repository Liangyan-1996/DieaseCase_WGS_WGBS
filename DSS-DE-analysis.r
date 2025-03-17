if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("DSS")
library(DSS)

control = read.table("control1.cov.txt",sep="\t",header=F,col.names=c("chr","start","end","mePCT","meCounts","nonmeCounts"))
patient = read.table("patient1.cov.txt",sep="\t",header=F,col.names=c("chr","start","end","mePCT","meCounts","nonmeCounts"))
head(control)
head(patient)

# test on chr1
# control = control[control$chr=="chr1",]
# patient = patient[patient$chr=="chr1",]

# generate dataframe in the format of chr,pos,N,X
control["N"] = control["meCounts"] + control["nonmeCounts"]
control["X"] = control["meCounts"]
control["pos"] = control["start"]
control1 = control[,c("chr","pos","N","X")]
patient["N"] = patient["meCounts"] + patient["nonmeCounts"]
patient["X"] = patient["meCounts"]
patient["pos"] = patient["start"]
patient1 = patient[,c("chr","pos","N","X")]

BSobj = makeBSseqData(list(control1, patient1), c("control","patient"))
dmlTest = DMLtest(BSobj, group1 = c("patient"), group2= c("control"), smoothing = TRUE)
# dmls = callDML(dmlTest, p.threshold = 0.001)
dmrs = callDMR(dmlTest, p.threshold = 0.01)
showOneDMR(dmrs[1,], BSobj)
save.image(file="DSS-DE-analysis.RData")

# ------------------------------------------
load("DSS-DE-analysis.RData")
write.table(dmrs,file="DSS-DE-analysis.tab",sep="\t",row.names=F,quote=F,col.names=T)
dmr_ITGB2 = read.table("DSS-DE-analysis.ITGB2.tab",sep="\t",header=F,col.names=colnames(dmrs))

for (i in seq(dim(dmr_ITGB2)[1])){
  id = paste0("ITGB2_",dmr_ITGB2[i,1],":",dmr_ITGB2[i,2],"-",dmr_ITGB2[i,3])
  pdf(paste0(id,".pdf"))
  showOneDMR(dmr_ITGB2[i,], BSobj,ext=1000)
  dev.off()
}