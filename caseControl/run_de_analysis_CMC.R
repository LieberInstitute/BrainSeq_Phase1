####
library(limma)
library(GenomicRanges)
library(sva)

## load data
load("/dcl01/lieber/ajaffe/PublicData/CMC/CMC_coverageMat_szControlEqtlDERs.rda")
load("/dcl01/lieber/ajaffe/PublicData/CMC/rpkmCounts_cmcDlpfc_szControl.rda")
load("/dcl01/lieber/ajaffe/PublicData/CMC/filteredTxData_CMC_szControl.rda")

## add ancestry
mds = read.table("/dcl01/lieber/ajaffe/PublicData/CMC/Genotypes/CMC_genotypes_imputed_common_cleaned.mds",
	header=TRUE, as.is=TRUE,row.names=1)
mds = mds[pd$Genotyping_Sample_ID,-(1:2)]
colnames(mds) = paste0("snpPC",1:10)
pd = cbind(pd,mds)

# fix age
pd$Age_of_Death[pd$Age_of_Death == "90+"] = 90
pd$Age_of_Death = as.numeric(pd$Age_of_Death)

# assignment rate
geneFn = paste0("/dcl01/lieber/ajaffe/PublicData/CMC/Counts/Gene/",
	pd$DLPFC_RNA_Sequencing_Sample_ID, "_Ensembl_v75_Genes.counts")
names(geneFn) = pd$DLPFC_RNA_Sequencing_Sample_ID
geneStatList = lapply(paste0(geneFn, ".summary"), 
	read.delim,row.names=1)
geneStats = do.call("cbind", geneStatList)
colnames(geneStats) = pd$DLPFC_RNA_Sequencing_Sample_ID
pd$totalAssignedGene = as.numeric(geneStats[1,] / colSums(geneStats))

##########################
## just take adults with some quality filters 
keepIndex = which(pd$totalAssignedGene > 0.3 & 
	pd$DLPFC_RNA_Sequencing_Percent_Aligned > 0.8 & 
	pd$DLPFC_RNA_isolation_RIN > 6 & 
	pd$Age < 80 & pd$snpPC3 > -0.01 & pd$snpPC5 < 0.01 & 
	pd$Ethnicity %in% c("African-American","Caucasian"))

## filter for people
pd = pd[keepIndex,]
jRpkm = jRpkm[,keepIndex]
exonRpkm = exonRpkm[,keepIndex]
geneRpkm = geneRpkm[,keepIndex]
tFpkmCmc = tFpkmCmc[,keepIndex]
regionMatCmc = regionMatCmc[,keepIndex]

### transform
yGene = as.matrix(log2(geneRpkm+1))
yExon = as.matrix(log2(exonRpkm+1))
yJxn = as.matrix(log2(jRpkm+1))
yTx = as.matrix(log2(tFpkmCmc+1))
yEr = as.matrix(log2(regionMatCmc+1))

###########################
## model based on EDA
mod = model.matrix(~ Dx + Age_of_Death + Institution + 
	DLPFC_RNA_isolation_RIN + DLPFC_RNA_Sequencing_Percent_Aligned + 
	Gender + Ethnicity + totalAssignedGene, data=pd)
rownames(mod) = pd$DLPFC_RNA_Sequencing_Sample_ID

## fits
fitGene = lmFit(yGene, mod)
ebGene = ebayes(fitGene)
fitExon = lmFit(yExon, mod)
ebExon = ebayes(fitExon)
fitJxn = lmFit(yJxn, mod)
ebJxn = ebayes(fitJxn)
fitTx = lmFit(yTx, mod)
ebTx = ebayes(fitTx)
fitEr = lmFit(yEr, mod)
ebEr = ebayes(fitEr)

############
# # adjust for degradation ###
degFiles = paste0("/users/ajaffe/Lieber/Projects/RNAseq/CMC/degrade/", 
	pd$DLPFC_RNA_Sequencing_Sample_ID, "_degradeStats.txt")
degCov = sapply(degFiles, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f10", x)), as.is=TRUE)$sum
})
colnames(degCov) = pd$DLPFC_RNA_Sequencing_Sample_ID
degCov  = degCov /100 # read length
bg = matrix(rep(pd$totalMapped/80e6), nc = nrow(pd), 
	nr = nrow(degCov), byrow=TRUE)
degCovAdj = degCov/bg
save(degCovAdj, file="rdas/degradation_mat_CMC_ribozero.rda")
load("rdas/degradation_mat_CMC_ribozero.rda")
degPca = prcomp(t(log2(degCovAdj+1)))

## how many PCs?
k = num.sv(log2(degCovAdj+1), mod)
qSVs = degPca$x[,1:k]

## fits
fitGeneQsva = lmFit(yGene, cbind(mod, qSVs))
ebGeneQsva = ebayes(fitGeneQsva)
fitExonQsva = lmFit(yExon, cbind(mod, qSVs))
ebExonQsva = ebayes(fitExonQsva)
fitJxnQsva = lmFit(yJxn, cbind(mod, qSVs))
ebJxnQsva = ebayes(fitJxnQsva)
fitTxQsva = lmFit(yTx, cbind(mod, qSVs))
ebTxQsva = ebayes(fitTxQsva)
fitErQsva = lmFit(yEr, cbind(mod, qSVs))
ebErQsva = ebayes(fitErQsva)

###########################
### pc adjustment #########
###########################

load("rdas/numPCs_caseControl_CMC.rda")

# ## gene
# geneIndex=order(rowSds(yGene),decreasing=TRUE)[1:20000]
# pcaGene = prcomp(t(yGene[geneIndex,]))
# kGene = num.sv(yGene[geneIndex,], mod)
# pcsGene = pcaGene$x[,1:kGene]
fitGenePca = lmFit(yGene, cbind(mod, pcsGene))
ebGenePca = ebayes(fitGenePca)

# ## exon
# exonIndex=order(rowSds(yExon),decreasing=TRUE)[1:50000]
# pcaExon = prcomp(t(yExon[exonIndex,]))
# kExon = num.sv(yExon[exonIndex,], mod)
# pcsExon = pcaExon$x[,1:kExon]
fitExonPca = lmFit(yExon, cbind(mod, pcsExon))
ebExonPca = ebayes(fitExonPca)

# ## transcript
# txIndex=order(rowSds(yTx),decreasing=TRUE)[1:50000]
# pcaTx = prcomp(t(yTx[txIndex,]))
# kTx = num.sv(yTx[txIndex,], mod)
# pcsTx = pcaTx$x[,1:kTx]
fitTxPca = lmFit(yTx, cbind(mod,pcsTx ))
ebTxPca = ebayes(fitTxPca)

# ## junction
# jxnIndex=order(rowSds(yJxn),decreasing=TRUE)[1:50000]
# pcaJxn = prcomp(t(yJxn[jxnIndex,]))
# kJxn = num.sv(yJxn[jxnIndex,], mod)
# pcsJxn = pcaJxn$x[,1:kJxn]
fitJxnPca = lmFit(yJxn, cbind(mod, pcsJxn ))
ebJxnPca = ebayes(fitJxnPca)

# ## er
# erIndex=order(rowSds(yEr),decreasing=TRUE)[1:50000]
# pcaEr = prcomp(t(yEr[erIndex,]))
# kEr = num.sv(yEr[erIndex,], mod)
# pcsEr = pcaEr$x[,1:kEr]
fitErPca = lmFit(yEr, cbind(mod, pcsEr))
ebErPca = ebayes(fitErPca)

###################
## write out

## gene
outGene = geneMap
outGene$meanExprsAdult = rowMeans(geneRpkm)
outGene$log2FC_adj = fitGene$coef[,2]
outGene$tstat_adj = ebGene$t[,2]
outGene$pval_adj = ebGene$p[,2]
outGene$log2FC_qsva = fitGeneQsva$coef[,2]
outGene$tstat_qsva = ebGeneQsva$t[,2]
outGene$pval_qsva = ebGeneQsva$p[,2]
outGene$log2FC_pca = fitGenePca$coef[,2]
outGene$tstat_pca = ebGenePca$t[,2]
outGene$pval_pca = ebGenePca$p[,2]
outGene = makeGRangesFromDataFrame(outGene, keep=TRUE)

### exon
outExon = exonMap
outExon$meanExprsAdult = rowMeans(exonRpkm)
outExon$log2FC_adj = fitExon$coef[,2]
outExon$tstat_adj = ebExon$t[,2]
outExon$pval_adj = ebExon$p[,2]
outExon$log2FC_qsva = fitExonQsva$coef[,2]
outExon$tstat_qsva = ebExonQsva$t[,2]
outExon$pval_qsva = ebExonQsva$p[,2]
outExon$log2FC_pca = fitExonPca$coef[,2]
outExon$tstat_pca = ebExonPca$t[,2]
outExon$pval_pca = ebExonPca$p[,2]
outExon = makeGRangesFromDataFrame(outExon, keep=TRUE)

### junction
outJxn = jMap
outJxn$meanExprsAdult = rowMeans(jRpkm)
outJxn$log2FC_adj = fitJxn$coef[,2]
outJxn$tstat_adj = ebJxn$t[,2]
outJxn$pval_adj = ebJxn$p[,2]
outJxn$log2FC_qsva = fitJxnQsva$coef[,2]
outJxn$tstat_qsva = ebJxnQsva$t[,2]
outJxn$pval_qsva = ebJxnQsva$p[,2]
outJxn$log2FC_pca = fitJxnPca$coef[,2]
outJxn$tstat_pca = ebJxnPca$t[,2]
outJxn$pval_pca = ebJxnPca$p[,2]

### transcript
outTx = tMap
outTx$meanExprsAdult = rowMeans(tFpkmCmc)
outTx$log2FC_adj = fitTx$coef[,2]
outTx$tstat_adj = ebTx$t[,2]
outTx$pval_adj = ebTx$p[,2]
outTx$log2FC_qsva = fitTxQsva$coef[,2]
outTx$tstat_qsva = ebTxQsva$t[,2]
outTx$pval_qsva = ebTxQsva$p[,2]
outTx$log2FC_pca = fitTxPca$coef[,2]
outTx$tstat_pca = ebTxPca$t[,2]
outTx$pval_pca = ebTxPca$p[,2]

## ERs
outEr = regionsCmc
outEr$meanExprsAdult = rowMeans(regionMatCmc)
outEr$log2FC_adj = fitEr$coef[,2]
outEr$tstat_adj = ebEr$t[,2]
outEr$pval_adj = ebEr$p[,2]
outEr$log2FC_qsva = fitErQsva$coef[,2]
outEr$tstat_qsva = ebErQsva$t[,2]
outEr$pval_qsva = ebErQsva$p[,2]
outEr$log2FC_pca = fitErPca$coef[,2]
outEr$tstat_pca = ebErPca$t[,2]
outEr$pval_pca = ebErPca$p[,2]

save(outGene, outExon, outJxn, outTx, outEr,
	file="rdas/DE_statistics_adjAndQsva_CMC.rda",
	compress=TRUE)
	
save(mod, pcsGene, pcsExon, pcsJxn, pcsTx, pcsEr,
	file="rdas/numPCs_caseControl_CMC.rda")


###############################
##### check degradation #####

#### read in degradation data
load("/users/ajaffe/Lieber/Projects/RNAseq/Degradation/rdas/degradeStats_byLibrary.rda")
degradeStatsGene = degradeStats[rownames(geneRpkm),]
load("/users/ajaffe/Lieber/Projects/RNAseq/Degradation/rdas/degradeStats_junctions_byLibrary.rda")
mmJxn = match(rownames(jRpkm), rownames(degradeStats)) # match
degradeStatsJxn = degradeStats[mmJxn[!is.na(mmJxn)],]

## gene
plot(ebGene$t[,2], degradeStatsGene$ribo_Tstat,pch=21,bg="grey")
plot(fitGene$coef[,2], degradeStatsGene$ribo_timeFC,pch=21,bg="grey")
table(sign(degradeStatsGene$ribo_Tstat)[ebGene$p[,2] < 1e-5], 
	sign(fitGene$coef[,2])[ebGene$p[,2] < 1e-5])
cor.test(ebGene$t[,2], degradeStatsGene$ribo_Tstat)
cor.test(fitGene$coef[,2], degradeStatsGene$ribo_timeFC)
