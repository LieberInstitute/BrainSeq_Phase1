####
library(limma)
library(GenomicRanges)
library(sva)
library(GOstats)

## load data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/transcript/transcript_data_filtered_n495.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda")

##########################
## just take adults with some quality filters 
table(pd$Dx[pd$Age > 17])

## filter
keepIndex = which(pd$Age > 17 & 
	pd$totalAssignedGene > 0.5 & 
	pd$mappingRate > 0.7 & pd$RIN > 6 & 
	pd$Age < 80 & pd$snpPC2 < 0.03)

## filter for people
pd = pd[keepIndex,]
jRpkm = jRpkm[,keepIndex]
exonRpkm = exonRpkm[,keepIndex]
geneRpkm = geneRpkm[,keepIndex]
tFpkm = tFpkm[,keepIndex]
regionMat = regionMat[,keepIndex]

table(pd$Dx[pd$Age > 17])

### transform
yGene = as.matrix(log2(geneRpkm+1))
yExon = as.matrix(log2(exonRpkm+1))
yJxn = as.matrix(log2(jRpkm+1))
yTx = as.matrix(log2(tFpkm+1))
yEr = as.matrix(log2(regionMat+1))

####################
### univariate #####
####################
modUniv = model.matrix(~pd$Dx)
fitGeneUniv = lmFit(yGene, modUniv)
ebGeneUniv = ebayes(fitGeneUniv)
sum(p.adjust(ebGeneUniv$p[,2], "fdr") < 0.05)
mean(p.adjust(ebGeneUniv$p[,2], "fdr") < 0.05)

###########################
## model based on EDA
mod = model.matrix(~ Dx + Age + Sex + 
	snpPC1 + snpPC5 + snpPC6 + snpPC9 + snpPC10 +
	mitoRate + RIN + totalAssignedGene, data=pd)
rownames(mod) = pd$RNum

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
# adjust ###

degFiles = paste0("degrade/DLPFC_PolyA_", 
	pd$RNum, "_degradeStats.txt")
degCov = sapply(degFiles, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f10", x)), 
		as.is=TRUE)$sum
})
colnames(degCov) = pd$RNum
degCov = degCov/100 # read length
bg = matrix(rep(pd$totalMapped/80e6), nc = nrow(pd), 
	nr = nrow(degCov),	byrow=TRUE)
degCovAdj = degCov/bg
save(degCovAdj, file="rdas/degradation_mat_LIBD_polyA.rda")
load("rdas/degradation_mat_LIBD_polyA.rda")
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
load("rdas/numPCs_caseControl_LIBD.rda")

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
###########################

## gene
outGene = geneMap
outGene$meanExprsAdult = rowMeans(geneRpkm)
outGene$isExp = outGene$meanExprsAdult > 0.1
outGene$code = "InEns"
outGene$log2FC_adj = fitGene$coef[,2]
outGene$tstat_adj = ebGene$t[,2]
outGene$pval_adj = ebGene$p[,2]
outGene$fdr_adj = NA
outGene$fdr_adj[outGene$isExp] = p.adjust(outGene$pval_adj[outGene$isExp],"fdr")
outGene$log2FC_qsva = fitGeneQsva$coef[,2]
outGene$tstat_qsva = ebGeneQsva$t[,2]
outGene$pval_qsva = ebGeneQsva$p[,2]
outGene$fdr_qsva = NA
outGene$fdr_qsva[outGene$isExp] = p.adjust(outGene$pval_qsva[outGene$isExp],"fdr")
outGene$log2FC_pca = fitGenePca$coef[,2]
outGene$tstat_pca = ebGenePca$t[,2]
outGene$pval_pca = ebGenePca$p[,2]
outGene$fdr_pca = NA
outGene$fdr_pca[outGene$isExp] = p.adjust(outGene$pval_pca[outGene$isExp],"fdr")
outGene = makeGRangesFromDataFrame(outGene, keep=TRUE)

### exon
outExon = exonMap
outExon$meanExprsAdult = rowMeans(exonRpkm)
outExon$isExp = outExon$meanExprsAdult > 0.2
outExon$code = "InEns"
outExon$log2FC_adj = fitExon$coef[,2]
outExon$tstat_adj = ebExon$t[,2]
outExon$pval_adj = ebExon$p[,2]
outExon$fdr_adj = NA
outExon$fdr_adj[outExon$isExp] = p.adjust(outExon$pval_adj[outExon$isExp],"fdr")
outExon$log2FC_qsva = fitExonQsva$coef[,2]
outExon$tstat_qsva = ebExonQsva$t[,2]
outExon$pval_qsva = ebExonQsva$p[,2]
outExon$fdr_qsva = NA
outExon$fdr_qsva[outExon$isExp] = p.adjust(outExon$pval_qsva[outExon$isExp],"fdr")
outExon$log2FC_pca = fitExonPca$coef[,2]
outExon$tstat_pca = ebExonPca$t[,2]
outExon$pval_pca = ebExonPca$p[,2]
outExon$fdr_pca = NA
outExon$fdr_pca[outExon$isExp] = p.adjust(outExon$pval_pca[outExon$isExp],"fdr")
outExon = makeGRangesFromDataFrame(outExon, keep=TRUE)

### junction
outJxn = jMap
outJxn$meanExprsAdult = rowMeans(jRpkm)
outJxn$isExp = outJxn$meanExprsAdult > 1 & outJxn$code != "Novel"
outJxn$log2FC_adj = fitJxn$coef[,2]
outJxn$tstat_adj = ebJxn$t[,2]
outJxn$pval_adj = ebJxn$p[,2]
outJxn$fdr_adj = NA
outJxn$fdr_adj[outJxn$isExp] = p.adjust(outJxn$pval_adj[outJxn$isExp],"fdr")
outJxn$log2FC_qsva = fitJxnQsva$coef[,2]
outJxn$tstat_qsva = ebJxnQsva$t[,2]
outJxn$pval_qsva = ebJxnQsva$p[,2]
outJxn$fdr_qsva = NA
outJxn$fdr_qsva[outJxn$isExp] = p.adjust(outJxn$pval_qsva[outJxn$isExp],"fdr")
outJxn$log2FC_pca = fitJxnPca$coef[,2]
outJxn$tstat_pca = ebJxnPca$t[,2]
outJxn$pval_pca = ebJxnPca$p[,2]
outJxn$fdr_pca = NA
outJxn$fdr_pca[outJxn$isExp] = p.adjust(outJxn$pval_pca[outJxn$isExp],"fdr")

### transcript
outTx = tMap
outTx$meanExprsAdult = rowMeans(tFpkm)
outTx$isExp = outTx$meanExprsAdult > 0.2

code = tMap$class_code
outTx$code = "Novel"
outTx$code[code == "="] = "InEns"
outTx$code[code == "j"] = "ExonSkip"
outTx$code[code %in% c("c", "e", "o", "p", "x") ] = "AltStartEnd"

outTx$log2FC_adj = fitTx$coef[,2]
outTx$tstat_adj = ebTx$t[,2]
outTx$pval_adj = ebTx$p[,2]
outTx$fdr_adj = NA
outTx$fdr_adj[outTx$isExp] = p.adjust(outTx$pval_adj[outTx$isExp],"fdr")
outTx$log2FC_qsva = fitTxQsva$coef[,2]
outTx$tstat_qsva = ebTxQsva$t[,2]
outTx$pval_qsva = ebTxQsva$p[,2]
outTx$fdr_qsva = NA
outTx$fdr_qsva[outTx$isExp] = p.adjust(outTx$pval_qsva[outTx$isExp],"fdr")
outTx$log2FC_pca = fitTxPca$coef[,2]
outTx$tstat_pca = ebTxPca$t[,2]
outTx$pval_pca = ebTxPca$p[,2]
outTx$fdr_pca = NA
outTx$fdr_pca[outTx$isExp] = p.adjust(outTx$pval_pca[outTx$isExp],"fdr")

## ERs
outEr = regions
outEr$meanExprsAdult = rowMeans(regionMat)
outEr$isExp = TRUE
code2 = regions$annoClass
outEr$code = "Novel"
outEr$code[code2 == "strictExonic"] = "InEns"
outEr$code[code2 %in% c("exonIntron", "extendUTR")] = "AltStartEnd"

outEr$log2FC_adj = fitEr$coef[,2]
outEr$tstat_adj = ebEr$t[,2]
outEr$pval_adj = ebEr$p[,2]
outEr$fdr_adj = p.adjust(outEr$pval_adj,"fdr")
outEr$log2FC_qsva = fitErQsva$coef[,2]
outEr$tstat_qsva = ebErQsva$t[,2]
outEr$pval_qsva = ebErQsva$p[,2]
outEr$fdr_qsva= p.adjust(outEr$pval_qsva,"fdr")
outEr$log2FC_pca = fitErPca$coef[,2]
outEr$tstat_pca = ebErPca$t[,2]
outEr$pval_pca = ebErPca$p[,2]
outEr$fdr_pca= p.adjust(outEr$pval_pca,"fdr")

save(outGene, outExon, outJxn, outTx,outEr, 
	file="rdas/DE_statistics_adjAndQsva.rda",
	compress=TRUE)
# save(mod, pcsGene, pcsExon, pcsJxn, pcsTx, pcsEr,
	# file="rdas/numPCs_caseControl_LIBD.rda")
	
######################

summary(lm(abs(outJxn$tstat_qsva) ~ outJxn$code))

###############################
##### check degradation #####

#### read in degradation data
load("/users/ajaffe/Lieber/Projects/RNAseq/Degradation/rdas/degradeStats_byLibrary.rda")
degradeStatsGene = degradeStats[names(outGene)[outGene$isExp],]
uniStats = data.frame(log2FC = fitGeneUniv$coef[,2],
	tstat = ebGeneUniv$t[,2], 
	pval = ebGeneUniv$t[,2])[names(outGene),]

## gene, tstatistics
pdf("plots/degradation_plots_geneLevel_tstats.pdf")
par(mar=c(5,6,3,2), cex.lab=2,cex.axis=2,cex.main=1.6)
plot(uniStats$tstat, degradeStatsGene$polyA_Tstat,
	xlab="SZ vs Control (Univariate)", 
	ylab= "Degradation Stat",
	main = "T-statistics", pch=21,bg="grey")
legend("topleft", paste0("r=",signif(cor(
	uniStats$tstat, degradeStatsGene$polyA_Tstat),3)),
	cex=1.4)
abline(h=0, v=c(-6,6))
table(sign(uniStats$tstat[abs(uniStats$tstat) > 6]),
	sign(degradeStatsGene$polyA_Tstat[abs(uniStats$tstat) > 6]),
	dnn = c("SZ","Degrade"))

plot(outGene$tstat_adj, degradeStatsGene$polyA_Tstat,
	xlab="SZ vs Control (Adjusted)", 
	ylab= "Degradation Stat",
	main = "T-statistics", pch=21,bg="grey")
abline(h=0, v=c(-4,4))
legend("topleft", paste0("r=",signif(cor(
	outGene$tstat_adj, degradeStatsGene$polyA_Tstat),3)),
	cex=1.4)
table(sign(outGene$tstat_adj[abs(outGene$tstat_adj) > 4]),
	sign(degradeStatsGene$polyA_Tstat[abs(outGene$tstat_adj) > 4]),
	dnn = c("SZ","Degrade"))

plot(outGene$tstat_qsva, degradeStatsGene$polyA_Tstat,
	xlab="SZ vs Control (Adj+Qual)", 
	ylab= "Degradation Stat",
	main = "T-statistics", pch=21,bg="grey")
abline(h=0, v=c(-4,4))
legend("topleft", paste0("r=",signif(cor(
	outGene$tstat_qsva, degradeStatsGene$polyA_Tstat),3)),
	cex=1.4)
dev.off()

##############33
## fold change
pdf("plots/degradation_plots_geneLevel_foldChange.pdf")
par(mar=c(5,6,3,2), cex.lab=2,cex.axis=2,cex.main=1.6)
plot(uniStats$log2FC, degradeStatsGene$polyA_timeFC,
	xlab="SZ vs Control (Univariate)", 
	ylab= "Degradation Stat",
	main = "Fold Changes", pch=21,bg="grey")
legend("topleft", paste0("r=",signif(cor(
	uniStats$log2FC, degradeStatsGene$polyA_timeFC),3)),
	cex=1.4)
abline(h=0, v=c(-0.25,0.25))

plot(outGene$log2FC_adj, degradeStatsGene$polyA_timeFC,
	xlab="SZ vs Control (Adjusted)", 
	ylab= "Degradation Stat",
	main = "Fold Changes", pch=21,bg="grey")
abline(h=0, v=c(-0.25,0.25))
legend("topleft", paste0("r=",signif(cor(
	outGene$log2FC_adj, degradeStatsGene$polyA_timeFC),3)),
	cex=1.4)

plot(outGene$log2FC_qsva, degradeStatsGene$polyA_timeFC,
	xlab="SZ vs Control (Adj+Qual)", 
	ylab= "Degradation Stat",
	main = "Fold Changes", pch=21,bg="grey")
abline(h=0, v=c(-0.25,0.25))
legend("topleft", paste0("r=",signif(cor(
	outGene$log2FC_qsva, degradeStatsGene$polyA_timeFC),3)),
	cex=1.4)
dev.off()


## "classic" model
pheno = read.csv("/users/ajaffe/Lieber/Projects/RNAseq/master_phenotype_list.csv",
	as.is=TRUE)
pheno$BrNum = paste0("Br", pheno$BRNum)
pd$pH = as.numeric(pheno$PH[match(pd$BrNum, pheno$BrNum)])

ind = which(!is.na(pd$pH))
mod2 = model.matrix(~ Dx + Age + Sex + 
	snpPC1 + snpPC5 + snpPC6 + snpPC9 + snpPC10 +
	RIN + PMI + pH , data=pd[ind,])
fitGene2 = lmFit(yGene[,ind], mod2)	
ebGene2 = ebayes(fitGene2)
adjStats = data.frame(log2FC = fitGene2$coef[,2],
	tstat = ebGene2$t[,2], 
	pval = ebGene2$t[,2])[names(outGene)[outGene$isExp],]
sum(p.adjust(adjStats$pval,"fdr") < 0.05)

pdf("plots/degradation_plots_geneLevel_adjModel.pdf")
par(mar=c(5,6,3,2), cex.lab=2,cex.axis=2,cex.main=1.6)
plot(adjStats$tstat, degradeStatsGene$polyA_Tstat,
	xlab="SZ vs Control (Age,Sex,Race,PMI,RIN,PH)", 
	ylab= "Degradation Stat",
	main = "T-statistics", pch=21,bg="grey")
legend("topleft", paste0("r=",signif(cor(
	adjStats$tstat, degradeStatsGene$polyA_Tstat),3)),
	cex=1.4)
abline(h=0, v=c(-6,6))


plot(adjStats$log2FC, degradeStatsGene$polyA_timeFC,
	xlab="SZ vs Control (Age,Sex,Race,PMI,RIN,PH)", 
	ylab= "Degradation Stat",
	main = "Fold Changes", pch=21,bg="grey")
abline(h=0, v=c(-0.25,0.25))
legend("topleft", paste0("r=",signif(cor(
	outGene$log2FC_qsva, degradeStatsGene$polyA_timeFC),3)),
	cex=1.4)
dev.off()