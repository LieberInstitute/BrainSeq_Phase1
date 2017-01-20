###
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

library(GenomicRanges)
library(limma)
library(sva)

# #### load data
# load("/dcs01/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")

# pdL = pd
# geneRpkmL = geneRpkm
# exonRpkmL = exonRpkm
# jRpkmL = jRpkm
# jMapL = jMap
# pdL$Dataset = "LIBD"

# ## filter LIBD to adult
# keepIndex = which(pdL$age > 16)
# pdL = pdL[keepIndex,]
# jRpkmL = jRpkmL[,keepIndex]
# exonRpkmL = exonRpkmL[,keepIndex]
# geneRpkmL = geneRpkmL[,keepIndex]

#### gene assignment rate
# a = read.delim("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/Counts/Ensembl_Genes_n746.counts.summary",
	# row.names=1)
# colnames(a) = ss(colnames(a), "_", 4)
# assignRate = as.numeric(a[1,] / colSums(a))
# pdL$totalAssignedGene = assignRate[match(pdL$RNum, colnames(a))]
# pdL$mitoRate = pdL$mitoMapped / (pdL$mitoMapped + pdL$totalMapped)

# #################
# ### CMC data #####
# load("/dcl01/lieber/ajaffe/PublicData/CMC/rpkmCounts_cmcDlpfc_szControl.rda")

# geneRpkmCmc = geneRpkm
# exonRpkmCmc = exonRpkm
# jRpkmCmc = jRpkm
# jMapCmc = jMap
# pdCmc = pd
# pdCmc$Dataset = "CMC"

# ##### process junctions #####
# ## drop out junctions not in LIBD samples
# dropIndex=which(rowSums(jRpkmL) == 0)
# jRpkmL = jRpkmL[-dropIndex,]
# jMapL = jMap[-dropIndex]

# ## number of junctions in common
# nam = intersect(names(jMapCmc) ,names(jMapL))

# # filter to common junctions
# jRpkmCmc = jRpkmCmc[nam,]
# jRpkmL = jRpkmL[nam,]
# jMap = jMapL[nam]

# ## filter expression to LIBD
# gIndex=which(rowMeans(geneRpkmL) > 0.1)
# geneRpkmL = geneRpkmL[gIndex,]
# geneRpkmCmc = geneRpkmCmc[gIndex,]
# geneMap = geneMap[gIndex,]

# eIndex = which(rowMeans(exonRpkmL) > 0.2)
# exonRpkmL = exonRpkmL[eIndex,]
# exonRpkmCmc = exonRpkmCmc[eIndex,]
# exonMap = exonMap[eIndex,]

# pdCmc$Age_of_Death[pdCmc$Age_of_Death == "90+"] = 90
# pdCmc$Age_of_Death = as.numeric(pdCmc$Age_of_Death)

# # number of reads assigned
# geneFn = paste0("/dcl01/lieber/ajaffe/PublicData/CMC/Counts/Gene/",
	# pdCmc$DLPFC_RNA_Sequencing_Sample_ID, "_Ensembl_v75_Genes.counts")
# names(geneFn) = pdCmc$DLPFC_RNA_Sequencing_Sample_ID
# geneStatList = lapply(paste0(geneFn, ".summary"), 
	# read.delim,row.names=1)
# geneStats = do.call("cbind", geneStatList)
# colnames(geneStats) = pdCmc$DLPFC_RNA_Sequencing_Sample_ID
# pdCmc$totalAssignedGene = as.numeric(geneStats[1,] / colSums(geneStats))

# save(pdL, geneRpkmL, exonRpkmL, jRpkmL, 
	# pdCmc, geneRpkmCmc, exonRpkmCmc, jRpkmCmc,
	# jMap, geneMap, exonMap, compress=TRUE,
	# file = "/dcs01/ajaffe/Brain/DLPFC_PolyA/szControl/data/CMC_LIBD_merged.rda")

#######################
## load back in data ##
load("/dcs01/ajaffe/Brain/DLPFC_PolyA/szControl/data/CMC_LIBD_merged.rda")
colnames(pdL)[18] = "mapRate"

## models
modL = model.matrix(~ Dx + age + Sex + snpPC1 + 
	snpPC2 + snpPC3 + mitoRate + mapRate + 
	totalAssignedGene, data=pdL)
fitGeneL = lmFit(log2(geneRpkmL+1), modL)
ebGeneL = ebayes(fitGeneL)

modCmc = model.matrix(~ Dx + Age_of_Death + Institution + 
	Gender + Ethnicity + totalAssignedGene, data=pdCmc)
fitGeneCmc = lmFit(log2(geneRpkmCmc+1), modCmc)
ebGeneCmc = ebayes(fitGeneCmc)

### young CMC
youngIndex = which(pdCmc$Age_of_Death < 65)
modCmcYoung = model.matrix(~ Dx + Age_of_Death + Institution + 
	Gender + Ethnicity + totalAssignedGene, data=pdCmc[youngIndex,])
fitGeneCmcYoung = lmFit(log2(geneRpkmCmc[,youngIndex]+1), modCmcYoung)
ebGeneCmcYoung = ebayes(fitGeneCmc)

compareMods = function(fit1,eb1, fit2,eb2,coi1 = 2,
	coi2=2, pCut = 0.05, ptype = "fdr") {
	
	cc = cor(fit1$coef[,coi1], fit2$coef[,coi2])
	signs = mean(sign(fit1$coef[,coi1]) == 
		sign(fit2$coef[,coi2]))

	p1 = p.adjust(eb1$p[,coi1], ptype)
	p2 = p.adjust(eb2$p[,coi2], ptype)
	n1 = sum(p1 < pCut)
	n2 = sum(p2 < pCut)
	
	signsPlusP = sum(sign(fit1$coef[,coi1]) == 
		sign(fit2$coef[,coi2]) & eb1$p[,coi1] < pCut & 
		eb2$p[,coi2]  < pCut) / (sum(eb1$p[,coi1] < pCut & 
		eb2$p[,coi2]  < pCut))
	signsPlusSig = sum(sign(fit1$coef[,coi1]) == 
		sign(fit2$coef[,coi2]) & p1 < pCut & 
		p2 < pCut) # / (sum(p1 < fdrCut & p2 < fdrCut))
	data.frame(pType = ptype, numSig1 = n1, 
		numSig2 = n2,corr = cc, concord = signs, 
		concordPlusMarg = signsPlusP, 
		concordPlusFdr = signsPlusSig)
}	
univStats = compareMods(fitGeneL, ebGeneL, fitGeneCmc, ebGeneCmc)
compareMods(fitGeneL, ebGeneL, fitGeneCmc, ebGeneCmc, ptype="bonf")
univStatsYoung = compareMods(fitGeneL, ebGeneL, fitGeneCmcYoung, ebGeneCmcYoung)
compareMods(fitGeneL, ebGeneL, fitGeneCmcYoung, ebGeneCmcYoung, ptype="bonf")

#### load degradation data
load("/users/ajaffe/Lieber/Projects/RNAseq/Degradation/rdas/brain_DIGs_3type.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/Degradation/rdas/degradeStats_byLibrary.rda")

# get qSVs, CMC
digIndexCmc = which(rownames(geneMap) %in% rownames(digsRibo))
kCmc = num.sv(log2(geneRpkmCmc[digIndexCmc,]+1), modCmc)
qsvCmc = prcomp(t(log2(geneRpkmCmc[digIndexCmc,]+1)))$x[,1:kCmc]
kCmcYoung = num.sv(log2(geneRpkmCmc[digIndexCmc,youngIndex]+1),
	modCmcYoung)
qsvCmcYoung = prcomp(t(log2(geneRpkmCmc[digIndexCmc,youngIndex]+1)))$x[,1:kCmcYoung]

# get qSVs, LIBD
digIndexL= which(rownames(geneMap) %in% rownames(digsPolyA))
kL = num.sv(log2(geneRpkmL[digIndexL,]+1), modL)
qsvL = prcomp(t(log2(geneRpkmL[digIndexL,]+1)))$x[,1:kL]

## full fits
fitGeneCmcQ = lmFit(log2(geneRpkmCmc+1), 
	cbind(modCmc, qsvCmc))
ebGeneCmcQ = ebayes(fitGeneCmcQ)
fitGeneCmcYoungQ = lmFit(log2(geneRpkmCmc[,youngIndex]+1), 
	cbind(modCmcYoung, qsvCmcYoung))
ebGeneCmcYoungQ = ebayes(fitGeneCmcYoungQ)

fitGeneLQ= lmFit(log2(geneRpkmL+1), 
	cbind(modL, qsvL))
ebGeneLQ = ebayes(fitGeneLQ)

qsvStats = compareMods(fitGeneLQ, ebGeneLQ, fitGeneCmcQ, ebGeneCmcQ)
compareMods(fitGeneLQ, ebGeneLQ, fitGeneCmcQ, ebGeneCmcQ,ptype="bonf")
qsvStatsYoung = compareMods(fitGeneLQ, ebGeneLQ, 
	fitGeneCmcYoungQ, ebGeneCmcYoungQ)
compareMods(fitGeneLQ, ebGeneLQ, 
	fitGeneCmcYoungQ, ebGeneCmcYoungQ,ptype="bonf")

sigIndexL = which(p.adjust(ebGeneLQ$p[,2], "bonf") < 0.05)
sigIndexL_fdr = which(p.adjust(ebGeneLQ$p[,2], "fdr") < 0.05)

cbind(fitGeneLQ$coef[sigIndexL,2], geneMap[sigIndexL,"Symbol"])
plot(fitGeneCmcYoungQ$coef[sigIndexL,2], 
	fitGeneLQ$coef[sigIndexL,2])
abline(v=0,h=0)
geneMap[names(which(ebGeneCmcYoungQ$p[sigIndexL,2] < 0.05 & 
	(sign(fitGeneCmcYoungQ$coef[sigIndexL,2]) ==
		sign(fitGeneLQ$coef[sigIndexL,2])))),]

plot(fitGeneCmcYoungQ$coef[sigIndexL_fdr,2], 
	fitGeneLQ$coef[sigIndexL_fdr,2])
abline(v=0,h=0)

table(ebGeneCmcYoungQ$p[sigIndexL_fdr,2] < 0.05 & 
	(sign(fitGeneCmcYoungQ$coef[sigIndexL_fdr,2]) ==
		sign(fitGeneLQ$coef[sigIndexL_fdr,2])))

### libd by race
rIndexesL = split(1:nrow(pdL), pdL$Race)
fitGeneRaceL = lapply(rIndexesL, function(ii) {
	 lmFit(log2(geneRpkmL[,ii]+1), cbind(modL, qsvL)[ii,])	
})
ebGeneRaceL = lapply(fitGeneRaceL,ebayes)

rIndexesCmc = split((1:nrow(pdCmc))[youngIndex], 
	pdCmc$Ethnicity[youngIndex])
rIndexesCmc = rIndexesCmc[c("Caucasian", "African-American")]
fitGeneRaceCmc = lapply(rIndexesCmc, function(ii) {
	 lmFit(log2(geneRpkmCmc[,ii]+1), cbind(modCmc[,-(7:10)], qsvCmc)[ii,])	
})
ebGeneRaceCmc = lapply(fitGeneRaceCmc,ebayes)

plot(sapply(fitGeneRaceL, function(x) x$coef[,2]))
plot(sapply(fitGeneRaceCmc, function(x) x$coef[,2]))

plot(fitGeneRaceL$CAUC$coef[,2], 
	fitGeneRaceCmc$Caucasian$coef[,2])
plot(fitGeneRaceL$AA$coef[,2], 
	fitGeneRaceCmc$"African-American"$coef[,2])
	
plot(sapply(ebGeneRaceL, function(x) x$t[,2]))
cor(sapply(ebGeneRaceL, function(x) x$t[,2]))
colSums(sapply(ebGeneRaceL, function(x) p.adjust(x$p[,2],"fdr")) < 0.05)


############
### explore plots	
plot(ebGeneLQ$t[,2], 
	degradeStats[rownames(geneMap), "polyA_Tstat"])
plot(ebGeneL$t[,2], 
	degradeStats[rownames(geneMap), "polyA_Tstat"])
plot(ebGeneCmcQ$t[,2], 
	degradeStats[rownames(geneMap), "ribo_Tstat"])
plot(ebGeneCmc$t[,2], 
	degradeStats[rownames(geneMap), "ribo_Tstat"])
	
fitRinL = lmFit(log2(geneRpkmL+1), model.matrix(~pdL$RIN))
plot(fitGeneLQ$coef[,2], fitRinL$coef[,2])
