######
library(minfi)
library(GenomicRanges)

#### load data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/snpData_LIBD_szControls_cleaned.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/transcript/transcript_data_filtered_n495.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")

### read in eQTLs to filter
pgcEqtl = read.csv("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/tables/suppTable_allPgcEqtlsToIndexSnps.csv",
	as.is=TRUE)

#################	
#### FILTER #####

## expression ###
geneRpkmSub = geneRpkm[rownames(geneRpkm) %in% pgcEqtl$Feature,]
geneMapSub = geneMap[rownames(geneMap) %in% pgcEqtl$Feature,]
exonRpkmSub = exonRpkm[rownames(exonRpkm) %in% pgcEqtl$Feature,]
exonMapSub = exonMap[rownames(exonMap) %in% pgcEqtl$Feature,]
jRpkmSub = jRpkm[rownames(jRpkm) %in% pgcEqtl$Feature,]
jMapSub = jMap[names(exonMap) %in% pgcEqtl$Feature]
tFpkmSub = tFpkm[rownames(tFpkm) %in% pgcEqtl$Feature,]
jMapSub = jMap[names(exonMap) %in% pgcEqtl$Feature]
regionMatSub = regionMat[rownames(regionMat) %in% pgcEqtl$Feature,]
regionsSub = regions[names(regions) %in% pgcEqtl$Feature]

exprsSub = rbind(geneRpkmSub, exonRpkmSub, jRpkmSub, 
	tFpkmSub, regionMatSub)
exprsSub = exprsSub[pgcEqtl$Feature,]

#### snp 
snpSub = snp[pgcEqtl$SNP,]
snpMapSub = snpMap[match(pgcEqtl$SNP, snpMap$SNP),]

#################################
######### project for PCs #######

aIndex = which(pd$Age > 13)
pd2 = pd[aIndex,]
geneRpkm2 = as.matrix(log2(geneRpkm[,aIndex]+1))
exonRpkm2 = as.matrix(log2(exonRpkm[,aIndex]+1))
jRpkm2 = as.matrix(log2(jRpkm[,aIndex]+1))
tFpkm2 = as.matrix(log2(tFpkm[,aIndex]+1))
regionMat2 = as.matrix(log2(regionMat[,aIndex]+1))

##### filter gene, exon, junction for PCA
gExpIndex=which(rowMeans(geneRpkm2) > 0.01 & 
	geneMap$Chr %in% paste0("chr", c(1:22, "X","Y")))
geneRpkm2 = geneRpkm2[gExpIndex,]
geneMap2 = geneMap[gExpIndex,]

eExpIndex=which(rowMeans(exonRpkm2) > 0.1 & 
	exonMap$Chr %in% paste0("chr", c(1:22,"X","Y")))
exonRpkm2 = exonRpkm2[eExpIndex,]
exonMap2 = exonMap[eExpIndex,]

jExpIndex=which(rowMeans(jRpkm2) > 0.2 & jMap$code != "Novel")
jRpkm2 = jRpkm2[jExpIndex,]
jMap2 = jMap[jExpIndex]

#### do PCA 
pcaGene = prcomp(t(as.matrix(geneRpkm2)))
pcaExon = prcomp(t(as.matrix(exonRpkm2)))
pcaJxn = prcomp(t(as.matrix(jRpkm2)))
pcaErs = prcomp(t(as.matrix(regionMat2)))
pcaTx = prcomp(t(as.matrix(tFpkm2)))

####################
##### project ######

# gene
geneRpkmScaled = scale(as.matrix(t(
		log2(geneRpkm[rownames(geneRpkm2),]+1))),
	pcaGene$center, pcaGene$scale) 
newGenePCs = geneRpkmScaled %*% pcaGene$rotation 
plot(newGenePCs[aIndex,1], pcaGene$x[,1]) 

# exon
exonRpkmScaled = scale(as.matrix(t(
		log2(exonRpkm[rownames(exonRpkm2),]+1))),
	pcaExon$center, pcaExon$scale) 
newExonPCs = exonRpkmScaled %*% pcaExon$rotation

# junction
jRpkmScaled = scale(as.matrix(t(
		log2(jRpkm[rownames(jRpkm2),]+1))),
	pcaJxn$center, pcaJxn$scale) 
newJxnPCs = jRpkmScaled %*% pcaJxn$rotation

# ER
regionMatScaled = scale(as.matrix(t(
		log2(regionMat[rownames(regionMat2),]+1))),
	pcaErs$center, pcaErs$scale) 
newErPCs = regionMatScaled %*% pcaErs$rotation

# Tx
tFpkmScaled = scale(as.matrix(t(
		log2(tFpkm[rownames(tFpkm2),]+1))),
	pcaTx$center, pcaTx$scale) 
newTxPCs = tFpkmScaled %*% pcaTx$rotation

################################
##### just PCs within fetal ####
fIndex = which(pd$Age < 0)
fetalExprsList = list(Gene = geneRpkm[gExpIndex,fIndex],
	Exon = exonRpkm[eExpIndex,fIndex],
	Transcript = tFpkm[,fIndex],
	Junction = jRpkm[jExpIndex, fIndex],
	ER = regionMat[,fIndex])
## xform
fetalExprsList = lapply(fetalExprsList, function(x) as.matrix(log2(x+1)))

## do PCA
fetalPcaList = parallel::mclapply(fetalExprsList, function(x) {
	prcomp(t(x))$x[,1:15]
}, mc.cores=5)

## save stuff ##
save(snpMapSub, snpSub, pd, geneRpkmSub, exonRpkmSub, 
	jRpkmSub, tFpkmSub, regionMatSub, newGenePCs, newExonPCs,
	newJxnPCs, newErPCs, newTxPCs, fetalPcaList, compress=TRUE,
	file = "rdas/pgcEqtl_indexSnps_subsetData.rda")
	
	
#######################
## bring in meth ######

# load data
load("/dcs01/lieber/ajaffe/Brain/DNAm/ECD2014/Mset_withReps.rda")
Mset = updateObject(Mset)

# filter probes on chrX,Y and containing snps at SBE and target CpG
Mset = addSnpInfo(Mset)
Mset = Mset[is.na(getSnpInfo(Mset)$CpG_rs) & 
	is.na(getSnpInfo(Mset)$SBE_rs),]

## filter people
keepIndex = which(pData(Mset)$bestQC & 
	pData(Mset)$Gender == pData(Mset)$predictedSex)
Mset = Mset[,keepIndex]





####################
## fetal effects ###
fetalIndex = which(pd$Age < 0)







