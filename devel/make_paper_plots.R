###
####
# R-devel
source("../eqtl_functions.R")

library(GenomicRanges)
library(limma)

# transcript
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/transcript/transcript_data_filtered_n495.rda")
# gene exon junction
load("/dcs01/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
# ERs
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda")

pd$ageGroup = cut(pd$Age, c(-1,0,1,10,20,50,100))

### add entrez ID
geneMapGR = makeGRangesFromDataFrame(geneMap,keep=TRUE)
jMap$EntrezID = geneMapGR$EntrezID[match(jMap$newGeneID,names(geneMapGR))]
tMap$EntrezID = geneMapGR$EntrezID[match(tMap$EnsemblGeneID,	names(geneMapGR))]

# add annotation to regions
dA = distanceToNearest(regions, geneMapGR)
regions$EnsemblGeneID = names(geneMapGR)[subjectHits(dA)]
regions$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
regions$EnsemblGeneID[mcols(dA)$distance > 100] = NA
regions$EntrezID[mcols(dA)$distance > 100] = NA

### filter for control
controlIndex = which(pd$Dx == "Control")
pd = pd[controlIndex,]
geneRpkm = geneRpkm[,controlIndex]
exonRpkm = exonRpkm[,controlIndex]
jRpkm = jRpkm[,controlIndex]
tFpkm = tFpkm[,controlIndex]
regionMat = regionMat[,controlIndex]

###### filter gene, exon, junction
exprsGeneIndex = which(rowMeans(geneRpkm) > 0.01)
exprsExonIndex = which(rowMeans(exonRpkm) > 0.1)
exonMap$coord = paste0("chr",exonMap$Chr, ":", exonMap$Start, "-",
	exonMap$End,"(",exonMap$Strand, ")")
exprsJxnIndex = which(rowMeans(jRpkm > 0.2) & jMap$code !="Novel")


###############################
### model, linear spline
fetal = ifelse(pd$Age < 0, 1,0)
birth = pd$Age
birth[birth < 0] = 0 # linear spline
infant = pd$Age - 1
infant[infant < 0] = 0 # linear spline
child = pd$Age - 10
child[child < 0] = 0 # linear spline
teen = pd$Age - 20
teen[teen < 0] = 0 # linear spline
adult = pd$Age - 50
adult[adult < 0] = 0 # linear spline

### adjust for race, assign, sex, and mito
mod = model.matrix(~Age + fetal + birth + infant +
	child + teen + adult + snpPC1 + snpPC2 + snpPC3 + 
	 Sex, data=pd)
	 
	
## load developmental stats
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/isoform_switch_devel_byFeature.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/devStats_controlSamples.rda")

########## #########
## PCA checks ###
geneRpkmClean = cleaningY(log2(geneRpkm+1), mod, P=8)

# PCA only on the consistent genes

allGenes = unique(unlist(lapply(statList, function(x) 
	x$EnsemblGeneID[which(x$p_bonf < 0.05)])))
allGenes = allGenes[!grepl("-", allGenes)] # drop fusions
allGenes = allGenes[!is.na(allGenes)] # universe
geneMat = sapply(statList, function(x) allGenes %in% x$EnsemblGeneID)
rownames(geneMat) = allGenes

rn = rownames(geneMat)[which(rowSums(geneMat) == 5)]
pcaCleanGene = prcomp(t(geneRpkmClean[rn,]))
pcaVars = getPcaVars(pcaCleanGene)[1:10]

pdf("plots/agePlotter_cleanPCs.pdf", h=3,w=4)
library(RColorBrewer)
for(i in 1:10) {
	agePlotter(pcaCleanGene$x[,i], pd$Age, mod, 
		ylab=paste0("PC", i, ": ", pcaVars[i],
			"% Var Expl"),	mainText="", smoothIt=FALSE)
}
dev.off()

####################################
#### dev plots ##################
####################################
xformData = list(Gene = log2(geneRpkm+1), Exon = log2(exonRpkm+1), 
	Transcript =  log2(tFpkm+1), Junction = log2(jRpkm+1), 
	DER = log2(regionMat+1))

## get features
downFeatures = lapply(statList, function(x) 
	names(x)[which(x$p_bonf < 0.05 & x$ageCorr < 0)])
upFeatures = lapply(statList, function(x) 
	names(x)[which(x$p_bonf < 0.05 & x$ageCorr > 0)])

## get values
downExprs = mapply(function(x, ind) x[ind,], 
	xformData, downFeatures,SIMPLIFY=FALSE) 
downExprs = do.call("rbind", downExprs)
downExprsScale = t(scale(t(downExprs), center=TRUE, scale=TRUE))

upExprs = mapply(function(x, ind) x[ind,], 
	xformData, upFeatures,SIMPLIFY=FALSE) 
upExprs = do.call("rbind", upExprs)
upExprsScale = t(scale(t(upExprs), center=TRUE, scale=TRUE))

## get means
downQuants = t(apply(downExprsScale,2, quantile, probs = c(0.5,0.25,0.75)))
upQuants = t(apply(upExprsScale,2, quantile, probs = c(0.5,0.25,0.75)))
allQuants = cbind(downQuants, upQuants)
colnames(allQuants) = paste0(rep(c("down_", "up_"),each=3), colnames(allQuants))

## smooth
allQuantsFit = apply(allQuants, 2, function(y)
	lm(y ~ mod[,2:8])$fitted)
oo = order(pd$Age)
aIndexes = split(oo, cut(pd$Age[oo], c(-1,0,10,100)))
fAges = seq(8,40,by=6)

pdf("plots/normalized_devFeatures.pdf",h=4,w=6)
## make plot
layout(matrix(c(1,1,1,1,1,2,2,2,3,3,3,3,3),nr = 1,byrow = TRUE))
## fetal
par(mar=c(5,6,2,2))
ii = aIndexes[[1]]
plot(allQuantsFit[ii,1] ~ pd$Age[ii], type="l",lwd=4,
	ylim=c(-2,2),bty="n", col="darkblue",ylab="",
		xaxt="n",xlab="",cex.axis=2)
axis(1,at=(fAges-40)/52,	labels = fAges, cex.axis=2)
polygon(c( pd$Age[ii],rev( pd$Age[ii])),
	c(allQuantsFit[ii,2],rev(allQuantsFit[ii,3])),
	col="lightblue",density = 25,border="darkblue")
lines(allQuantsFit[ii,4] ~  pd$Age[ii], type="l",lwd=4,
	col="darkorange")
polygon(c( pd$Age[ii],rev( pd$Age[ii])),
	c(allQuantsFit[ii,5],rev(allQuantsFit[ii,6])),
	col="orange",density = 25,border="darkorange")	

### infant
par(mar=c(5,0,2,0))
for(i in 2:3) {
	ii = aIndexes[[i]]
	plot(allQuantsFit[ii,1] ~ pd$Age[ii], type="l",lwd=4,
		ylim=c(-2,2),bty="n", cex.axis=2,
		col="darkblue",ylab="",xlab="",yaxt="n")
	polygon(c( pd$Age[ii],rev( pd$Age[ii])),
		c(allQuantsFit[ii,2],rev(allQuantsFit[ii,3])),
		col="lightblue",density = 25,border="darkblue")
	lines(allQuantsFit[ii,4] ~  pd$Age[ii], type="l",lwd=4,
		col="darkorange")
	polygon(c( pd$Age[ii],rev( pd$Age[ii])),
		c(allQuantsFit[ii,5],rev(allQuantsFit[ii,6])),
		col="orange",density = 25,border="darkorange")	
}
dev.off()

################
#### big pca ###
pcaDevFeatures = prcomp(t(rbind(downExprs,upExprs)))
pcaVarsFeatures = getPcaVars(pcaDevFeatures)[1:10]

pdf("plots/agePlotter_cleanPCs_devFeatures.pdf", h=3,w=4)
library(RColorBrewer)
for(i in 1:10) {
	agePlotter(pcaDevFeatures$x[,i], pd$Age, mod, 
		ylab=paste0("PC", i, ": ", pcaVarsFeatures[i],
			"% Var Expl"),	mainText="", smoothIt=FALSE)
}
dev.off()

#######################################
###### all features in switches #######
switchStats = mapply(function(x,y) {
	x[x$EnsemblGeneID %in% rownames(y)]
}, statList[-1], switchList, SIMPLIFY=FALSE)

## get features
downFeatures = lapply(statList, function(x) 
	names(x)[which(x$p_bonf < 0.05 & x$ageCorr < 0)])
upFeatures = lapply(switchStats, function(x) 
	names(x)[which(x$p_bonf < 0.05 & x$ageCorr > 0)])

## get values
downExprs = mapply(function(x, ind) x[ind,], 
	xformData[-1], downFeatures,SIMPLIFY=FALSE) 
downExprs = do.call("rbind", downExprs)
downExprsScale = t(scale(t(downExprs), center=TRUE, scale=TRUE))

upExprs = mapply(function(x, ind) x[ind,], 
	xformData[-1], upFeatures,SIMPLIFY=FALSE) 
upExprs = do.call("rbind", upExprs)
upExprsScale = t(scale(t(upExprs), center=TRUE, scale=TRUE))

## get means
downQuants = t(apply(downExprsScale,2, quantile, probs = c(0.5,0.25,0.75)))
upQuants = t(apply(upExprsScale,2, quantile, probs = c(0.5,0.25,0.75)))
allQuants = cbind(downQuants, upQuants)
colnames(allQuants) = paste0(rep(c("down_", "up_"),each=3), colnames(allQuants))

## smooth
allQuantsFit = apply(allQuants, 2, function(y)
	lm(y ~ mod[,2:8])$fitted)
oo = order(pd$Age)
aIndexes = split(oo, cut(pd$Age[oo], c(-1,0,10,100)))
fAges = seq(8,40,by=6)

pdf("plots/normalized_isoSwitchFeatures.pdf",h=4,w=6)
## make plot
layout(matrix(c(1,1,1,1,1,2,2,2,3,3,3,3,3),nr = 1,byrow = TRUE))
## fetal
par(mar=c(5,6,2,2))
ii = aIndexes[[1]]
plot(allQuantsFit[ii,1] ~ pd$Age[ii], type="l",lwd=4,
	ylim=c(-2,2),bty="n", col="darkblue",ylab="",
		xaxt="n",xlab="",cex.axis=2)
axis(1,at=(fAges-40)/52,	labels = fAges, cex.axis=2)
polygon(c( pd$Age[ii],rev( pd$Age[ii])),
	c(allQuantsFit[ii,2],rev(allQuantsFit[ii,3])),
	col="lightblue",density = 25,border="darkblue")
lines(allQuantsFit[ii,4] ~  pd$Age[ii], type="l",lwd=4,
	col="darkorange")
polygon(c( pd$Age[ii],rev( pd$Age[ii])),
	c(allQuantsFit[ii,5],rev(allQuantsFit[ii,6])),
	col="orange",density = 25,border="darkorange")	

### infant
par(mar=c(5,0,2,0))
for(i in 2:3) {
	ii = aIndexes[[i]]
	plot(allQuantsFit[ii,1] ~ pd$Age[ii], type="l",lwd=4,
		ylim=c(-2,2),bty="n", cex.axis=2,
		col="darkblue",ylab="",xlab="",yaxt="n")
	polygon(c( pd$Age[ii],rev( pd$Age[ii])),
		c(allQuantsFit[ii,2],rev(allQuantsFit[ii,3])),
		col="lightblue",density = 25,border="darkblue")
	lines(allQuantsFit[ii,4] ~  pd$Age[ii], type="l",lwd=4,
		col="darkorange")
	polygon(c( pd$Age[ii],rev( pd$Age[ii])),
		c(allQuantsFit[ii,5],rev(allQuantsFit[ii,6])),
		col="orange",density = 25,border="darkorange")	
}
dev.off()