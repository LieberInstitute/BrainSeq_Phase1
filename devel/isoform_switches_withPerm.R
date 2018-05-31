###
source("../eqtl_functions.R")

library(GenomicRanges)
library(clusterProfiler)
library(limma)
library(parallel)

##### summary statistics ####
load("rdas/devStats_controlSamples.rda")

#################################
###### Create isoform switches ##########

# drop genes
rangeList = lapply(statList[-1], function(x) {
	cat(".")
	xSig = x[which(x$p_bonf < 0.05),]
	xList = split(xSig, factor(xSig$EnsemblGeneID,
		levels = unique(x$EnsemblGeneID)))
	xList = xList[lengths(xList) > 0]
	theRange = t(sapply(split(xSig$ageCorr,factor(xSig$EnsemblGeneID,
		levels = unique(xSig$EnsemblGeneID))), range))
	theRange = as.data.frame(theRange)
	colnames(theRange) = c("negCorr","posCorr") 
	
	## get min and max feature
	mins = xSig[order(xSig$ageCorr)]
	mins = mins[!duplicated(mins$EnsemblGeneID) & 
		!is.na(mins$EnsemblGeneID)]
	theRange$minFeature = names(mins)[match(rownames(theRange),
		mins$EnsemblGeneID)] 
	maxs = xSig[order(xSig$ageCorr,decreasing=TRUE)]
	maxs = maxs[!duplicated(maxs$EnsemblGeneID) & 
		!is.na(maxs$EnsemblGeneID)]
	theRange$maxFeature = names(maxs)[match(rownames(theRange),
		maxs$EnsemblGeneID)] 
		
	## other metrics
	theRange$numFeatures = table(x$EnsemblGeneID)[rownames(theRange)]
	theRange$numSigFeatures = lengths(xList)
	theRange$Symbol = x$Symbol[match(rownames(theRange), x$EnsemblGeneID)]
	theRange$EntrezID = x$EntrezID[match(rownames(theRange), x$EnsemblGeneID)]
	
	return(theRange)
})

### significant switches
switchList = lapply(rangeList, function(x) {
	x$corDiff = x$posCorr - x$negCorr
	x[which(x$negCorr < 0 & x$posCorr > 0 ),]
})
sapply(switchList, nrow)

## permutation ###
set.seed(45232)
B=100
nullCorList = mclapply(statList[-1], function(x) { 
	
	x$factorGene = factor(x$EnsemblGeneID,
		levels = unique(x$EnsemblGeneID))
	## permute gene
	nullCorList = lapply(1:B, function(i) {
		if(i %% 20 == 0) cat(".")
		ageCorrStar = sample(x$ageCorr)
		dat = data.frame(t(sapply(split(ageCorrStar, 
			x$factorGene),range)))
		names(dat) = c("negCorr", "posCorr")
		dat
	})
	nullCor = do.call("rbind", nullCorList)
})

sigSwitchList = mapply(function(obs, null) {
	cat(".")
	N = nrow(null) # number of null
	G = nrow(null)/B 

	## make null stats
	nullDiff = null[,2] - null[,1] # corr diff
	
	## must switch sign
	nullDiff[sign(null[,1]) == sign(null[,2])] = NA
	nullDiff[is.na(nullDiff)] = 0

	## reshape
	nullMat = matrix(nullDiff, nrow = length(G),
		ncol = B, byrow=FALSE)
	
	obs$empPvalue = edge.pvalue(obs$corDiff, nullMat, pool=TRUE)
	sig = obs[obs$empPvalue < 0.05,]
	
	sig = sig[order(sig$corDiff, decreasing=TRUE),]
	
	return(sig)
}, switchList, nullCorList, SIMPLIFY=FALSE)
	
# # save
save(sigSwitchList,	file="rdas/isoform_switch_devel_byFeature_withPerm.rda")

##################
## brainspan replication
load("rdas/signif_dev_stats_plusBrainSpan.rda")

sigSwitchList = mapply(function(x,y) {
	cat(".")
	x$negCorr_span = y$ageCorr_span[match(x$minFeature, names(y))]
	x$posCorr_span = y$ageCorr_span[match(x$maxFeature, names(y))]
	x$corDiff_span = x$posCorr_span - x$negCorr_span
	return(x)
}, x = sigSwitchList, y=statListSig[-1], SIMPLIFY=FALSE)

## how many genes are stil significan
sapply(switchList, nrow)
sapply(sigSwitchList, nrow)

sapply(switchList, function(x) round(100*quantile(x$posCorr),2))
sapply(switchList, function(x) round(100*quantile(x$negCorr),2))
sapply(sigSwitchList, function(x) round(100*quantile(x$posCorr),2))
sapply(sigSwitchList, function(x) round(100*quantile(x$negCorr),2))

pdf("plots/isoswitch_featureByDataset.pdf")
par(mar=c(5,6,3,2), cex.axis=2,cex.lab=2,cex.main=1.8)
for(i in seq(along=sigSwitchList)) {
	plot(sigSwitchList[[i]]$negCorr,
		sigSwitchList[[i]]$negCorr_span,
		xlim = c(-1,1), ylim = c(-1,1),
		pch = 21, bg="blue",main = names(sigSwitchList)[i],
		xlab = "LIBD: Age Corr", ylab="BrainSpan: Age Corr")
	points(sigSwitchList[[i]]$posCorr,
		sigSwitchList[[i]]$posCorr_span,
		xlim = c(-1,1), ylim = c(-1,1),
		pch = 21, bg="red")
	abline(v=0,h=0)
}
dev.off()

pdf("plots/isoswitch_deltaCorrelation_byFeature.pdf")
par(mar=c(5,6,3,2), cex.axis=2,cex.lab=2,cex.main=1.8)
for(i in seq(along=sigSwitchList)) {
	plot(sigSwitchList[[i]]$corDiff,
		sigSwitchList[[i]]$corDiff_span,
		main = names(sigSwitchList)[i],
		xlim = c(-2,2), ylim = c(-2,2),pch=21,bg="grey",
		xlab = "LIBD: Delta Corr", ylab="BrainSpan: Delta Corr")
	abline(v=0,h=0)
}
dev.off()


## check replication rates by cut
corCut= seq(0,1.5,by=0.01)
checkCor = t(sapply(corCut, function(cc) {
	sapply(sigSwitchList, function(x) {
		xx = x[x$corDiff > cc,]
		length(which(xx$negCorr_span < 0 & 
			xx$posCorr_span > 0 ))/nrow(xx)
	})
}))
matplot(corCut, checkCor)

############################################
#### gene ontology on genes that switch ####
load("rdas/isoform_switch_devel_byFeature_withPerm.rda")

geneSwitchList = lapply(sigSwitchList, rownames)
geneSwitchList = lapply(geneSwitchList, function(x) 
	x[!grepl("-",x) & !is.na(x)]) # non fusion 

## venn diagram of the IDs by the 5 features
allGenesSwitch = unique(unlist(geneSwitchList))
geneMatSwitch = sapply(geneSwitchList, function(x) allGenesSwitch %in% x)
rownames(geneMatSwitch) = allGenesSwitch

pdf("plots/venn_geneIDs_devChanges_withSwitch.pdf",h=4.5,w=4.5)
vennDiagram(vennCounts(geneMatSwitch))
dev.off()

### write CSV
geneSwitchAll = do.call("rbind", sigSwitchList)
geneSwitchAll$Type = ss(rownames(geneSwitchAll), "\\.")
geneSwitchAll$EnsemblID = ss(rownames(geneSwitchAll), "\\.",2)
rownames(geneSwitchAll)= NULL
geneSwitchAll = geneSwitchAll[,c(11,10,1:9)]
write.csv(geneSwitchAll, file="tables/suppTable5_isoformSwitches.csv",
	row.names=FALSE)
	
########################################
## get entrez id
entrezBgList = lapply(statList[-1], function(x) {
	o = x$EntrezID[!is.na(x$p_bonf)]
	unique(o[!is.na(o)])
})
lengths(entrezBgList)

entrezGeneSwitchList = lapply(sigSwitchList, function(x) {
	unique(x$EntrezID[!is.na(x$EntrezID)])
})
lengths(entrezGeneSwitchList)

## also just regulated genes
entrezBgList_dev = lapply(statList[-1], function(x) {
	o = x$EntrezID[which(x$p_bonf < 0.05)]
	unique(o[!is.na(o)])
})
lengths(entrezBgList_dev)
(lengths(entrezBgList)-lengths(entrezBgList_dev))/lengths(entrezBgList)


############################
### kegg on switches #######

## kegg on dev reg
keggListSwitch_dev = mapply(function(g, bg) {
	ht=enrichKEGG(as.character(g), 
		organism="human", pvalueCutoff=1, 
		universe= as.character(bg),minGSSize=5,
		pAdjustMethod="none", qvalueCutoff=1)
	as.data.frame(ht) 
}, entrezGeneSwitchList, entrezBgList_dev, SIMPLIFY=FALSE)

keggSwitchMat = do.call("rbind", lapply(keggListSwitch_dev,
	function(x) {
		x$SetSize = as.integer(ss(as.character(x$BgRatio), "/", 1))
		x[,c("ID", "Description","SetSize")]}))
keggSwitchMat = keggSwitchMat[!duplicated(keggSwitchMat$ID),]
rownames(keggSwitchMat) = keggSwitchMat$ID

keggSwitchMat2 = do.call("cbind", lapply(keggListSwitch_dev, function(x) 
		x[match(keggSwitchMat$ID,x$ID),c("pvalue", "qvalue")]))
rownames(keggSwitchMat2) = keggSwitchMat$ID
keggSwitchMat = cbind(keggSwitchMat, keggSwitchMat2)
keggSwitchMat$Type = "KEGG"

## numbers	
colSums(keggSwitchMat[,grep("qvalue", colnames(keggSwitchMat))] < 0.05,
	na.rm=TRUE)
table(rowSums(keggSwitchMat[,grep("qvalue", colnames(keggSwitchMat))] < 0.05,
	na.rm=TRUE))

###################################
#### gene ontology on switches ####
###################################

## development background
goListSwitch_MF_dev = mapply(function(g, bg) {
	ht=enrichGO(as.character(g), 
		OrgDb="org.Hs.eg.db", pvalueCutoff=1, 
		universe= as.character(bg),minGSSize=5,
		pAdjustMethod="none", qvalueCutoff=1)
	as.data.frame(ht) 
}, entrezGeneSwitchList, entrezBgList_dev, SIMPLIFY=FALSE)

goListSwitch_BP_dev = mapply(function(g, bg) {
	ht=enrichGO(as.character(g), ont = "BP", 
		OrgDb="org.Hs.eg.db", pvalueCutoff=1, 
		universe= as.character(bg),minGSSize=5,
		pAdjustMethod="none", qvalueCutoff=1, readable=TRUE)
	as.data.frame(ht) 
}, entrezGeneSwitchList, entrezBgList_dev, SIMPLIFY=FALSE)

################
### combine ####

## make matrix
goSwitchMat_MF = do.call("rbind", lapply(goListSwitch_MF_dev,
	function(x) {
		x$SetSize = as.integer(ss(as.character(x$BgRatio), "/", 1))
		x[,c("ID", "Description","SetSize")]}))
goSwitchMat_MF = goSwitchMat_MF[!duplicated(goSwitchMat_MF$ID),]
rownames(goSwitchMat_MF) = goSwitchMat_MF$ID

goSwitchMat_BP = do.call("rbind", lapply(goListSwitch_BP_dev,
		function(x) {
		x$SetSize = as.integer(ss(as.character(x$BgRatio), "/", 1))
		x[,c("ID", "Description","SetSize")]}))
goSwitchMat_BP = goSwitchMat_BP[!duplicated(goSwitchMat_BP$ID),]
rownames(goSwitchMat_BP) = goSwitchMat_BP$ID

goSwitchMat_MF2 = do.call("cbind", lapply(goListSwitch_MF_dev, function(x) 
		x[match(goSwitchMat_MF$ID,x$ID),c("pvalue", "qvalue")]))
rownames(goSwitchMat_MF2) = goSwitchMat_MF$ID
goSwitchMat_MF = cbind(goSwitchMat_MF, goSwitchMat_MF2)

goSwitchMat_BP2 = do.call("cbind", lapply(goListSwitch_BP_dev, function(x) 
		x[match(goSwitchMat_BP$ID,x$ID),c("pvalue", "qvalue")]))
rownames(goSwitchMat_BP2) = goSwitchMat_BP$ID
goSwitchMat_BP = cbind(goSwitchMat_BP, goSwitchMat_BP2)

## merge again
goSwitchMat = rbind(goSwitchMat_BP, goSwitchMat_MF)
goSwitchMat$Type = rep(c("BP", "MF"), 
	times = c(nrow(goSwitchMat_BP), nrow(goSwitchMat_MF)))

colSums(goSwitchMat[,grep("qvalue", colnames(keggSwitchMat))] < 0.05,
	na.rm=TRUE)
table(rowSums(goSwitchMat[,grep("qvalue", colnames(keggSwitchMat))] < 0.05,
	na.rm=TRUE))

##################################
### filter to dev and q-value
geneSetSwitch_dev = rbind(goSwitchMat[,c(1:3,12,
	grep("qvalue",names(goSwitchMat)))],
		keggSwitchMat[,c(1:3,12,
			grep("qvalue",names(goSwitchMat)))])
geneSetSwitch_dev = geneSetSwitch_dev[order(rowMeans(geneSetSwitch_dev[,5:8])),]
geneSetSwitch_dev = geneSetSwitch_dev[geneSetSwitch_dev$SetSize < 5000,]

colSums(geneSetSwitch_dev[,grep("qvalue",
	names(geneSetSwitch_dev))] < 0.05,na.rm=TRUE)
	
write.csv(geneSetSwitch_dev, file="tables/geneSets_isoformSwitches_withPerm.csv",
	row.names=FALSE)
	
## filter to only significant
numSig = rowSums(geneSetSwitch_dev[,grep("qvalue",
	names(geneSetSwitch_dev))] < 0.05,na.rm=TRUE)
sigGeneSetSwitch_dev = geneSetSwitch_dev[numSig > 0,]

write.csv(sigGeneSetSwitch_dev, row.names=FALSE,
	file="tables/geneSets_isoformSwitches_onlySigInOne_withPerm.csv")

### make figure for paper
keggMatSig = sigGeneSetSwitch_dev[sigGeneSetSwitch_dev$Type == "KEGG",]
	
library(lattice)
qMat = -log10(keggMatSig[,5:8])
rownames(qMat) = keggMatSig$Description
colnames(qMat) = ss(colnames(qMat),"\\.")

## reorder to old figure
oldSwitch = read.csv("tables/geneSets_isoformSwitches_onlySigInOne.csv",as.is=TRUE)
keggMatSigOld = oldSwitch[oldSwitch$Type == "KEGG",]
qMat = qMat[keggMatSigOld$Description,]

pdf("plots/geneSet_isoSwitch_heatmap_KEGG_withPerm.pdf",	useDingbats=FALSE,w=10)
theSeq = seq(0,7,by=0.1) 
my.col <- colorRampPalette(c("white","darkblue"))(length(theSeq))
print(levelplot(t(as.matrix(qMat)), at = theSeq,pretty=TRUE, 
	col.regions = my.col, scales=list(y=list(cex=1.5), 
		x=list(rot=90, cex=1.5)),aspect="fill",
	ylab = "", xlab = ""))
dev.off()

## summary statistics
colSums(sigGeneSetSwitch_dev[,grep("qvalue", 
	names(sigGeneSetSwitch_dev))] < 0.05, na.rm=TRUE)
table(rowSums(sigGeneSetSwitch_dev[,grep("qvalue", 
	names(sigGeneSetSwitch_dev))] < 0.05, na.rm=TRUE))
