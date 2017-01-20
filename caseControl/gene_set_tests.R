####
library(GenomicRanges)
library(GOstats)
library(clusterProfiler)

source("../eqtl_functions.R")

## load libd data
load("rdas/expressed_de_features.rda")

#####################
## gene set overlap #
#####################

### number significant
# sigStats = endoapply(outStatsExprs, function(x) x[x$fdr_qsva < 0.1])
sigStats = endoapply(outStatsExprs, 
	function(x) x[which(x$fdr_qsva < 0.1 & 
		sign(x$CMC_log2FC_qsva) == sign(x$log2FC_qsva) & 
		x$CMC_pval_qsva < 0.05)])
elementLengths(sigStats)
		# sigStatsSort = endoapply(outStatsExprs, function(x) x[order(x$pval_qsva)])

sigStatsAll = unlist(sigStats)
sigStatsAll$Feature = ss(names(sigStatsAll), "\\.")
sigStatsAll$Sign = ifelse(sigStatsAll$log2FC_qsva> 0,"szUp", "szDown")
sigStatsBySign = split(sigStatsAll, paste0(sigStatsAll$Feature,
	":", sigStatsAll$Sign))

### GO ####
bgGeneList = lapply(outStatsExprs[ss(names(sigStatsBySign), ":")], function(x) {
	unique(x$EntrezID[!is.na(x$EntrezID)])
})

## enr sets
geneList_fdr10 = lapply(sigStatsBySign, function(x) {
	unique(x$EntrezID[!is.na(x$EntrezID)])
})
elementLengths(geneList_fdr10)

#############################
## by feature #####

#### KEGG
keggList_fdr10 = mapply(function(g, bg) {
	ht=enrichKEGG(as.character(g), 
		organism="human", pvalueCutoff=1, 
		universe= as.character(bg),minGSSize=5,
		pAdjustMethod="BH", qvalueCutoff=1,readable=TRUE)
	summary(ht) 
}, geneList_fdr10, bgGeneList, SIMPLIFY=FALSE)

#### try the different GOs
goList_fdr10_MF = mapply(function(g, bg) {
	ht=enrichGO(as.character(g), 
		organism="human", pvalueCutoff=1, 
		universe= as.character(bg),minGSSize=5,
		pAdjustMethod="none", qvalueCutoff=1, readable=TRUE)
	summary(ht) 
}, geneList_fdr10, bgGeneList, SIMPLIFY=FALSE)
# lapply(goList_fdr10_MF, function(x) x[x$qvalue < 0.05,])

goList_fdr10_BP = mapply(function(g, bg) {
	ht=enrichGO(as.character(g), ont = "BP", 
		organism="human", pvalueCutoff=1, 
		universe= as.character(bg),minGSSize=5,
		pAdjustMethod="none", qvalueCutoff=1, readable=TRUE)
	summary(ht) 
}, geneList_fdr10, bgGeneList, SIMPLIFY=FALSE)

#################
##### merge #####
#################

### kegg, just one set
keggMat = do.call("rbind", lapply(keggList_fdr10[1:6],
	function(x) {
		x$SetSize = as.integer(ss(as.character(x$BgRatio), "/", 1))
		x[,c("ID", "Description","SetSize")]}))
keggMat = keggMat[!duplicated(keggMat$ID),]
rownames(keggMat) = keggMat$keggMat

keggMat2 = do.call("cbind", lapply(keggList_fdr10[1:6], function(x) 
		x[match(keggMat$ID,x$ID),c("pvalue", "qvalue")]))
keggMat = cbind(keggMat,keggMat2)
keggMat$"Junction:szDown.pvalue" = NA
keggMat$"Junction:szDown.qvalue" = NA
keggMat$"Junction:szUp.pvalue" = NA
keggMat$"Junction:szUp.qvalue" = NA
keggMat$Type = "KEGG"

###  go
goMat_MF = do.call("rbind", lapply(goList_fdr10_MF,
	function(x) {
		x$SetSize = as.integer(ss(as.character(x$BgRatio), "/", 1))
		x[,c("ID", "Description","SetSize")]}))
goMat_MF = goMat_MF[!duplicated(goMat_MF$ID),]
rownames(goMat_MF) = goMat_MF$ID

goMat_BP = do.call("rbind", lapply(goList_fdr10_BP,
		function(x) {
		x$SetSize = as.integer(ss(as.character(x$BgRatio), "/", 1))
		x[,c("ID", "Description","SetSize")]}))
goMat_BP = goMat_BP[!duplicated(goMat_BP$ID),]
rownames(goMat_BP) = goMat_BP$ID

goMat_MF2 = do.call("cbind", lapply(goList_fdr10_MF, function(x) 
		x[match(goMat_MF$ID,x$ID),c("pvalue", "qvalue")]))
rownames(goMat_MF2) = goMat_MF$ID
goMat_MF = cbind(goMat_MF, goMat_MF2)

goMat_BP2 = do.call("cbind", lapply(goList_fdr10_BP, function(x) 
		x[match(goMat_BP$ID,x$ID),c("pvalue", "qvalue")]))
rownames(goMat_BP2) = goMat_BP$ID
goMat_BP = cbind(goMat_BP, goMat_BP2)

## merge again
goMat = rbind(goMat_BP, goMat_MF)
goMat$Type = rep(c("BP", "MF"), 
	times = c(nrow(goMat_BP), nrow(goMat_MF)))
geneSetMat = rbind(keggMat, goMat)
rownames(geneSetMat) = geneSetMat$ID

geneSetMat = geneSetMat[,c(1:3,20,grep("qvalue", colnames(geneSetMat)))]
geneSetMat = geneSetMat[geneSetMat$SetSize < 5000,]
geneSetMat = geneSetMat[rowSums(geneSetMat[,5:12] < 0.05, na.rm=TRUE) > 0,]


## filter
sigMat =  geneSetMat[rowSums(geneSetMat[,5:12] < 0.05, na.rm=TRUE) > 1,]
	
qMat = -log10(sigMat[,c(5,7,9)])
rownames(qMat) = sigMat$Description
colnames(qMat) = ss(ss(colnames(qMat),"\\."),":")

library(lattice)
pdf("plots/geneSet_caseControl_heatmap_replication.pdf",h=6,w=10)
theSeq = seq(0,3,by=0.1) 
my.col <- colorRampPalette(c("white","darkblue"))(length(theSeq))
print(levelplot(t(as.matrix(qMat)), at = theSeq,pretty=TRUE, 
	col.regions = my.col, scales=list(y=list(cex=1.5), 
		x=list(rot=90, cex=2)),aspect="fill",
	ylab = "", xlab = ""))
dev.off()

geneSetMat = geneSetMat[order(apply(geneSetMat[,5:12],1,min,na.rm=TRUE)),]
write.csv(geneSetMat, file="tables/geneSetEnrichments_caseControl.csv"
#############################
### joint analysis ########
# geneEither_fdr10 = lapply(sigStats, function(x) 
	# x$EntrezID[!is.na(x$EntrezID)])


compareKegg = compareCluster(geneList_fdr10, 
	fun = "enrichKEGG",qvalueCutoff = 0.05, pvalueCutoff = 0.1)
plot(compareKegg,colorBy="qvalue",include = TRUE)

compareGoMf = compareCluster(geneList_fdr10, 
	fun = "enrichGO",ont = "MF", 
	qvalueCutoff = 0.05, pvalueCutoff = 0.1)
plot(compareGoMf,colorBy="qvalue")

compareGoBp = compareCluster(geneList_fdr10, 
	fun = "enrichGO",ont = "BP", 
	qvalueCutoff = 0.05, pvalueCutoff = 0.1)
plot(compareGoBp,colorBy="qvalue")


## trick into 1 list
library(lattice)

tmp  = rbind(compareKegg@compareClusterResult,
	compareGoMf@compareClusterResult,compareGoBp@compareClusterResult)
tmp$Type = rep(c("KEGG", "MF", "BP"), c(nrow(compareKegg@compareClusterResult),
	nrow(compareGoMf@compareClusterResult),nrow(compareGoBp@compareClusterResult)))
tmp$PlotValue = -log10(tmp$pvalue)
tmp$SetSize = as.numeric(ss(as.character(tmp$BgRatio), "/"))

## filter
tab = table(tmp$ID)
tab = tab[tab>1]
tmp = tmp[(tmp$Type == "KEGG" | tmp$ID %in% names(tab)),]
tmp = tmp[tmp$SetSize < 6000,]
tmp$Cluster=  factor(tmp$Cluster, 
	levels = c("Gene:szDown","Exon:szDown", "ER:szDown",
		"Gene:szUp", "ER:szUp",  "Exon:szUp"))
tmp = tmp[order(tmp$Type, tmp$PlotValue,decreasing=TRUE),]		
tmp$Description = factor(tmp$Description, levels = unique(tmp$Description))

pdf("plots/geneSet_caseControl_heatmap_replication.pdf",w=8)
theSeq = seq(0,7,by=0.1) 
my.col <- colorRampPalette(c("white","darkblue"))(length(theSeq))
print(levelplot(PlotValue ~ Cluster + Description, 
	data= tmp, at = theSeq,pretty=TRUE, 
	col.regions = my.col, scales=list(y=list(cex=1.5), 
		x=list(rot=90, cex=1.5)),
	ylab = "", xlab = ""))
dev.off()


