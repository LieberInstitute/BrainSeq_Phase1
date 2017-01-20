###
library(GenomicRanges)

## load developmental stats
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/isoform_switch_devel_byFeature.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/devStats_controlSamples.rda")

## load case-control stats
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/caseControl/rdas/DE_statistics_adjAndQsva.rda")

### add entrez ID
outJxn$EntrezID = outGene$EntrezID[match(outJxn$newGeneID,	names(outGene))]
outTx$EntrezID = outGene$EntrezID[match(outTx$EnsemblGeneID, names(outGene))]

dA = distanceToNearest(outEr, outGene)
outEr$EntrezID = outGene$EntrezID[subjectHits(dA)]
outEr$EntrezID[mcols(dA)$distance > 100] = NA

### merge together into one GRList
## only some columns
n = c("EntrezID", "meanExprsAdult", "isExp", 
	"log2FC_adj", "tstat_adj", "pval_adj", "fdr_adj", 
	"log2FC_qsva" , "tstat_qsva","pval_qsva", "fdr_qsva",
	"log2FC_pca" , "tstat_pca","pval_pca", "fdr_pca")

## merge
dxStats = c(Gene = outGene[,n], Exon = outExon[,n],
	Transcript = outTx[,n],Junction = outJxn[,n],  ER = outEr[,n])
dxStats = GRangesList(dxStats)

########################################################
##### enrichment among developmentally regulated gene ##
########################################################
dxStats = GRangesList(mapply(function(x,y) {
	cat(".")
	y$devRegulated = x$p_bonf < 0.05
	return(y)
}, statList, dxStats))

## expressed only
dxStatsExprs = endoapply(dxStats, function(x) x[x$isExp])

################################
### add info about directionality
dxStats2 = GRangesList(mapply(function(x,y) {
	cat(".")
	y$negCor =names(y) %in% names(x)[x$p_bonf < 0.05 & x$ageCorr < 0]
	y$posCor =names(y) %in% names(x)[x$p_bonf < 0.05 & x$ageCorr > 0]
	return(y)
}, statList, dxStats))

sapply(dxStats2, function(x) table(x$negCor))
sapply(dxStats2, function(x) table(x$posCor))

#### filter to expression
dxStatsExprs = endoapply(dxStats2, function(x) x[x$isExp])
sapply(dxStatsExprs, function(x) table(x$negCor))
sapply(dxStatsExprs, function(x) table(x$posCor))

allDxStatsExprs = unlist(dxStatsExprs)

### directionality
xx = lapply(dxStatsExprs, function(x) {
	tmp = x$posCor
	tmp[x$negCor == 1] = -1
	table(tmp, sign(x$log2FC_qsva))
})

xx[[1]][,2]/xx[[1]][,1]

#### overall
summary(lm(allDxStatsExprs$tstat_qsva ~ allDxStatsExprs$negCor))$coef
summary(lm(allDxStatsExprs$tstat_qsva ~ allDxStatsExprs$posCor))$coef

## gene plots
pdf("plots/densityPlots_devReg_byDir_geneLevel.pdf",w=5,h=4)
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2)
plot(density(dxStatsExprs$Gene$tstat_qsva[dxStatsExprs$Gene$negCor==0]),
	col="grey",lwd=2,main="",xlab="SZ vs Control",xlim=c(-7,7))
lines(density(dxStatsExprs$Gene$tstat_qsva[dxStatsExprs$Gene$negCor==1]),
	col="darkblue",lwd=3)
plot(density(dxStatsExprs$Gene$tstat_qsva[dxStatsExprs$Gene$posCor==0]),
	col="grey",lwd=2,main="",xlab="SZ vs Control",xlim=c(-7,7))
lines(density(dxStatsExprs$Gene$tstat_qsva[dxStatsExprs$Gene$posCor==1]),
	col="darkorange",lwd=3)
dev.off()

### find overlaps by features
oStatsQual = lapply(dxStatsExprs, function(x) {
	o = rbind(summary(lm(x$tstat_qsva ~ x$negCor))$coef[2,c(1,4)],
		summary(lm(x$tstat_qsva ~ x$posCor))$coef[2,c(1,4)])
	rownames(o) = c("fetal", "postnatal")
	return(o)
})
oStatsQual

### find overlaps
oStatsAdj = lapply(dxStatsExprs, function(x) {
	o = rbind(summary(lm(x$tstat_adj ~ x$negCor))$coef[2,c(1,4)],
		summary(lm(x$tstat_adj ~ x$posCor))$coef[2,c(1,4)])
	rownames(o) = c("fetal", "postnatal")
	return(o)
})
oStatsAdj

### find overlaps
oStatsPca= lapply(dxStatsExprs, function(x) {
	o = rbind(summary(lm(x$tstat_pca ~ x$negCor))$coef[2,c(1,4)],
		summary(lm(x$tstat_pca ~ x$posCor))$coef[2,c(1,4)])
	rownames(o) = c("fetal", "postnatal")
	return(o)
})
oStatsPca

###### PLOTS ########

## all plots
pdf("plots/densityPlots_devReg_byDir.pdf",w=10,h=4)
par(mar=c(5,6,2,2), mfrow = c(1,2), cex.axis=2,cex.lab=2)
for(i in seq(along=dxStatsExprs)) {
	x = dxStatsExprs[[i]]
	plot(density(x$tstat_qsva[x$negCor==0]),
		col="grey",lwd=2,main=names(dxStatsExprs)[i],
		xlab="SZ vs Control",xlim=c(-7,7))
	lines(density(x$tstat_qsva[x$negCor==1]),
		col="darkblue",lwd=3)
	legend("topright", paste0("p=", signif(oStatsQual[[i]][1,2],3)),
		bty="n", cex=1.5)
		
	plot(density(x$tstat_qsva[x$posCor==0]),
		col="grey",lwd=2,main=names(dxStatsExprs)[i],
		xlab="SZ vs Control",xlim=c(-7,7))
	lines(density(x$tstat_qsva[x$posCor==1]),
		col="darkorange",lwd=3)
	legend("topright", paste0("p=", signif(oStatsQual[[i]][2,2],3)),
		bty="n", cex=1.5)
}
dev.off()
