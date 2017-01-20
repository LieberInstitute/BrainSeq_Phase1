###
library(GenomicRanges)
source("../eqtl_functions.R")
library(GenomicFeatures)
library(org.Hs.eg.db)
library(biomaRt)
library(clusterProfiler)
library(RColorBrewer)

### load output
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/caseControl/rdas/expressed_de_features.rda")

###### which significant
pvalCuts = c(0.05, 0.01, 0.005, 0.001, 1e-4, 1e-5,1e-6)
names(pvalCuts) = paste0("p<", pvalCuts)

#### adjustment models
concordAdj = sapply(outStatsExprs, function(y) {
	sapply(pvalCuts, function(x) {
		sum(y$pval_adj < x & 
			sign(y$log2FC_adj) == 
				sign(y$CMC_log2FC_adj),
			na.rm=TRUE) / sum(y$pval_adj < x)
	})
})

concordAndSigAdj = sapply(outStatsExprs, function(y) {
	sapply(pvalCuts, function(x) {
		sum(y$pval_adj < x & 
			sign(y$log2FC_adj) == 
				sign(y$CMC_log2FC_adj) &
			y$CMC_pval_adj < 0.05,
		na.rm=TRUE) / sum(y$pval_adj < x)
	})
})
numAdj = sapply(outStatsExprs, function(y) {
	sapply(pvalCuts, function(x) sum(y$pval_adj < x))})

## qSVA models
concordQual = sapply(outStatsExprs, function(y) {
	sapply(pvalCuts, function(x) {
		sum(y$pval_qsva < x & 
			sign(y$log2FC_qsva) == 
				sign(y$CMC_log2FC_qsva),
			na.rm=TRUE) / sum(y$pval_qsva < x)
	})
})
concordAndSigQual = sapply(outStatsExprs, function(y) {
	sapply(pvalCuts, function(x) {
		sum(y$pval_qsva < x & 
			sign(y$log2FC_qsva) == 
				sign(y$CMC_log2FC_qsva) &
			y$CMC_pval_qsva < 0.05,
		na.rm=TRUE) / sum(y$pval_qsva < x)
	})
})
numQual = sapply(outStatsExprs, function(y) {
	sapply(pvalCuts, function(x) sum(y$pval_qsva < x))})

## PCA models
concordPca = sapply(outStatsExprs, function(y) {
	sapply(pvalCuts, function(x) {
		sum(y$pval_pca < x & 
			sign(y$log2FC_pca) == 
				sign(y$CMC_log2FC_pca),
			na.rm=TRUE) / sum(y$pval_pca < x)
	})
})
concordAndSigPca = sapply(outStatsExprs, function(y) {
	sapply(pvalCuts, function(x) {
		sum(y$pval_pca < x & 
			sign(y$log2FC_pca) == 
				sign(y$CMC_log2FC_pca) &
			y$CMC_pval_pca < 0.05,
		na.rm=TRUE) / sum(y$pval_pca < x)
	})
})
numPca = sapply(outStatsExprs, function(y) {
	sapply(pvalCuts, function(x) sum(y$pval_pca < x))})

###### combine
library(reshape2)
concordStats = rbind(melt(concordAdj), melt(concordAndSigAdj),
	melt(concordQual), melt(concordAndSigQual),
	melt(concordPca), melt(concordAndSigPca))
colnames(concordStats) = c("Cutoff", "Feature", "Concordance")

concordStats$Model = rep(c("Adjust", "Qual","PCA"), each= 70)
concordStats$Model = factor(concordStats$Model, 
	c("Adjust", "Qual","PCA"))
concordStats$Type = rep(c("Dir", "Dir+Sig"), each= 35)
concordStats$numSig = c(rep(unlist(numAdj), 2), 
	rep(unlist(numQual),2),rep(unlist(numPca),2))
concordStats$ptcex = log10(concordStats$numSig)

concordStatsList = split(concordStats, concordStats$Feature)

### replication
pdf("plots/replicationRates_byFeature_withPca.pdf",w=6)
for(i in seq(along=concordStatsList)) {
	x = concordStatsList[[i]]
	x = x[x$Type == "Dir+Sig",]
	x$g = paste0(x$Model, ":", x$Type)
	par(mar=c(8,6,2,2))
	palette(brewer.pal(4, "Dark2"))
	plot(as.numeric(x$Cutoff), x$Concordance, 
		pch= 21, bg = x$Model,
		cex = x$ptcex, xaxt="n", ylim = c(0,0.7),
		main = names(concordStatsList)[i],
		xlab="", ylab="Replication Rate",
		cex.axis=2, cex.lab=2, las = 3, cex.main = 1.8)
	axis(1,1:7, levels(x$Cutoff),cex.axis=1.8,las=3)
	gIndexes = split(seq(along=x$g), x$g)
	for(j in seq(along=gIndexes)) {
		lines( x$Concordance ~ as.numeric(x$Cutoff), 
			col = as.numeric(factor(x$Model)), 
			lty = 1, lwd = 3, subset=gIndexes[[j]])
	}
	legend("topleft", c("Adj", "Qual","PCA"), col = 1:3, 
		pch = 15,cex=1.8,nc=3)
	abline(h=0.05*0.5, lty=2, col="black",lwd=2)
}
dev.off()

## directional consistency alone
pdf("plots/directionalRates_byFeature_withPca.pdf",w=6)
for(i in seq(along=concordStatsList)) {
	x = concordStatsList[[i]]
	x = x[x$Type == "Dir",]
	x$g = paste0(x$Model, ":", x$Type)
	par(mar=c(8,6,2,2))
	palette(brewer.pal(4, "Dark2"))
	plot(as.numeric(x$Cutoff), x$Concordance, 
		pch= 21, bg = x$Model,
		cex = x$ptcex, xaxt="n", ylim = c(0.4,1),
		main = names(concordStatsList)[i],
		xlab="", ylab="Consistency Rate",
		cex.axis=2, cex.lab=2, las = 3, cex.main = 1.8)
	axis(1,1:7, levels(x$Cutoff),cex.axis=1.8,las=3)
	gIndexes = split(seq(along=x$g), x$g)
	for(j in seq(along=gIndexes)) {
		lines( x$Concordance ~ as.numeric(x$Cutoff), 
			col = as.numeric(factor(x$Model)), 
			lty = 1, lwd = 3, subset=gIndexes[[j]])
	}
	legend("topleft", c("Adj", "Qual","PCA"), col = 1:3, 
		pch = 15,cex=1.8,nc=3)
	abline(h=0.5, lty=2, col="black",lwd=2)
}
dev.off()

###############################	
###### check correlations #####
corChecks = data.frame(adjDir = sapply(outStatsExprs[-4], function(x) {
	cor(x$log2FC_adj, x$CMC_log2FC_adj,use="comp")
	}), adjDirAndSig = sapply(outStatsExprs[-4], function(x) {
	cor(x$log2FC_adj[x$fdr_adj < 0.05], 
		x$CMC_log2FC_adj[x$fdr_adj < 0.05],	use="comp")
	}), qualDir = sapply(outStatsExprs[-4], function(x) {
	cor(x$log2FC_qsva, x$CMC_log2FC_qsva,use="comp")
	}), qualDirAndSig = sapply(outStatsExprs[-4], function(x) {
	cor(x$log2FC_qsva[x$fdr_qsva < 0.1], 
		x$CMC_log2FC_qsva[x$fdr_qsva < 0.1 ],use="comp")
	}), pcaDir = sapply(outStatsExprs[-4], function(x) {
	cor(x$log2FC_pca, x$CMC_log2FC_pca,use="comp")
	}), pcaDirAndSig = sapply(outStatsExprs[-4], function(x) {
		cor(x$log2FC_pca[x$fdr_pca < 0.1], 
		x$CMC_log2FC_pca[x$fdr_pca < 0.1 ],use="comp")
	}))
	
corChecks

## read in GWAS regions
pgcFinal = read.delim("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/tables/pgc2_loci.txt", as.is=TRUE)
pgcFinal$chr = ss(pgcFinal$Position..hg19.,":")
tmp = ss(pgcFinal$Position..hg19.,":",2)
pgcFinal$start = as.numeric(ss(tmp, "-"))
pgcFinal$end = as.numeric(ss(tmp, "-",2))
gr108 = GRanges(pgcFinal$chr, IRanges(pgcFinal$start,pgcFinal$end))

### overall
allStatsExprs = unlist(outStatsExprs)
allStatsExprs$Type = ss(names(allStatsExprs), "\\.")
allStatsExprs$inPGC = countOverlaps(allStatsExprs, gr108) > 0
table(allStatsExprs$inPGC)
table(allStatsExprs$inPGC, allStatsExprs$Type)

pvalQual = summary(lm(allStatsExprs$tstat_qsva ~ allStatsExprs$inPGC + 
		allStatsExprs$meanExprsAdult))$coef[2,c(1,4)]
pvalPca = summary(lm(allStatsExprs$tstat_pca ~ allStatsExprs$inPGC + 
		allStatsExprs$meanExprsAdult))$coef[2,c(1,4)]
pvalAdj = summary(lm(allStatsExprs$tstat_adj ~ allStatsExprs$inPGC + 
		allStatsExprs$meanExprsAdult))$coef[2,c(1,4)]


#######################		
## gene set test ######
sigStats = endoapply(outStatsExprs, 
	function(x) x[which(x$fdr_pca < 0.1 & 
		sign(x$CMC_log2FC_pca) == sign(x$log2FC_pca) & 
		x$CMC_pval_pca < 0.05)])
	
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

### compare
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