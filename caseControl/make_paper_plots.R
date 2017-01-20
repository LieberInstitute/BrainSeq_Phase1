####################
library(GenomicRanges)
library(RColorBrewer)

source("../eqtl_functions.R")

##### load output ####
load("rdas/expressed_de_features.rda")

### number significant
sapply(outStatsExprs, function(x) sum(x$fdr_qsva < 0.05))
sapply(outStatsExprs, function(x) sum(x$fdr_qsva < 0.1))

######################################
### use FDR < 10% and replicated #####
sigStatsExprs = endoapply(outStatsExprs, 
	function(x) x[which(x$fdr_qsva < 0.1 & 
		sign(x$CMC_log2FC_qsva) == sign(x$log2FC_qsva) & 
		x$CMC_pval_qsva < 0.05)])
elementLengths(sigStatsExprs)


#### numbers
sigStats = unlist(sigStatsExprs)
length(sigStats)

sigStats$Type = ss(names(sigStats), "\\.")
table(sigStats$Type)
length(unique(sigStats$Symbol))
length(unique(sigStats$EnsemblGeneID))

## tab of gene by features
tabId_plusRep = table(factor(sigStats$EnsemblGeneID, 
	levels = unique(sigStats$EnsemblGeneID)), 
		factor(sigStats$Type,
		levels = c("Gene","Exon","Junction", "ER"))) > 0
colSums(tabId_plusRep)
colMeans(tabId_plusRep)
table(rowSums(tabId_plusRep))

### which genes
tabSym_plusRep = table(factor(sigStats$Symbol, 
	levels = unique(sigStats$Symbol)), 
		factor(sigStats$Type,
		levels = c("Gene","Exon","Junction", "ER"))) > 0
colnames(tabSym_plusRep) = c("G", "E", "J", "R")
tabSymOut_plusRep = data.frame(Symbol = rownames(tabSym_plusRep),
	Features = apply(tabSym_plusRep, 1, function(x) {
	paste(colnames(tabSym_plusRep)[x], collapse=",")
	}), stringsAsFactors=FALSE)
tabSymOut_plusRep$Dir = sapply(tabSymOut_plusRep$Symbol, function(x) {
	tmp = sigStats[sigStats$Symbol == x]
	unique(sign(tmp$log2FC_qsva))
})
tabSymOut_plusRep = tabSymOut_plusRep[
	rownames(tabSymOut_plusRep) != "",]
tabSymOut_plusRep$Dir = ifelse(tabSymOut_plusRep$Dir == 1, "+","-")
tabSymOut_plusRep = tabSymOut_plusRep[order(nchar(tabSymOut_plusRep$Features),decreasing=TRUE),]

write.csv(tabSymOut_plusRep, row.names=FALSE, 
	file="tables/suppTable9_consistentGenes_caseControls_allHits.csv")
	

### which genes by ID
tabID_plusRep = table(factor(sigStats$EnsemblGeneID, 
	levels = unique(sigStats$EnsemblGeneID)), 
		factor(sigStats$Type,
		levels = c("Gene","Exon","Junction", "ER"))) > 0
colnames(tabID_plusRep) = c("G", "E", "J", "R")
sum(tabID_plusRep[,1] & rowSums(tabID_plusRep[,2:4]) == 0)
table(rowSums(tabID_plusRep))

tabIDOut_plusRep = data.frame(EnsemblGeneID = rownames(tabID_plusRep),
	Features = apply(tabID_plusRep, 1, function(x) {
	paste(colnames(tabID_plusRep)[x], collapse=",")
	}), stringsAsFactors=FALSE)
	
tabIDOut_plusRep$Dir = sapply(tabIDOut_plusRep$EnsemblGeneID, function(x) {
	tmp = sigStats[which(sigStats$EnsemblGeneID == x)]
	unique(sign(tmp$log2FC_qsva))
})
tabIDOut_plusRep = tabIDOut_plusRep[
	rownames(tabIDOut_plusRep) != "",]
tabIDOut_plusRep$Dir = ifelse(tabIDOut_plusRep$Dir == 1, "+","-")
tabIDOut_plusRep = tabIDOut_plusRep[order(
	nchar(tabIDOut_plusRep$Features),decreasing=TRUE),]

table(tabIDOut_plusRep$Features)
table(nchar(tabIDOut_plusRep$Features))

###########################
#### make some plots ######

## which genes
sigStats[abs(sigStats$log2FC_qsva) > 0.4]

####################################
#### histogram of effect sizes #####
hist(sigStats$log2FC_qsva, breaks=50)
breaks = seq(-0.4, 1.1,by=0.02)

breakCounts = sapply(sigStatsExprs[-4], function(x) 
	hist(x$log2FC_qsva, breaks=breaks, plot=FALSE)$counts)

pdf("plots/histogram_of_Dx_effect_sizes.pdf", h = 4,w=5)
palette(brewer.pal(4, "Set1"))
par(mar=c(5,6,2,2))
barplot(t(breakCounts), col = 1:4,
	names.arg=breaks[-1], space=0, las=1,xaxt="n",
	cex.axis=1.6, cex.lab=1.8, xlim = c(1,40),
	xlab="Fold Change, SZ vs Control",ylab="Count")
axis(1, at = c(5,13,21,28,34,40,45.5,50.3,55),
	seq(0.8,1.6,by=0.1),cex.axis=1.4)
legend("topright", names(sigStatsExprs)[-4],
	col = 1:4, pch = 15, cex=1.3,bty="n")
dev.off()

##############################
####### scatterplot of fold changes
pdf("plots/scatterplot_of_Dx_effect_sizes.pdf")
palette(brewer.pal(4, "Set1"))
par(mar=c(5,6,4,2))
plot(sigStats$log2FC_qsva, sigStats$CMC_log2FC_qsva,
	pch = 21, cex=1.4, bg = factor(sigStats$Type, 
		levels = names(sigStatsExprs[-4])),
	cex.axis=2, cex.lab=2,  cex.main = 2,
	main = "Fold Changes (SZ vs Control)",
	xlim = c(-0.4, 0.4), ylim = c(-0.4,0.4),
	axes = FALSE, 	ylab="Replication (CMC)", xlab="Discovery (LIBD)")
axis(1,at = log2(seq(0.7,1.3, by=0.1)), seq(0.7,1.3, by=0.1),
	cex.axis = 2)
axis(2,at = log2(seq(0.7,1.3, by=0.1)), seq(0.7,1.3, by=0.1),
	cex.axis=2)
abline(0,1,lty=2,lwd=2)
abline(h=0,v=0, lwd=2)
dev.off()

####
### all FDR < 10%
sigStatsDiscovery = endoapply(outStatsExprs, 
	function(x) x[which(x$fdr_qsva < 0.1)])
	sigStats = unlist(sigStatsExprs)
sigStatsDiscovery = unlist(sigStatsDiscovery)
sigStatsDiscovery$Type = ss(names(sigStatsDiscovery), "\\.")	


tt = table(sign(sigStatsDiscovery$log2FC_qsva), 
	sign(sigStatsDiscovery$CMC_log2FC_qsva),
	dnn = c("Disc","Rep"))

pdf("plots/scatterplot_of_Dx_effect_sizes_allFDR10.pdf")
palette(brewer.pal(4, "Set1"))
par(mar=c(5,6,4,2))
plot(sigStatsDiscovery$log2FC_qsva, 
	sigStatsDiscovery$CMC_log2FC_qsva,
	pch = 21, cex=1.4, bg = factor(sigStatsDiscovery$Type, 
		levels = names(sigStatsExprs)[c(1:3,5,4)]),
	cex.axis=2, cex.lab=2,  cex.main = 2,
	main = "Fold Changes (SZ vs Control)",
	xlim = c(-0.4, 0.4), ylim = c(-0.4,0.4),
	axes = FALSE, 	ylab="Replication (CMC)", xlab="Discovery (LIBD)")
axis(1,at = log2(seq(0.7,1.3, by=0.1)), seq(0.7,1.3, by=0.1),
	cex.axis = 2)
axis(2,at = log2(seq(0.7,1.3, by=0.1)), seq(0.7,1.3, by=0.1),
	cex.axis=2)
abline(0,1,lty=2,lwd=2)
abline(h=0,v=0, lwd=2)
text(x=c(-0.3,-0.3,0.3,0.3), y = c(-0.3,0.3,0.3,-0.3),
	paste0("N=", c(tt[1,1], tt[1,2], tt[2,2],tt[2,1])),cex=1.5)
dev.off()


#####################
#### PGC Overlap ####
#####################

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
pvalAdj = summary(lm(allStatsExprs$tstat_adj ~ allStatsExprs$inPGC + 
		allStatsExprs$meanExprsAdult))$coef[2,c(1,4)]
pvalQualCmc = summary(lm(allStatsExprs$CMC_tstat_qsva ~ allStatsExprs$inPGC))$coef[2,c(1,4)]
pvalAdjCmc = summary(lm(allStatsExprs$CMC_tstat_adj ~ allStatsExprs$inPGC))$coef[2,c(1,4)]

pdf("plots/densityPlots_pgcRegionOverlaps_Quality.pdf",h=4,w=5)
par(mar=c(5,6,2,2))
plot(density(allStatsExprs$tstat_qsva[!allStatsExprs$inPGC]), 
		xlim = c(-7,7), lwd=2,cex.axis=1.8,
		cex.lab=1.8, cex.main = 1.8,  main = "Overall",
		xlab = "T-statistic (SZ vs Control)")
lines(density(allStatsExprs$tstat_qsva[allStatsExprs$inPGC]),
	col="red",lwd=2)
legend("topright", paste0("p=", signif(pvalQual[2], 3)),
	cex=1.3, bty="n")
legend("topleft", "inPGC", col="red", 
	pch = 15, bty="n",cex=1.3)
dev.off()

#############################
##### by feature ###########
############################

### find overlaps to PGC
outStatsExprs = endoapply(outStatsExprs, function(x) {
	x$inPGC = countOverlaps(x, gr108) > 0
	x
})
outStatsExprs = endoapply(outStatsExprs, function(x) {
	x$sigAndRep = x$fdr_qsva < 0.1 & 
		sign(x$CMC_log2FC_qsva) == sign(x$log2FC_qsva) & 
		x$CMC_pval_qsva < 0.05
	x
})

### metrics
oStatsQual = t(sapply(outStatsExprs, function(x) {
	summary(lm(x$tstat_qsva ~ x$inPGC + x$meanExprsAdult))$coef[2,c(1,4)]
}))
oStatsQual = data.frame(numExprs = sapply(outStatsExprs,length),
	numInPGC = sapply(outStatsExprs, function(x) sum(x$inPGC)),
	oStatsQual)
colnames(oStatsQual)[3:4] =paste0("Qual_", c("Change","Pvalue"))

oStatsAdj = t(sapply(outStatsExprs, function(x) {
	summary(lm(x$tstat_adj ~ x$inPGC + x$meanExprsAdult))$coef[2,c(1,4)]
}))
colnames(oStatsAdj) =  paste0("Adj_", c("Change","Pvalue"))
oStats = cbind(oStatsQual, oStatsAdj)
oStats
write.csv(oStats, file="tables/suppTable_pgcOverlap_stats.csv")

### density plots ###
pdf("plots/densityPlots_pgcRegionOverlaps_Quality_byFeature.pdf",h=4,w=5)
par(mar=c(5,6,2,2))
for(i in seq(along=outStatsExprs)) {
	x = outStatsExprs[[i]]
	plot(density(x$tstat_qsva[!x$inPGC]), 
		xlim = c(-7,7), lwd=2,cex.axis=1.8,cex.lab=1.8,
		cex.main = 1.8, main = names(outStatsExprs)[i],
		xlab = "T-statistic (SZ vs Control)")
	lines(density(x$tstat_qsva[x$inPGC]), col="red",lwd=2)
	legend("topright", paste0("p=", signif(oStats[i,4], 3)),
		cex=1.3, bty="n")
	legend("topleft", "inPGC", col="red", 
		pch = 15, bty="n",cex=1.3)
}
dev.off()

### table
pgcTabList = lapply(outStatsExprs, function(x) {
	table(x$sigAndRep, x$inPGC,dnn=c("Sig+Rep","inPGC"))
})
pgcTabList
lapply(outStatsExprs, function(x) x[x$sigAndRep & x$inPGC])

### check with CMC
t(sapply(outStatsExprs, function(x) {
	summary(lm(x$CMC_tstat_qsva ~ x$inPGC))$coef[2,c(1,4)]
}))
t(sapply(outStatsExprs, function(x) {
	summary(lm(x$CMC_tstat_adj ~ x$inPGC))$coef[2,c(1,4)]
}))

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
	}))
corChecks

###### which significant
pvalCuts = c(0.05, 0.01, 0.005, 0.001, 1e-4, 1e-5)
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

###### combine
library(reshape2)
concordStats = rbind(melt(concordAdj), melt(concordAndSigAdj),
	melt(concordQual), melt(concordAndSigQual))
colnames(concordStats) = c("Cutoff", "Feature", "Concordance")
concordStats$Model = rep(c("Adjust", "Qual"), each= 60)
concordStats$Type = rep(c("Dir", "Dir+Sig"), each= 30)
concordStats$numSig = c(rep(unlist(numAdj), 2), rep(unlist(numQual),2))
concordStats$ptcex = log10(concordStats$numSig)
concordStatsList = split(concordStats, concordStats$Feature)

pdf("plots/concordanceRates_byFeature.pdf", h= 6, w = 10)
for(i in seq(along=concordStatsList)) {
	x = concordStatsList[[i]]
	x$g = paste0(x$Model, ":", x$Type)
	par(mar=c(6,6,2,2))
	palette(brewer.pal(4, "Dark2"))
	plot(as.numeric(x$Cutoff), x$Concordance, 
		pch= 21, bg = as.numeric(factor(x$Model)),
		cex = x$ptcex, xaxt="n", ylim = c(0,1),
		main = names(concordStatsList)[i],
		xlab="P-value Threshold", ylab="Concordance Rate",
		cex.axis=2, cex.lab=2, las = 3, cex.main = 1.8)
	axis(1,1:6, levels(x$Cutoff),cex.axis=1.8)
	gIndexes = split(seq(along=x$g), x$g)
	for(j in seq(along=gIndexes)) {
		lines( x$Concordance ~ as.numeric(x$Cutoff), 
			col = as.numeric(factor(x$Model)), 
			lty = c(1,2,1,2)[j], lwd = 3, subset=gIndexes[[j]])
	}
	legend("topleft", c("Adj", "Qual"), col = 1:2, 
		pch = 15,cex=1.8,nc=2)
	legend("topright", c("Dir", "Dir+Sig"), nc = 2,
		col = "black", lty = 1:2,cex=1.8,lwd=4)
	abline(h=c(0.5, 0.05*0.5), lty=1:2, col="black")
}
dev.off()

#########################
### make paper plots ####