####
library(GenomicRanges)

source("../../eqtl_functions.R")

## data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/cleaned_eqtl_data_subset_n412.rda")
snpMapSub$chrpos= paste0("chr", snpMapSub$CHR, ":", snpMapSub$POS)

#################
#### metrics ####
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/allEqtlsRep.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/sigEqtlsRep.rda")
sigEqtlsRep$Type = factor(as.character(sigEqtlsRep$Type),
	levels = c("Gene", "Exon", "Transcript", "Junction", "ER"))
sigEqtlsRep$snpChrPos = snpMapSub$chrpos[match(sigEqtlsRep$SNP, snpMapSub$SNP)]
sigEqtlsRep$countedAllele = snpMapSub$COUNTED[match(sigEqtlsRep$SNP, snpMapSub$SNP)]
sigEqtlsRep$refAllele = snpMapSub$ALT[match(sigEqtlsRep$SNP, snpMapSub$SNP)]

length(unique(allEqtlsRep$SNP))
length(unique(allEqtlsRep$SNP[allEqtlsRep$bonf < 0.05]))

### LD independence ###
pruned = read.table("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/plink/LIBD_Brain_DLPFC_szControls_imputed_indep.prune.in")
allEqtlsRep$ldIndep = allEqtlsRep$SNP %in% pruned$V1

length(unique(allEqtlsRep$SNP[allEqtlsRep$ldIndep]))
length(unique(allEqtlsRep$SNP[
	allEqtlsRep$bonf < 0.05 & allEqtlsRep$ldIndep]))

######################
# make list
eqtlList = split(sigEqtlsRep, sigEqtlsRep$Type)
	
### FDR filters ####
numFdr = sapply(eqtlList, function(x) sum(x$FDR < 0.01))	
numFdr

# bonf
numBonf = sapply(eqtlList, function(x) sum(x$bonf < 0.05))	
numBonf

## significant ones at bonf
sigList = lapply(eqtlList, 
	function(x) x[which(x$bonf < 0.05),])

# p-value cutoff for BONF < 5%
pCutFdr = sapply(eqtlList, function(x) max(x$pvalue))
pCutFdr
pCutBonf = sapply(sigList, function(x) max(x$pvalue))
pCutBonf

# number of unique genes
nGeneFdrID = sapply(eqtlList, function(x) length(unique(x$EnsemblID)))
nGeneFdrSym = sapply(eqtlList, function(x) length(unique(x$Symbol)))
nGeneBonfID = sapply(sigList, function(x) length(unique(x$EnsemblID)))
nGeneBonfSym = sapply(sigList, function(x) length(unique(x$Symbol)))

## effect size
effectFdr = t(sapply(eqtlList, function(x) 2^quantile(abs(x$beta))[c(3,2,4)]))
effectFdr = apply(signif(effectFdr,3), 1, function(x) 
	paste0(x[1], " (",	x[2],"-",x[3],")"))
effectBonf = t(sapply(sigList, function(x) 2^quantile(abs(x$beta))[c(3,2,4)]))
effectBonf = apply(signif(effectBonf,3), 1, function(x) 
	paste0(x[1], " (",	x[2],"-",x[3],")"))


pdf("plots/effectSize_boxplots.pdf", h=5.5,w=11)
par(mar=c(3,6,2,2), cex.axis=2,cex.lab=2, cex.main=2)
boxplot(abs(sigEqtlsRep$beta) ~ sigEqtlsRep$Type, yaxt = "n",
	ylim = c(0,2),ylab="Effect Size")	
axis(2, at = 0:4, 2^(0:4))
dev.off()


summary(lm(abs(sigEqtlsRep$beta) ~ sigEqtlsRep$Type))

# dist
distFdr = t(sapply(eqtlList, function(x) quantile(x$snpDistToFeature)[c(3,2,4)]))
distFdr = apply(distFdr, 1, function(x) 
	paste0(x[1], "bp (",	signif(x[2]/1000,2),"kb-",signif(x[3]/1000,2),"kb)"))
distBonf = t(sapply(sigList, function(x) quantile(x$snpDistToFeature)[c(3,2,4)]))
distBonf = apply(distBonf, 1, function(x) 
	paste0(x[1], "bp (",	signif(x[2]/1000,2),"kb-",signif(x[3]/1000,2),"kb)"))
hist(sigList$Gene$snpDistToFeature,breaks=100)

######################
### concordance ######

# snpByTypeFdr = table(allEqtlsRep$SNP, allEqtlsRep$Type)
allEqtlsRep$snpByGene = paste0(allEqtlsRep$SNP, ".", allEqtlsRep$EnsemblID)
snpByTypeBonf = table(allEqtlsRep$snpByGene[allEqtlsRep$bonf < 0.05],
	as.character(allEqtlsRep$Type[allEqtlsRep$bonf < 0.05]))

# how much gene specificity?
snpSpecTable = table(table(ss(rownames(snpByTypeBonf), "\\.")))

# by feature
snpSpecTableByFeature = apply(snpByTypeBonf, 2, function(x) {
	cat(".")
	ii = which(x > 0)
	table(table(ss(rownames(snpByTypeBonf)[ii], "\\.")))
})
sapply(snpSpecTableByFeature, function(x) 100*prop.table(x)[1])

### ld indep SNPs
snpByTypeBonfLd = table(allEqtlsRep$snpByGene[
		allEqtlsRep$bonf < 0.05 & allEqtlsRep$ldIndep],
	as.character(allEqtlsRep$Type[
		allEqtlsRep$bonf < 0.05& allEqtlsRep$ldIndep]))
snpSpecTableLd = table(table(ss(rownames(snpByTypeBonfLd), "\\.",2)))
100*prop.table(snpSpecTableLd)[1:5]


########################
### tx specifity ###
numTxFdr = sapply(eqtlList, function(x) 
	c(nrow(x), sum(x$NumTxGene > 1), sum(x$NumTxEqtl <= 1 & x$NumTxGene > 1)))
numTxBonf = sapply(sigList, function(x) 
	c(nrow(x), sum(x$NumTxGene > 1), sum(x$NumTxEqtl <= 1 & x$NumTxGene > 1)))
numTxFdr[2,]/numTxFdr[1,]
numTxBonf[2,]/numTxBonf[1,]
txUniqueFdr = numTxFdr[3,]/numTxFdr[2,]
txUniqueBonf = numTxBonf[3,]/numTxBonf[2,]
	
numTxFdr2 = sapply(eqtlList, function(x) 
	c(nrow(x), sum(x$NumTxGene > 2), sum(x$NumTxEqtl <= 2 & x$NumTxGene > 2)))
numTxBonf2 = sapply(sigList, function(x) 
	c(nrow(x), sum(x$NumTxGene > 2), sum(x$NumTxEqtl <= 2 & x$NumTxGene > 2)))

	
########novel
novelFdr = sapply(eqtlList, function(x) mean(x$Class != "InEns"))
sapply(eqtlList, function(x) table(x$Class))
sapply(eqtlList, function(x) sum(x$Class != "InEns"))

novelBonf = sapply(sigList, function(x) mean(x$Class != "InEns"))
novelBonf
sapply(sigList, function(x) sum(x$Class != "InEns"))
sapply(sigList, function(x) table(x$Class))

####################################
##### make table to write out ####
outTableBonf = data.frame(numBonf, pCutBonf, 
	nGeneBonfID, nGeneBonfSym, effectBonf, 
	distBonf, txUniqueBonf, novelBonf)
write.csv(outTableBonf, file="tables/bonf_eqtl_figure3a.csv")

outTableFdr = data.frame(numFdr, pCutFdr, 
	nGeneFdrID, nGeneFdrSym, effectFdr, 
	distFdr, txUniqueFdr, novelFdr)
write.csv(outTableFdr, file="tables/fdr_eqtl_suppTable.csv")


#################### 
## all genes ####
allGenesSym = unique(unlist(lapply(sigList, "[[", "Symbol")))
allGenesSym = allGenesSym[!is.na(allGenesSym) & !grepl("-", allGenesSym)]
length(allGenesSym)

allGenesID = unique(unlist(lapply(sigList, "[[", "EnsemblID")))
allGenesID = allGenesID[!is.na(allGenesID) & !grepl("-", allGenesID)]
length(allGenesID)

allGenesFdrID = unique(unlist(lapply(eqtlList, "[[", "EnsemblID")))
allGenesFdrID = allGenesFdrID[!is.na(allGenesFdrID) & 
	!grepl("-", allGenesFdrID)]
length(allGenesFdrID)
allGenesFdrSym = unique(unlist(lapply(eqtlList, "[[", "Symbol")))
allGenesFdrSym = allGenesFdrSym[!is.na(allGenesFdrSym) & 
	!grepl("-", allGenesFdrSym)]
length(allGenesFdrSym)

convergeMatSym = sapply(sigList, function(x) allGenesSym %in% x$Symbol)
colSums(convergeMatSym)
table(rowSums(convergeMatSym))
convergeMatID = sapply(sigList, function(x) allGenesID %in% x$EnsemblID)
colSums(convergeMatID)
table(rowSums(convergeMatID))

convergeMatSymFdr = sapply(eqtlList, function(x) allGenesFdrSym %in% x$Symbol)
colSums(convergeMatSymFdr)
table(rowSums(convergeMatSymFdr))
convergeMatIdFdr = sapply(eqtlList, function(x) allGenesFdrID %in% x$EnsemblID)
colSums(convergeMatIdFdr)
table(rowSums(convergeMatIdFdr))

library(gplots)
pdf("plots/vennDiagram_annotated.pdf")
par(mar=c(1,1,1,1))
venn(as.data.frame(convergeMatID))
venn(as.data.frame(convergeMatID[,-3])) # no tx
venn(as.data.frame(convergeMatID[,-5])) # no der
venn(as.data.frame(convergeMatID[,c(1,2,4)])) # j,e,g
dev.off()

###############################
#### effect size histograms ###

pdf("plots/histogram_distanceToFeature_bestAndAll.pdf",h=6,w=8)
par(mar=c(6,6,2,2), cex.axis=2,cex.lab=2, cex.main=2)
hist(sigEqtlsRep$snpDistToFeature,breaks=100, main="",
	col = "grey", xlab="",xaxt="n",yaxt="n")
axis(1, at = seq(-5e5,5e5,by=1e5), 
	seq(-500,500,by=100), las=3)
axis(2, at = seq(0,25000, by=5000), seq(0,25,5))
hist(allEqtlsRep$snpDistToFeature,breaks=100, main="",
	col = "grey", xlab="",xaxt="n")
axis(1, at = seq(-5e5,5e5,by=1e5), 
	seq(-500,500,by=100), las=3)
dev.off()

pdf("plots/CAUC_vs_AA_posthoc_repEqtls.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2, cex.main=2)
for(i in seq(along=eqtlList)) {
	hist(
	plot(AA_tstat ~ CAUC_tstat, data = repList[[i]],
		pch = 21, bg="grey", xlab="CAUC", ylab="AA",
		main = paste0(names(repList)[i], " (T-statistics)"))
	abline(v=0,h=0,lty=2)
	abline(0,1,col="blue",lty=2,lwd=2)
}	
dev.off()

### write out tables
outEqtl = as.data.frame(sigEqtlsRep)
outEqtl$WhichTx = sapply(outEqtl$WhichTx, paste, collapse=";")
colnames(outEqtl)[c(1,17)] = c("BestSNP","snpLoc")
outEqtl = outEqtl[,c(1:12,14:19,13)]
rownames(outEqtl) = NULL
outEqtl$eqtlIndex = paste0("eQTL_", 1:nrow(outEqtl))
write.csv(outEqtl, row.names=FALSE,
	file=gzfile("tables/LIBD_eQTLs_sigAndRep_FDR_allFeatures.csv.gz"))
	
### make all boxplots
cleanExprs = rbind(cleanGeneSub, cleanExonSub, cleanJxnSub,
	cleanTxSub, cleanErSub)

theSnps = snpSub[outEqtl$BestSNP,]
theSnpMap = snpMapSub[match(outEqtl$BestSNP, snpMapSub$SNP),]

plotFns = paste0("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/boxplots/",
	outEqtl$eqtlIndex, ".png")
mainTxt = paste0(outEqtl$EnsemblID, " (", outEqtl$Symbol, 
	") - ",	outEqtl$Class)
ll = ifelse(outEqtl$beta> 0, "topleft", "topright")
	
for(i in seq(along=plotFns)) {
	if(i %% 500 == 0) cat(".")
	png(plotFns[i], width=600,height=600, type="cairo-png")
	par(mar=c(5,6,4,2), cex.axis=2, cex.lab=2, cex.main=1.6)
	boxplot(cleanExprs[i,] ~ theSnps[i,], outline=FALSE,
		ylim= range(cleanExprs[i,]), 
		ylab=paste(outEqtl$Feature[i], "(log2, Adjusted)"),
		xlab=paste0(outEqtl$snpRsNum[i], ": # of ",
			outEqtl$countedAllele[i]), main = mainTxt[i])
	points(cleanExprs[i,] ~ jitter(theSnps[i,]+1,amount=0.1),
		pch = 21, bg="grey",cex=1.4)
	legend(ll[i], paste0("p=",signif(outEqtl$pvalue[i], 3)),cex=1.8)
	dev.off()
}
