####
library(GenomicRanges)

source("../../eqtl_functions.R")

#############################
### replication stats  ######
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/discovery_eqtls_forRepStats.rda")

## fix exprs levels
colnames(sigGene)[27] = "meanExprs"
colnames(sigExon)[28] = "meanExprs"
colnames(sigTrans)[31] = "meanExprs"
colnames(sigJxn)[41] = "meanExprs"
colnames(sigDer)[22] = "meanExprs"

## FDR 
discoveryList = list(Gene = sigGene, Exon = sigExon, Transcript = sigTrans,
	Junction = sigJxn, ER = sigDer)

repTableFdr = data.frame(fdr01 = sapply(discoveryList, nrow), 
		inCMC = sapply(discoveryList, function(x) 
			sum(!is.na(x$cmcDirAndP01))), 
		repCMC = sapply(discoveryList, function(x) 
			sum(x$cmcDirAndP01,na.rm=TRUE)), 
		inGtex = sapply(discoveryList, function(x) 
			sum(!is.na(x$gtexDir))), 
		dirGtex = sapply(discoveryList, function(x) 
			sum(x$gtexDir,na.rm=TRUE)),
		inCmcAndGtex = sapply(discoveryList, function(x) 
			sum(!is.na(x$gtexDir) & !is.na(x$cmcDirAndP01))), 
		repCmcAndGtex = sapply(discoveryList, function(x) 
			sum(x$gtexDir & x$cmcDirAndP01, na.rm=TRUE)))
repTableFdr$Expected = round(repTableFdr$fdr01*(0.5*0.01)*0.5)
repTableFdr = t(repTableFdr)
write.csv(repTableFdr, file="tables/eqtl_replication_fdr.csv")

## bonferroni			
discoveryListBonf = lapply(discoveryList, function(x) x[x$bonf < 0.05,])
repTableBonf = data.frame(fdr01 = sapply(discoveryListBonf, nrow), 
		inCMC = sapply(discoveryListBonf, function(x) 
			sum(!is.na(x$cmcDirAndP01))), 
		repCMC = sapply(discoveryListBonf, function(x) 
			sum(x$cmcDirAndP01,na.rm=TRUE)), 
		inGtex = sapply(discoveryListBonf, function(x) 
			sum(!is.na(x$gtexDir))), 
		dirGtex = sapply(discoveryListBonf, function(x) 
			sum(x$gtexDir,na.rm=TRUE)),
		inCmcAndGtex = sapply(discoveryListBonf, function(x) 
			sum(!is.na(x$gtexDir) & !is.na(x$cmcDirAndP01))), 
		repCmcAndGtex = sapply(discoveryListBonf, function(x) 
			sum(x$gtexDir & x$cmcDirAndP01, na.rm=TRUE)))
repTableBonf$Expected = round(repTableBonf$fdr01*(0.5*0.01)*0.5)
repTableBonf = t(repTableBonf)
write.csv(repTableBonf, file="tables/eqtl_replication_bonf.csv")

########		
### what determines replication??
discoveryList = lapply(discoveryList, function(x) {
	x$cmcRep = ifelse(is.na(x$cmcDirAndP01), 1, 
		ifelse(!x$cmcDirAndP01, 2, 3))
	x$gtexRep = ifelse(is.na(x$gtexDir), 1, 
		ifelse(!x$gtexDir, 2, 3))
	x
})

## presence? 
lapply(discoveryList, function(dat) {
	ii = which(dat$cmcRep %in% 1:2)
	dat2 = dat[ii,]
	dat2$cmcRep = dat2$cmcRep-1
	summary(glm(cmcRep ~ abs(statistic) + log2(meanExprs+1) +
		inSampleMAF, data = dat2, family="binomial"))$coef
})

# replication?
lapply(discoveryList, function(dat) {
	ii = which(dat$cmcRep %in% 2:3)
	dat2 = dat[ii,]
	dat2$cmcRep = dat2$cmcRep-2
	summary(glm(cmcRep ~ abs(statistic) + log2(meanExprs+1) + 
		inSampleMAF, data = dat2, family="binomial"))$coef
})

#####################
## filter for rep
repList = lapply(discoveryList, function(x) 
	x[which(x$cmcDirAndP01 & x$gtexDir),])
repListBonf = lapply(repList, function(x) x[x$bonf < 0.05,])

### posthoc race and age
sapply(repList, function(x) cor(x$AA_tstat, x$CAUC_tstat,use="comp"))
lapply(repList, function(x) table(x$AA_pval < 1e-3, 
	x$CAUC_pval <1e-3,dnn = c("AA", "CAUC")))
sapply(repListBonf, function(x) cor(x$AA_tstat, x$CAUC_tstat,use="comp"))
lapply(repListBonf, function(x) table(x$AA_pval < 1e-3, 
	x$CAUC_pval <1e-3,dnn = c("AA", "CAUC")))

### plots ####
pdf("plots/CAUC_vs_AA_posthoc_repEqtls.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2, cex.main=2)
for(i in seq(along=repList)) {
	plot(AA_tstat ~ CAUC_tstat, data = repList[[i]],
		pch = 21, bg="grey", xlab="CAUC", ylab="AA",
		main = paste0(names(repList)[i], " (T-statistics)"))
	abline(v=0,h=0,lty=2)
	abline(0,1,col="blue",lty=2,lwd=2)
}	
dev.off()
pdf("plots/CAUC_vs_AA_posthoc_repEqtlsBonf.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2, cex.main=2)
for(i in seq(along=repListBonf)) {
	plot(AA_tstat ~ CAUC_tstat, data = repListBonf[[i]],
		pch = 21, bg="grey", xlab="CAUC", ylab="AA",
		main = paste0(names(repListBonf)[i], " (T-statistics)"))
	abline(v=0,h=0,lty=2)
	abline(0,1,col="blue",lty=2,lwd=2)
}	
dev.off()

## fetal
sapply(repList, function(x) cor(x$statistic,
	x$fetal_tstat,use="comp"))
	
pdf("plots/fetal_vs_adult_posthoc_repEqtls.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2, cex.main=2)
for(i in seq(along=repList)) {
	plot(statistic ~ fetal_tstat, data = repList[[i]],
		pch = 21, bg="grey", xlab="Fetal", 
		ylab="Adult (Discovery)",
		main = paste0(names(repList)[i], " (T-statistics)"))
	abline(v=0,h=0,lty=2)
	abline(0,1,col="blue",lty=2,lwd=2)
}	
dev.off()

pdf("plots/fetal_vs_adult_posthoc_repEqtlsBonf.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2, cex.main=2)
for(i in seq(along=repListBonf)) {
	plot(statistic ~ fetal_tstat, data = repListBonf[[i]],
		pch = 21, bg="grey", xlab="Fetal", 
		ylab="Adult (Discovery)",
		main = paste0(names(repListBonf)[i], " (T-statistics)"))
	abline(v=0,h=0,lty=2)
	abline(0,1,col="blue",lty=2,lwd=2)
}	
dev.off()

## by MAF
absTstatByMaf = lapply(repListBonf, function(x) 
	x[,c("inSampleMAF", "statistic","pvalue","bonf")])
absTstatByMaf = lapply(absTstatByMaf, function(x) {
	x$mafCut = cut(x$inSampleMAF, seq(0.05,0.5, by=0.05))
	x 
})

absT = do.call("rbind", absTstatByMaf)
absT$Type = ss(rownames(absT), "\\.")
absT$pvalue[absT$pvalue < 1e-40] = 1e-40
boxplot(-log10(pvalue) ~ mafCut, data=absT)
boxplot(-log10(pvalue) ~ mafCut, data=absT,
	subset = bonf < 0.05)
summary(lm(-log10(pvalue) ~ as.numeric(mafCut) + Type, data=absT))
