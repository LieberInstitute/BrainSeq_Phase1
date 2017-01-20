###########
library(GenomicRanges)
library(MatrixEQTL)

source("../../eqtl_functions.R")

### load data
load("rdas/pgcEqtl_indexSnps_subsetData.rda")

### load table from PGC
pgcEqtlsFinal = read.csv("tables/suppTable_allPgcEqtlsToIndexSnps.csv",
	as.is=TRUE)
length(unique(pgcEqtlsFinal$pgcFinalRank))
length(unique(pgcEqtlsFinal$pgcFinalRank[pgcEqtlsFinal$bonf < 0.05]))

length(unique(pgcEqtlsFinal$EnsemblID))
length(unique(pgcEqtlsFinal$EnsemblID[pgcEqtlsFinal$bonf < 0.05]))
length(unique(pgcEqtlsFinal$Symbol))
length(unique(pgcEqtlsFinal$Symbol[pgcEqtlsFinal$bonf < 0.05]))

### fetal effects ###
exprsList = list(Gene = geneRpkmSub, Exon = exonRpkmSub,
	Transcript = tFpkmSub, Junction = jRpkmSub,	ER = regionMatSub)
# pcList = list(Gene = newGenePCs, Exon = newExonPCs,
	# Transcript = newTxPCs, Junction = newJxnPCs, ER = newErPCs)

## fetal samples
fIndex= which(pd$Age < 0)
# modList = lapply(pcList, function(x) {
modList = lapply(fetalPcaList, function(x) {
	m = model.matrix(~snpPC1 + snpPC2 + snpPC3 + Sex,
		data=pd[fIndex,])
	cbind(m, x[,1:15])
})
fetalExprsList = lapply(exprsList, function(x) x[,fIndex])

## snp data
snpFetal = snpSub[pgcEqtlsFinal$SNP,fIndex]

## eqtls
eqtlList = split(pgcEqtlsFinal, pgcEqtlsFinal$Type)
eqtlList = eqtlList[names(modList)]

################
## get stats ###
fetalStatList = mapply(function(eqtl,y, mod) {
	yy = as.matrix(log2(y[eqtl$Feature,]+1))
	s = as.matrix(snpFetal[eqtl$SNP,])
	out = matrix(NA, nrow = nrow(s), nc = 3)
	for(i in 1:nrow(yy)) {
		f = lm(yy[i,] ~ s[i,] + mod - 1)
		out[i,] = summary(f)$coef[1,-2]
	}
	colnames(out) = paste0("fetal_", c("beta","tstat","pvalue"))
	cbind(eqtl, out)
}, eqtlList, fetalExprsList, modList,SIMPLIFY=FALSE)
pgcEqtlsFetal = do.call("rbind",fetalStatList) 
pgcEqtlsFetal = pgcEqtlsFetal[order(pgcEqtlsFetal$pvalue),]
rownames(pgcEqtlsFetal) = NULL

#############
### plot ####

## FDR < 1%
plot(beta ~ fetal_beta, data=pgcEqtlsFetal)
table(pgcEqtlsFetal$fetal_pvalue < 0.05)
mean(pgcEqtlsFetal$fetal_pvalue < 0.05)
table(sign(pgcEqtlsFetal$fetal_beta) == sign(pgcEqtlsFetal$beta))
prop.test(table(sign(pgcEqtlsFetal$fetal_beta) == sign(pgcEqtlsFetal$beta)))$p.value
mean(sign(pgcEqtlsFetal$fetal_beta) == sign(pgcEqtlsFetal$beta))

mean(pgcEqtlsFetal$fetal_pvalue < 0.05 & 
	sign(pgcEqtlsFetal$fetal_beta) == sign(pgcEqtlsFetal$beta))
table(pgcEqtlsFetal$fetal_pvalue < 0.05,  
	sign(pgcEqtlsFetal$fetal_beta) == sign(pgcEqtlsFetal$beta),
	dnn = c("P<0.05", "Dir"))

###################
## load developmental stats
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/isoform_switch_devel_byFeature.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/devStats_controlSamples.rda")

statListEqtl = endoapply(statList, function(x) 
	x[names(x) %in% pgcEqtlsFetal$Feature])
devStatsEqtl = unlist(statListEqtl)
names(devStatsEqtl) = ss(names(devStatsEqtl), "\\.", 2)
devStatsEqtl = devStatsEqtl[pgcEqtlsFetal$Feature]
pgcEqtlsFetal$ageCorr = devStatsEqtl$ageCorr
pgcEqtlsFetal$devPbonf = devStatsEqtl$p_bonf

# check devel
table(pgcEqtlsFetal$devPbonf < 0.05)
table(sign(pgcEqtlsFetal$ageCorr[which(pgcEqtlsFetal$devPbonf < 0.05)]))
table(sign(pgcEqtlsFetal$ageCorr), pgcEqtlsFetal$fetal_pvalue < 0.05)


## switches?
switchDat = do.call("rbind", switchList)
pgcEqtlsFetal$fetalIsoSwitch = pgcEqtlsFetal$EnsemblID %in% 
	ss(rownames(switchDat), "\\.",2)
table(pgcEqtlsFetal$fetalIsoSwitch)
table(pgcEqtlsFetal$fetalIsoSwitch[!duplicated(pgcEqtlsFetal$EnsemblID)])

## all expressed genes
bgGene = unique(unlist(lapply(statList, function(x) x$EnsemblGeneID[!is.na(x$p_bonf)])))
bgGene = bgGene[!is.na(bgGene) & !grepl("-", bgGene)]
switchGene = unique(unlist(lapply(switchList, rownames)))

mat = matrix(c(57,77,9148,33764), nr= 2, byrow=2)
mat[1,1]*mat[2,2]/mat[1,2]/mat[2,1]

write.csv(pgcEqtlsFetal, row.names=FALSE,
	"tables/suppTable_allPgcEqtlsToIndexSnps_withFetal.csv")

	##############
## BONF < 5% #
pgcEqtlsFetal = read.csv("tables/suppTable_allPgcEqtlsToIndexSnps_withFetal.csv",
	as.is=TRUE)
pgcEqtlsFetalBonf = pgcEqtlsFetal[pgcEqtlsFetal$bonf < 0.05,]
plot(beta ~ fetal_beta, data=pgcEqtlsFetalBonf)
mean(pgcEqtlsFetalBonf$fetal_pvalue < 0.05)
mean(sign(pgcEqtlsFetalBonf$fetal_beta) == sign(pgcEqtlsFetalBonf$beta))

mean(pgcEqtlsFetalBonf$fetal_pvalue < 0.05 & 
	sign(pgcEqtlsFetalBonf$fetal_beta) == sign(pgcEqtlsFetalBonf$beta))
table(pgcEqtlsFetalBonf$fetal_pvalue < 0.05,  
	sign(pgcEqtlsFetalBonf$fetal_beta) == sign(pgcEqtlsFetalBonf$beta),
	dnn = c("P<0.05", "Dir"))

table(pgcEqtlsFetalBonf$devPbonf < 0.05)
table(sign(pgcEqtlsFetalBonf$ageCorr[which(pgcEqtlsFetalBonf$devPbonf < 0.05)]))
table(sign(pgcEqtlsFetalBonf$ageCorr), pgcEqtlsFetalBonf$fetal_pvalue < 0.05)

table(pgcEqtlsFetalBonf$fetalIsoSwitch)
table(pgcEqtlsFetalBonf$fetalIsoSwitch[!duplicated(pgcEqtlsFetalBonf$EnsemblID)])



###################	
## check vs pgc ###
pgc = readxl::read_excel("nature13595-s3.xlsx",sheet=2)[,1:12]
colnames(pgc)=gsub(" ", "_", colnames(pgc))
# pgc = pgc[pgc$Index_SNPs %in% pgcEqtlsFetal$snpRsNum ,]
table(pgc$eQTL_gene %in% pgcEqtlsFetal$Symbol)
pgc[pgc$eQTL_gene %in% pgcEqtlsFetal$Symbol,]

table(pgc$eQTL_gene[pgc[,12] > 0.6] %in% pgcEqtlsFetal$Symbol)
