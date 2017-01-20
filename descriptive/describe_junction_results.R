#######
# R-devel
source("../eqtl_functions.R")

library(GenomicRanges)
library(limma)
library(parallel)

load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
pd$ageGroup = cut(pd$Age, c(-1,0,1,10,20,50,100),
	labels = c("Fetal","Infant","Child","Teen","Adult","50+"))

quantile(pd$totalMapped/1e6)

length(jMap)

## junctions
exprsJxnIndex = which(rowSums(jRpkm > 0) > 1) # 1987089 are unique to a single person
length(exprsJxnIndex)
jMap$geneIndex = match(jMap$ensemblGeneID, rownames(geneMap))
jMap$geneIndexAny = match(jMap$newGeneID, rownames(geneMap))

## make factor
jMap$code = factor(jMap$code,
	levels=c("InEns", "ExonSkip","AltStartEnd", "Novel"))
jMap$meanExprs = rowMeans(jRpkm)	

sum(jMap$meanExprs < 1)
mean(jMap$meanExprs < 1)

##############
# junction stats ####
### which comparisons?

juncStats = junctionStats(jRpkm, jMap)
##################
### other datasets
jMapLibd = jMap

### GTEX
load("/users/ajaffe/Lieber/Projects/RNAseq/GTEX/rdas/gtex_brain_jMap.rda")
mmGtex = match(names(jMapLibd), names(jMapGtexBrain))
jMapLibd$inGTEX = !is.na(mmGtex)

#### GEUVADIS
load("/dcs01/ajaffe/GEUVADIS/Counts/rpkmCounts_GEUVADIS.rda")
mmGeuv = match(names(jMapLibd), names(jMap))
jMapLibd$inGeuv = !is.na(mmGeuv)

### metrics for the paper
table(jMapLibd$code)
tab = with(jMapLibd, table(inGTEX, inGeuv, code))
round(prop.table(tab, 3),3)*100
tab

table(jMapLibd$code[jMapLibd$meanExprs > 1])
tab1 = with(jMapLibd[jMapLibd$meanExprs > 1], table(inGTEX, inGeuv, code))
round(prop.table(tab1, 3),3)*100
tab1

table(jMapLibd$code[jMapLibd$meanExprs > 5])
tab5 = with(jMapLibd[jMapLibd$meanExprs > 5], table(inGTEX, inGeuv, code))
round(prop.table(tab5, 3),3)*100
round(prop.table(tab5, 3),4)*100
tab5

tabFusion = with(jMapLibd[jMapLibd$meanExprs > 1], 
	table(inGTEX, isFusion))
round(prop.table(tabFusion, 2),3)*100
tab1

#### PGC regions?? #####

## read in GWAS
pgcFinal = read.delim("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/tables/pgc2_loci.txt", as.is=TRUE)
pgcFinal$chr = ss(pgcFinal$Position..hg19.,":")
tmp = ss(pgcFinal$Position..hg19.,":",2)
pgcFinal$start = as.numeric(ss(tmp, "-"))
pgcFinal$end = as.numeric(ss(tmp, "-",2))
gr108 = GRanges(pgcFinal$chr, IRanges(pgcFinal$start,pgcFinal$end))

# overlap junctions
oo = findOverlaps(jMap, gr108)
jMap$inPGC2 = FALSE
jMap$inPGC2[queryHits(oo)] = TRUE

## are novel junctions more likely to be represented?
thecuts = c(0,0.1,0.5,1,2,5)
fitList = lapply(thecuts, function(x) {
	cat(".")
	summary(glm(factor(jMap$inPGC2) ~ jMap$code +
		jMap$meanExprs + jMap$numTx, family=binomial,
			subset=jMap$meanExprs > x))
})
statList = lapply(fitList, function(x) cbind(exp(x$coef[,1]), x$coef[,4]))
statOverlap = do.call("cbind", statList)
colnames(statOverlap) = paste0(rep(c("OR", "pval"), 
	length(thecuts)),"_rp80M>", 	rep(thecuts, each=2))
rownames(statOverlap) = gsub("jMap$", "", rownames(statOverlap), fixed=TRUE)
write.csv(statOverlap, file="tables/PGC_enrichment_of_novel_isoforms_exprsCutoffs.csv")

# numbers for paper
range(statOverlap["codeExonSkip",grep("OR", colnames(statOverlap))])
range(statOverlap["codeExonSkip",grep("pval", colnames(statOverlap))])
range(statOverlap["codeAltStartEnd",grep("OR", colnames(statOverlap))])
range(statOverlap["codeAltStartEnd",grep("pval", colnames(statOverlap))])

#### other disorders
xx=load("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/gwasResults_lifted.rda")

## get genes in each region, all expressed first
gwasList = split(gwasLift, gwasLift$Dx)
ooOther_all = findOverlaps(gwasList, jMap)
gOther_all = split(subjectHits(ooOther_all), 
	names(gwasList)[queryHits(ooOther_all)])
names(gOther_all) = names(gwasList)	

otherMat= matrix(FALSE, nc = 3, nrow = length(jMap),
	dimnames = list(names(jMap), names(gwasList)))
for(i in seq(along=gOther_all)) otherMat[gOther_all[[i]], i] = TRUE

fitListOther = lapply(as.data.frame(otherMat), function(z) {
	tmp = lapply(thecuts, function(x) {
		cat(".")
		summary(glm(factor(z) ~ jMap$code +
			jMap$meanExprs + jMap$numTx, family=binomial,
				subset=jMap$meanExprs > x))
		})
	statList = lapply(tmp, function(x) cbind(exp(x$coef[,1]), x$coef[,4]))
	statOverlap = do.call("cbind", statList)
	colnames(statOverlap) = paste0(rep(c("OR", "pval"), length(thecuts)),
		"_rp80M>", 	rep(thecuts, each=2))
	statOverlap
})
lapply(fitListOther, function(x) x[2:3,])

save(fitList, fitListOther, file="rdas/gwas_enrichment_newJxns.rda")

## second analysis within moderate expression
fit01 = glm(factor(jMap$inPGC2 ) ~ jMap$code +
	jMap$meanExprs + jMap$numTx, family=binomial,
		subset = jMap$meanExprs >0.1)
summary(fit01)$coef

fit1 = glm(factor(jMap$inPGC2 ) ~ jMap$code +
	jMap$meanExprs + jMap$numTx, family=binomial,
		subset = jMap$meanExprs > 1)
summary(fit1)$coef

## MA plots
M = as.numeric(log2(jMapLibd$meanExprs[!is.na(mmGtex)]+1) - 
	log2(jMapGtex$meanExprs[mmGtex[!is.na(mmGtex)]]+1))
A = as.numeric((log2(jMapLibd$meanExprs[!is.na(mmGtex)]+1) +
	log2(jMapGtex$meanExprs[mmGtex[!is.na(mmGtex)]]+1))/2)
table(M < 1, A > 1, dnn = c("Within2", "Exprs"))

## by type
Indexes= split(1:sum(!is.na(mmGtex)), jMapLibd$code[!is.na(mmGtex)])
for(i in seq(along=Indexes)) {
	cat(".")
	png(paste0("plots/MA_plots_", names(Indexes)[i],".png"),
		height = 500, width=500)
	par(mar=c(5,6,3,2))
	palette(brewer.pal(8,"Dark2"))
	ii = Indexes[[i]]
	plot(M ~ A, subset=ii, pch=21,bg=i,cex=0.7,
		main = names(Indexes)[i], 
		xlab="Mean Expression (log2)",
		ylab="LIBD - GTEX log2FC",  cex.axis=2, 
		cex.lab=2,cex.main=2)
	dev.off()
}
lapply(Indexes, function(ii) table(M[ii] < 1,
	A[ii] > 1, dnn = c("Within2", "Exprs")))
	
## DLPFC only
Md = as.numeric(log2(jMapLibd$meanExprs[!is.na(mmGtex)]+1) - 
	log2(jMapGtex$meanExprsDlpfc[mmGtex[!is.na(mmGtex)]]+1))
Ad = as.numeric((log2(jMapLibd$meanExprs[!is.na(mmGtex)]+1) +
	log2(jMapGtex$meanExprsDlpfc[mmGtex[!is.na(mmGtex)]]+1))/2)
table(Md < 1, Ad > 1, dnn = c("Within2", "Exprs"))

Indexes= split(1:sum(!is.na(mmGtex)), jMapLibd$code[!is.na(mmGtex)])
for(i in seq(along=Indexes)) {
	cat(".")
	png(paste0("plots/MA_plots_", names(Indexes)[i],"_DLPFC.png"),
		height = 500, width=500)
	par(mar=c(5,6,3,2))
	palette(brewer.pal(8,"Dark2"))
	ii = Indexes[[i]]
	plot(Md ~ Ad, subset=ii, pch=21,bg=i,cex=0.7,
		main = names(Indexes)[i], 
		xlab="Mean Expression (log2)",
		ylab="LIBD - GTEX DLPFC log2FC",  cex.axis=2, 
		cex.lab=2,cex.main=2)
	dev.off()
}
lapply(Indexes, function(ii) table(M[ii] < 1,
	A[ii] > 1, dnn = c("Within2", "Exprs")))

### write out replication rates
repJxnOut = as.data.frame(mcols(jMapLibd)[,c("code", "inGeuv", "inGtex")])
rownames(repJxnOut) = names(jMapLibd)	
save(repJxnOut, file="rdas/junction_replication_GtexGeuvadis.rda")

	
##### OLD #####

# cortecon
load("/dcs01/ajaffe/LIBD_StemCell/Data/processed_RPKMs_Cortecon.rda")
pdCort = pd ; jRpkmCort = jRpkm ; jMapCort = jMap

mmCort = match(names(jMapLibd), names(jMapCort))
jMapLibd$inCortecon = !is.na(mmCort)
table(jMapLibd$inCortecon, jMapLibd$code, jMapLibd$meanExprs > 1)

## meissner
load("/users/ajaffe/Lieber/Projects/RNAseq/LibdStem/meissner/rdas/junctionCounts_meissner.rda")

mmMeis = match(names(jMapLibd), names(juncCounts$anno))
jMapLibd$inMeissner = !is.na(mmMeis)
table(jMapLibd$inMeissner, jMapLibd$code, jMapLibd$meanExprs> 1)

## body map
xx=load("/users/ajaffe/Lieber/Projects/RNAseq/AS3MT/rdas/junctionCounts_bodyMap.rda")

mmBody = match(names(jMapLibd), names(juncCounts$anno))
jMapLibd$inBodyMap = !is.na(mmBody)
table(jMapLibd$inBodyMap, jMapLibd$code, jMapLibd$meanExprs > 1)

## all 3
table(jMapLibd$inBodyMap | jMapLibd$inCortecon | jMapLibd$inMeissner,
	jMapLibd$code)
prop.table(table(jMapLibd$inBodyMap | jMapLibd$inCortecon | jMapLibd$inMeissner,
	jMapLibd$code),2)
table(jMapLibd$inBodyMap | jMapLibd$inCortecon | jMapLibd$inMeissner,
	jMapLibd$code, jMapLibd$meanExprs  > 1)
prop.table(table(jMapLibd$inBodyMap | jMapLibd$inCortecon | jMapLibd$inMeissner,
	jMapLibd$code),2)
table(jMapLibd$inBodyMap | jMapLibd$inCortecon | jMapLibd$inMeissner,
	jMapLibd$code, jMapLibd$meanExprs  > 5)

## all 3
table(jMapLibd$inBodyMap & jMapLibd$inCortecon & jMapLibd$inMeissner,
	jMapLibd$code)
prop.table(table(jMapLibd$inBodyMap & jMapLibd$inCortecon & jMapLibd$inMeissner,
	jMapLibd$code),2)
table(jMapLibd$inBodyMap & jMapLibd$inCortecon & jMapLibd$inMeissner,
	jMapLibd$code, jMapLibd$meanExprs  > 1)
table(jMapLibd$inBodyMap & jMapLibd$inCortecon & jMapLibd$inMeissner,
	jMapLibd$code, jMapLibd$meanExprs  > 5)
	
## what junctions???
checkTrans = jMapLibd[which((jMapLibd$inBodyMap | 
	jMapLibd$inCortecon | jMapLibd$inMeissner) & 
		jMapLibd$meanExprs > 5 & jMapLibd$code == "NovelTrans")]
length(unique(checkTrans$newGeneID))			

tab = table(checkTrans$newGeneSymbol)
sort(tab,decreasing=TRUE)[1:50]

checkJxn = jMapLibd[which((jMapLibd$inBodyMap | 
	jMapLibd$inCortecon | jMapLibd$inMeissner) & 
		jMapLibd$meanExprs > 5 & jMapLibd$code == "NovelJxn")]
length(intersect(unique(checkTrans$newGeneID),
	unique(checkJxn$newGeneID)))
tab2 = table(checkJxn$newGeneSymbol)	
sort(tab2,decreasing=TRUE)[1:50]
