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
mod0 = model.matrix(~ snpPC1 + snpPC2 + snpPC3 + 
	 Sex, data=pd)

#### test fit
# agePlotter(as.numeric(geneRpkm[1,]), pd$Age, mod, 
	# ylab=rownames(geneRpkm)[1], mainText="")

##### xform
xformData = list(Gene = log2(geneRpkm+1), Exon = log2(exonRpkm+1), 
	Junction = log2(jRpkm+1), Transcript =  log2(tFpkm+1), 
	DER = log2(regionMat+1))

##############
### model ####

statList = lapply(xformData, function(y) {
	cat(".")
	fit = lmFit(y, mod)
	eb = ebayes(fit)
	fit0 = lmFit(y, mod0)
	ff = getF(fit,fit0, y)
	yResid = cleaningY(y, mod, P=8)
	ff$ageCorr = cor(t(yResid), pd$Age) # for directionality switch
	ff$meanRpkm = 2^rowMeans(y)-1
	return(ff)
})

## add multiple testing corrected p-values
exprsIndexList = list(Gene=exprsGeneIndex, Exon = exprsExonIndex,
	Junction= exprsJxnIndex, Transcript = 1:nrow(xformData$Transcript), 
	DER = 1:nrow(xformData$DER))
	
# number features expressed
sapply(exprsIndexList, length)

## add bonferroni
statList = mapply(function(x, n) {
	cat(".")
	x$p_bonf = NA
	x$p_bonf[n] = p.adjust(x$f_pval[n], "bonf")
	return(x)
}, statList, exprsIndexList, SIMPLIFY=FALSE)

## just fetal versus postnatal
modFetal = model.matrix(~ifelse(pd$Age < 0, 1, 0))
statList = mapply(function(y, dat) {
	cat(".")
	fit = lmFit(y, modFetal)
	dat$log2FC_fetal = fit$coef[,2]
	eb = ebayes(fit)
	dat$pval_fetal = eb$p[,2]
	return(dat)
}, xformData, statList,SIMPLIFY=FALSE)

#########################
#### save stats #########

### add gene symbols
outGene = makeGRangesFromDataFrame(	
	cbind(geneMap, statList$Gene), keep=TRUE)
outGene$EnsemblGeneID = names(outGene)
outGene$code = "InEns"
outGene$meanControl = rowMeans(geneRpkm)

outExon = makeGRangesFromDataFrame(	
	cbind(exonMap, statList$Exon), keep=TRUE)
outExon$code = "InEns"
outExon$meanControl = rowMeans(exonRpkm)

outTx = tMap 
mcols(outTx) = cbind(mcols(tMap), statList$Transcript)
outTx$meanControl = rowMeans(tFpkm)
code = tMap$class_code
outTx$code = "Novel"
outTx$code[code == "="] = "InEns"
outTx$code[code == "j"] = "ExonSkip"
outTx$code[code %in% c("c", "e", "o", "p", "x") ] = "AltStartEnd"

outJxn = jMap
mcols(outJxn) = cbind(mcols(jMap), statList$Junction)
outJxn$meanControl = rowMeans(jRpkm)

outEr = regions
mcols(outEr) = cbind(mcols(regions), statList$DER)
outEr$meanControl = rowMeans(regionMat) 
code2 = regions$annoClass
outEr$code = "Novel"
outEr$code[code2 == "strictExonic"] = "InEns"
outEr$code[code2 %in% c("exonIntron", "extendUTR")] = "AltStartEnd"

colnames(mcols(outExon))[1] = "EnsemblGeneID"
colnames(mcols(outTx))[3] = "Symbol"
colnames(mcols(regions))[7] = "Symbol"
colnames(mcols(outJxn))[12:13] = c("EnsemblGeneID","Symbol")

n = c("Symbol", "EnsemblGeneID", "EntrezID", "meanControl",
	colnames(statList$Gene))
statListGR = GRangesList(c(Gene = outGene[,n], Exon = outExon[,n],
	Transcript = outTx[,n], Junction = outJxn[,n],
	DER = outEr[,n]))
save(statList, file="rdas/devStats_controlSamples.rda")

##########################
### metrics/output #######
##########################
load("rdas/devStats_controlSamples.rda")

# number features significant
sapply(statList, function(x) sum(!is.na(x$p_bonf)))
sapply(statList, function(x) sum(x$p_bonf < 0.05, na.rm=TRUE))
sapply(statList, function(x) mean(x$p_bonf < 0.05, na.rm=TRUE))
sapply(statList, function(x) max(x$f_pval[which(x$p_bonf < 0.05)]))

# number of unique genes, by ID
sapply(statList, function(x) {
	g = x$EnsemblGeneID[which(x$p_bonf < 0.05)]
	g = g[!is.na(g) & !grepl("-", g)]
	length(unique(g))
})

# number of unique genes, by symbol
sapply(statList, function(x) {
	g = x$Symbol[which(x$p_bonf < 0.05)]
	g = g[!is.na(g) & !grepl("-", g)]
	length(unique(g))
})

# mean expression 
t(sapply(statList, function(x) {
	quantile(x$meanControl[which(x$p_bonf < 0.05)],
		prob = c(0.5,0.25,0.75))
}))

# mean abs correlation 
t(sapply(statList, function(x) {
	quantile(abs(x$ageCorr[which(x$p_bonf < 0.05)]),
		prob = c(0.5,0.25,0.75))
}))

# novel
lapply(statList, function(x) {
	prop.table(table(x$code[which(x$p_bonf < 0.05)]))
})
lapply(statList, function(x) {
	table(x$code[which(x$p_bonf < 0.05)])
})

##############################################
## venn diagram of the IDs by the 5 features
allGenes = unique(unlist(lapply(statList, function(x) 
	x$EnsemblGeneID[which(x$p_bonf < 0.05)])))
allGenes = allGenes[!grepl("-", allGenes)] # drop fusions
allGenes = allGenes[!is.na(allGenes)] # universe

geneMat = sapply(statList, function(x) allGenes %in% x$EnsemblGeneID)
rownames(geneMat) = allGenes

pdf("plots/venn_geneIDs_devChanges.pdf")
vennDiagram(vennCounts(geneMat))
dev.off()

write.csv(geneMat, file="tables/suppTable4_devChanges_collapsedByGene.csv")

table(rowSums(geneMat))

## by symbol 
allGenesSym = unique(unlist(lapply(statList, function(x) 
	x$Symbol[which(x$p_bonf < 0.05)])))
allGenesSym = allGenesSym[!grepl("-", allGenesSym)] # drop fusions
allGenesSym = allGenesSym[!is.na(allGenesSym)] # universe
length(unique(allGenesSym))

##### novel transcription results #####
allStatsExprs = endoapply(statList, function(x) 
	x[!is.na(x$p_bonf)])
allStatsExprs = unlist(allStatsExprs)
allStatsExprs$Type = ss(names(allStatsExprs), "\\.")
names(allStatsExprs) =  ss(names(allStatsExprs), "\\.", 2)
table(allStatsExprs$Type, allStatsExprs$code)

### Write out as CSV for paper
exprsStatsDf = as.data.frame(allStatsExprs)
colnames(exprsStatsDf)[1] = "chr"
write.csv(exprsStatsDf, file=gzfile("tables/suppTable_develChanges_allExprsFeatures.csv.gz"))

# filter to significant
sigStatsExprs = allStatsExprs[allStatsExprs$p_bonf < 0.05]
length(sigStatsNovel)
### Write out as CSV for paper
sigStatsDf = as.data.frame(sigStatsExprs)
colnames(sigStatsDf)[1] = "chr"
write.csv(sigStatsDf, file=gzfile("tables/suppTable3_develChanges_sigExprsFeatures.csv.gz"))


sigStatsNovel = sigStatsExprs[sigStatsExprs$code %in% 
	c("ExonSkip", "AltStartEnd")]
table(	sigStatsNovel$Type, sigStatsNovel$code)
	
# by ID
novelTabId = t(table(sigStatsNovel$Type, 
	sigStatsNovel$EnsemblGeneID))
novelTabId = novelTabId[rownames(novelTabId)!="",]
novelTabId = novelTabId[!grepl("-", rownames(novelTabId)),]
dim(novelTabId)
table(rowSums(novelTabId > 0))
# by symbol
novelTabSym = t(table(sigStatsNovel$Type, sigStatsNovel$Symbol))
novelTabSym = novelTabSym[rownames(novelTabSym)!="",]
novelTabSym = novelTabSym[rownames(novelTabSym)!="",]
dim(novelTabSym)
table(rowSums(novelTabSym > 0))

## numbers for paper
table(sigStatsNovel$Type)
100*prop.table(table(sigStatsNovel$Type))
table(sigStatsNovel$Type, sigStatsNovel$code)

