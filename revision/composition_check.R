###

## load expression data

## load data
load("/dcl01/lieber/ajaffe/lab/libd_stem_timecourse/deconvolution/brainSeq_compEsts_viaSingleCell.rda")

## load pheno
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/phenotype_annotated_szControlEqtl_DLPFC.rda")

## merge
mm = match(pd$RNum, rownames(dlpfcPropEsts))
dlpfcPropEsts = dlpfcPropEsts[mm[!is.na(mm)],]
pd = pd[!is.na(mm),]
dlpfcPropEsts = as.data.frame(dlpfcPropEsts)

## cut age
pd$ageGroup = cut(pd$Age, c(-1,0,1,10,20,50,100))

# make plots
boxplot(dlpfcPropEsts$NPC ~ pd$ageGroup)
boxplot(dlpfcPropEsts$Neurons ~ pd$ageGroup)

## load DNAm data
pheno = read.csv("/users/ajaffe/Lieber/Projects/450k/ECD2014/final_pheno_data.csv",
	as.is=TRUE)
dnamCounts = pheno[match(pd$BrNum,pheno$BrNum),23:27]

plot(dlpfcPropEsts$Neurons, dnamCounts$NeuN_pos)
plot(dlpfcPropEsts$NPC, dnamCounts$NPC)

####################################
### cell type specificity ###########
 
### retain specific single cell samples
load("/dcl01/lieber/ajaffe/PublicData/Brain/Darmanis_Quake/rdas/rpkmCounts_darmanisSingleCell.rda")

library(limma)
library(genefilter)

## drop hydrid 
hybridIndex=which(pd$Cell_type == "Hybrid")
pd = pd[-hybridIndex,]
geneRpkm = geneRpkm[,-hybridIndex]
pd$Cell_type = factor(pd$Cell_type)

## read in isoform shift genes
load("../devel/rdas/isoform_switch_devel_byFeature.rda")
load("../devel/rdas/devStats_controlSamples.rda")

geneStat = statList$Gene
geneStat = geneStat[!is.na(geneStat$p_bonf),]

## get F-stat for cell type
ff = rowFtests(log2(geneRpkm[names(geneStat),]+1), pd$Cell_type)

## check overall
plot(-log10(geneStat$f_pval), -log10(ff$p.value))
cor(-log10(geneStat$f_pval), -log10(ff$p.value),use="comp")
plot(log10(geneStat$fstat), log10(ff$statistic))

## switch genes?
switchGenes = lapply(switchList,rownames)
sapply(switchGenes, function(x) {
	z = rep(0,length(geneStat))
	z[names(geneStat) %in% x] = 1
	summary(lm(log10(ff$statistic) ~ z + log2(geneStat$meanRpkm+1)))$coef[2,]
})
	