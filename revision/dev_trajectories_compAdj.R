####
# R-devel
source("../eqtl_functions.R")

library(GenomicRanges)
library(limma)

# gene exon junction
load("/dcs01/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
ageGroup = cut(pd$Age, c(-1,0,1,10,20,50,100))
geneRpkmBrain = geneRpkm # just gene
pdBrain = pd
rm(exonRpkm, jRpkm, jMap, exonMap)

### add entrez ID
geneMapGR = makeGRangesFromDataFrame(geneMap,keep=TRUE)

### filter for control
controlIndex = which(pdBrain$Dx == "Control")
pdBrain = pdBrain[controlIndex,]
geneRpkmBrain = geneRpkmBrain[,controlIndex]

###### filter gene, exon, junction
exprsGeneIndex = which(rowMeans(geneRpkmBrain) > 0.01)

###############################
### model, linear spline
fetal = ifelse(pdBrain$Age < 0, 1,0)
birth = pdBrain$Age
birth[birth < 0] = 0 # linear spline
infant = pdBrain$Age - 1
infant[infant < 0] = 0 # linear spline
child = pdBrain$Age - 10
child[child < 0] = 0 # linear spline
teen = pdBrain$Age - 20
teen[teen < 0] = 0 # linear spline
adult = pdBrain$Age - 50
adult[adult < 0] = 0 # linear spline

### adjust for race, assign, sex, and mito
mod = model.matrix(~Age + fetal + birth + infant +
	child + teen + adult + snpPC1 + snpPC2 + snpPC3 + 
	 Sex, data=pdBrain)
mod0 = model.matrix(~ snpPC1 + snpPC2 + snpPC3 + 
	 Sex, data=pdBrain)

## get composition estimates
### retain specific single cell samples
load("/dcl01/lieber/ajaffe/PublicData/Brain/Darmanis_Quake/rdas/rpkmCounts_darmanisSingleCell.rda")
geneRpkmSingle = geneRpkm
pdSingle = pd

## drop hydrid 
hybridIndex=which(pdSingle$Cell_type == "Hybrid")
pdSingle = pdSingle[-hybridIndex,]
geneRpkmSingle = geneRpkmSingle[,-hybridIndex]
pdSingle$Cell_type = factor(pdSingle$Cell_type)
geneExprsSingle = log2(geneRpkmSingle+1)

library(genefilter)
tIndexes <- splitit(pdSingle$Cell_type)
tstatList <- lapply(tIndexes, function(i) {
	x <- rep(0, ncol(geneExprsSingle))
	x[i] <- 1
    return(rowttests(geneExprsSingle, factor(x)))
})
numProbes=20
probeList <- lapply(tstatList, function(x) {
	y <- x[x[, "p.value"] < 1e-08, ]
	yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
	yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
	 c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
})

trainingProbes <- unique(unlist(probeList))
geneProfiles <- geneExprsSingle[trainingProbes, ]
pMeans <- colMeans(geneProfiles)
names(pMeans) <- pdSingle$Cell_type

form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pdSingle$Cell_type),
	collapse = "+")))
phenoDF <- as.data.frame(model.matrix(~pdSingle$Cell_type - 1))
colnames(phenoDF) <- sub("^pdSingle\\$Cell_type", "", colnames(phenoDF))
tmp <- minfi:::validationCellType(Y = geneProfiles, 
		pheno = phenoDF, modelFix = form)
coefEsts <- tmp$coefEsts
save(coefEsts, file = "BrainSeq_Phase1_DLPFC_darmanis_hg19-based_cellComp.rda")
write.csv(coefEsts, file = "BrainSeq_Phase1_DLPFC_darmanis_hg19-based_cellComp.csv")

## estimate cell type
geneExprsBrain = log2(geneRpkmBrain+1)
compEsts = minfi:::projectCellType(geneExprsBrain[rownames(coefEsts),], coefEsts)
compEsts = as.data.frame(compEsts)
write.csv(compEsts, file="BrainSeq_Phase1_DLPFC_darmanis_hg19_estimatedComp.csv")
## 
compEsts = read.csv("BrainSeq_Phase1_DLPFC_darmanis_hg19_estimatedComp.csv", row.names=1)

## plot each comp versus age
nms = c("Adult: Astrocytes", "Adult: Endothelial", 
	"Fetal: Quiescent Neurons", "Fetal: Replicating Neurons",
	"Adult: Microglia", "Adult: Neurons", "Adult: Oligodendrocytes",
	"Adult: OPCs")
pdf("cellComp_overAge.pdf", h=4,w=5.5)
for(i in 1:ncol(compEsts)) {
	agePlotter(compEsts[,i], pdBrain$Age, mod, 
		smoothIt=FALSE, ylab="RNA proportion", 	mainText=nms[i],
		ylim = c(0,1), ageLabel=FALSE)
}
dev.off()

####
# do PCA of expression
pca = prcomp(t(geneExprsBrain))

signif(cor(pca$x[,1:15], compEsts),3)

pdf("pc1_vs_FetalQuiescent_Comp.pdf")
palette(brewer.pal(3,"Set1"))
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
plot(1-compEsts$Fetal_replicating, pca$x[,1],
	xlab = "(1 - Fetal Replicating Comp)", ylab="PC1",
		pch = 21, bg=2,cex=1.3)
plot(1-compEsts$Fetal_replicating[pdBrain$Age < 0], 
	pca$x[pdBrain$Age < 0,1],
	xlab = "(1 - Fetal Quiescent Comp)", ylab="PC1",
		pch = 21, bg=2, cex=1.3)
plot(1-compEsts$Fetal_quiescent, pca$x[,1],
	xlab = "(1 - Fetal Quiescent Comp)", ylab="PC1",
		pch = 21, bg=2,cex=1.3)
plot(1-compEsts$Fetal_quiescent[pdBrain$Age < 0], 
	pca$x[pdBrain$Age < 0,1],
	xlab = "(1 - Fetal Quiescent Comp)", ylab="PC1",
		pch = 21, bg=2, cex=1.3)
plot(1-(rowSums(compEsts[,c("Fetal_quiescent", "Fetal_replicating")])),
	pca$x[,1],	xlab = "(1 - Fetal Cell Proportion)", ylab="PC1",
		pch = 21, bg=2,cex=1.3)
plot(1-(rowSums(compEsts[pdBrain$Age < 0,c("Fetal_quiescent", "Fetal_replicating")])),
	pca$x[pdBrain$Age < 0,1], xlab = "(1 - Fetal Cell Proportion)", ylab="PC1",
		pch = 21, bg=2, cex=1.3)
dev.off()

plot(compEsts$Fetal_quiescent ~ log(pdBrain$Age+1))
plot(pca$x[,1] ~ compEsts$Fetal_quiescent)

plot(pca$x[pdBrain$Age < 0,1] ~ compEsts$Fetal_quiescent[pdBrain$Age < 0])
plot(pca$x[pdBrain$Age < 0,1] ~ compEsts$Fetal_replicating[pdBrain$Age < 0])

plot(pca$x[,1] ~ compEsts$Neurons)

#### test fit
# agePlotter(as.numeric(geneRpkm[1,]), pdBrain$Age, mod, 
	# ylab=rownames(geneRpkm)[1], mainText="")

##############
### model ####
modComp = cbind(mod, compEsts)
mod0Comp = cbind(mod0, compEsts)

fitComp = lmFit(geneExprsBrain, modComp)
fit0Comp = lmFit(geneExprsBrain, mod0Comp)
ffComp = getF(fitComp,fit0Comp, geneExprsBrain)
ffComp$p_bonf = NA
ffComp$p_bonf[exprsGeneIndex] = p.adjust(ffComp$f_pval[exprsGeneIndex], "bonf")
sum(ffComp$p_bonf < 0.05, na.rm=TRUE)

## compare to previous
fit = lmFit(geneExprsBrain, mod)
fit0 = lmFit(geneExprsBrain, mod0)
ff = getF(fit,fit0, geneExprsBrain)
ff$p_bonf = NA
ff$p_bonf[exprsGeneIndex] = p.adjust(ff$f_pval[exprsGeneIndex], "bonf")
sum(ff$p_bonf < 0.05, na.rm=TRUE)

plot(log2(ff$fstat[exprsGeneIndex]), 
	log2(ffComp$fstat[exprsGeneIndex]),
	pch = 21, bg="grey", xlim = c(-5,10), ylim=c(-5,10),
	xlab="Original Model", ylab="Comp-Adj Model")
abline(0,1,col="blue", lty=2)

table(ff$p_bonf < 0.05, ffComp$p_bonf < 0.05,
	dnn = c("Original","Comp-Adj"))

### just cell type
fit_onlyComp = lmFit(geneExprsBrain, cbind(mod0, compEsts))
fit0_onlyComp = lmFit(geneExprsBrain, mod0)
ff_onlyComp = getF(fit_onlyComp,fit0_onlyComp, geneExprsBrain)
ff_onlyComp$p_bonf = NA
ff_onlyComp$p_bonf[exprsGeneIndex] = p.adjust(ff_onlyComp$f_pval[exprsGeneIndex], "bonf")
mean(ff_onlyComp$p_bonf < 0.05, na.rm=TRUE)
table(ff$p_bonf < 0.05, ff_onlyComp$p_bonf < 0.05,
	dnn = c("Original","OnlyComp"))
prop.table(table(ff$p_bonf < 0.05, ff_onlyComp$p_bonf < 0.05,
	dnn = c("Original","OnlyComp")),1)

pdf("suppFigure_compVsSpline_Fstat_models.pdf")
par(mar=c(5,6,3,2), cex.axis=1.7,cex.lab=1.7,cex.main = 1.6)
plot(-log10(ff$f_pval[exprsGeneIndex]), 
	-log10(ff_onlyComp$f_pval[exprsGeneIndex]),
	pch = 21, bg="grey", xlim = c(0,300), ylim=c(0,300),
	xlab="Age Spline Model", ylab="Composition Model",
	main = "-log10 p-values")
abline(0,1,col="blue", lty=1,lwd=3)
dev.off()
### DNAm estimates
pheno_DNAm = read.csv("/users/ajaffe/Lieber/Projects/450k/ECD2014/final_pheno_data.csv",
	as.is=TRUE)	
pheno_DNAm = pheno_DNAm[match(pdBrain$BrNum, pheno_DNAm$BrNum),]

pdf("compCompare_RNAvsDNA.pdf")
par(mar=c(5,6,3,2), cex.axis=1.6,cex.lab=1.6,cex.main=1.5)
plot(compEsts$Fetal_quiescent,pheno_DNAm$NPC,
	pch = 21, bg=ifelse(pheno_DNAm$Age < 0, 1,2),
	xlab = "RNA-based (Fetal Quiescent)",
	ylab= "DNAm-based (NPC)", main = "All Controls") 
plot(pca$x[pdBrain$Age<0,1],pheno_DNAm$NPC[pdBrain$Age<0],
	pch = 21, bg=ifelse(pheno_DNAm$Age[pdBrain$Age<0] < 0, 1,2),
	xlab = "PC1 of Expression",
	ylab= "DNAm-based (NPC)", main = "Prenatal Controls") 

plot(compEsts$Neurons,pheno_DNAm$NeuN_pos,
	pch = 21, bg=ifelse(pheno_DNAm$Age < 0, 1,2),
	xlab = "RNA-based (Adult Neurons)",
	ylab= "DNAm-based (Adult Neurons)",
	main = "All Controls") 
plot(compEsts$Neurons[pheno_DNAm$Age>18],
	pheno_DNAm$NeuN_pos[pheno_DNAm$Age>18],
	pch = 21, bg=2,main = "Controls > 18 yrs",
	xlab = "RNA-based (Adult Neurons)",
	ylab= "DNAm-based (Adult Neurons)")
dev.off()

cor.test(compEsts$Fetal_quiescent,pheno_DNAm$NPC)
cor.test(compEsts$Fetal_quiescent[pheno_DNAm$Age<0],
	pheno_DNAm$NPC[pheno_DNAm$Age<0])
cor.test(compEsts$Fetal_quiescent[pheno_DNAm$Age>0],
	pheno_DNAm$NPC[pheno_DNAm$Age>0])

cor.test(compEsts$Neurons,pheno_DNAm$NeuN_pos)
cor.test(compEsts$Neurons[pheno_DNAm$Age<0],
	pheno_DNAm$NeuN_pos[pheno_DNAm$Age<0])
cor.test(compEsts$Neurons[pheno_DNAm$Age>0],
	pheno_DNAm$NeuN_pos[pheno_DNAm$Age>0])

summary(lm(compEsts$Neurons ~ pheno_DNAm$NeuN_pos,
	subset = pheno_DNAm$Age>18))