##
library(limma)

load("../caseControl/rdas/geneLevel_LIBD_qSVsAndMod.rda")
load("../caseControl/rdas/degradation_mat_LIBD_polyA.rda")
identical(pd$RNum, colnames(degCovAdj)) # TRUE

fit_libd = lmFit(log2(degCovAdj+1), mod)
eb_libd = eBayes(fit_libd)
de_libd = topTable(eb_libd, coef=2, number = nrow(degCovAdj))
hist(de_libd$P.Value)
sum(de_libd$adj.P.Val < 0.05)

### read in pd
load("../phenotype_annotated_szControlEqtl_DLPFC.rda")
controlIndex = which(pd$Dx == "Control")
pdControl = pd[controlIndex,]

degCovAdj = read.csv("/users/ajaffe/Lieber/Projects/RNAseq/Consortium/PhaseI/degradationMat_DLPFC_polyA_Phase1.csv",
	as.is=TRUE, row.names=1)
degCovAdj = degCovAdj[,pdControl$RNum]

###############################
### model, linear spline
fetal = ifelse(pdControl$Age < 0, 1,0)
birth = pdControl$Age
birth[birth < 0] = 0 # linear spline
infant = pdControl$Age - 1
infant[infant < 0] = 0 # linear spline
child = pdControl$Age - 10
child[child < 0] = 0 # linear spline
teen = pdControl$Age - 20
teen[teen < 0] = 0 # linear spline
adult = pdControl$Age - 50
adult[adult < 0] = 0 # linear spline

### adjust for race, assign, sex, and mito
mod = model.matrix(~Age + fetal + birth + infant +
	child + teen + adult + snpPC1 + snpPC2 + snpPC3 + 
	 Sex, data=pdControl)
mod0 = model.matrix(~ snpPC1 + snpPC2 + snpPC3 + 
	 Sex, data=pdControl)
	
fit_dev = lmFit(log2(degCovAdj+1), mod)
eb_dev = eBayes(fit_dev)
de_dev = topTable(eb_dev, coef=2:8, number = nrow(degCovAdj))
sum(de_dev$adj.P.Val < 0.05)

summary(lm(RIN ~ Age, data=pdControl))
t.test(RIN ~ Age < 0, data=pdControl)
t.test(mitoRate ~ Age < 0, data=pdControl)$p.value


# gene assignment rate
library(jaffelab)
summ = read.delim("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/Counts/Ensembl_Exons_n746.counts.summary",
	row.names=1)
colnames(summ) = ss(colnames(summ), "_", 4)
summ = summ[,pdControl$RNum]
pdControl$assignRate = as.numeric(summ[1,] / colSums(summ))
t.test(assignRate ~ Age < 0, data=pdControl)$p.value

# check age ~ chrM
cor.test(pdControl$mitoRate[pd$Age < 0], pdControl$RIN[pd$Age < 0])$p.value
cor.test(pdControl$mitoRate[pd$Age > 0], pdControl$RIN[pd$Age > 0])$p.value