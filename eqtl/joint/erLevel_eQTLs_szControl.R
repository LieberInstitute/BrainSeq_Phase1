###
source("/users/ajaffe/Lieber/lieber_functions_aj.R")

library(MatrixEQTL)
library(GenomicRanges)
library(sva)
library(derfinder)

#### load data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/snpData_LIBD_szControls_cleaned.rda")
pd$DxNum = as.numeric(factor(pd$Dx))-1

##### log xform
y = log2(regionMat + 1)

# filter for Age
aIndex= which(pd$Age > 13)
pd2= pd[aIndex,]
snp2 = as.matrix(snp[,aIndex])
y2 = y[,aIndex]
mod = model.matrix(~snpPC1 + snpPC2 + snpPC3 + Sex,data=pd2)

# ###################
# # do cis eQTLS, maybe 250kb? 1M
# pcaDer = prcomp(t(y2))
# nsvDer = num.sv(y2, mod, vfilter=1e5)
# pcsDer = pcaDer$x[,1:nsvDer]

# covsDer = SlicedData$new(t(cbind(mod[,-1], pcsDer, pd2$DxNum)))
# exprsDer = SlicedData$new(y2)
# exprsDer$ResliceCombined(sliceSize = 1000)

# theSnps = SlicedData$new(snp2)
# theSnps$ResliceCombined(sliceSize = 50000)

# ### JUST CIS
# snpspos = snpMap[,c("SNP","CHR","POS")]
# snpspos$CHR = paste0("chr",snpspos$CHR)
# colnames(snpspos) = c("name","chr","pos")

# posDer = as.data.frame(regions)[,1:3]
# colnames(posDer)[1] = "chr"
# posDer$name = names(regions)
# posDer = posDer[,c(4,1:3)]

# meDerCis = Matrix_eQTL_main(snps=theSnps, gene = exprsDer, 
	# cvrt = covsDer, output_file_name.cis =  ".txt" ,
	# pvOutputThreshold.cis = 0.001, pvOutputThreshold=0,
	# snpspos = snpspos, genepos = posDer, 
	# useModel = modelLINEAR,	cisDist=5e5,
	# pvalue.hist = 100,min.pv.by.genesnp = TRUE)

# save(meDerCis, nsvDer, pcsDer, file="rdas/der_eqtl_szControl_cisOnly.rda")

#### annotate
load("rdas/der_eqtl_szControl_cisOnly.rda")
eqtl = meDerCis$cis$eqtl

# anntoate
eqtl$gene = as.character(eqtl$gene)
eqtl$snps = as.character(eqtl$snps)
colnames(eqtl)[2] = "der"
eqtl = eqtl[eqtl$FDR < 0.01,]

### number of effective tests
pruned = read.table("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/plink/LIBD_Brain_DLPFC_szControls_imputed_indep.prune.in")
snpMap$isLdIndep = snpMap$SNP %in% pruned$V1
snpMapGR = GRanges(paste0("chr", snpMap$CHR), IRanges(snpMap$POS, width=1))

nTests = length(findOverlaps(snpMapGR[snpMap$isLdIndep], 
	regions, maxgap=5e5))
eqtl$bonf = p.adjust(eqtl$pvalue, "bonferroni", n= nTests)
length(unique(eqtl$der[eqtl$bonf < 0.05])) # 21448

# snp coordinates
m = match(eqtl$snps, snpMap$SNP)
eqtl$snp_chr = snpMap$CHR[m]
eqtl$snp_chr = paste0("chr",eqtl$snp_chr)
eqtl$snp_chr = gsub("chr23","chrX", eqtl$snp_chr)
eqtl$snp_pos = snpMap$POS[m]
eqtl$snpRsNum = snpMap$name[m]

# take best eqtl per der
sig = eqtl[!duplicated(eqtl$der),]
tt = table(eqtl$der)
sig$numSnps = tt[sig$der]

# other snp info
snprange = t(sapply(split(eqtl$snp_pos, eqtl$der), range))
snprange= snprange[as.character(sig$der),]
sig$startSnpLD = snprange[,1]
sig$endSnpLD = snprange[,2]
sig$snpLDLength = snprange[,2]-snprange[,1]

# other snp annotation
m = match(sig$snps, snpMap$SNP)
sig$snpCounted = snpMap$COUNTED[m]
sig$snpAlt = snpMap$ALT[m]
sig$inSampleMAF = rowSums(snp2[m,],na.rm=TRUE)/
	(2*rowSums(!is.na(snp2[m,])))
## annotate
g = as.data.frame(regions)[sig$der,c(1:4,6,12:17)]
colnames(g)[c(1,5)] = c("chr", "meanCov")
colnames(g) = paste0("exprs_", colnames(g))
sig = cbind(sig,g)
rownames(sig) = NULL

## cis or trans
sig$cisOrTrans = ifelse(sig$snp_chr == sig$exprs_chr,"cis","trans")
start = sig$exprs_start
start[which(sig$exprs_strand == "-")] = sig$exprs_end[
	which(sig$exprs_strand == "-")]
sig$cisDistToDer  = sig$snp_pos - start
sig$cisDistToDer[sig$cisOrTrans=="trans"]=NA

### add race sensitivity analyses
rIndexes=splitit(pd2$Race)
raceGene = mclapply(rIndexes, function(ii) {
	ssnp = as.matrix(snp2[sig$snps,ii])
	yyGene = as.matrix(y2[sig$der,ii])
	out = matrix(nrow = nrow(ssnp), nc = 3)
	for(j in 1:nrow(ssnp)) {
		if(j %% 5000 == 0) cat(".")
		out[j,] = summary(lm(yyGene[j,]~ssnp[j,] + cbind(
			mod[ii,-1],pcsDer[ii,], pd2$DxNum[ii])))$coef[2,c(1,3,4)]
	}
	colnames(out) = c("slope","tstat","pval")
	cat("\n")
	return(out)
},mc.cores=2)
xx = do.call("cbind", raceGene)
colnames(xx) = paste0(rep(names(rIndexes), each=3), "_", colnames(xx))
sigDer = cbind(sig, xx)

mafRace = sapply(rIndexes, function(ii){
	ssnp = as.matrix(snp2[sig$snps,ii])
	rowSums(ssnp, na.rm=TRUE)/(2*rowSums(!is.na(ssnp)))
})
colnames(mafRace) = paste0(colnames(mafRace),"_inSampleMAF")
sigDer = cbind(sigDer, mafRace)

#### add fetal effect
fIndex=which(pd$Age < 0)
pdFetal= pd[fIndex,]
pcaFetal = prcomp(t(log2(regionMat[,fIndex]+1)))
modFetal = model.matrix(~ pdFetal$snpPC1+
	pdFetal$snpPC2 + pdFetal$snpPC3+ pdFetal$Sex)
nsvFetal = num.sv(log2(regionMat[,fIndex]+1), modFetal)

snpFetal = as.matrix(snp[match(sigDer$snps,rownames(snp)),fIndex])
derExprsRpkmFetal = as.matrix(log2(regionMat[sigDer$der,fIndex]+1))
outFetal = matrix(nrow = nrow(sigDer), nc = 3)

for(j in 1:nrow(sigDer)) {
	if(j %% 5000 == 0) cat(".")
	outFetal[j,] = summary(lm(derExprsRpkmFetal[j,]~snpFetal[j,] + 
		cbind(modFetal[,-1], pcaFetal$x[,1:nsvFetal])))$coef[2,c(1,3,4)]
}
colnames(outFetal) = paste0("fetal_", c("slope","tstat","pval"))
sigDer = cbind(sigDer, outFetal)
sigDer$fetal_MAF = rowSums(snpFetal, na.rm=TRUE)/(2*rowSums(!is.na(snpFetal)))

#### add postnatal effect
pIndex=which(pd$Age > 0 & pd$Age < 13)
pdPost= pd[pIndex,]
pcaPost = prcomp(t(log2(regionMat[,pIndex]+1)))
modPost = model.matrix(~ pdPost$snpPC1+
	pdPost$snpPC2 + pdPost$snpPC3+ pdPost$Sex)
nsvPost = num.sv(log2(regionMat[,pIndex]+1), modPost)

snpPost = as.matrix(snp[match(sigDer$snps,rownames(snp)),pIndex])
derExprsRpkmPost = as.matrix(log2(regionMat[sigDer$der,pIndex]+1))
outPost = matrix(nrow = nrow(sigDer), nc = 3)

for(j in 1:nrow(sigDer)) {
	if(j %% 5000 == 0) cat(".")
	outPost[j,] = summary(lm(derExprsRpkmPost[j,]~snpPost[j,] + 
		cbind(modPost[,-1], pcaPost$x[,1:nsvPost])))$coef[2,c(1,3,4)]
}
colnames(outPost) = paste0("postnatal_", c("slope","tstat","pval"))
sigDer = cbind(sigDer, outPost)
sigDer$postnatal_MAF = rowSums(snpPost, na.rm=TRUE)/(2*rowSums(!is.na(snpPost)))

## save output
save(sigDer, file="rdas/annotated_der_eqtl_szControl_cisOnly.rda")
