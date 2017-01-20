###
# qsub -t 1-23 -cwd -l mem_free=50G,h_vmem=60G,h_stack=256M -b y R CMD BATCH --no-save transcriptFpkm_eQTLs_adults.R

source("/users/ajaffe/Lieber/lieber_functions_aj.R")

library(MatrixEQTL)
library(GenomicRanges)
library(sva)
library(ballgown)
library(GenomicFeatures)

#### load data
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/phenotype_annotated_szControlEqtl_DLPFC.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/transcript/transcript_data_filtered_n495.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/snpData_LIBD_szControls_cleaned.rda")
pd$DxNum = as.numeric(factor(pd$Dx))-1

# filter for Age
aIndex= which(pd$Age > 13)
pd2= pd[aIndex,]
snp2 = as.matrix(snp[,aIndex] )
tFpkm2 = as.matrix(log2(tFpkm[,aIndex]+1))

mod = model.matrix(~snpPC1 + snpPC2 + snpPC3 + Sex,data=pd2)

# ####################
# ##### transcript
# pcaTrans = prcomp(t(tFpkm2))
# nsvTrans = num.sv(tFpkm2, mod)
# pcsTrans = pcaTrans$x[,1:nsvTrans]

# covsTrans = SlicedData$new(t(cbind(mod[,-1], pcsTrans, pd2$DxNum)))
# exprsTrans = SlicedData$new(as.matrix(tFpkm2))
# exprsTrans$ResliceCombined(sliceSize = 5000)

# # JUST CIS
# theSnps = SlicedData$new(snp2)
# theSnps$ResliceCombined(sliceSize = 50000)

# snpspos = snpMap[,c("SNP","CHR","POS")]
# snpspos$CHR = paste0("chr",snpspos$CHR)
# colnames(snpspos) = c("name","chr","pos")

# posTrans = as.data.frame(tMap)[,1:3]
# posTrans$name = rownames(tFpkm2)
# posTrans = posTrans[,c(4,1:3)]
# names(posTrans)[2] = "chr"

# meTransCis = Matrix_eQTL_main(snps=theSnps, gene = exprsTrans, 
	# cvrt = covsTrans, output_file_name.cis =  ".txt" ,
	# pvOutputThreshold.cis = 0.001,  pvOutputThreshold=0,
	# snpspos = snpspos, genepos = posTrans, 
	# useModel = modelLINEAR,	cisDist=5e5,
	# pvalue.hist = 100,min.pv.by.genesnp = TRUE)
# save(meTransCis,nsvTrans, pcsTrans, file="rdas/transcript_eqtl_szControl_cisOnly.rda")

#### load back in
load("rdas/transcript_eqtl_szControl_cisOnly.rda")
eqtl = meTransCis$cis$eqtl

# annotate
eqtl$gene = as.character(eqtl$gene)
eqtl$snps = as.character(eqtl$snps)
colnames(eqtl)[2] = "tx"
eqtl = eqtl[order(eqtl$pvalue),]

# take best eqtl per transcript
eqtl = eqtl[eqtl$FDR < 0.01,]

###### number of effective tests
pruned = read.table("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/plink/LIBD_Brain_DLPFC_szControls_imputed_indep.prune.in")
snpMap$isLdIndep = snpMap$SNP %in% pruned$V1
snpMapGR = GRanges(paste0("chr", snpMap$CHR), IRanges(snpMap$POS, width=1))

nTests = length(findOverlaps(snpMapGR[snpMap$isLdIndep], 
	tMap, maxgap=5e5))
eqtl$bonf = p.adjust(eqtl$pvalue, "bonferroni", n= nTests)
length(unique(eqtl$tx[eqtl$bonf < 0.05])) # 7018

# clean up
sig = eqtl[!duplicated(eqtl$tx),]
tt = table(eqtl$tx)
sig$numSnps = tt[sig$tx]

# snp coordinates
m = match(eqtl$snps, snpMap$SNP)
eqtl$snp_chr = snpMap$CHR[m]
eqtl$snp_chr = paste0("chr",eqtl$snp_chr)
eqtl$snp_chr = gsub("chr23","chrX", eqtl$snp_chr)
eqtl$snp_pos = snpMap$POS[m]
eqtl$snpRsNum = snpMap$name[m]

# take best eqtl per tx
sig = eqtl[!duplicated(eqtl$tx),]
tt = table(eqtl$tx)
sig$numSnps = as.numeric(tt[sig$tx])

# other snp info
snprange = t(sapply(split(eqtl$snp_pos, eqtl$tx), range))
snprange= snprange[as.character(sig$tx),]
sig$startSnpLD = snprange[,1]
sig$endSnpLD = snprange[,2]
sig$snpLDLength = snprange[,2]-snprange[,1]

# other annotation
m = match(sig$snps, snpMap$SNP)
sig$snpCounted = snpMap$COUNTED[m]
sig$snpAlt = snpMap$ALT[m]
sig$inSampleMAF = rowSums(snp2[m,],na.rm=TRUE)/
	(2*rowSums(!is.na(snp2[m,])))
	
# annotate gene
g = as.data.frame(tMap[match(sig$tx, names(tMap))])
colnames(g)[1] = "chr"
colnames(g) = paste0("exprs_", colnames(g))
sig = cbind(sig,g)
rownames(sig) = NULL

## cis or trans
tss = ifelse(sig$exprs_strand =="+",
	sig$exprs_start, sig$exprs_end)
sig$cisOrTrans = ifelse(sig$snp_chr == sig$exprs_chr & 
	abs(sig$snp_pos - tss) < 1e6,	"cis","trans")
sig$distToTss  = sig$snp_pos - tss
sig$distToTss[sig$cisOrTrans=="trans"] = NA

sig$meanFPKM = 2^rowMeans(tFpkm2)[sig$tx]-1

### add race sensitivity analyses
rIndexes=splitit(pd2$Race)
raceTrans = mclapply(rIndexes, function(ii) {
	ssnp = as.matrix(snp2[sig$snps,ii])
	yyTx = as.matrix(tFpkm2[sig$tx,ii])
	out = matrix(nrow = nrow(ssnp), nc = 3)
	for(j in 1:nrow(ssnp)) {
		if(j %% 5000 == 0) cat(".")
		out[j,] = summary(lm(yyTx[j,]~ssnp[j,] + cbind(
			mod[ii,-1], pcsTrans[ii,])))$coef[2,c(1,3,4)]
	}
	colnames(out) = c("slope","tstat","pval")
	cat("\n")
	return(out)
},mc.cores=2)
xx = do.call("cbind", raceTrans)
colnames(xx) = paste0(rep(names(rIndexes), each=3), "_", colnames(xx))
sigTrans = cbind(sig, xx)

# race MAFs
mafRace = sapply(rIndexes, function(ii){
	ssnp = as.matrix(snp2[sig$snps,ii])
	rowSums(ssnp, na.rm=TRUE)/(2*rowSums(!is.na(ssnp)))
})
colnames(mafRace) = paste0(colnames(mafRace),"_inSampleMAF")
sigTrans = cbind(sigTrans, mafRace)

#### add fetal effect
fIndex=which(pd$Age < 0)
pdFetal= pd[fIndex,]
pcaFetal = prcomp(t(log2(tFpkm[,fIndex]+1)))
modFetal = model.matrix(~ pdFetal$snpPC1+
	pdFetal$snpPC2 + pdFetal$snpPC3+pdFetal$Sex)
nsvFetal = num.sv(log2(tFpkm[,fIndex]+1), modFetal)

snpFetal = as.matrix(snp[match(sigTrans$snps,rownames(snp)),fIndex])
tFpkmFetal = as.matrix(log2(tFpkm[sigTrans$tx,fIndex]+1))
outFetal = matrix(nrow = nrow(sigTrans), nc = 3)

for(j in 1:nrow(sigTrans)) {
	if(j %% 5000 == 0) cat(".")
	outFetal[j,] = summary(lm(tFpkmFetal[j,]~snpFetal[j,] + 
		cbind(modFetal[,-1], pcaFetal$x[,1:nsvFetal])))$coef[2,c(1,3,4)]
}
colnames(outFetal) = paste0("fetal_", c("slope","tstat","pval"))
sigTrans = cbind(sigTrans, outFetal)
sigTrans$fetal_MAF = rowSums(snpFetal, na.rm=TRUE)/(2*rowSums(!is.na(snpFetal)))

#### add postnatal effect
pIndex=which(pd$Age > 0 & pd$Age < 13)
pdPost= pd[pIndex,]
rownames(tFpkm) = tn
pcaPost = prcomp(t(log2(tFpkm[,pIndex]+1)))
modPost = model.matrix(~ pdPost$snpPC1+
	pdPost$snpPC2 + pdPost$snpPC3+pdPost$Sex)
nsvPost = num.sv(log2(tFpkm[,pIndex]+1), modPost)

snpPost = as.matrix(snp[match(sigTrans$snps,rownames(snp)),pIndex])
tFpkmPost = as.matrix(log2(tFpkm[sigTrans$tx,pIndex]+1))
outPost = matrix(nrow = nrow(sigTrans), nc = 3)

for(j in 1:nrow(sigTrans)) {
	if(j %% 5000 == 0) cat(".")
	outPost[j,] = summary(lm(tFpkmPost[j,]~snpPost[j,] + 
		cbind(modPost[,-1], pcaPost$x[,1:nsvPost])))$coef[2,c(1,3,4)]
}
colnames(outPost) = paste0("postnatal_", c("slope","tstat","pval"))
sigTrans = cbind(sigTrans, outPost)
sigTrans$postnatal_MAF = rowSums(snpPost, na.rm=TRUE)/(2*rowSums(!is.na(snpPost)))


save(sigTrans, file="rdas/annotated_transcript_eqtl_szControl_cisOnly.rda")
