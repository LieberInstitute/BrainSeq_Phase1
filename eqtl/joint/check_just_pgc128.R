###
###
source("/users/ajaffe/Lieber/lieber_functions_aj.R")

library(MatrixEQTL)
library(GenomicRanges)
library(sva)

#### load data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/snpData_LIBD_szControls_cleaned.rda")
pd$DxNum = as.numeric(factor(pd$Dx))-1

# filter for Age
aIndex= which(pd$Age > 13)
pd2= pd[aIndex,]
snp2 = as.matrix(snp[,aIndex])
geneRpkm2 = as.matrix(log2(geneRpkm[,aIndex]+1))
mod = model.matrix(~snpPC1 + snpPC2 + snpPC3 + Sex,data=pd2)

### filter rpkm
expIndex=which(rowMeans(geneRpkm2) > 0.01 & 
	geneMap$Chr %in% paste0("chr", c(1:22, "X","Y")))
geneRpkm2 = geneRpkm2[expIndex,]
geneMap = geneMap[expIndex,]

# for PCs
load("rdas/gene_eqtl_szControl_cisOnly.rda")

# model
covsGene = SlicedData$new(t(cbind(mod[,-1], pcsGene, pd2$DxNum)))
exprsGene = SlicedData$new(geneRpkm2)
exprsGene$ResliceCombined(sliceSize = 1000)


## drop to pgc128
pgcSig = read.delim("/users/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/pgc/pgc2_128loci.txt")
load("/users/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/rdas/granges_pgc2.rda")
pgcSigGR = pgcFullGR[pgcSig$Index_SNP]
pgcSig$chrpos = paste0(seqnames(pgcSigGR), ":", start(pgcSigGR))

## filter
snpMap$chrpos = paste0("chr", snpMap$CHR, ":", snpMap$POS)
sIndex= which(snpMap$chrpos %in% pgcSig$chrpos)
snpMapSub = snpMap[sIndex,]
snpSub = snp2[sIndex,]

theSnps = SlicedData$new(as.matrix(snpSub))

# ## JUST CIS
snpspos = snpMapSub[,c("SNP","CHR","POS")]
snpspos$CHR = paste0("chr",snpspos$CHR)
colnames(snpspos) = c("name","chr","pos")

posGene = geneMap[,1:3]
posGene$name = rownames(geneMap)
posGene = posGene[,c(4,1:3)]

meGeneCis = Matrix_eQTL_main(snps=theSnps, gene = exprsGene, 
	cvrt = covsGene, output_file_name.cis =  ".txt" ,
	pvOutputThreshold.cis = 1, pvOutputThreshold=0,
	snpspos = snpspos, genepos = posGene, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)

eqtl = meGeneCis$cis$eqtl
eqtl$gene = as.character(eqtl$gene)
eqtl$snps = as.character(eqtl$snps)

# snp coordinates
m = match(eqtl$snps, snpMap$SNP)
eqtl$snp_chr = snpMap$CHR[m]
eqtl$snp_chr = paste0("chr",eqtl$snp_chr)
eqtl$snp_chr = gsub("chr23","chrX", eqtl$snp_chr)
eqtl$snp_pos = snpMap$POS[m]
eqtl$snpRsNum = snpMap$name[m]

eqtl$snpCounted = snpMap$COUNTED[m]
eqtl$snpAlt = snpMap$ALT[m]
eqtl$inSampleMAF = rowSums(snp2[m,],na.rm=TRUE)/
	(2*rowSums(!is.na(snp2[m,])))

# annotate gene
g = geneMap[match(eqtl$gene, rownames(geneMap)),]
colnames(g) = tolower(colnames(g))
colnames(g) = paste0("exprs_", colnames(g))
eqtl = cbind(eqtl,g)
rownames(eqtl) = NULL

## cis or trans
eqtl$cisOrTrans = ifelse(eqtl$snp_chr == eqtl$exprs_chr,"cis","trans")
start = eqtl$exprs_start
start[which(eqtl$exprs_strand == "-")] = eqtl$exprs_end[
	which(eqtl$exprs_strand == "-")]
eqtl$cisDistToTss  = eqtl$snp_pos - start
eqtl$cisDistToTss[eqtl$cisOrTrans=="trans"]=NA
eqtl$meanRPKM = 2^rowMeans(geneRpkm2)[eqtl$gene]-1
save(eqtl, file = "rdas/gene_eqtl_szControl_cisOnly_pgcOnly_annotated.rda")