######
library(GenomicRanges)
library(matrixStats)
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

## load SNP data for map
load("/dcs01/ajaffe/Brain/DLPFC_PolyA/controlEqtl/data/snpData_LIBD_controls_cleaned.rda")
snpMap$CHR[snpMap$CHR == "23"] = "X"
snpMap$chrpos = paste0("chr", snpMap$CHR, ":", snpMap$POS)

## load eqtls
xx=load("rdas/annotated_der_eqtl_szControl_cisOnly.rda")
xx=load("rdas/annotated_exon_eqtl_szControl_cisOnly.rda")
xx=load("rdas/annotated_gene_eqtl_szControl_cisOnly.rda")
xx=load("rdas/annotated_junction_eqtl_szControl_cisOnly.rda")
xx=load("rdas/annotated_transcript_eqtl_szControl_cisOnly.rda")

sigGene$chrpos = paste0(sigGene$snp_chr, ":", sigGene$snp_pos)
sigExon$chrpos = paste0(sigExon$snp_chr, ":", sigExon$snp_pos)
sigJxn$chrpos = paste0(sigJxn$snp_chr, ":", sigJxn$snp_pos)
sigDer$chrpos = paste0(sigDer$snp_chr, ":", sigDer$snp_pos)
sigTrans$chrpos = paste0(sigTrans$snp_chr, ":", sigTrans$snp_pos)

# raw
load("rdas/gene_eqtl_szControl_cisOnly.rda")
geneEqtl = meGeneCis$cis$eqtls
geneEqtl = geneEqtl[geneEqtl$FDR < 0.01,]
geneEqtl$gene = as.character(geneEqtl$gene)
geneEqtl$snps = as.character(geneEqtl$snps)
geneEqtl$snp_chrpos = snpMap$chrpos[match(as.character(geneEqtl$snps), snpMap$SNP)]
geneEqtl$snp_rsNumber = snpMap$name[match(as.character(geneEqtl$snps), snpMap$SNP)]
geneEqtl$rank = match(as.character(geneEqtl$gene), sigGene$gene)

load("rdas/exon_eqtl_szControl_cisOnly.rda")
exonEqtl = meExonCis$cis$eqtls
exonEqtl = exonEqtl[exonEqtl$FDR < 0.01,]
exonEqtl$gene = as.character(exonEqtl$gene)
exonEqtl$snps = as.character(exonEqtl$snps)
exonEqtl$snp_chrpos = snpMap$chrpos[match(as.character(exonEqtl$snps), snpMap$SNP)]
exonEqtl$snp_rsNumber = snpMap$name[match(as.character(exonEqtl$snps), snpMap$SNP)]
exonEqtl$rank = match(as.character(exonEqtl$gene), sigExon$exon)

load("rdas/junction_eqtl_szControl_cisOnly.rda")
jxnEqtl = meJxnCis$cis$eqtls
jxnEqtl = jxnEqtl[jxnEqtl$FDR < 0.01,]
jxnEqtl$gene = as.character(jxnEqtl$gene)
jxnEqtl$snps = as.character(jxnEqtl$snps)
jxnEqtl$snp_chrpos = snpMap$chrpos[match(as.character(jxnEqtl$snps), snpMap$SNP)]
jxnEqtl$snp_rsNumber = snpMap$name[match(as.character(jxnEqtl$snps), snpMap$SNP)]
jxnEqtl$rank = match(as.character(jxnEqtl$gene), sigJxn$jxn)

load("rdas/der_eqtl_szControl_cisOnly.rda")
derEqtl = meDerCis$cis$eqtls
derEqtl = derEqtl[derEqtl$FDR < 0.01,]
derEqtl$gene = as.character(derEqtl$gene)
derEqtl$snps = as.character(derEqtl$snps)
derEqtl$snp_chrpos = snpMap$chrpos[match(as.character(derEqtl$snps), snpMap$SNP)]
derEqtl$snp_rsNumber = snpMap$name[match(as.character(derEqtl$snps), snpMap$SNP)]
derEqtl$rank = match(as.character(derEqtl$gene), sigDer$der)

load("rdas/transcript_eqtl_szControl_cisOnly.rda")
transEqtl = meTransCis$cis$eqtls
transEqtl = transEqtl[transEqtl$FDR < 0.01,]
transEqtl$gene = as.character(transEqtl$gene)
transEqtl$snps = as.character(transEqtl$snps)
transEqtl$snp_chrpos = snpMap$chrpos[match(as.character(transEqtl$snps), snpMap$SNP)]
transEqtl$snp_rsNumber = snpMap$name[match(as.character(transEqtl$snps), snpMap$SNP)]
transEqtl$rank = match(as.character(transEqtl$gene), sigTrans$tx)

####

###############################
#### snps from CMC #########

snpMapCmc = read.delim("/dcl01/lieber/ajaffe/PublicData/CMC/Genotypes/CMC_genotypes_imputed_common_cleaned.bim",
	as.is=TRUE, header=FALSE)
colnames(snpMapCmc)[-3] = colnames(snpMap)[1:5]
snpMapCmc$CHR[snpMapCmc$CHR == "23"] = "X"
snpMapCmc$chrpos = paste0("chr", snpMapCmc$CHR, ":", snpMapCmc$POS)

## match up to eqtls
geneEqtl$inCMC = match(geneEqtl$snp_chrpos, snpMapCmc$chrpos)

mean(!is.na(geneEqtl$inCMC))
mean(!is.na(geneEqtl$inCMC[!duplicated(geneEqtl$rank)]))

geneEqtl2 = geneEqtl[!is.na(geneEqtl$inCMC),]
geneEqtl2$snpCounted = snpMap$COUNTED[
	match(as.character(geneEqtl2$snps), snpMap$SNP)]

exonEqtl$inCMC = match(exonEqtl$snp_chrpos, snpMapCmc$chrpos)
exonEqtl2 = exonEqtl[!is.na(exonEqtl$inCMC),]
exonEqtl2$snpCounted = snpMap$COUNTED[
	match(as.character(exonEqtl2$snps), snpMap$SNP)]

jxnEqtl$inCMC = match(jxnEqtl$snp_chrpos, snpMapCmc$chrpos)
jxnEqtl2 = jxnEqtl[!is.na(jxnEqtl$inCMC),]
jxnEqtl2$snpCounted = snpMap$COUNTED[
	match(as.character(jxnEqtl2$snps), snpMap$SNP)]

derEqtl$inCMC = match(derEqtl$snp_chrpos, snpMapCmc$chrpos)
derEqtl2 = derEqtl[!is.na(derEqtl$inCMC),]
derEqtl2$snpCounted = snpMap$COUNTED[
	match(as.character(derEqtl2$snps), snpMap$SNP)]
	
transEqtl$inCMC = match(transEqtl$snp_chrpos, snpMapCmc$chrpos)
transEqtl2 = transEqtl[!is.na(transEqtl$inCMC),]
transEqtl2$snpCounted = snpMap$COUNTED[
	match(as.character(transEqtl2$snps), snpMap$SNP)]

### cat out
snpIndex=unique(c(geneEqtl2$inCMC, exonEqtl2$inCMC,
	jxnEqtl2$inCMC, derEqtl2$inCMC, transEqtl2$inCMC))
variantsToPull = snpMapCmc$SNP[snpIndex]
cat(variantsToPull, file="snps_to_extract_CMC.txt", sep="\n")

# run plink
system("plink --bfile /dcl01/lieber/ajaffe/PublicData/CMC/Genotypes/CMC_genotypes_imputed_common_cleaned --extract snps_to_extract_CMC.txt --recode A-transpose --out CMC_SNPs_LIBDeqtls")

##### read in SNPs
genotypes = read.delim("CMC_SNPs_LIBDeqtls.traw",as.is=TRUE)
snpMapSub = genotypes[,1:6]
snpMapSub$CHR[snpMapSub$CHR == "23"] = "X"
snpMapSub$chrpos = paste0("chr", snpMapSub$CHR, ":", snpMapSub$POS)
snpSub = genotypes[,-(1:6)]
rownames(snpMapSub) = rownames(snpSub) = snpMapSub$SNP

### add MDS
mds = read.table("/dcl01/lieber/ajaffe/PublicData/CMC/Genotypes/CMC_genotypes_imputed_common_cleaned.mds",
	header=TRUE, as.is=TRUE)

########	
##### load expression data
load("/dcl01/lieber/ajaffe/PublicData/CMC/CMC_coverageMat_szControlEqtlDERs.rda")
load("/dcl01/lieber/ajaffe/PublicData/CMC/filteredTxData_CMC_szControl.rda")
load("/dcl01/lieber/ajaffe/PublicData/CMC/rpkmCounts_cmcDlpfc_szControl.rda")

## just controls
pd$snpColumn = match(paste0(pd$Genotyping_Sample_ID, "_", 
	pd$Genotyping_Sample_ID), colnames(snpSub))

pdSub = pd
mds = mds[match(pdSub$Genotyping_Sample_ID, mds$FID),]
snpSub = snpSub[,match(paste0(pd$Genotyping_Sample_ID, "_", 
	pd$Genotyping_Sample_ID), colnames(snpSub))]
colnames(snpSub) = pdSub$Individual_ID

# ##################################
# ## subset to eqtl features #######
geneRpkmSub = geneRpkm[unique(geneEqtl2$gene),]
geneMapSub = geneMap[unique(geneEqtl2$gene),]

exonRpkmSub = exonRpkm[unique(exonEqtl2$gene),]
exonMapSub = exonMap[unique(exonEqtl2$gene),]

mmJxn = match(unique(jxnEqtl2$gene), rownames(jRpkm))
jRpkmSub = jRpkm[mmJxn[!is.na(mmJxn)],]
jMapSub = jMap[mmJxn[!is.na(mmJxn)]]

regionMatSub = regionMatCmc[unique(derEqtl2$gene),]
covMapSub = regionsCmc[unique(derEqtl2$gene)]

tFpkmSub = tFpkmCmc[unique(transEqtl2$gene),]
tMapSub = tMap[unique(transEqtl2$gene)]

###############
## do PCA #####
ooGene = order(rowSds(log2(geneRpkm + 1)),decreasing=TRUE)[1:15000]
pcaGene = prcomp(t(log2(geneRpkm[ooGene,]+1)))

ooExon = order(rowSds(log2(exonRpkm + 1)),decreasing=TRUE)[1:50000]
pcaExon = prcomp(t(log2(exonRpkm[ooExon,]+1)))

ooJxn = order(rowSds(log2(as.matrix(jRpkm) + 1)),decreasing=TRUE)[1:50000]
pcaJxn = prcomp(t(log2(jRpkm[ooJxn,]+1)))

ooDer = order(rowSds(log2(as.matrix(regionMatCmc) + 1)),decreasing=TRUE)[1:50000]
pcaDer = prcomp(t(log2(as.matrix(regionMatCmc[ooDer,]) + 1)))

ooTrans = order(rowSds(log2(as.matrix(tFpkmCmc) + 1)),decreasing=TRUE)[1:50000]
pcaTrans = prcomp(t(log2(as.matrix(tFpkmCmc[ooTrans,]) + 1)))

pcList = list(gene = pcaGene$x[,1:10], exon = pcaExon$x[,1:10],
	transcript = pcaTrans$x[,1:10],
	junction = pcaJxn$x[,1:10], der = pcaDer$x[,1:10])

save(pdSub, mds, snpSub, snpMapSub, geneMapSub, geneRpkmSub, 
	exonRpkmSub, exonMapSub, jRpkmSub, jMapSub, 
	tFpkmSub, tMapSub, 	pcList,  regionMatSub, covMapSub, 
	geneEqtl2, exonEqtl2, jxnEqtl2,transEqtl2,derEqtl2,
	file="/dcl01/lieber/ajaffe/PublicData/CMC/CMC_brainEqtl_subsets.rda",compress=TRUE)
