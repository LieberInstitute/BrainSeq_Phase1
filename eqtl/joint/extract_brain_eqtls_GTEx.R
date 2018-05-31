###
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

###############################
#### snps from GTEX #########

# #### MAF filter
# bfile = "/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Genotypes/Merged/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed"
# system(paste("plink --bfile", bfile, "--maf 0.05 --hwe 0.000001 --geno 0.1 --make-bed",
	# "--out /dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Genotypes/Merged/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf05_geno10_hwe1e6"))

# # get MDS components
# bfile = "/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Genotypes/Merged/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf05_geno10_hwe1e6"
# system(paste("plink --bfile", bfile,
	# " --indep 100 10 1.25 --geno 0.1 --out", bfile))

# ## MDS
# system(paste0("plink --bfile ", bfile, 
	# " --cluster --mds-plot 10 --extract ",
	# bfile, ".prune.in --out ", bfile))

### lighter MAF filter
# bfile = "/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Genotypes/Merged/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed"
# system(paste("plink --bfile", bfile, "--maf 0.005 --hwe 0.000001 --geno 0.1 --make-bed",
	# "--out /dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Genotypes/Merged/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6"))
snpMapGtex = read.delim("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Genotypes/Merged/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6.bim",
	as.is=TRUE, header=FALSE)
colnames(snpMapGtex)[-3] = colnames(snpMap)[1:5]
snpMapGtex$CHR[snpMapGtex$CHR == "23"] = "X"
snpMapGtex$chrpos = paste0("chr", snpMapGtex$CHR, ":", snpMapGtex$POS)

## match up to eqtls
geneEqtl$inGTEX = match(geneEqtl$snp_chrpos, snpMapGtex$chrpos)
geneEqtl2 = geneEqtl[!is.na(geneEqtl$inGTEX),]
geneEqtl2$snpCounted = snpMap$COUNTED[
	match(as.character(geneEqtl2$snps), snpMap$SNP)]

exonEqtl$inGTEX = match(exonEqtl$snp_chrpos, snpMapGtex$chrpos)
exonEqtl2 = exonEqtl[!is.na(exonEqtl$inGTEX),]
exonEqtl2$snpCounted = snpMap$COUNTED[
	match(as.character(exonEqtl2$snps), snpMap$SNP)]

jxnEqtl$inGTEX = match(jxnEqtl$snp_chrpos, snpMapGtex$chrpos)
jxnEqtl2 = jxnEqtl[!is.na(jxnEqtl$inGTEX),]
jxnEqtl2$snpCounted = snpMap$COUNTED[
	match(as.character(jxnEqtl2$snps), snpMap$SNP)]

derEqtl$inGTEX = match(derEqtl$snp_chrpos, snpMapGtex$chrpos)
derEqtl2 = derEqtl[!is.na(derEqtl$inGTEX),]
derEqtl2$snpCounted = snpMap$COUNTED[
	match(as.character(derEqtl2$snps), snpMap$SNP)]

transEqtl$inGTEX = match(transEqtl$snp_chrpos, snpMapGtex$chrpos)
transEqtl2 = transEqtl[!is.na(transEqtl$inGTEX),]
transEqtl2$snpCounted = snpMap$COUNTED[
	match(as.character(transEqtl2$snps), snpMap$SNP)]

snpIndex=unique(c(geneEqtl2$inGTEX, exonEqtl2$inGTEX, 
	jxnEqtl2$inGTEX, derEqtl2$inGTEX, transEqtl2$inGTEX))
variantsToPull = snpMapGtex$SNP[snpIndex]
cat(variantsToPull, file="snps_to_extract.txt", sep="\n")

## run plink
system("plink --bfile /dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Genotypes/Merged/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf005_geno10_hwe1e6 --extract snps_to_extract.txt --recode A-transpose --out GTEX_SNPs_LIBDeqtls")

##### read in SNPs
genotypes = read.delim("GTEX_SNPs_LIBDeqtls.traw",as.is=TRUE)
snpMapSub = genotypes[,1:6]
snpMapSub$CHR[snpMapSub$CHR == "23"] = "X"
snpMapSub$chrpos = paste0("chr", snpMapSub$CHR, ":", snpMapSub$POS)
snpSub = genotypes[,-(1:6)]
colnames(snpSub) = ss(colnames(snpSub), "\\.",2)
rownames(snpMapSub) = rownames(snpSub) = snpMapSub$SNP

### add MDS
mds = read.table("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Genotypes/Merged/GTEX_Brain_Illumina_Omni5M_Omni2pt5M_imputed_maf05_geno10_hwe1e6.mds",
	header=TRUE, as.is=TRUE)

########	
##### load expression data
xx=load("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Counts/rpkmCounts_brainGTEX.rda")
pdBrain$snpColumn = match(ss(pdBrain$SUBJID, "-", 2), colnames(snpSub))
pdBrainSub = pdBrain[!is.na(pdBrain$snpColumn),]

# ##################################
# ## subset to eqtl features #######
geneRpkmGeno = geneRpkm[,!is.na(pdBrain$snpColumn)]
geneRpkmSub = geneRpkmGeno[unique(geneEqtl2$gene),]
geneMapSub = geneMap[unique(geneEqtl2$gene),]

exonRpkmGeno = exonRpkm[,!is.na(pdBrain$snpColumn)]
exonRpkmSub = exonRpkmGeno[unique(exonEqtl2$gene),]
exonMapSub = exonMap[unique(exonEqtl2$gene),]

rownames(countsM)=names(jMap)
jRpkmLogical = DataFrame(sapply(countsM, function(x) x > 0))
jRpkmLogical = as.data.frame(jRpkmLogical)
jxnIndex = which(rowSums(jRpkmLogical) >= 5) 
countsM = countsM[jxnIndex,]
jMap = jMap[jxnIndex]

mmJxn = match(unique(jxnEqtl2$gene), rownames(countsM))

jRpkmGeno = as.matrix(as.data.frame(countsM[,!is.na(pdBrain$snpColumn)]))
jRpkmSub = jRpkmGeno[mmJxn[!is.na(mmJxn)],]
jMapSub = jMap[mmJxn[!is.na(mmJxn)]]

## ERs
xx=load("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/coverage/GTEX_coverageMat_szControlEqtlDERs.rda")
covMatGeno = covMatGtex[,pdBrainSub$sra_accession]
covMatSub = covMatGeno[unique(derEqtl2$gene),]
covMapSub = regionsGtex[unique(derEqtl2$gene)]

## Txs
xx = load("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/filteredTxData_GTEX_Brain_szControl.rda")
colnames(tFpkmGtex) = ss(colnames(tFpkmGtex),"\\.",2)
tFpkmGeno = tFpkmGtex[,pdBrainSub$sra_accession]
tFpkmSub = tFpkmGeno[unique(transEqtl2$gene),]
tMapSub = tMap[unique(transEqtl2$gene)]


###############
## do PCA #####
pcaGene = prcomp(t(log2(geneRpkmGeno+1)))

ooExon = order(rowSds(log2(exonRpkmGeno + 1)),decreasing=TRUE)[1:100000]
pcaExon = prcomp(t(log2(exonRpkmGeno[ooExon,]+1)))

ooJxn = order(rowSds(jRpkmGeno),decreasing=TRUE)[1:100000]
pcaJxn = prcomp(t(log2(jRpkmGeno[ooJxn,]+1)))

ooDer = order(rowSds(log2(as.matrix(covMatGeno) + 1)),decreasing=TRUE)[1:100000]
pcaDer = prcomp(t(log2(as.matrix(covMatGeno[ooDer,]) + 1)))

ooTx = order(rowSds(log2(as.matrix(tFpkmGeno) + 1)),decreasing=TRUE)[1:100000]
pcaTx = prcomp(t(log2(as.matrix(tFpkmGeno[ooTx,]) + 1)))

pcList = list(gene = pcaGene$x[,1:10], exon = pcaExon$x[,1:10],
	junction = pcaJxn$x[,1:10], transcript = pcaTx$x[,1:10],
	der = pcaDer$x[,1:10])

save(pdBrainSub, mds, snpSub, snpMapSub, geneMapSub, geneRpkmSub, 
	exonRpkmSub, exonMapSub, jRpkmSub, jMapSub, 
	covMatSub, covMapSub, tFpkmSub, tMapSub, pcList, 
	geneEqtl2, exonEqtl2, jxnEqtl2,derEqtl2,transEqtl2,
	file="/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/gtex_brainEqtl_subsets.rda")
