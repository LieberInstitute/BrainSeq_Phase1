###
source("/users/ajaffe/Lieber/lieber_functions_aj.R")

library(MatrixEQTL)
library(GenomicRanges)
library(sva)

#### load data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/snpData_LIBD_szControls_cleaned.rda")

# filter for Age
aIndex= which(pd$Age > 13 & pd$Dx == "Control")
pd2= pd[aIndex,]
snp2 = as.matrix(snp[,aIndex])
geneRpkm2 = as.matrix(log2(geneRpkm[,aIndex]+1))
mod = model.matrix(~snpPC1 + snpPC2 + snpPC3 + Sex,data=pd2)

### filter rpkm
expIndex=which(rowMeans(geneRpkm2) > 0.01 & 
	geneMap$Chr %in% paste0("chr", c(1:22, "X","Y")))
geneRpkm2 = geneRpkm2[expIndex,]
geneMap = geneMap[expIndex,]
	
# # ####################
# # ## do cis eQTLS, maybe 500kb
# pcaGene = prcomp(t(geneRpkm2))
# nsvGene = num.sv(geneRpkm2, mod)
# pcsGene = pcaGene$x[,1:nsvGene]

# covsGene = SlicedData$new(t(cbind(mod[,-1], pcsGene)))
# exprsGene = SlicedData$new(geneRpkm2)
# exprsGene$ResliceCombined(sliceSize = 1000)

# theSnps = SlicedData$new(snp2)
# theSnps$ResliceCombined(sliceSize = 50000)

# # ## JUST CIS
# snpspos = snpMap[,c("SNP","CHR","POS")]
# snpspos$CHR = paste0("chr",snpspos$CHR)
# colnames(snpspos) = c("name","chr","pos")

# posGene = geneMap[,1:3]
# posGene$name = rownames(geneMap)
# posGene = posGene[,c(4,1:3)]

# meGeneCis = Matrix_eQTL_main(snps=theSnps, gene = exprsGene, 
	# cvrt = covsGene, output_file_name.cis =  ".txt" ,
	# pvOutputThreshold.cis = 0.001, pvOutputThreshold=0,
	# snpspos = snpspos, genepos = posGene, 
	# useModel = modelLINEAR,	cisDist=5e5,
	# pvalue.hist = 100,min.pv.by.genesnp = TRUE)

# save(meGeneCis, nsvGene, pcsGene, file="rdas/gene_eqtl_control_13plus_cisOnly.rda")

# #### annotate
load("rdas/gene_eqtl_control_13plus_cisOnly.rda")
eqtl = meGeneCis$cis$eqtl
eqtl$gene = as.character(eqtl$gene)
eqtl$snps = as.character(eqtl$snps)

### number of effective tests
snpMapGR = GRanges(paste0("chr", snpMap$CHR), IRanges(snpMap$POS, width=1))
geneMapGR = makeGRangesFromDataFrame(geneMap, keep=TRUE)
pruned = read.table("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/plink/LIBD_Brain_DLPFC_szControls_imputed_indep.prune.in")
snpMap$isLdIndep = snpMap$SNP %in% pruned$V1

nTests = length(findOverlaps(snpMapGR[snpMap$isLdIndep], 
	geneMapGR, maxgap=5e5))
eqtl$bonf = p.adjust(eqtl$pvalue, "bonferroni", n= nTests)
length(unique(eqtl$gene[eqtl$bonf < 0.05])) # 4301

# snp coordinates
m = match(eqtl$snps, snpMap$SNP)
eqtl$snp_chr = snpMap$CHR[m]
eqtl$snp_chr = paste0("chr",eqtl$snp_chr)
eqtl$snp_chr = gsub("chr23","chrX", eqtl$snp_chr)
eqtl$snp_pos = snpMap$POS[m]
eqtl$snpRsNum = snpMap$name[m]

# take best eqtl per gene
eqtl = eqtl[eqtl$FDR < 0.01,]
sig = eqtl[!duplicated(eqtl$gene),]
tt = table(eqtl$gene)
sig$numSnps = tt[sig$gene]

# other snp info
snprange = t(sapply(split(eqtl$snp_pos, eqtl$gene), range))
snprange= snprange[sig$gene,]
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
g = geneMap[match(sig$gene, rownames(geneMap)),]
colnames(g) = tolower(colnames(g))
colnames(g) = paste0("exprs_", colnames(g))
sig = cbind(sig,g)
rownames(sig) = NULL

## cis or trans
sig$cisOrTrans = ifelse(sig$snp_chr == sig$exprs_chr,"cis","trans")
start = sig$exprs_start
start[which(sig$exprs_strand == "-")] = sig$exprs_end[
	which(sig$exprs_strand == "-")]
sig$cisDistToTss  = sig$snp_pos - start
sig$cisDistToTss[sig$cisOrTrans=="trans"]=NA
sig$meanRPKM = 2^rowMeans(geneRpkm2)[sig$gene]-1

### add race sensitivity analyses
rIndexes=splitit(pd2$Race)
raceGene = mclapply(rIndexes, function(ii) {
	ssnp = as.matrix(snp2[sig$snps,ii])
	yyGene = as.matrix(geneRpkm2[sig$gene,ii])
	out = matrix(nrow = nrow(ssnp), nc = 3)
	for(j in 1:nrow(ssnp)) {
		if(j %% 5000 == 0) cat(".")
		out[j,] = summary(lm(yyGene[j,]~ssnp[j,] + cbind(
			mod[ii,-1], pcsGene[ii,])))$coef[2,c(1,3,4)]
	}
	colnames(out) = c("slope","tstat","pval")
	cat("\n")
	return(out)
},mc.cores=2)
xx = do.call("cbind", raceGene)
colnames(xx) = paste0(rep(names(rIndexes), each=3), "_", colnames(xx))
sigGene = cbind(sig, xx)

mafRace = sapply(rIndexes, function(ii){
	ssnp = as.matrix(snp2[sig$snps,ii])
	rowSums(ssnp, na.rm=TRUE)/(2*rowSums(!is.na(ssnp)))
})
colnames(mafRace) = paste0(colnames(mafRace),"_inSampleMAF")
sigGene = cbind(sigGene, mafRace)

save(sigGene, file="rdas/annotated_gene_eqtl_control_13plus_cisOnly.rda")
