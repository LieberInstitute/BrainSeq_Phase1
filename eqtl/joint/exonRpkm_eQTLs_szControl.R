##
# qsub -t 1-23 -cwd -l mem_free=50G,h_vmem=55G,h_stack=256M -b y R CMD BATCH --no-save exonRpkm_eQTLs_adults.R
source("/users/ajaffe/Lieber/lieber_functions_aj.R")

library(sva)
library(MatrixEQTL)
library(GenomicRanges)

#### load data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/snpData_LIBD_szControls_cleaned.rda")
pd$DxNum = as.numeric(factor(pd$Dx))-1

# filter for Age
aIndex= which(pd$Age > 13)
pd2= pd[aIndex,]
snp2 = as.matrix(snp[,aIndex])
exonRpkm2 = as.matrix(log2(exonRpkm[,aIndex]+1))

mod = model.matrix(~snpPC1 + snpPC2 + snpPC3 + Sex,data=pd2)

### filter rpkm
expIndex=which(rowMeans(exonRpkm2) > 0.1 & 
	exonMap$Chr %in% paste0("chr", c(1:22,"X","Y")))
exonRpkm2 = exonRpkm2[expIndex,]
exonMap = exonMap[expIndex,]
	
# #########
# #### exon
# pcaExon = prcomp(t(exonRpkm2))
# nsvExon = num.sv(exonRpkm2, mod, vfilter=1e5) # loaded below
# pcsExon = pcaExon$x[,1:nsvExon]

# covsExon = SlicedData$new(t(cbind(mod[,-1],pcsExon,pd$DxNum)))
# exprsExon = SlicedData$new(as.matrix(exonRpkm2))
# exprsExon$ResliceCombined(sliceSize = 5000)

# ### just CIS
# theSnps = SlicedData$new(snp2)
# theSnps$ResliceCombined(sliceSize = 50000)
# snpspos = snpMap[,c("SNP","CHR","POS")]
# snpspos$CHR = paste0("chr",snpspos$CHR)
# colnames(snpspos) = c("name","chr","pos")

# posExon = exonMap[,2:4]
# posExon$name = rownames(exonMap)
# posExon = posExon[,c(4,1:3)]

# meExonCis = Matrix_eQTL_main(snps=theSnps, gene = exprsExon, 
	# cvrt = covsExon, output_file_name.cis =  ".txt" ,
	# pvOutputThreshold.cis = 0.0001,  pvOutputThreshold=0,
	# snpspos = snpspos, genepos = posExon, 
	# useModel = modelLINEAR,	cisDist=5e5,
	# pvalue.hist = 100,min.pv.by.genesnp = TRUE)

# save(meExonCis, nsvExon, pcsExon, file="rdas/exon_eqtl_szControl_cisOnly.rda")

#### load back in
load("rdas/exon_eqtl_szControl_cisOnly.rda")
eqtl = meExonCis$cis$eqtl

## number of effective tests
snpMapGR = GRanges(paste0("chr", snpMap$CHR), IRanges(snpMap$POS, width=1))
exonMapGR = makeGRangesFromDataFrame(exonMap, keep=TRUE)
pruned = read.table("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/plink/LIBD_Brain_DLPFC_szControls_imputed_indep.prune.in")
snpMap$isLdIndep = snpMap$SNP %in% pruned$V1

nTests = length(findOverlaps(snpMapGR[snpMap$isLdIndep], 
	exonMapGR, maxgap=5e5))
eqtl$bonf = p.adjust(eqtl$pvalue, "bonferroni", n= nTests)
length(unique(eqtl$gene[eqtl$bonf < 0.05])) # 27553

eqtl = eqtl[eqtl$FDR < 0.01,] # FDR < 1%

# annotate
eqtl$gene = as.character(eqtl$gene)
eqtl$snps = as.character(eqtl$snps)
colnames(eqtl)[2] = "exon"
eqtl = eqtl[order(eqtl$pvalue),]

# snp coordinates
m = match(eqtl$snps, snpMap$SNP)
eqtl$snp_chr = snpMap$CHR[m]
eqtl$snp_chr = paste0("chr",eqtl$snp_chr)
eqtl$snp_chr = gsub("chr23","chrX", eqtl$snp_chr)
eqtl$snp_pos = snpMap$POS[m]
eqtl$snpRsNum = snpMap$name[m]

# take best eqtl per exon
sig = eqtl[!duplicated(eqtl$exon),]
tt = table(eqtl$exon)
sig$numSnps = tt[sig$exon]

# other snp info
snprange = t(sapply(split(eqtl$snp_pos, eqtl$exon), range))
snprange= snprange[as.character(sig$exon),]
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
g = exonMap[match(sig$exon, rownames(exonMap)),]
colnames(g) = tolower(colnames(g))
colnames(g) = paste0("exprs_", colnames(g))
sig = cbind(sig,g)
rownames(sig) = NULL

## cis or trans
tss = ifelse(sig$exprs_strand =="+",
	sig$exprs_start, sig$exprs_end)
sig$cisOrTrans = ifelse(sig$snp_chr == sig$exprs_chr & 
	abs(sig$snp_pos - tss) < 1e6,	"cis","trans")
sig$distToExon  = sig$snp_pos - tss
sig$distToExon[sig$cisOrTrans=="trans"] = NA

sig$meanRPKM = 2^rowMeans(exonRpkm2)[sig$exon]-1

### add race sensitivity analyses
rIndexes=splitit(pd2$Race)
raceExon = mclapply(rIndexes, function(ii) {
	ssnp = as.matrix(snp2[sig$snps,ii])
	yyExon = as.matrix(exonRpkm2[sig$exon,ii])
	out = matrix(nrow = nrow(ssnp), nc = 3)
	for(j in 1:nrow(ssnp)) {
		if(j %% 5000 == 0) cat(".")
		out[j,] = summary(lm(yyExon[j,]~ssnp[j,] + cbind(
			mod[ii,-1], pcsExon[ii,])))$coef[2,c(1,3,4)]
	}
	colnames(out) = c("slope","tstat","pval")
	cat("\n")
	return(out)
},mc.cores=2)
xx = do.call("cbind", raceExon)
colnames(xx) = paste0(rep(names(rIndexes), each=3), "_", colnames(xx))
sigExon = cbind(sig, xx)

# race MAFs
mafRace = sapply(rIndexes, function(ii){
	ssnp = as.matrix(snp2[sig$snps,ii])
	rowSums(ssnp, na.rm=TRUE)/(2*rowSums(!is.na(ssnp)))
})
colnames(mafRace) = paste0(colnames(mafRace),"_inSampleMAF")
sigExon = cbind(sigExon, mafRace)

#### add fetal effect
fIndex=which(pd$Age < 0)
pdFetal= pd[fIndex,]
pcaFetal = prcomp(t(log2(exonRpkm[expIndex,fIndex]+1)))
modFetal = model.matrix(~ pdFetal$snpPC1+
	pdFetal$snpPC2 + pdFetal$snpPC3 + pdFetal$Sex)
nsvFetal = num.sv(log2(exonRpkm[expIndex,fIndex]+1), modFetal)

snpFetal = as.matrix(snp[match(sigExon$snps,rownames(snp)),fIndex])
exonRpkmFetal = as.matrix(log2(exonRpkm[sigExon$exon,fIndex]+1))
outFetal = matrix(nrow = nrow(sigExon), nc = 3)

for(j in 1:nrow(sigExon)) {
	if(j %% 5000 == 0) cat(".")
	outFetal[j,] = summary(lm(exonRpkmFetal[j,]~snpFetal[j,] + 
		cbind(modFetal[,-1], pcaFetal$x[,1:nsvFetal])))$coef[2,c(1,3,4)]
}
colnames(outFetal) = paste0("fetal_", c("slope","tstat","pval"))
sigExon = cbind(sigExon, outFetal)
sigExon$fetal_MAF = rowSums(snpFetal, na.rm=TRUE)/(2*rowSums(!is.na(snpFetal)))

save(sigExon, file="rdas/annotated_exon_eqtl_szControl_cisOnly.rda")
