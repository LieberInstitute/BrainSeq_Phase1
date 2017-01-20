###
# qsub -pe local 2 -cwd -l mem_free=70G,h_vmem=80G,h_stack=256M -b y R CMD BATCH --no-save junctionRpm_eQTLs_adults.R

source("/users/ajaffe/Lieber/lieber_functions_aj.R")

library(MatrixEQTL)
library(GenomicRanges)
library(sva)

#### load data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/snpData_LIBD_szControls_cleaned.rda")

# filter for age
aIndex= which(pd$Age > 13 & pd$Dx == "Control")
pd2= pd[aIndex,]
snp2 = as.matrix(snp[,aIndex] )
jRpkm2 = as.matrix(log2(jRpkm[,aIndex]+1))

mod = model.matrix(~snpPC1 + snpPC2 + snpPC3 + Sex,data=pd2)

### filter rp80m
expIndex=which(rowMeans(jRpkm2) > 0.2 & jMap$code != "Novel")
jRpkm2 = jRpkm2[expIndex,]
jMap = jMap[expIndex]

# # ####################
# # # ##### junction
# pcaJxn = prcomp(t(jRpkm2))
# nsvJxn = num.sv(jRpkm2, mod,vfilter=1e5)
# pcsJxn = pcaJxn$x[,1:nsvJxn]

# covsJxn = SlicedData$new(t(cbind(mod[,-1], pcsJxn)))
# exprsJxn = SlicedData$new(as.matrix(jRpkm2))
# exprsJxn$ResliceCombined(sliceSize = 5000)

# ## JUST CIS
# posJxn = as.data.frame(jMap)[,1:3]
# posJxn$name = names(jMap)
# posJxn = posJxn[,c(4,1:3)]
# names(posJxn)[2] = "chr"

# theSnps = SlicedData$new(snp2)
# theSnps$ResliceCombined(sliceSize = 50000)

# snpspos = snpMap[,c("SNP","CHR","POS")]
# snpspos$CHR = paste0("chr",snpspos$CHR)
# colnames(snpspos) = c("name","chr","pos")

# meJxnCis = Matrix_eQTL_main(snps=theSnps, gene = exprsJxn, 
	# cvrt = covsJxn, output_file_name.cis =  ".txt" ,
	# pvOutputThreshold.cis = 0.0001,  pvOutputThreshold=0,
	# snpspos = snpspos, genepos = posJxn, 
	# useModel = modelLINEAR,	cisDist=5e5,
	# pvalue.hist = 100,min.pv.by.genesnp = TRUE)
# save(meJxnCis, nsvJxn, pcsJxn, file="rdas/junction_eqtl_control_13plus_cisOnly.rda")

#### load back in
load("rdas/junction_eqtl_control_13plus_cisOnly.rda")
eqtl = meJxnCis$cis$eqtl
eqtl = eqtl[eqtl$FDR < 0.01,] # FDR < 1%

## number of effective tests
snpMapGR = GRanges(paste0("chr", snpMap$CHR), IRanges(snpMap$POS, width=1))
pruned = read.table("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/plink/LIBD_Brain_DLPFC_szControls_imputed_indep.prune.in")
snpMap$isLdIndep = snpMap$SNP %in% pruned$V1

nTests = length(findOverlaps(snpMapGR[snpMap$isLdIndep], 
	jMap, maxgap=5e5))
eqtl$bonf = p.adjust(eqtl$pvalue, "bonferroni", n= nTests)
length(unique(eqtl$gene[eqtl$bonf < 0.05])) # 11073

# annotate
eqtl$gene = as.character(eqtl$gene)
eqtl$snps = as.character(eqtl$snps)
colnames(eqtl)[2] = "jxn"

# snp coordinates
m = match(eqtl$snps, snpMap$SNP)
eqtl$snp_chr = snpMap$CHR[m]
eqtl$snp_chr = paste0("chr",eqtl$snp_chr)
eqtl$snp_chr = gsub("chr23","chrX", eqtl$snp_chr)
eqtl$snp_pos = snpMap$POS[m]
eqtl$snpRsNum = snpMap$name[m]

# take best eqtl per jxn
sig = eqtl[!duplicated(eqtl$jxn),]
tt = table(eqtl$jxn)
sig$numSnps = tt[sig$jxn]

# other snp info
snprange = t(sapply(split(eqtl$snp_pos, eqtl$jxn), range))
snprange= snprange[as.character(sig$jxn),]
sig$startSnpLD = snprange[,1]
sig$endSnpLD = snprange[,2]
sig$snpLDLength = snprange[,2]-snprange[,1]

# other annotation
m = match(sig$snps, snpMap$SNP)
sig$snpCounted = snpMap$COUNTED[m]
sig$snpAlt = snpMap$ALT[m]
sig$inSampleMAF = rowSums(snp2[m,],na.rm=TRUE)/
	(2*rowSums(!is.na(snp2[m,])))
	
# annotate junction
library(GenomicRanges)
g = as.data.frame(jMap[match(sig$jxn, names(jMap)),])
colnames(g)[1]="chr"
colnames(g) = paste0("exprs_", colnames(g))
g$exprs_numTx[g$exprs_numTx==0] = NA
g$exprs_ensemblTx = sapply(g$exprs_ensemblTx, paste, collapse=";")
sig = cbind(sig,g)
rownames(sig) = NULL

## cis or trans
start = ifelse(sig$exprs_strand =="+",
	sig$exprs_start, sig$exprs_end)
sig$cisOrTrans = ifelse(sig$snp_chr == sig$exprs_chr & 
	abs(sig$snp_pos - start) < 1e6,	"cis","trans")
sig$distToStart  = sig$snp_pos - start
sig$distToStart[sig$cisOrTrans=="trans"] = NA

sig$meanRPKM = 2^rowMeans(jRpkm2)[sig$jxn]-1

### add race sensitivity analyses
rIndexes=splitit(pd2$Race)
raceJxn = mclapply(rIndexes, function(ii) {
	ssnp = as.matrix(snp2[sig$snps,ii])
	yyJxn = as.matrix(jRpkm2[sig$jxn,ii])
	out = matrix(nrow = nrow(ssnp), nc = 3)
	for(j in 1:nrow(ssnp)) {
		if(j %% 5000 == 0) cat(".")
		out[j,] = summary(lm(yyJxn[j,]~ssnp[j,] + cbind(
			mod[ii,-1], pcsJxn[ii,])))$coef[2,c(1,3,4)]
	}
	colnames(out) = c("slope","tstat","pval")
	cat("\n")
	return(out)
},mc.cores=2)
xx = do.call("cbind", raceJxn)
colnames(xx) = paste0(rep(names(rIndexes), each=3), "_", colnames(xx))
sigJxn = cbind(sig, xx)

mafRace = sapply(rIndexes, function(ii){
	ssnp = as.matrix(snp2[sig$snps,ii])
	rowSums(ssnp, na.rm=TRUE)/(2*rowSums(!is.na(ssnp)))
})
colnames(mafRace) = paste0(colnames(mafRace),"_inSampleMAF")
sigJxn = cbind(sigJxn, mafRace)


save(sigJxn, file="rdas/annotated_junction_eqtl_control_13plus_cisOnly.rda")
