###
library(GenomicRanges)
library(matrixStats)
library(MatrixEQTL)
library(jaffelab)
library(parallel)

##########################
### load GTEx data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/gtex_brainEqtl_subsets.rda")

### get indices to split by region ####
rIndexes=split(1:nrow(pdBrainSub), pdBrainSub$SMTSD)
lengths(rIndexes)
names(rIndexes) = gsub(" ", "_", ss(ss(names(rIndexes), 
	"Brain - ", 2), " \\("))
	
########################### 
## analysis by region #####
###########################

inds = seq(along=rIndexes)
names(inds) = names(rIndexes)

## gene level
geneEqtlGtex = mclapply(inds, function(i) {
	ii = rIndexes[[i]]
	## model
	mod = cbind(as.matrix(mds[pdBrainSub$snpColumn[ii],4:8]), pcList$gene[ii,])
	colnames(mod) = paste0(rep(c("snpPC", "exprsPC"), times=c(5,10)), c(1:5,1:10))
	covs = SlicedData$new(t(mod))

	#####################
	#### gene eqtls #####
	exprsGene = SlicedData$new(as.matrix(log2(geneRpkmSub[,ii]+1)))
	exprsGene$ResliceCombined(sliceSize = 1000)

	posGene = geneMapSub[,1:3]
	posGene$name = rownames(geneMapSub)
	posGene = posGene[,c(4,1:3)]

	## make SNP objects
	sIndex = which(snpMapSub$chrpos %in% geneEqtl2$snp_chrpos)
	theSnpsGene = SlicedData$new(as.matrix(snpSub)[sIndex,pdBrainSub$snpColumn[ii]])
	theSnpsGene$ResliceCombined(sliceSize = 50000)

	snpsposGene = snpMapSub[sIndex,c("SNP","CHR","POS")]
	snpsposGene$CHR = paste0("chr",snpsposGene$CHR)
	colnames(snpsposGene) = c("name","chr","pos")

	## run qtl
	meGeneReplication = Matrix_eQTL_main(snps=theSnpsGene, 
		gene = exprsGene, cvrt = covs, 
		output_file_name.cis =  paste0("/dcl01/ajaffe/data/scratch/gene/",
				names(rIndexes)[i], ".txt"),
		pvOutputThreshold.cis = 1, pvOutputThreshold=0,
		snpspos = snpsposGene, genepos = posGene, 
		useModel = modelLINEAR,	cisDist=5e5)	
	## merge back into LIBD
	geneEqtlRep = meGeneReplication$cis$eqtls
	geneEqtlRep$snp_chrpos = snpMapSub$chrpos[
		match(as.character(geneEqtlRep$snps), snpMapSub$SNP)]
	geneEqtlRep$gtexSnpCounted= snpMapSub$COUNTED[
		match(as.character(geneEqtlRep$snps), snpMapSub$SNP)]

	geneEqtlRep$pairID = paste0(geneEqtlRep$snp_chrpos,
		"_", as.character(geneEqtlRep$gene))
	geneEqtl2$pairID = paste0(geneEqtl2$snp_chrpos, "_", geneEqtl2$gene)
	geneEqtlRep = geneEqtlRep[match(geneEqtl2$pairID, geneEqtlRep$pairID),c(3:6,8)]

	### flip those as needed
	toFlipGene = which(geneEqtl2$snpCounted != geneEqtlRep$gtexSnpCounted)
	geneEqtlRep[toFlipGene,c(1,4)] = -1*geneEqtlRep[toFlipGene,c(1,4)]
	geneEqtlRep$gtexSnpCounted = NULL
	return(geneEqtlRep)
},mc.cores=6)



#####################
#### exon eqtls #####
#####################
# exonMapSub$Chr = paste0("chr", exonMapSub$Chr)

exonEqtlGtex = mclapply(inds, function(i) {
	ii = rIndexes[[i]]
		
	## model
	mod = cbind(as.matrix(mds[pdBrainSub$snpColumn[ii],4:8]), pcList$exon[ii,])
	colnames(mod) = paste0(rep(c("snpPC", "exprsPC"), times=c(5,10)), c(1:5,1:10))
	covs = SlicedData$new(t(mod))

	#####################
	#### gene eqtls #####
	exprsExon = SlicedData$new(as.matrix(log2(exonRpkmSub[,ii]+1)))
	exprsExon$ResliceCombined(sliceSize = 1000)

	posExon = exonMapSub[,2:4]
	posExon$name = rownames(exonMapSub)
	posExon = posExon[,c(4,1:3)]

	## make SNP objects
	sIndex = which(snpMapSub$chrpos %in% exonEqtl2$snp_chrpos)
	theSnpsExon = SlicedData$new(as.matrix(snpSub)[sIndex,pdBrainSub$snpColumn[ii]])
	theSnpsExon$ResliceCombined(sliceSize = 50000)

	snpsposExon = snpMapSub[sIndex,c("SNP","CHR","POS")]
	snpsposExon$CHR = paste0("chr",snpsposExon$CHR)
	colnames(snpsposExon) = c("name","chr","pos")

	## run qtl
	meExonReplication = Matrix_eQTL_main(snps=theSnpsExon, 
		gene = exprsExon, cvrt = covs, 
		output_file_name.cis =  paste0("/dcl01/ajaffe/data/scratch/exon/",
				names(rIndexes)[i], ".txt"),
			pvOutputThreshold.cis = 1, pvOutputThreshold=0,
		snpspos = snpsposExon, genepos = posExon, 
		useModel = modelLINEAR,	cisDist=5e5)

	## merge back into LIBD
	exonEqtlRep = meExonReplication$cis$eqtls
	exonEqtlRep$snp_chrpos = snpMapSub$chrpos[
		match(as.character(exonEqtlRep$snps), snpMapSub$SNP)]
	exonEqtlRep$gtexSnpCounted= snpMapSub$COUNTED[
		match(as.character(exonEqtlRep$snps), snpMapSub$SNP)]

	exonEqtlRep$pairID = paste0(exonEqtlRep$snp_chrpos,
		"_", as.character(exonEqtlRep$gene))
	exonEqtl2$pairID = paste0(exonEqtl2$snp_chrpos, "_", exonEqtl2$gene)
	exonEqtlRep = exonEqtlRep[match(exonEqtl2$pairID, exonEqtlRep$pairID),c(3:6,8)]

	### flip those as needed
	toFlipExon = which(exonEqtl2$snpCounted != exonEqtlRep$gtexSnpCounted)
	exonEqtlRep[toFlipExon,c(1,4)] = -1*exonEqtlRep[toFlipExon,c(1,4)]
	exonEqtlRep$gtexSnpCounted = NULL
	return(exonEqtlRep)
},mc.cores=6)


####################################
######### junction eqtls ############

jxnEqtlGtex = mclapply(inds, function(i) {
	ii = rIndexes[[i]]
		
	## model
	mod = cbind(as.matrix(mds[pdBrainSub$snpColumn[ii],4:8]), pcList$junction[ii,])
	colnames(mod) = paste0(rep(c("snpPC", "exprsPC"), times=c(5,10)), c(1:5,1:10))
	covs = SlicedData$new(t(mod))

	#####################
	#### gene eqtls #####
	exprsJxn = SlicedData$new(as.matrix(log2(jRpkmSub[,ii]+1)))
	exprsJxn$ResliceCombined(sliceSize = 5000)

	posJxn = as.data.frame(jMapSub)[,1:3]
	posJxn$name = names(jMapSub)
	posJxn = posJxn[,c(4,1:3)]
	
	## make SNP objects
	sIndex = which(snpMapSub$chrpos %in% jxnEqtl2$snp_chrpos)
	theSnpsJxn = SlicedData$new(as.matrix(snpSub)[sIndex,pdBrainSub$snpColumn[ii]])
	theSnpsJxn$ResliceCombined(sliceSize = 50000)

	snpsposJxn = snpMapSub[sIndex,c("SNP","CHR","POS")]
	snpsposJxn$CHR = paste0("chr",snpsposJxn$CHR)
	colnames(snpsposJxn) = c("name","chr","pos")

	## run qtl
	meJxnReplication = Matrix_eQTL_main(snps=theSnpsJxn, 
		gene = exprsJxn, cvrt = covs, 
		output_file_name.cis =  paste0("/dcl01/ajaffe/data/scratch/junction/",
				names(rIndexes)[i], ".txt"),
		pvOutputThreshold.cis = 1, pvOutputThreshold=0,
		snpspos = snpsposJxn, genepos = posJxn, 
		useModel = modelLINEAR,	cisDist=5e5,
		pvalue.hist = 100,min.pv.by.genesnp = TRUE)

	## merge back into LIBD
	jxnEqtlRep = meJxnReplication$cis$eqtls
	jxnEqtlRep$snp_chrpos = snpMapSub$chrpos[
		match(as.character(jxnEqtlRep$snps), snpMapSub$SNP)]
	jxnEqtlRep$gtexSnpCounted= snpMapSub$COUNTED[
		match(as.character(jxnEqtlRep$snps), snpMapSub$SNP)]

	jxnEqtlRep$pairID = paste0(jxnEqtlRep$snp_chrpos,
		"_", as.character(jxnEqtlRep$gene))
	jxnEqtl2$pairID = paste0(jxnEqtl2$snp_chrpos, "_", jxnEqtl2$gene)
	jxnEqtlRep = jxnEqtlRep[match(jxnEqtl2$pairID, jxnEqtlRep$pairID),c(3:6,8)]

	### flip those as needed
	toFlipJxn = which(jxnEqtl2$snpCounted != jxnEqtlRep$gtexSnpCounted)
	jxnEqtlRep[toFlipJxn,c(1,4)] = -1*jxnEqtlRep[toFlipJxn,c(1,4)]
	jxnEqtlRep$gtexSnpCounted = NULL
	return(jxnEqtlRep)
},mc.cores=6)


################
## DERs ##
derEqtl2$pairID = paste0(derEqtl2$snp_chrpos, "_", derEqtl2$gene)

derEqtlGtex = mclapply(inds, function(i) {
	ii = rIndexes[[i]]
		
	## model
	mod = cbind(as.matrix(mds[pdBrainSub$snpColumn[ii],4:8]), pcList$der[ii,])
	colnames(mod) = paste0(rep(c("snpPC", "exprsPC"), times=c(5,10)), c(1:5,1:10))
	covs = SlicedData$new(t(mod))

	#####################
	#### der eqtls #####
	exprsDer = SlicedData$new(as.matrix(log2(covMatSub[,ii]+1)))
	exprsDer$ResliceCombined(sliceSize = 5000)

	posDer = as.data.frame(covMapSub)[,1:3]
	posDer$name = names(covMapSub)
	posDer = posDer[,c(4,1:3)]
	
	## make SNP objects
	sIndex = which(snpMapSub$chrpos %in% derEqtl2$snp_chrpos)
	theSnpsDer = SlicedData$new(as.matrix(snpSub)[sIndex,pdBrainSub$snpColumn[ii]])
	theSnpsDer$ResliceCombined(sliceSize = 50000)

	snpsposDer = snpMapSub[sIndex,c("SNP","CHR","POS")]
	snpsposDer$CHR = paste0("chr",snpsposDer$CHR)
	colnames(snpsposDer) = c("name","chr","pos")

	## run qtl
	meDerReplication = Matrix_eQTL_main(snps=theSnpsDer, 
		gene = exprsDer, cvrt = covs, 
		output_file_name.cis =  paste0("/dcl01/ajaffe/data/scratch/er/",
				names(rIndexes)[i], ".txt"),
		pvOutputThreshold.cis = 1, pvOutputThreshold=0,
		snpspos = snpsposDer, genepos = posDer, 
		useModel = modelLINEAR,	cisDist=5e5)
	derEqtlRep = meDerReplication$cis$eqtls

	## extract same eQTLs
	derEqtlRep$snp_chrpos = snpMapSub$chrpos[
		match(as.character(derEqtlRep$snps), snpMapSub$SNP)]
	derEqtlRep$gtexSnpCounted= snpMapSub$COUNTED[
		match(as.character(derEqtlRep$snps), snpMapSub$SNP)]

	derEqtlRep$pairID = paste0(derEqtlRep$snp_chrpos,
		"_", as.character(derEqtlRep$gene))
	derEqtlRep = derEqtlRep[match(derEqtl2$pairID, derEqtlRep$pairID),c(3:6,8)]

	### flip those as needed
	toFlipDer = which(derEqtl2$snpCounted != derEqtlRep$gtexSnpCounted)
	derEqtlRep[toFlipDer,c(1,4)] = -1*derEqtlRep[toFlipDer,c(1,4)]
	derEqtlRep$gtexSnpCounted = NULL
	return(derEqtlRep)
	## merge back into LIBD
}, mc.cores=6)

	


## Txs ##
txEqtlGtex = mclapply(inds, function(i) {
	ii = rIndexes[[i]]
		
	## model
	mod = cbind(as.matrix(mds[pdBrainSub$snpColumn[ii],4:8]), pcList$transcript[ii,])
	colnames(mod) = paste0(rep(c("snpPC", "exprsPC"), times=c(5,10)), c(1:5,1:10))
	covs = SlicedData$new(t(mod))

	#####################
	#### tx eqtls #####
	exprsTrans = SlicedData$new(as.matrix(log2(tFpkmSub[,ii]+1)))
	exprsTrans$ResliceCombined(sliceSize = 5000)

	posTrans = as.data.frame(tMapSub)[,1:3]
	posTrans$name = names(tMapSub)
	posTrans = posTrans[,c(4,1:3)]

	## make SNP objects
	sIndex = which(snpMapSub$chrpos %in% transEqtl2$snp_chrpos)
	theSnpsTrans = SlicedData$new(as.matrix(snpSub)[sIndex,pdBrainSub$snpColumn[ii]])
	theSnpsTrans$ResliceCombined(sliceSize = 50000)

	snpsposTrans = snpMapSub[sIndex,c("SNP","CHR","POS")]
	snpsposTrans$CHR = paste0("chr",snpsposTrans$CHR)
	colnames(snpsposTrans) = c("name","chr","pos")

	## run qtl
	meTransReplication = Matrix_eQTL_main(snps=theSnpsTrans, 
		gene = exprsTrans, cvrt = covs, 
		output_file_name.cis =  paste0("/dcl01/ajaffe/data/scratch/transcript/",
				names(rIndexes)[i], ".txt"),
		pvOutputThreshold.cis = 1, pvOutputThreshold=0,
		snpspos = snpsposTrans, genepos = posTrans, 
		useModel = modelLINEAR,	cisDist=5e5,
		pvalue.hist = 100,min.pv.by.genesnp = TRUE)

	## merge back into LIBD
	transEqtlRep = meTransReplication$cis$eqtls
	transEqtlRep$snp_chrpos = snpMapSub$chrpos[
		match(as.character(transEqtlRep$snps), snpMapSub$SNP)]
	transEqtlRep$gtexSnpCounted= snpMapSub$COUNTED[
		match(as.character(transEqtlRep$snps), snpMapSub$SNP)]

	transEqtlRep$pairID = paste0(transEqtlRep$snp_chrpos,
		"_", as.character(transEqtlRep$gene))
	transEqtl2$pairID = paste0(transEqtl2$snp_chrpos, "_", transEqtl2$gene)
	transEqtlRep = transEqtlRep[match(transEqtl2$pairID, transEqtlRep$pairID),c(3:6,8)]

	### flip those as needed
	toFlipTrans = which(transEqtl2$snpCounted != transEqtlRep$gtexSnpCounted)
	transEqtlRep[toFlipTrans,c(1,4)] = -1*transEqtlRep[toFlipTrans,c(1,4)]
	transEqtlRep$gtexSnpCounted = NULL
	return(transEqtlRep)
},mc.cores=6)

## save
save(geneEqtlGtex, exonEqtlGtex, jxnEqtlGtex,
	txEqtlGtex,  derEqtlGtex, compress=TRUE,
	file = "/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/gtex_replication_stats.rda")

#### combine
geneEqtlGtexRep = do.call("cbind", geneEqtlGtex)
geneEqtlGtexRep$pairID= paste0(geneEqtl2$snps, ".", geneEqtl2$gene)
exonEqtlGtexRep = do.call("cbind", exonEqtlGtex)
exonEqtlGtexRep$pairID= paste0(exonEqtl2$snps, ".", exonEqtl2$gene)
jxnEqtlGtexRep = do.call("cbind", jxnEqtlGtex)
jxnEqtlGtexRep$pairID= paste0(jxnEqtl2$snps, ".", jxnEqtl2$gene)

derEqtlGtexRep = do.call("cbind", derEqtlGtex)
derEqtlGtexRep$pairID= paste0(derEqtl2$snps, ".", derEqtl2$gene)
txEqtlGtexRep = do.call("cbind", txEqtlGtex)
txEqtlGtexRep$pairID= paste0(transEqtl2$snps, ".", transEqtl2$gene)

## add pair IDs as labels
rownames(geneEqtlGtexRep) = paste0(geneEqtl2$snps, ".", geneEqtl2$gene)
rownames(exonEqtlGtexRep) = paste0(exonEqtl2$snps, ".", exonEqtl2$gene)
rownames(jxnEqtlGtexRep) = paste0(jxnEqtl2$snps, ".", jxnEqtl2$gene)
rownames(txEqtlGtexRep) = paste0(transEqtl2$snps, ".", transEqtl2$gene)
rownames(derEqtlGtexRep) = paste0(derEqtl2$snps, ".", derEqtl2$gene)

jxnEqtlRep = jEqtlRep # rename

## combine
allEqtlRepGtex = rbind(geneEqtlGtexRep, exonEqtlGtexRep, 
	jxnEqtlGtexRep, txEqtlGtexRep, derEqtlGtexRep)

## save
save(allEqtlRepGtex, compress=TRUE,
	file = "/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/gtex_brainEqtl_replication_stats.rda")
	
#####	
# check stats

geneStats = as.data.frame(geneStats)
geneStats$tested = rowSums(is.na(geneStats)) < 2
geneStats$direction = sign(geneStats$beta) ==sign(geneEqtl2$beta)
geneStats$dirAndSig = sign(geneStats$beta) ==sign(geneEqtl2$beta) & 
	geneStats$pval < 0.01
geneStats$dirAndSig5 = sign(geneStats$beta) ==sign(geneEqtl2$beta) & 
	geneStats$pval < 1e-5

exonStats = as.data.frame(exonStats)
exonStats$tested = rowSums(is.na(exonStats)) < 2
exonStats$direction = sign(exonStats$beta) ==sign(exonEqtl2$beta)
exonStats$dirAndSig = sign(exonStats$beta) ==sign(exonEqtl2$beta) & 
	exonStats$pval < 0.01
exonStats$dirAndSig5 = sign(exonStats$beta) ==sign(exonEqtl2$beta) & 
	exonStats$pval < 1e-5

jxnStats = as.data.frame(jxnStats)
jxnStats$tested = rowSums(is.na(jxnStats)) < 2
jxnStats$direction = sign(jxnStats$beta) ==sign(jxnEqtl2_filter$beta)
jxnStats$dirAndSig = sign(jxnStats$beta) ==sign(jxnEqtl2_filter$beta) &
	jxnStats$pval < 0.01
jxnStats$dirAndSig5 = sign(jxnStats$beta) ==sign(jxnEqtl2_filter$beta) &
	jxnStats$pval < 1e-5

derStats = as.data.frame(derStats)
derStats$tested = rowSums(is.na(derStats)) < 2
derStats$direction = sign(derStats$beta) ==sign(derEqtl2$beta)
derStats$dirAndSig = sign(derStats$beta) ==sign(derEqtl2$beta) & 
	derStats$pval < 0.01
derStats$dirAndSig5 = sign(derStats$beta) ==sign(derEqtl2$beta) & 
	derStats$pval < 1e-5

transStats = as.data.frame(transStats)
transStats$tested = rowSums(is.na(transStats)) < 2
transStats$direction = sign(transStats$beta) ==sign(transEqtl2$beta)
transStats$dirAndSig = sign(transStats$beta) ==sign(transEqtl2$beta) & 
	transStats$pval < 0.01
transStats$dirAndSig5 = sign(transStats$beta) ==sign(transEqtl2$beta) & 
	transStats$pval < 1e-5

################
## load eqtls for means and metrics

load("rdas/annotated_junction_eqtl_szControl_cisOnly.rda")
load("rdas/annotated_gene_eqtl_szControl_cisOnly.rda")
load("rdas/annotated_exon_eqtl_szControl_cisOnly.rda")
load("rdas/annotated_der_eqtl_szControl_cisOnly.rda")
load("rdas/annotated_transcript_eqtl_szControl_cisOnly.rda")

sigGene2 = sigGene[match(rownames(geneStats), sigGene$gene),]
sigExon2 = sigExon[match(rownames(exonStats), sigExon$exon),]
sigJxn2 = sigJxn[match(rownames(jxnStats), sigJxn$jxn),]
sigDer2 = sigDer[match(rownames(derStats), sigDer$der),]
sigTrans2 = sigTrans[match(rownames(transStats), sigTrans$tx),]

geneStats$meanRpkm = sigGene2$meanRPKM
exonStats$meanRpkm = sigExon2$meanRPKM
jxnStats$meanRp80m = sigJxn2$meanRPKM
derStats$meanRp80m = sigDer2$exprs_meanCov
transStats$meanFpkm = sigTrans2$exprs_meanCov


save(geneStats, exonStats, jxnStats, derStats,transStats,
	file="rdas/CMC_replication_stats_overall.rda")
	
	
#############################
### NUMBERS FOR THE PAPER ###

load("rdas/CMC_replication_stats_overall.rda")

##### gene
## gene replication ##
	
gInd = sigGene2$bonf < 0.05
colSums(geneStats[,3:6], na.rm=TRUE)
colSums(geneStats[gInd,3:6], na.rm=TRUE)

plot(-log10(geneStats[,2]), -log10(geneEqtl2$pval))
plot(geneStats[,1], geneEqtl2$beta,xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
plot(-log10(geneStats[gInd,2]), -log10(geneEqtl2$pval[gInd]))
plot(geneStats[gInd,1], geneEqtl2$beta[gInd],xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))

##########
#### exon?
eInd = sigExon2$bonf < 0.05
colSums(exonStats[,3:6], na.rm=TRUE)
colSums(exonStats[eInd,3:6], na.rm=TRUE)

plot(-log10(exonStats[,2]), -log10(exonEqtl2$pval))
plot(exonStats[,1], exonEqtl2$beta,xlim=c(-1,1),ylim=c(-1,1))
plot(-log10(exonStats[eInd,2]), -log10(exonEqtl2$pval[eInd]))
plot(exonStats[eInd,1], exonEqtl2$beta[eInd],xlim=c(-1,1),ylim=c(-1,1))

##########
#### transcript?
tInd = sigTrans2$bonf < 0.05
colSums(transStats[,3:6], na.rm=TRUE)
colSums(transStats[tInd,3:6], na.rm=TRUE)

plot(-log10(exonStats[,2]), -log10(exonEqtl2$pval))
plot(exonStats[,1], exonEqtl2$beta,xlim=c(-1,1),ylim=c(-1,1))
plot(-log10(exonStats[eInd,2]), -log10(exonEqtl2$pval[eInd]))
plot(exonStats[eInd,1], exonEqtl2$beta[eInd],xlim=c(-1,1),ylim=c(-1,1))

## junction
jInd = sigJxn2$bonf < 0.05
colSums(jxnStats[,3:6], na.rm=TRUE)
colSums(jxnStats[jInd,3:6], na.rm=TRUE)

## unannoated
jInd2 = sigJxn2$bonf < 0.05 & sigJxn2$exprs_code != "InEns"
colSums(jxnStats[jInd2,3:6], na.rm=TRUE)

plot(-log10(jxnStats[,2]), -log10(jxnEqtl2_filter$pval))
plot(jxnStats[,1], jxnEqtl2_filter$beta,xlim=c(-1,1),ylim=c(-1,1))
plot(-log10(jxnStats[jInd,2]), -log10(jxnEqtl2_filter$pval[jInd]))
plot(jxnStats[jInd,1], jxnEqtl2_filter$beta[jInd],xlim=c(-1,1),ylim=c(-1,1))

## der
dInd = sigDer2$bonf < 0.05
colSums(derStats[,3:6], na.rm=TRUE)
colSums(derStats[dInd,3:6], na.rm=TRUE)

table(derStats$tested)
mean(derStats$direction,na.rm=TRUE)
mean(derStats$dirAndSig,na.rm=TRUE)
mean(derStats$dirAndSig5,na.rm=TRUE)

table(derStats$tested[dInd])
mean(derStats$direction[dInd],na.rm=TRUE)
mean(derStats$dirAndSig[dInd],na.rm=TRUE)
mean(derStats$dirAndSig5[dInd],na.rm=TRUE)

plot(-log10(derStats[,2]), -log10(derEqtl2$pval))
plot(derStats[,1], derEqtl2$beta,xlim=c(-1,1),ylim=c(-1,1))
plot(-log10(derStats[dInd,2]), -log10(derEqtl2$pval[dInd]))
plot(derStats[dInd,1], derEqtl2$beta[dInd],xlim=c(-1,1),ylim=c(-1,1))

## trans
tInd = sigTrans2$bonf < 0.05
table(transStats$tested)
mean(transStats$direction,na.rm=TRUE)
mean(transStats$dirAndSig,na.rm=TRUE)
mean(transStats$dirAndSig5,na.rm=TRUE)

table(transStats$tested[tInd])
mean(transStats$direction[tInd],na.rm=TRUE)
mean(transStats$dirAndSig[tInd],na.rm=TRUE)
mean(transStats$dirAndSig5[tInd],na.rm=TRUE)

plot(-log10(transStats[,2]), -log10(transEqtl2$pval))
plot(transStats[,1], transEqtl2$beta,xlim=c(-1,1),ylim=c(-1,1))
plot(-log10(transStats[dInd,2]), -log10(transEqtl2$pval[dInd]))
plot(transStats[dInd,1], transEqtl2$beta[dInd],xlim=c(-1,1),ylim=c(-1,1))


#### check factors for predicting replication
geneStats$negLog10pval = -log10(geneStats$pval)
geneStats$Symbol = sigGene2$exprs_symbol
fGene = summary(glm(dirAndSig ~ abs(beta) + meanRpkm + 
	negLog10pval + (Symbol != ""), data=geneStats))
fGene2 = summary(glm(dirAndSig ~ abs(beta) + meanRpkm + 
	negLog10pval + (Symbol != ""), data=geneStats,subset=gInd))

exonStats$negLog10pval = -log10(exonStats$pval)
exonStats$Symbol = sigExon2$exprs_symbol
fExon = summary(glm(dirAndSig ~ abs(beta) + meanRpkm + 
	negLog10pval + (Symbol != ""), data=exonStats))
fExon2 = summary(glm(dirAndSig ~ abs(beta) + meanRpkm + 
	negLog10pval + (Symbol != ""), data=exonStats,subset=eInd))

jxnStats$negLog10pval = -log10(jxnStats$pval)
jxnStats$Symbol = sigJxn2$exprs_newGeneSymbol
fJxn = summary(glm(dirAndSig ~ abs(beta) + meanRp80m + 
	negLog10pval + !is.na(Symbol), data=jxnStats))
fJxn2 = summary(glm(dirAndSig ~ abs(beta) + meanRp80m + 
	negLog10pval + !is.na(Symbol), data=jxnStats,subset=jInd))

derStats$negLog10pval = -log10(derStats$pval)
derStats$Symbol = sigDer2$exprs_symbol
fDer = summary(glm(dirAndSig ~ abs(beta) + meanRp80m + 
	negLog10pval + (Symbol != ""), data=derStats))
fDer2 = summary(glm(dirAndSig ~ abs(beta) + meanRp80m + 
	negLog10pval + (Symbol != ""), data=derStats,subset=dInd))

repStats = data.frame(gene_OR = exp(fGene$coef[,1]),
	gene_pval = fGene$coef[,4], exon_OR = exp(fExon$coef[,1]),
	exon_pval = fExon$coef[,4], jxn_OR = exp(fJxn$coef[,1]),
	jxn_pval = fJxn$coef[,4], der_OR = exp(fDer$coef[,1]),
	der_pval = fDer$coef[,4])
repStatsBonf = data.frame(gene_OR = exp(fGene2$coef[,1]),
	gene_pval = fGene2$coef[,4], exon_OR = exp(fExon2$coef[,1]),
	exon_pval = fExon2$coef[,4], jxn_OR = exp(fJxn2$coef[,1]),
	jxn_pval = fJxn2$coef[,4], der_OR = exp(fDer2$coef[,1]),
	der_pval = fDer2$coef[,4])
rownames(repStats) =rownames(repStatsBonf) = c("Intercept", 
	"Abs(beta)", "ExprsLevel","-log10(p)","hasSymbol")
repStats
repStatsBonf