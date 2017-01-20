###
library(GenomicRanges)
library(matrixStats)

ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

##########################
### load CMC data
load("rdas/CMC_brainEqtl_subsets.rda")

#####################
#### gene eqtls #####
snpGene = as.matrix(snpSub[match(geneEqtl2$snp_chrpos, 
	snpMapSub$chrpos),pdSub$snpColumn])
snpGeneMap = snpMapSub[match(geneEqtl2$snp_chrpos, snpMapSub$chrpos),]

## try overall data, pooling across regions
geneStats = matrix(NA,  nrow = nrow(snpGene), ncol=2,
			dimname = list(rownames(geneRpkmSub), c("beta","pval")))
g = log2(geneRpkmSub+1)
identical(rownames(g), geneEqtl2$gene) # TRUE
mmds = as.matrix(mds[match(pdSub$Genotyping_Sample_ID, mds$FID),4:6])
pcs = pcList$gene[,1:10]

usnps = apply(snpGene, 1, table)
useIndex=as.numeric(which(sapply(usnps,length)>1 &	
	rowSums(g == 0) < ncol(g)*.75))

for(j in useIndex) {
	if(j %% 500 == 0) cat(".")
	f=summary(lm(g[j,] ~ snpGene[j,] + mmds + pcs))
	geneStats[j,] = f$coef[2,c(1,4)]
}

### flip those as needed
toFlipGene = which(geneEqtl2$snpCounted != snpGeneMap$COUNTED)
geneStats[toFlipGene,1] = -1*geneStats[toFlipGene,1]
	
#####################
#### exon eqtls #####
snpExon = as.matrix(snpSub[match(exonEqtl2$snp_chrpos, 
	snpMapSub$chrpos),pdSub$snpColumn])
snpExonMap = snpMapSub[match(exonEqtl2$snp_chrpos, snpMapSub$chrpos),]

## try overall data, pooling across regions
exonStats = matrix(NA,  nrow = nrow(snpExon), ncol=2,
			dimname = list(rownames(exonRpkmSub), c("beta","pval")))
g = log2(exonRpkmSub+1)
pcs = pcList$exon[,1:10]

usnps = apply(snpExon, 1, table)
useIndex=which(sapply(usnps,length)>1 &	rowSums(g == 0) < ncol(g)*.75)

for(j in useIndex) {
	if(j %% 500 == 0) cat(".")
	f=summary(lm(g[j,] ~ snpExon[j,] + mmds + pcs))
	exonStats[j,] = f$coef[2,c(1,4)]
}

### flip those as needed
toFlipExon = which(exonEqtl2$snpCounted != snpExonMap$COUNTED)
exonStats[toFlipExon,1] = -1*exonStats[toFlipExon,1]

####################################
######### junction eqtls ############

snpJxn = as.matrix(snpSub[match(jxnEqtl2$snp_chrpos, 
	snpMapSub$chrpos),pdSub$snpColumn])
snpJxnMap = snpMapSub[match(jxnEqtl2$snp_chrpos, snpMapSub$chrpos),]

## some junctions werent present
jInd = match(rownames(jRpkmSub), jxnEqtl2$gene)
snpJxn = snpJxn[jInd,]
snpJxnMap = snpJxnMap[jInd,]
jxnEqtl2_filter = jxnEqtl2[jInd,]

## try overall data, pooling across regions
jxnStats = matrix(NA,  nrow = nrow(jRpkmSub), ncol=2,
			dimname = list(rownames(jRpkmSub), c("beta","pval")))
g = as.matrix(log2(jRpkmSub+1))
pcs = pcList$junction[,1:10]

usnps = apply(snpJxn, 1, table)
useIndex=which(sapply(usnps,length)>1 &	rowSums(g == 0) < ncol(g)*.75)

for(j in useIndex) {
	if(j %% 500 == 0) cat(".")
	f=summary(lm(g[j,] ~ snpJxn[j,] + mmds + pcs))
	jxnStats[j,] = f$coef[2,c(1,4)]
}

### flip those as needed
toFlipJxn = which(jxnEqtl2_filter$snpCounted != snpJxnMap$COUNTED)
jxnStats[toFlipJxn,1] = -1*jxnStats[toFlipJxn,1]

################
## DERs ##
snpDer = as.matrix(snpSub[match(derEqtl2$snp_chrpos, 
	snpMapSub$chrpos),pdSub$snpColumn])
snpDerMap = snpMapSub[match(derEqtl2$snp_chrpos, snpMapSub$chrpos),]

## try overall data, pooling across regions
derStats = matrix(NA,  nrow = nrow(snpDer), ncol=2,
			dimname = list(rownames(regionMatSub), c("beta","pval")))

g = log2(regionMatSub+1)
pcs = pcList$der[,1:10]

usnps = apply(snpDer, 1, table)
useIndex=which(sapply(usnps,length)>1 &	rowSums(g == 0) < ncol(g)*.75)

for(j in useIndex) {
	if(j %% 500 == 0) cat(".")
	f=summary(lm(g[j,] ~ snpDer[j,] + mmds + pcs))
	derStats[j,] = f$coef[2,c(1,4)]
}

### flip those as needed
toFlipDer = which(derEqtl2$snpCounted != snpDerMap$COUNTED)
derStats[toFlipDer,1] = -1*derStats[toFlipDer,1]

################
## Txs ##
snpTrans = as.matrix(snpSub[match(transEqtl2$snp_chrpos, 
	snpMapSub$chrpos),pdSub$snpColumn])
snpTransMap = snpMapSub[match(transEqtl2$snp_chrpos, snpMapSub$chrpos),]

## try overall data, pooling across regions
transStats = matrix(NA,  nrow = nrow(snpTrans), ncol=2,
			dimname = list(rownames(tFpkmSub), c("beta","pval")))

g = log2(tFpkmSub+1)
pcs = pcList$transcript[,1:10]

usnps = apply(snpTrans, 1, table)
useIndex=which(sapply(usnps,length)>1 &	rowSums(g == 0) < ncol(g)*.75)

for(j in useIndex) {
	if(j %% 500 == 0) cat(".")
	f=summary(lm(g[j,] ~ snpTrans[j,] + mmds + pcs))
	transStats[j,] = f$coef[2,c(1,4)]
}

### flip those as needed
toFlipTrans = which(transEqtl2$snpCounted != snpTransMap$COUNTED)
transStats[toFlipTrans,1] = -1*transStats[toFlipTrans,1]

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