#######
library(GenomicRanges)
source("../../eqtl_functions.R")


### load data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/cleaned_eqtl_data_subset_n412.rda")
snpMapSub$chrpos= paste0("chr", snpMapSub$CHR, ":", snpMapSub$POS)

### load eQTL results
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/allEqtlsRep.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/sigEqtlsRep.rda")
allEqtlsRep$chrpos = snpMapSub$chrpos[match(allEqtlsRep$SNP, snpMapSub$SNP)]

### load SZ GWAS clumps ###
pgc =read.delim("/users/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/pgc/daner_PGC_SCZ52_0513a.gz.p4.clump.areator.sorted.1mhc.txt",
	as.is=TRUE, header=TRUE)
tmp = unlist(strsplit(pgc$LD.friends.0.1..p0.001,","))
tmp2 = strsplit(pgc$LD.friends.0.1..p0.001,",")

theSnps = data.frame(name = c(pgc$SNP, ss(tmp, "\\(")),
	R2 = as.numeric(c(rep(1,nrow(pgc)), ss(ss(tmp, "\\(",2),"/"))),
	Dist = as.numeric(c(rep(0,nrow(pgc)), 
		gsub(")", "", ss(ss(tmp, "\\(",2),"/",2),fixed=TRUE))))

# map back
theSnps$hitIndex=c(1:nrow(pgc), rep(seq(along=tmp2), times=sapply(tmp2,length)))
theSnps$Genes = pgc$genes.6.50kb.dist2index.[theSnps$hitIndex]
theSnps$pvalue = pgc$P[theSnps$hitIndex]

# drop those that don't map, and add coordinate info
library(GenomicRanges)
load("/users/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/rdas/granges_pgc2.rda")
pgcFullGR$riskAllele = ifelse(pgcFullGR$or > 1, 
	pgcFullGR$A1, pgcFullGR$A2) # risk allele

## match up
mm = match(theSnps$name, names(pgcFullGR))
theSnps = theSnps[!is.na(mm),]
pgcStats = pgcFullGR[mm[!is.na(mm)]]
theSnps$chrpos = paste0(seqnames(pgcStats), ":", start(pgcStats))
theSnps$riskAllele = pgcStats$riskAllele

## drop those not with R^2 > 0.6
theSnps = theSnps[which(theSnps$R2 > 0.6),]

######################
### pgc index SNPs ###
pgcSig = read.delim("/users/ajaffe/Lieber/Projects/RNAseq/theControlEqtl/pgc/pgc2_128loci.txt")
pgcSig = pgcSig[order(pgcSig$Rank),]
mm = match(pgcSig$Index_SNP, names(pgcFullGR))
pgcStatsSig = pgcFullGR[mm]
pgcSig$chrpos = paste0(seqnames(pgcStatsSig), 
	":", start(pgcStatsSig))
theSnps$finalHitIndex = match(theSnps$chrpos, pgcSig$chrpos)
	
#########################
## how many are in our data to begin with?
# load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/snpData_LIBD_szControls_cleaned.rda")
# snpMap$chrpos= paste0("chr", snpMap$CHR, ":", snpMap$POS)
# save(snpMap, file="rdas/snpMap_LIBD_maf05.rda")
load("rdas/snpMap_LIBD_maf05.rda")

theSnps$inLibd = theSnps$chrpos %in% snpMap$chrpos  | 
	theSnps$name %in% snpMap$name
table(theSnps$inLibd)
theSnps = theSnps[theSnps$inLibd,] # filter

#### filter our eQTLs
pgcEqtls = allEqtlsRep[allEqtlsRep$chrpos %in% theSnps$chrpos,]

## match up alleles
pgcEqtls$snpCountedAllele = snpMapSub$COUNTED[match(pgcEqtls$SNP, snpMapSub$SNP)]
pgcEqtls$snpRefAllele = snpMapSub$ALT[match(pgcEqtls$SNP, snpMapSub$SNP)]

## add stats
pgcEqtls$pgcRow = match(pgcEqtls$chrpos, theSnps$chrpos)
pgcEqtls$pgcDiscRank = theSnps$hitIndex[pgcEqtls$pgcRow]
pgcEqtls$pgcFinalRank = theSnps$finalHitIndex[pgcEqtls$pgcRow]
pgcEqtls$riskAllele = theSnps$riskAllele[pgcEqtls$pgcRow]

#################
## flip slopes  #
pgcEqtls$flipSlopes = 0
pgcEqtls$flipSlopes[pgcEqtls$snpCountedAllele == pgcEqtls$riskAllele] = 1
pgcEqtls$flipSlopes[pgcEqtls$snpRefAllele == pgcEqtls$riskAllele] = -1

## CNVs
pgcEqtls$flipSlopes[nchar(pgcEqtls$snpCountedAllele) > nchar(pgcEqtls$snpRefAllele) &
	grepl("I", pgcEqtls$riskAllele)] = 1
pgcEqtls$flipSlopes[nchar(pgcEqtls$snpCountedAllele) < nchar(pgcEqtls$snpRefAllele) &
	grepl("I", pgcEqtls$riskAllele)] = -1
pgcEqtls$flipSlopes[nchar(pgcEqtls$snpCountedAllele) < nchar(pgcEqtls$snpRefAllele) &
	grepl("D", pgcEqtls$riskAllele)] = 1
pgcEqtls$flipSlopes[nchar(pgcEqtls$snpCountedAllele) > nchar(pgcEqtls$snpRefAllele) &
	grepl("D", pgcEqtls$riskAllele)] = -1
table(pgcEqtls$flipSlopes)

## look manually, 21 SNPs
as.data.frame(pgcEqtls[!duplicated(pgcEqtls$SNP) & pgcEqtls$flipSlopes==0,])
pgcEqtls = pgcEqtls[pgcEqtls$flipSlopes != 0,] # and drop

## flip
pgcEqtls$riskDir = ifelse(sign(pgcEqtls$beta*
	pgcEqtls$flipSlopes)==1,"UP","DOWN")

##############
### bring in case-control diffs
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/caseControl/rdas/all_de_features.rda")
outStatsEqtl = endoapply(outStats, function(x) x[names(x) %in% pgcEqtls$Feature])
dxStats = unlist(outStatsEqtl)
dxStats$Feature = ss(names(dxStats), "\\.",2)
dxStatsFilter = dxStats[match(pgcEqtls$Feature, dxStats$Feature)]

## add stats
pgcEqtls$dxDir = ifelse(sign(dxStatsFilter$log2FC_qsva) == 1, "UP", "DOWN")
pgcEqtls$dxPval =dxStatsFilter$pval_qsva
save(pgcEqtls, file = "rdas/pgcEqtls.rda")

#########################
### metrics for paper ###

## check table overall
dxVsRiskTab = table(pgcEqtls$dxDir,pgcEqtls$riskDir,dnn = c("dx","risk"))
getOR(dxVsRiskTab)
chisq.test(dxVsRiskTab)

## only take the index SNPs
pgcEqtlsFinal = pgcEqtls[!is.na(pgcEqtls$pgcFinalRank),]

## write out
firstSuppTable =as.data.frame(pgcEqtlsFinal)
firstSuppTable$WhichTx = sapply(firstSuppTable$WhichTx, paste, collapse=";")
firstSuppTable = firstSuppTable[,c(1:19, 22:27)]
colnames(firstSuppTable)[17] = "snpLoc"
write.csv(firstSuppTable, file="tables/suppTable_allPgcEqtlsToIndexSnps.csv",
	row.names=FALSE)
	
## check zhu
gg = c("SF3B1", "PCCB", "ENDOG", "ZDHHC5", "SNX19", "C12orf24",
	"ABCB9", "BAG5", "PSMA4","NMB", "NMRAL1", "KCTD13", "DUS2L",
	"C17ORF39", "GATAD2A" ,"IRF3")
gg[which(gg%in% pgcEqtlsFinal$Symbol)]
# how many loci?
length(unique(pgcEqtlsFinal$pgcFinalRank))
length(unique(pgcEqtlsFinal$pgcFinalRank[pgcEqtlsFinal$bonf < 0.05]))

## p-value
dxVsRiskTabFinal = table(pgcEqtlsFinal$dxDir,pgcEqtlsFinal$riskDir,dnn = c("dx","risk"))
getOR(dxVsRiskTabFinal)
chisq.test(dxVsRiskTabFinal)

# check bonf
dxVsRiskTabFinalBonf = table(pgcEqtlsFinal$dxDir[pgcEqtlsFinal$bonf < 0.05],
	pgcEqtlsFinal$riskDir[pgcEqtlsFinal$bonf < 0.05],dnn = c("dx","risk"))
getOR(dxVsRiskTabFinalBonf)
chisq.test(dxVsRiskTabFinalBonf)

#######################################
## only keep directional consistency
pgcEqtlsFinalDir = pgcEqtlsFinal[pgcEqtlsFinal$riskDir == pgcEqtlsFinal$dxDir,]
pgcEqtlsFinalDir$Type =  factor(as.character(pgcEqtlsFinalDir$Type),
		levels = c("Gene", "Exon", 	"Junction", "Transcript", "ER"), 
		labels = c("G","E","J","T","R"))
		
length(unique(pgcEqtlsFinalDir$pgcFinalRank))
nrow(pgcEqtlsFinalDir)
length(unique(pgcEqtlsFinalDir$EnsemblID))

length(unique(pgcEqtlsFinalDir$pgcFinalRank[pgcEqtlsFinalDir$bonf < 0.05]))
nrow(pgcEqtlsFinalDir[pgcEqtlsFinalDir$bonf < 0.05 ,])
length(unique(pgcEqtlsFinalDir$EnsemblID[pgcEqtlsFinalDir$bonf < 0.05]))

#### collapse by locus
pgcSigList = split(pgcEqtlsFinalDir, pgcEqtlsFinalDir$pgcFinalRank)
pgcSigList = lapply(pgcSigList, function(x) {
	x$pvalRatio = -log10(x$pvalue)/-log10(x$pvalue[1])
	x
})
pgcSigList = lapply(pgcSigList, function(x) 
	x[x$pvalRatio > 0.5,c(1:19,22,23,25:27)])

## write table out
pgcSigFiltered = do.call("rbind", pgcSigList)
nrow(pgcSigFiltered)

pgcSigFilteredDf = as.data.frame(pgcSigFiltered)
pgcSigFilteredDf$WhichTx = sapply(pgcSigFilteredDf$WhichTx, paste, collapse=";")
colnames(pgcSigFilteredDf)[17] = "snpLoc"

write.csv(pgcSigFilteredDf, row.names=FALSE,
	file = "tables/supp_table_pgcTop128_dirConsFiltered.csv")

#### collapse by locus
pgcListSym = lapply(pgcSigList, function(x) {
	x[which(!duplicated(x[,c("Symbol", "Type")]) & x$Symbol != ""),] 
})

table(sapply(pgcListSym, function(x) 
	length(unique(x$Symbol[!is.na(x$Symbol)]))))
table(sapply(pgcListSym, function(x) 
	length(unique(x$EnsemblID[!is.na(x$EnsemblID) & !is.na(x$Symbol)]))))
	
outList = lapply(pgcListSym, function(x) {
	risk = x$riskDir[!duplicated(x$Symbol)]
	
	types  = sapply(split(as.character(x$Type), 
		factor(x$Symbol,unique(x$Symbol))), 
			paste, collapse=",")
	pv = x$pvalue[!duplicated(x$Symbol)]

	sym = names(types) 

	x$newTx = x$WhichTx
	x$newTx[x$Type == "G"] = ""
	tx = sapply(split(x$newTx, 
		factor(x$Symbol,unique(x$Symbol))), function(y) 
			paste(unlist(y[y!=""]), collapse=","))
	
	data.frame(Direction = risk, Genes = sym,
		Types = types, BestPval = pv, Tx = tx,
			stringsAsFactors = FALSE)
})

outList = outList[sapply(outList,ncol) ==5]

outTable = do.call("rbind", outList)
outTable$pgcRank = ss(rownames(outTable), "\\.")

outTable = outTable[,c(6,1:5)]

write.csv(outTable, file="tables/pgcEqtlTable_forPaper_withTx.csv",
	row.names=FALSE)

tab = elementLengths(strsplit(outTable$Tx, ","))
sum(tab==1)/sum(tab>0)
###################
####### plot ######
###################

# best eqtl per gene/locus
bestEqtl = pgcSigFilteredDf[
	!duplicated(pgcSigFilteredDf[,c("pgcFinalRank","EnsemblID")]),] 

### make all boxplots
cleanExprs = rbind(cleanGeneSub, cleanExonSub, cleanJxnSub, # these are log2
	cleanTxSub, cleanErSub)
cleanExprs = cleanExprs[bestEqtl$Feature,]

theSnps = snpSub[bestEqtl$SNP,]
theSnpMap = snpMapSub[match(bestEqtl$SNP, snpMapSub$SNP),]

bestEqtl$fullType = factor(as.character(bestEqtl$Type), 
	levels = levels(bestEqtl$Type), 
	labels = c("Gene","Exon", "Junction", "Transcript", "Region"))

mainTxt = paste0("PGC Rank: ", bestEqtl$pgcFinalRank, " (",
	bestEqtl$snpRsNum, ")\n",
	bestEqtl$fullType, " - ", bestEqtl$Symbol, 
	" - ",	bestEqtl$Class)
ll = ifelse(bestEqtl$beta> 0, "topleft", "topright")
	
pdf("plots/eqtl_gwas_boxplots.pdf",useDingbats=FALSE)
par(mar=c(3,6,6,2), cex.axis=2, cex.lab=2, cex.main=1.7)
for(i in 1:nrow(bestEqtl)) {
	s = as.numeric(theSnps[i,])
	s = gsub(0, paste0(bestEqtl$snpRefAllele[i],
		bestEqtl$snpRefAllele[i]), s)
	s = gsub(1, paste0(bestEqtl$snpRefAllele[i], 
		bestEqtl$snpCountedAllele[i]), s)
	s = gsub(2, paste0(bestEqtl$snpCountedAllele[i], 
		bestEqtl$snpCountedAllele[i]), s)
	s = factor(s, levels = c(paste0(bestEqtl$snpRefAllele[i],
		bestEqtl$snpRefAllele[i]),paste0(bestEqtl$snpRefAllele[i], 
		bestEqtl$snpCountedAllele[i]), paste0(bestEqtl$snpCountedAllele[i], 
		bestEqtl$snpCountedAllele[i])))

	boxplot(cleanExprs[i,] ~ s, outline=FALSE,
		ylim= range(cleanExprs[i,]), 
		ylab= "Adjusted log2 Expression",
		xlab="", main = mainTxt[i])
	points(cleanExprs[i,] ~ jitter(as.numeric(s),amount=0.1),
		pch = 21, bg="grey",cex=1.4)
	legend(ll[i], paste0("p=",signif(bestEqtl$pvalue[i],
		3)),cex=1.8)
}
dev.off()

#### CASE CONTROL DIFFS
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/caseControl/rdas/degradation_mat_LIBD_polyA.rda")
rownames(pdSub) = pdSub$RNum
pdSub2 = pdSub[colnames(degCovAdj),]
allExprs = rbind(geneRpkmSub, exonRpkmSub, jRpkmSub,
	tFpkmSub, regionMatSub)
allExprs = allExprs[bestEqtl$Feature,colnames(degCovAdj)]

mod = model.matrix(~ Dx + Age + Sex + 
	snpPC1 + snpPC5 + snpPC6 + snpPC9 + snpPC10 +
	mitoRate + RIN + totalAssignedGene, data=pdSub2)
rownames(mod) = pdSub2$RNum

qSVs = prcomp(t(log2(degCovAdj+1)))$x[,1:12]

cleanExprsDx = cleaningY(allExprs, cbind(mod,qSVs), P=2)

library(RColorBrewer)
pdf("plots/eqtl_gwas_caseControl.pdf",useDingbats=FALSE)
par(mar=c(3,6,6,2), cex.axis=2, cex.lab=2, cex.main=1.7)
palette(brewer.pal(6,"Set1"))
for(i in 1:nrow(bestEqtl)) {
	boxplot(cleanExprs[i,] ~ pdSub$Dx, outline=FALSE,
		ylim= range(cleanExprs[i,]), 
		ylab= "Adjusted log2 Expression",
		xlab="", main = mainTxt[i])
	points(cleanExprs[i,] ~ jitter(as.numeric(factor(pdSub$Dx)),amount=0.1),
		pch = 21, bg=factor(pdSub$Dx),cex=1.4)
}
dev.off()


############## other analyses ##############
library(GenomicRanges)

## isoform switch + devel stats ##
pgcSigFilteredDf = read.csv("tables/supp_table_pgcTop128_dirConsFiltered.csv",
	as.is=TRUE)
	
## devel stats
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/isoform_switch_devel_byFeature.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/devStats_controlSamples.rda")
statListPgc = endoapply(statList, function(x) x[names(x) %in% pgcSigFilteredDf$Feature])
devStatsPgc = unlist(statListPgc)
names(devStatsPgc) = ss(names(devStatsPgc),"\\.",2)
devStatsPgc = devStatsPgc[pgcSigFilteredDf$Feature]

## tabulate
mean(devStatsPgc$p_bonf < 0.05) # 70%
sum(sapply(statList, function(x) sum(x$p_bonf < 0.05,na.rm=TRUE)))/sum(elementLengths(statList)) # 16%
devTable = table(pgcSigFilteredDf$riskDir[devStatsPgc$p_bonf < 0.05], 
	ifelse(devStatsPgc$ageCorr[devStatsPgc$p_bonf < 0.05] < 0, "Fetal","Postnatal"))
chisq.test(devTable)

## check isoswitch
isoGenes = unique(unlist(lapply(switchList, rownames)))
isoGenes = isoGenes[!is.na(isoGenes) & !grepl("-", isoGenes)]
table(pgcSigFilteredDf$EnsemblID %in% isoGenes)

## whats bg?
bgEns = unique(unlist(lapply(statList,
	function(x) x$EnsemblGeneID[which(x$p_bonf < 0.05)])))
bgEns = bgEns[!is.na(bgEns) & !grepl("-", bgEns)]

topleft = sum(pgcSigFilteredDf$EnsemblID %in% isoGenes)
topright = sum(!pgcSigFilteredDf$EnsemblID %in% isoGenes)
bottomleft = length(isoGenes) - topleft
bottomright = length(bgEns) - topleft - topright - bottomleft
isoTab= matrix(c(topleft, topright, bottomleft, bottomright),
	nrow = 2, byrow = TRUE,dimnames=list(c("pgcEqtl","noPgcEql"), 
		c("isoSwitch", "noSwitch")))
chisq.test(isoTab)$p.value
getOR(isoTab)

####################################
#### transcript plot for GRanges DER
library(GenomicFeatures)
ensTxDb = loadDb("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/ensembl_v75_txdb.sqlite")
seqlevels(ensTxDb,force=TRUE) = c(1:22,"X","Y","MT")
seqlevels(ensTxDb) = paste0("chr", c(1:22,"X","Y","MT"))

# DER for figure
gr = GRanges("chr2", IRanges(198350929,198376874))

pdf("plots/HSPD1_structure.pdf",h=3.5,w=7)
par(mar=c(5,5,1,1))
plotTranscripts(gr, ensTxDb)
text(x=(198375211+198375272)/2, y=1, "*", font=2,cex=2)
dev.off()

# jxn for figure
gr2 = GRanges("chr7", IRanges(111016890,111152010))

pdf("plots/IMMP2L_structure.pdf",h=3.5,w=7)
par(mar=c(5,5,1,1))
plotTranscripts(gr2, ensTxDb)
rect(xleft = 111057399, xright = 111127293,
	ytop = 0.5, ybottom = -0.5,col="black")
dev.off()