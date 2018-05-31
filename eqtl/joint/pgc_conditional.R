##
library(jaffelab)
library(GenomicRanges)

## load data
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/cleaned_eqtl_data_n412_allStats.rda")

## load pgc hits
load("rdas/PGC_SZ_hits_allEqtl_subset_fdr.rda")

## fix some ERs
pgcEqtls$Symbol[is.na(pgcEqtls$EnsemblID)] = NA # some ERs

## only top hits
pgcEqtlsFinal = pgcEqtls[which(!is.na(pgcEqtls$pgcFinalRank)),]

## filter expression
yExprs = rbind(cleanGeneSub[rownames(cleanGeneSub) %in% pgcEqtlsFinal$Feature,],
	cleanExonSub[rownames(cleanExonSub) %in% pgcEqtlsFinal$Feature,],
	cleanJxnSub[rownames(cleanJxnSub) %in% pgcEqtlsFinal$Feature,],
	cleanTxSub[rownames(cleanTxSub) %in% pgcEqtlsFinal$Feature,],
	cleanErSub[rownames(cleanErSub) %in% pgcEqtlsFinal$Feature,])
yExprs = yExprs[match(pgcEqtlsFinal$Feature, rownames(yExprs)),]

s = snpSub[match(pgcEqtlsFinal$SNP, rownames(snpSub)),]

# ##### analyses by region, post hoc on clean to compare
# statPost = t(sapply(1:nrow(pgcEqtlsFinal), function(i) {
	# summary(lm(yExprs[i,] ~ s[i,]))$coef[2,-2]
# }))
# statPost = as.data.frame(statPost)
# colnames(statPost) = c("beta", "statistic", "pvalue")

# plot(statPost$statistic ~ pgcEqtlsFinal$statistic) # looks okay

# ###### split by region
# rIndexes = splitit(pgcEqtlsFinal$pgcFinalRank)
# sigIndexes = sapply(rIndexes, function(ii) {
	# cat(".")
	# best = ii[1]
	# rest = ii[-1]
	# theSnp = s[best,]
	# while(length(rest) > 1) { 
		# if(length(best) == 1) z = yExprs[best,] else z = t(yExprs[best,])
		# coefs = as.data.frame(t(sapply(rest, function(i) summary(lm(yExprs[i,] ~ 
			# theSnp + z))$coef[2,-2])))
		# colnames(coefs) = c("beta", "statistic", "pvalue")
		
		# if(sum(coefs$pvalue < 0.05) == 0) {
			# return(best) 
		# } else {
			# best = c(best, rest[which(coefs$pvalue < 0.05)[1]])
			# rest = rest[which(coefs$pvalue < 0.05)[-1]]
		# }
	# }
	
	# return(best)
# })
# indeptIndex = unlist(sigIndexes)

# pgcEqtlsFinal$condIndept = 0
# pgcEqtlsFinal$condIndept[indeptIndex] = 1

# ########################
# ## check results out ###
# ########################

# ### bring in case-control
# load("/users/ajaffe/Lieber/Projects/RNAseq/SzControl_DE_paper/rdas/all_de_features.rda")
# outStatsSub = unlist(endoapply(outStats, function(x) x[names(x) %in% pgcEqtlsFinal$Feature]))
# outStatsSub$Type = ss(names(outStatsSub), "\\.",1)
# outStatsSub$Feature = ss(names(outStatsSub), "\\.",2)
# outStatsMatch = outStatsSub[match(pgcEqtlsFinal$Feature, outStatsSub$Feature)]

# pgcEqtlsFinal$SzDir_qSVA = ifelse(outStatsMatch$log2FC_qsva > 0, "UP", "DOWN")
# pgcEqtlsFinal$SzPval_qSVA = outStatsMatch$pval_qsva

# chisq.test(table(pgcEqtlsFinal$SzDir_qSVA, pgcEqtlsFinal$riskDir))
# chisq.test(table(pgcEqtlsFinal$SzDir_qSVA[pgcEqtlsFinal$GTEx_MetaPval < 1e-8], 
	# pgcEqtlsFinal$riskDir[pgcEqtlsFinal$GTEx_MetaPval < 1e-8]))
# chisq.test(table(pgcEqtlsFinal$SzDir_qSVA[pgcEqtlsFinal$bonf < 0.05], 
	# pgcEqtlsFinal$riskDir[pgcEqtlsFinal$bonf < 0.05]))

# chisq.test(table(pgcEqtlsFinal$SzDir_qSVA[pgcEqtlsFinal$SzPval_qSVA < 0.1], 
	# pgcEqtlsFinal$riskDir[pgcEqtlsFinal$SzPval_qSVA < 0.1]))	
	
# ### clean up
# pgcEqtlsFinal$WhichTx[pgcEqtlsFinal$Type == "Gene"] = NA
# rownames(pgcEqtlsFinal) = NULL
# pgcEqtlsFinal$flipSlopes = NULL
# pgcEqtlsFinal$pgcRow = NULL

# ## write table
# pgcOut = as.data.frame(pgcEqtlsFinal)
# pgcOut$WhichTx = sapply(pgcOut$WhichTx, paste, collapse=";")
# pgcOut = pgcOut[order(pgcOut$pgcFinalRank, pgcOut$pvalue),]
# write.csv(pgcOut, file = "tables/pgcEqtls_fdr01_replStats_withIndept.csv",
	# row.names=FALSE)

# save(pgcEqtlsFinal, file = "rdas/PGC_SZ_indexSnps_annotated_conditionalEqtls.rda")

#### metrics
load("rdas/PGC_SZ_indexSnps_annotated_conditionalEqtls.rda")

pgcEqtlsIndept = pgcEqtlsFinal[pgcEqtlsFinal$condIndept == 1,] 

## fix as3mt
pgcEqtlsIndept$EnsemblID[pgcEqtlsIndept$Symbol == "C10orf32-ASMT"] = "ENSG00000270316"
pgcEqtlsIndept$EnsemblID[pgcEqtlsIndept$Symbol == "AS3MT-C10orf32-ASMT"] = "ENSG00000270316"
pgcEqtlsIndept$Symbol[pgcEqtlsIndept$Symbol == "C10orf32-ASMT"] = "AS3MT"
pgcEqtlsIndept$Symbol[pgcEqtlsIndept$Symbol == "AS3MT-C10orf32-ASMT"] = "AS3MT"

tableCompare = sapply(list(FDR = pgcEqtlsFinal, Indept= pgcEqtlsIndept), function(x) {
	ttNovel = table(x$pgcFinalRank, x$Class)	## novelty
	ttType = table(x$pgcFinalRank, x$Type)	## type

	data.frame(numGwas = length(unique(x$pgcFinalRank)),
		numGwasNoGene = sum(ttType[,"Gene"] == 0),
		numGwasNoTxOrGene = sum(rowSums(ttType[,c("Gene", "Transcript")]) == 0),
		numGwasUn = sum(rowSums(ttNovel[,-3]) > 0),
		numGwasOnlyUn = sum(ttNovel[,"InEns"] == 0),
		numGwasTx = sum(lengths(sapply(split(x$WhichTx, x$pgcFinalRank), 
			function(y) unlist(unique(y)))) <= 1))
})
tableCompare

## summary stats
dim(pgcEqtlsIndept)
table(pgcEqtlsIndept$Type)
length(unique(pgcEqtlsIndept$EnsemblID[!grepl("-", pgcEqtlsIndept$EnsemblID) & 
	!is.na(pgcEqtlsIndept$EnsemblID) ]))
length(unique(pgcEqtlsIndept$Symbol[!grepl("-", pgcEqtlsIndept$Symbol) & 
	!is.na(pgcEqtlsIndept$Symbol) ]))
	
# number genes
table(lengths(sapply(split(pgcEqtlsFinal, pgcEqtlsFinal$pgcFinalRank), function(x) 
	unique(x$EnsemblID[!is.na(x$EnsemblID) & !grepl("-", x$EnsemblID)]))))
table(lengths(sapply(split(pgcEqtlsIndept, pgcEqtlsIndept$pgcFinalRank), function(x) 
	unique(x$EnsemblID[!is.na(x$EnsemblID) & !grepl("-", x$EnsemblID)]))))


locusList = split(pgcEqtlsIndept, pgcEqtlsIndept$pgcFinalRank)

## number of genes
table(lengths(sapply(locusList, function(x) 
	unique(x$EnsemblID[!is.na(x$EnsemblID)]))))


# num tx
t(sapply(locusList, function(x) 
	c(length(unique(x$EnsemblID[!is.na(x$EnsemblID)])), 
		length(unique(unlist(x$WhichTx[!is.na(x$WhichTx) & x$WhichTx!=""]))))))
		
checkInd0 = which(lengths(sapply(locusList, function(x) 
	unique(x$EnsemblID[!is.na(x$EnsemblID)]))) == 0)
as.list(locusList[names(checkInd0)])

checkInd1 = which(lengths(sapply(locusList, function(x) 
	unique(x$EnsemblID[!is.na(x$EnsemblID)]))) == 1)
as.list(locusList[names(checkInd1)])

checkInd2 = which(lengths(sapply(locusList, function(x) 
	unique(x$EnsemblID[!is.na(x$EnsemblID)]))) == 2)
as.list(locusList[names(checkInd2)])

## one and two gene loci
oneOrTwo = sapply(locusList, function(x) 
	length(unique(x$EnsemblID)) %in% 1:2)
tmp = sapply(locusList[oneOrTwo], function(x) unique(x$EnsemblID))
outTable = data.frame(locus = rep(names(tmp), lengths(tmp)), EnsemblID = unlist(tmp),
	stringsAsFactors=FALSE)
outTable$Gene = pgcEqtlsIndept$Symbol[match(outTable$EnsemblID, pgcEqtlsIndept$EnsemblID)]
outTable$Gene[is.na(outTable$EnsemblID)] = "Intergenic"
outTable$EnsemblID[is.na(outTable$EnsemblID)] = "Intergenic"
outTable$Gene[outTable$Gene == ""] = c("ZSCAN26", "AC068831.1", "AC011816.1", "AC117382.2",
	"SNORA77", "AC090568.2", "AP003049.1")
outTable$SNP = pgcEqtlsIndept$snpRsNum[match(outTable$locus, pgcEqtlsIndept$pgcFinalRank)]
outTable = outTable[,c(1,4,3,2)]

write.csv(outTable, file="tables/table2_eqtl_pgc_oneOrTwo.csv",row.names=FALSE)
	
###############
# make plots ##
###############


####################################
#### transcript plot for GRanges DER
library(GenomicFeatures)
library(rtracklayer)

library("EnsDb.Hsapiens.v75")
ensTxDb = EnsDb.Hsapiens.v75
seqlevelsStyle(ensTxDb) <- "UCSC"


plotTranscripts = function(gr, txByExon) {
	require(derfinderPlot)
	require(GenomicFeatures)
	ov <- findOverlaps(gr, txByExon)
	txList <- split(txByExon[subjectHits(ov)], queryHits(ov))
	poly.data <- lapply(txList, derfinderPlot:::.plotData)[[1]]
	yrange <- range(poly.data$exon$y)
	xrange <- c(start(gr),end(gr))
    plot(0,0, type="n", xlim = xrange , ylim = yrange + c(-0.75, 0.75),
		yaxt="n", ylab="", xlab=as.character(seqnames(gr)),
		cex.axis = 1.5, cex.lab =1.5)
	derfinderPlot:::.plotPolygon(poly.data$exon, 'blue')
    derfinderPlot:::.plotPolygon(poly.data$intron, 'lightblue')
	yTx <- unique(yrange / abs(yrange))
	if(length(yTx) > 1) {
		axis(2, yTx, c('-', '+')[c(-1, 1) %in% yTx], 
			tick = FALSE, las = 1, cex.axis = 3)
		abline(h = 0, lty = 3)
	}
}
# TranscriptDb = makeTxDbFromGFF("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/Counts/Ensembl/Homo_sapiens.GRCh37.75.gtf",
	# format="gtf")
# saveDb(TranscriptDb, file="/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/ensembl_v75_txdb_update.sqlite")
ensTxDb = loadDb("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/ensembl_v75_txdb_update.sqlite")
ensTxDb = keepStandardChromosomes(ensTxDb, species = "Homo_sapiens")
ensTxDb = renameSeqlevels(ensTxDb, paste0("chr", seqlevels(ensTxDb)))

## extract tx
txByExon <- exonsBy(ensTxDb)

# DER for figure, ZSCAN23
pdf("plots/ZSCAN23_structure.pdf",h=3.5,w=7)
gr = GRanges("chr6", IRanges(28378778,28412830))
par(mar=c(1,5,1,1))
layout(matrix(c(1,1, 2,2,2), nc= 1))
meanCov = import("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/hub/hg19/Adult_coverage.bw",
	which = gr)
plot(log2(meanCov$score+1) ~ start(meanCov),type="l",xaxt="n", xlab="",
	ylab="Mean Cov",cex.axis=1.4,cex.lab=1.5)
abline(v=c(as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$snpRsNum=="rs1233578"], ":",2)[1],"-",1)),
	as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$snpRsNum=="rs1233578"],"-",2)[1], "\\("))),lty=2)
par(mar=c(5,5,1,1))
plotTranscripts(gr, txByExon)
abline(v=c(as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$snpRsNum=="rs1233578"], ":",2)[1],"-",1)),
	as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$snpRsNum=="rs1233578"],"-",2)[1], "\\("))),lty=2)
dev.off()


# jxn for figure
# chr5:137,837,170-138,284,761
gr2 = GRanges("chr5", IRanges(137837170,138280000))

pdf("plots/CTNNA1_structure.pdf",h=3.5,w=7)
par(mar=c(5,5,1,1))
plotTranscripts(gr2, txByExon)
rect(xleft = 137946779, xright = 137988752,
	ytop = 1, ybottom = -1,col="red")
dev.off()

# der for figure
# chr3:180,730,454-181,525,311
gr = GRanges("chr3", IRanges(180730454,181525311))
meanCov = import("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/hub/hg19/Adult_coverage.bw",
	which = gr)

pdf("plots/SOX2OT_structure.pdf",h=3.5,w=7)
par(mar=c(1,5,1,1))
layout(matrix(c(1,1, 2,2,2), nc= 1))
plot(log2(meanCov$score+1) ~ end(meanCov),type="l",xaxt="n", xlab="",
	ylab="Mean Cov",cex.axis=1.4,cex.lab=1.5, xlim = c(start(gr), end(gr)))
abline(v=c(as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$snpRsNum=="rs9841616"], ":",2)[1],"-",1)),
	as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$snpRsNum=="rs9841616"],"-",2)[1], "\\("))),lty=2)
par(mar=c(5,5,1,1))
plotTranscripts(gr, txByExon)
abline(v=c(as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$snpRsNum=="rs9841616"], ":",2)[1],"-",1)),
	as.numeric(ss(ss(pgcEqtlsIndept$Coordinates[pgcEqtlsIndept$snpRsNum=="rs9841616"],"-",2)[1], "\\("))),lty=2)
dev.off()


#######

	
## check zhu
gg = c("SF3B1", "PCCB", "ENDOG", "ZDHHC5", "SNX19", "C12orf24",
	"ABCB9", "BAG5", "PSMA4","NMB", "NMRAL1", "KCTD13", "DUS2L",
	"C17ORF39", "GATAD2A" ,"IRF3")
gg[which(gg%in% pgcEqtlsFinal$Symbol)]
gg[which( ! gg%in% pgcEqtlsFinal$Symbol)]


## devel stats
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/isoform_switch_devel_byFeature.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/devStats_controlSamples.rda")
statListPgc = endoapply(statList, function(x) x[names(x) %in% pgcEqtlsFinal$Feature])

devStatsPgc = unlist(statListPgc)
names(devStatsPgc) = ss(names(devStatsPgc),"\\.",2)
devStatsPgc = devStatsPgc[pgcEqtlsFinal$Feature]

## tabulate
mean(devStatsPgc$p_bonf < 0.05, na.rm=TRUE) # 58.7%
sum(sapply(statList, function(x) sum(x$p_bonf < 0.05,na.rm=TRUE)))/sum(lengths(statList)) # 16%
devTable = table(pgcEqtlsFinal$riskDir[devStatsPgc$p_bonf < 0.05], 
	ifelse(devStatsPgc$ageCorr[devStatsPgc$p_bonf < 0.05] < 0, "Fetal","Postnatal"))
chisq.test(devTable)
getOR(devTable)

#####
## check isoswitch
isoGenes = unique(unlist(lapply(switchList, rownames)))
isoGenes = isoGenes[!is.na(isoGenes) & !grepl("-", isoGenes)]
table(unique(pgcEqtlsFinal$EnsemblID) %in% isoGenes)

## whats bg?
bgEns = unique(unlist(lapply(statList,
	function(x) x$EnsemblGeneID[which(x$p_bonf < 0.05)])))
bgEns = bgEns[!is.na(bgEns) & !grepl("-", bgEns)]

# any
topleft = sum(unique(pgcEqtlsFinal$EnsemblID) %in% isoGenes)
topright = sum(! unique(pgcEqtlsFinal$EnsemblID) %in% isoGenes)
bottomleft = length(unique(isoGenes)) - topleft
bottomright = length(bgEns) - topleft - topright - bottomleft
isoTab= matrix(c(topleft, topright, bottomleft, bottomright),
	nrow = 2, byrow = TRUE,dimnames=list(c("pgcEqtl","noPgcEql"), 
		c("isoSwitch", "noSwitch")))
chisq.test(isoTab)$p.value
getOR(isoTab)

# cond indep
topleft = sum(unique(pgcEqtlsIndept$EnsemblID) %in% isoGenes)
topright = sum(! unique(pgcEqtlsIndept$EnsemblID) %in% isoGenes)
bottomleft = length(unique(isoGenes)) - topleft
bottomright = length(bgEns) - topleft - topright - bottomleft
isoTab_cond= matrix(c(topleft, topright, bottomleft, bottomright),
	nrow = 2, byrow = TRUE,dimnames=list(c("pgcEqtl","noPgcEql"), 
		c("isoSwitch", "noSwitch")))
chisq.test(isoTab_cond)$p.value
getOR(isoTab_cond)
prop.table(isoTab_cond,1)
