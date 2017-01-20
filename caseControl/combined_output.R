###
library(GenomicRanges)
library(RColorBrewer)

source("../eqtl_functions.R")

###########################
##### load data ##########
##########################

##### load CMC data
load("rdas/DE_statistics_adjAndQsva_CMC.rda")
outGeneCmc = outGene
outExonCmc = outExon
outJxnCmc = outJxn
outTxCmc = outTx
outErCmc = outEr

## load libd data
load("rdas/DE_statistics_adjAndQsva.rda")

###########################
########## add additional annotation ######
###########################

### add entrez ID
outJxn$EntrezID = outGene$EntrezID[match(outJxn$newGeneID, names(outGene))]
colnames(mcols(outJxn))[13] = "Symbol"
outTx$EntrezID = outGene$EntrezID[match(outTx$EnsemblGeneID, names(outGene))]
colnames(mcols(outTx))[3] = "Symbol"

dA = distanceToNearest(outEr, outGene)
outEr$EntrezID = outGene$EntrezID[subjectHits(dA)]
outEr$EntrezID[mcols(dA)$distance > 100] = NA
outEr$EnsemblGeneID[queryHits(dA) ] = names(outGene)[subjectHits(dA)]
outEr$EnsemblGeneID[mcols(dA)$distance > 100] = NA
colnames(mcols(outEr))[7] = "Symbol"

## clean up labels
colnames(mcols(outExon))[1] = "EnsemblGeneID"
colnames(mcols(outTx))[3] = "Symbol"
colnames(mcols(outEr))[7] = "Symbol"
colnames(mcols(outJxn))[12:13] = c("EnsemblGeneID","Symbol")
outGene$EnsemblGeneID = names(outGene)

############################
## add cmc stats ###########
#############################

## gene
outGene$CMC_log2FC_adj = outGeneCmc$log2FC_adj
outGene$CMC_tstat_adj = outGeneCmc$tstat_adj
outGene$CMC_pval_adj = outGeneCmc$pval_adj
outGene$CMC_log2FC_qsva = outGeneCmc$log2FC_qsva
outGene$CMC_tstat_qsva = outGeneCmc$tstat_qsva
outGene$CMC_pval_qsva = outGeneCmc$pval_qsva
outGene$CMC_log2FC_pca = outGeneCmc$log2FC_pca
outGene$CMC_tstat_pca = outGeneCmc$tstat_pca
outGene$CMC_pval_pca = outGeneCmc$pval_pca

# exon
outExon$CMC_log2FC_adj = outExonCmc$log2FC_adj
outExon$CMC_tstat_adj = outExonCmc$tstat_adj
outExon$CMC_pval_adj = outExonCmc$pval_adj
outExon$CMC_log2FC_qsva = outExonCmc$log2FC_qsva
outExon$CMC_tstat_qsva = outExonCmc$tstat_qsva
outExon$CMC_pval_qsva = outExonCmc$pval_qsva
outExon$CMC_log2FC_pca = outExonCmc$log2FC_pca
outExon$CMC_tstat_pca = outExonCmc$tstat_pca
outExon$CMC_pval_pca = outExonCmc$pval_pca

# junction
matchJxn = findOverlaps(outJxn, outJxnCmc, 
	ignore.strand=TRUE, type="equal")
outJxn$CMC_log2FC_adj = outJxn$CMC_pval_adj = NA
outJxn$CMC_log2FC_adj[queryHits(matchJxn)] = 
	outJxnCmc$log2FC_adj[subjectHits(matchJxn)]
outJxn$CMC_tstat_adj[queryHits(matchJxn)] = 
	outJxnCmc$tstat_adj[subjectHits(matchJxn)]
outJxn$CMC_pval_adj[queryHits(matchJxn)] = 
	outJxnCmc$pval_adj[subjectHits(matchJxn)]

outJxn$CMC_log2FC_qsva = outJxn$CMC_pval_qsva = NA
outJxn$CMC_log2FC_qsva[queryHits(matchJxn)] = 
	outJxnCmc$log2FC_qsva[subjectHits(matchJxn)]
outJxn$CMC_pval_qsva[queryHits(matchJxn)] = 
	outJxnCmc$pval_qsva[subjectHits(matchJxn)]
outJxn$CMC_tstat_qsva[queryHits(matchJxn)] = 
	outJxnCmc$tstat_qsva[subjectHits(matchJxn)]
outJxn$CMC_log2FC_pca = outJxn$CMC_pval_pca = NA
outJxn$CMC_log2FC_pca[queryHits(matchJxn)] = 
	outJxnCmc$log2FC_pca[subjectHits(matchJxn)]
outJxn$CMC_pval_pca[queryHits(matchJxn)] = 
	outJxnCmc$pval_pca[subjectHits(matchJxn)]
outJxn$CMC_tstat_pca[queryHits(matchJxn)] = 
	outJxnCmc$tstat_pca[subjectHits(matchJxn)]

## transcript
outTx$CMC_log2FC_adj = outTxCmc$log2FC_adj
outTx$CMC_tstat_adj = outTxCmc$tstat_adj
outTx$CMC_pval_adj = outTxCmc$pval_adj
outTx$CMC_log2FC_qsva = outTxCmc$log2FC_qsva
outTx$CMC_tstat_qsva = outTxCmc$tstat_qsva
outTx$CMC_pval_qsva = outTxCmc$pval_qsva
outTx$CMC_log2FC_pca = outTxCmc$log2FC_pca
outTx$CMC_tstat_pca = outTxCmc$tstat_pca
outTx$CMC_pval_pca = outTxCmc$pval_pca

## expressed region
outEr$CMC_log2FC_adj = outErCmc$log2FC_adj
outEr$CMC_tstat_adj = outErCmc$tstat_adj
outEr$CMC_pval_adj = outErCmc$pval_adj
outEr$CMC_log2FC_qsva = outErCmc$log2FC_qsva
outEr$CMC_tstat_qsva = outErCmc$tstat_qsva
outEr$CMC_pval_qsva = outErCmc$pval_qsva
outEr$CMC_log2FC_pca = outErCmc$log2FC_pca
outEr$CMC_tstat_pca = outErCmc$tstat_pca
outEr$CMC_pval_pca = outErCmc$pval_pca

### merge
n = c("Symbol", "EntrezID", "EnsemblGeneID", "meanExprsAdult", "isExp", 
	"code", "log2FC_adj", "pval_adj", "tstat_adj", "fdr_adj", 
	"log2FC_qsva" , "pval_qsva", "tstat_qsva", "fdr_qsva",
	"log2FC_pca" , "pval_pca", "tstat_pca", "fdr_pca",
	"CMC_log2FC_adj", "CMC_tstat_adj", "CMC_pval_adj", 
	"CMC_log2FC_qsva", "CMC_tstat_qsva", "CMC_pval_qsva",
	"CMC_log2FC_pca", "CMC_tstat_pca", "CMC_pval_pca")

outStats = c(Gene = outGene[,n], Exon = outExon[,n],
	Junction = outJxn[,n], Transcript = outTx[,n], ER = outEr[,n])
outStats = GRangesList(outStats)

save(outStats, file="rdas/all_de_features.rda")

### only those expressed
outStatsExprs = endoapply(outStats, function(x) x[x$isExp])
sapply(outStatsExprs,length)
save(outStatsExprs, file="rdas/expressed_de_features.rda")

#################
## write table ##
outTable = unlist(outStatsExprs)
outTable$Type =  ss(names(outTable), "\\.", 1)
outTable$Feature =  ss(names(outTable), "\\.", 2)

outTable = as.data.frame(outTable)
outTable$Coordinates = paste0(outTable$seqnames, 
	":", outTable$start, "-", outTable$end, "(",
	outTable$strand, ")")
outTable$SignifAndRep = outTable$fdr_qsva < 0.1 & 
		sign(outTable$CMC_log2FC_qsva) == 
			sign(outTable$log2FC_qsva) & 	
				outTable$CMC_pval_qsva < 0.05

outTable = outTable[,c(34, 33, 6:8, 11, 9, 35, 36, 12:32)]
rownames(outTable) = NULL
names(outTable)[6] = "AnnoClass"
outTable = outTable[,!grepl("tstat", colnames(outTable))]

write.table(outTable, row.names=FALSE, sep = "\t",
	file=gzfile("tables/suppTable_allExprs_DEwithRep.tsv.gz"))
	
## filter to signif
outTableSig = outTable[outTable$SignifAndRep,]
outTableSig$SignifAndRep = NULL
outTableSig = outTableSig[order(outTableSig$pval_qsva),]
write.table(outTableSig, row.names=FALSE, sep = "\t",
	file=gzfile("tables/suppTable_allExprs_sigDEwithRep.tsv.gz"))
