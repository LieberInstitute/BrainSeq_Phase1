####
###
library(GenomicRanges)
library(limma)
source("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl_functions.R")

#################################
#### LOAD DATA ##################
#################################

# load SNPs for annotation info
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/snpData_LIBD_szControls_cleaned.rda")
snpMap$CHR[snpMap$CHR == "23"] = "X"


# load expression annotation
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
exonMap$coord = with(exonMap, paste0(Chr, ":", Start, "-", End, "(", Strand, ")"))

# ERs + Tx
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/transcript/transcript_data_filtered_n495.rda")

oo = findOverlaps(regions, makeGRangesFromDataFrame(geneMap))
regions$EnsemblID = NA
regions$EnsemblID[queryHits(oo) ] = rownames(geneMap)[subjectHits(oo)]

################################
## bring back in eqtls #########
################################

###### GENE
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/annotated_gene_eqtl_szControl_cisOnly.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/gene_eqtl_szControl_cisOnly.rda")
geneEqtl = meGeneCis$cis$eqtls
geneEqtl = geneEqtl[geneEqtl$FDR < 0.01,]
geneEqtl$gene = as.character(geneEqtl$gene)
geneEqtl$snps = as.character(geneEqtl$snps)
geneEqtl$bonf = p.adjust(geneEqtl$pvalue, "bonf", 
	n = (sigGene$bonf[1]/sigGene$pvalue[1]))

###### EXON
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/annotated_exon_eqtl_szControl_cisOnly.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/exon_eqtl_szControl_cisOnly.rda")
exonEqtl = meExonCis$cis$eqtls
exonEqtl = exonEqtl[exonEqtl$FDR < 0.01,]
exonEqtl$gene = as.character(exonEqtl$gene)
exonEqtl$snps = as.character(exonEqtl$snps)
exonEqtl$bonf = p.adjust(exonEqtl$pvalue, "bonf", 
	n = (sigExon$bonf[1]/sigExon$pvalue[1]))

######### TX
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/annotated_transcript_eqtl_szControl_cisOnly.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/transcript_eqtl_szControl_cisOnly.rda")
transcriptEqtl = meTransCis$cis$eqtls
transcriptEqtl = transcriptEqtl[transcriptEqtl$FDR < 0.01,]
transcriptEqtl$gene = as.character(transcriptEqtl$gene)
transcriptEqtl$snps = as.character(transcriptEqtl$snps)
transcriptEqtl$bonf = p.adjust(transcriptEqtl$pvalue, "bonf", 
	n = (sigTrans$bonf[1]/sigTrans$pvalue[1]))
	
###### JUNCTION
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/annotated_junction_eqtl_szControl_cisOnly.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/junction_eqtl_szControl_cisOnly.rda")
junctionEqtl = meJxnCis$cis$eqtls
junctionEqtl = junctionEqtl[junctionEqtl$FDR < 0.01,]
junctionEqtl$gene = as.character(junctionEqtl$gene)
junctionEqtl$snps = as.character(junctionEqtl$snps)
junctionEqtl$bonf = p.adjust(junctionEqtl$pvalue, "bonf", 
	n = (sigJxn$bonf[1]/sigJxn$pvalue[1]))

###### DER
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/annotated_der_eqtl_szControl_cisOnly.rda")
load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/der_eqtl_szControl_cisOnly.rda")
derEqtl = meDerCis$cis$eqtls
derEqtl = derEqtl[derEqtl$FDR < 0.01,]
derEqtl$gene = as.character(derEqtl$gene)
derEqtl$snps = as.character(derEqtl$snps)
derEqtl$bonf = p.adjust(derEqtl$pvalue, "bonf", 
	n = (sigDer$bonf[1]/sigDer$pvalue[1]))

###########################
### add gene and transcript info
###############################

geneEqtl$Symbol = geneMap$Symbol[match(geneEqtl$gene, rownames(geneMap))]
geneEqtl$EnsemblID = as.character(geneEqtl$gene)
geneEqtl$Class = "InEns"

exonEqtl$Symbol = exonMap$Symbol[match(exonEqtl$gene, rownames(exonMap))]
exonEqtl$EnsemblID = exonMap$Geneid[match(exonEqtl$gene, rownames(exonMap))]
exonEqtl$Class = "InEns"

transcriptEqtl$Symbol = tMap$gene_name[match(transcriptEqtl$gene, names(tMap))]
transcriptEqtl$EnsemblID = tMap$EnsemblGeneID[match(transcriptEqtl$gene, names(tMap))]
code = tMap$class_code[match(transcriptEqtl$gene, names(tMap))]
transcriptEqtl$Class = "Novel"
transcriptEqtl$Class[code == "="] = "InEns"
transcriptEqtl$Class[code == "j"] = "ExonSkip"
transcriptEqtl$Class[code %in% c("c", "e", "o", "p", "x") ] = "AltStartEnd"

junctionEqtl$Symbol = jMap$newGeneSymbol[match(junctionEqtl$gene, names(jMap))]
junctionEqtl$EnsemblID = jMap$newGeneID[match(junctionEqtl$gene, names(jMap))]
junctionEqtl$Class = jMap$code[match(junctionEqtl$gene, names(jMap))]

derEqtl$Symbol = regions$nearestSymbol[match(derEqtl$gene, names(regions))]
derEqtl$EnsemblID = regions$EnsemblID[match(derEqtl$gene, names(regions))]
code2 = regions$annoClass[match(derEqtl$gene, names(regions))]
derEqtl$Class = "Novel"
derEqtl$Class[code2 == "strictExonic"] = "InEns"
derEqtl$Class[code2 %in% c("exonIntron", "extendUTR")] = "AltStartEnd"

###################
### try merging ###
###################

allEqtls = rbind(geneEqtl, exonEqtl, 
	transcriptEqtl, junctionEqtl, derEqtl)
allEqtls$snpRsNum = snpMap$name[match(allEqtls$snps, snpMap$SNP)]
allEqtls$Type = rep(c("Gene","Exon", "Transcript", "Junction", "ER"),
	c(nrow(geneEqtl), nrow(exonEqtl), nrow(transcriptEqtl), 
		nrow(junctionEqtl), nrow(derEqtl)))

## make smaller 
allEqtls = DataFrame(allEqtls)
allEqtls$Class = Rle(allEqtls$Class)
allEqtls$Type = Rle(allEqtls$Type)
allEqtls$Type = Rle(allEqtls$Type)
allEqtls = allEqtls[order(allEqtls$pvalue),]
names(allEqtls)[1:2] = c("SNP", "Feature")

#############################
### get transcript info #####
xx=load("/users/ajaffe/Lieber/Projects/RNAseq/ensembl_V75_feature_to_Tx.rda")

## add annotation from tx
theTxs = CharacterList(as.list(tMap$nearest_ref))
names(theTxs) = names(tMap)
allTx = c(allTx, theTxs)

mmTx = match(allEqtls$Feature, names(allTx))
tx = CharacterList(vector("list", nrow(allEqtls)))
tx[!is.na(mmTx)] = allTx[mmTx[!is.na(mmTx)]]
allEqtls$NumTxEqtl = elementNROWS(tx)
allEqtls$WhichTx = tx

mmTxGene = match(allEqtls$EnsemblID, names(allTx))
txGene = CharacterList(vector("list", nrow(allEqtls)))
txGene[!is.na(mmTxGene)] = allTx[mmTxGene[!is.na(mmTxGene)]]
allEqtls$NumTxGene = elementNROWS(txGene)

allEqtls$NumTxEqtl = Rle(allEqtls$NumTxEqtl)
allEqtls$NumTxGene = Rle(allEqtls$NumTxGene)

## add coordinates
geneMap$coord = with(geneMap, paste0(Chr, ":", Start, "-", End, "(", Strand, ")"))
tMap$coord = paste0(seqnames(tMap), ":", start(tMap), "-", end(tMap),
	"(", strand(tMap), ")")
regions$coord = paste0(seqnames(regions), ":", start(regions), 
	"-", end(regions),	"(", strand(regions), ")")

coords = c(geneMap$coord, exonMap$coord, tMap$coord, 
	names(jMap), regions$coord)
names(coords) = c(rownames(geneMap), rownames(exonMap),
		names(tMap), names(jMap), names(regions))
allEqtls$Coordinates = coords[match(allEqtls$Feature, names(coords))]

### add distance
starts = c(geneMap$Start, exonMap$Start, start(tMap), 
	start(jMap), start(regions))
ends = c(geneMap$End, exonMap$End, end(tMap), 
	end(jMap), end(regions))
names(starts) = names(ends) = c(rownames(geneMap), rownames(exonMap),
		names(tMap), names(jMap), names(regions))
snpPos = snpMap$POS[match(allEqtls$SNP, snpMap$SNP)]	
tmp = cbind(snpPos - starts[match(allEqtls$Feature, names(coords))],
		snpPos - ends[match(allEqtls$Feature, names(coords))])
mins = matrixStats::rowMins(abs(tmp))
allEqtls$snpDistToFeature = sign(rowMeans(tmp))*mins

rownames(allEqtls) = paste0(allEqtls$SNP, ".", allEqtls$Feature)
allEqtls = allEqtls[order(allEqtls$pvalue),]

#### CMC replication #####
load("/dcl01/lieber/ajaffe/PublicData/CMC/replication_eqtl_stats.rda")
mm = match(rownames(allEqtls), rownames(allEqtlRep))

allEqtls$CMC_statistic = allEqtlRep$statistic[mm]
allEqtls$CMC_pvalue = allEqtlRep$pvalue[mm]
allEqtls$CMC_beta = allEqtlRep$beta[mm]
save(allEqtls, compress=TRUE,
	file="/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/allEqtls.rda")
	
#### GTEx replication ####
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/gtex_brainEqtl_replication_stats.rda")
mmGtex = match(rownames(allEqtls), rownames(allEqtlRepGtex))
allEqtlRepGtex = allEqtlRepGtex[mmGtex,]
rownames(allEqtlRepGtex) = rownames(allEqtls)
gtexTstats = allEqtlRepGtex[,grep("statistic", colnames(allEqtlRepGtex))]
colnames(gtexTstats) = paste0("GTEx_", colnames(gtexTstats))
allEqtls = cbind(allEqtls, gtexTstats)

save(allEqtls, compress=TRUE,
	file="/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/allEqtls_withGtex.rda")
	
#########################
######### DATABASE ######
#########################


# filter for Age
pd$usedInEqtlDiscovery = pd$Age > 13

geneExprs = as.matrix(log2(geneRpkm+1))
exonExprs = as.matrix(log2(exonRpkm+1))
jxnExprs = as.matrix(log2(jRpkm+1))
txExprs = as.matrix(log2(tFpkm+1))
regionExprs = as.matrix(log2(regionMat+1))

## residualize
mod = model.matrix(~snpPC1 + snpPC2 + snpPC3 + Sex + Dx,
	data=pd[pd$usedInEqtlDiscovery,])

cleanGeneSub = cleaningY(geneExprs[rownames(geneRpkm) %in% allEqtls$Feature,
	pd$usedInEqtlDiscovery], cbind(mod, pcsGene),P=1)
geneMapSub = geneMap[rownames(geneRpkm) %in% allEqtls$Feature,]

cleanExonSub = cleaningY(exonExprs[rownames(exonRpkm) %in% allEqtls$Feature,
	pd$usedInEqtlDiscovery], cbind(mod, pcsExon),P=1)
exonMapSub = exonMap[rownames(exonRpkm) %in% allEqtls$Feature,]

cleanJxnSub = cleaningY(jxnExprs[rownames(jRpkm) %in% allEqtls$Feature,
	pd$usedInEqtlDiscovery], cbind(mod, pcsJxn),P=1)
jMapSub = jMap[rownames(jRpkm) %in% allEqtls$Feature,]

cleanTxSub = cleaningY(txExprs[rownames(tFpkm) %in% allEqtls$Feature,
	pd$usedInEqtlDiscovery], cbind(mod, pcsTrans),P=1)
tMapSub = tMap[rownames(tFpkm) %in% allEqtls$Feature,]

cleanErSub = cleaningY(regionExprs[rownames(regionMat) %in% allEqtls$Feature,
	pd$usedInEqtlDiscovery], cbind(mod, pcsDer),P=1)
regionsSub = regions[rownames(regionMat) %in% allEqtls$Feature,]

snpSub = as.matrix(snp[rownames(snp) %in% allEqtls$SNP,pd$usedInEqtlDiscovery])
snpMapSub = snpMap[rownames(snp) %in% allEqtls$SNP,]
pdSub = pd[pd$usedInEqtlDiscovery,]

save(pdSub, snpSub, snpMapSub, cleanGeneSub, geneMapSub,
	cleanExonSub, exonMapSub,cleanJxnSub, jMapSub,
	cleanTxSub, tMapSub, cleanErSub, regionsSub, compress=TRUE,
	file = "/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/cleaned_eqtl_data_n412_allStats.rda")
	
		
#################################
### export database files #######
eqtlStatsOut = allEqtls[,c(1:7,18)]
eqtlStatsOut = as.data.frame(eqtlStatsOut)
write.table(eqtlStatsOut, row.names=FALSE, sep = "\t",
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/AllEqtls/eqtl_stats_allStats.txt.gz"))

## expression anno
exprsAnnoOut = allEqtls[,c(2,9:11,13:17)]
exprsAnnoOut = exprsAnnoOut[!duplicated(exprsAnnoOut$Feature),]
exprsAnnoOut$chr = ss(exprsAnnoOut$Coordinates, ":")
exprsAnnoOut$start = as.integer(ss(ss(exprsAnnoOut$Coordinates, ":",2),"-"))
exprsAnnoOut$end = as.integer(ss(ss(ss(exprsAnnoOut$Coordinates, ":",2),"-",2),"\\("))
exprsAnnoOut = as.data.frame(exprsAnnoOut)
exprsAnnoOut$WhichTx = sapply(exprsAnnoOut$WhichTx, paste, collapse=";")
write.table(exprsAnnoOut, row.names=FALSE, sep = "\t",
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/AllEqtls/eqtl_exprsAnno_allStats.txt.gz"))

## SNP anno	
snpAnnoOut = snpMap
colnames(snpAnnoOut)[c(1,3,4:5)] = c("chr","pos","Counted","Ref")
snpAnnoOut$chr = paste0("chr", snpAnnoOut$chr)
snpAnnoOut = snpAnnoOut[,c(2,6,1,3:5,7)]

snpAnnoOut = as.data.frame(snpAnnoOut)
write.table(snpAnnoOut, row.names=FALSE, sep = "\t",
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/AllEqtls/snpAnno_allConsidered.txt.gz"))

## expression data
cleanExprs = rbind(cleanGeneSub, cleanExonSub, cleanJxnSub,	cleanTxSub, cleanErSub)
write.table(cleanExprs, sep = "\t",
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/AllEqtls/eqtl_log2cleanExprsData_allStats.txt.gz"))

## snp data
write.table(snpSub, sep = "\t",
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/AllEqtls/eqtl_snpData_allStats.txt.gz"))

## replication stats
repStatsOut = allEqtls[,c(1:2, 19:34)]
write.table(repStatsOut, sep = "\t",row.names=FALSE,
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/AllEqtls/replication_stats_allStats.txt.gz"))

	
################
#### OLD  ######
################


allEqtlRep = allEqtlRep[rownames(allEqtls),]
colnames(allEqtlRep) = paste0("CMC_", colnames(allEqtlRep))
allEqtls = cbind(allEqtls, allEqtlRep[,1:4])


load("/users/ajaffe/Lieber/Projects/RNAseq/DLPFC_eQTL_paper/joint/rdas/CMC_replication_stats_overall.rda")

sigGene$cmcTested = geneStats$tested[match(sigGene$gene, rownames(geneStats))]
sigGene$cmcBeta = geneStats$beta[match(sigGene$gene, rownames(geneStats))]
sigGene$cmcPvalue = geneStats$pval[match(sigGene$gene, rownames(geneStats))]

sigExon$cmcTested = exonStats$tested[match(sigExon$exon, rownames(exonStats))]
sigExon$cmcBeta = exonStats$beta[match(sigExon$exon, rownames(exonStats))]
sigExon$cmcPvalue = exonStats$pval[match(sigExon$exon, rownames(exonStats))]

sigTrans$cmcTested = transStats$tested[match(sigTrans$tx, rownames(transStats))]
sigTrans$cmcBeta = transStats$beta[match(sigTrans$tx, rownames(transStats))]
sigTrans$cmcPvalue = transStats$pval[match(sigTrans$tx, rownames(transStats))]

sigJxn$cmcTested = jxnStats$tested[match(sigJxn$jxn, rownames(jxnStats))]
sigJxn$cmcBeta = jxnStats$beta[match(sigJxn$jxn, rownames(jxnStats))]
sigJxn$cmcPvalue = jxnStats$pval[match(sigJxn$jxn, rownames(jxnStats))]

sigDer$cmcTested = derStats$tested[match(sigDer$der, rownames(derStats))]
sigDer$cmcBeta = derStats$beta[match(sigDer$der, rownames(derStats))]
sigDer$cmcPvalue = derStats$pval[match(sigDer$der, rownames(derStats))]
	
### add 
sigGene$cmcDir = geneStats$direction[match(sigGene$gene, rownames(geneStats))]
sigExon$cmcDir = exonStats$direction[match(sigExon$exon, rownames(exonStats))]
sigJxn$cmcDir= jxnStats$direction[match(sigJxn$jxn, rownames(jxnStats))]
sigTrans$cmcDir = transStats$direction[match(sigTrans$tx, rownames(transStats))]
sigDer$cmcDir = derStats$direction[match(sigDer$der, rownames(derStats))]
	
### add 
sigGene$cmcDirAndP01 = geneStats$dirAndSig[match(sigGene$gene, rownames(geneStats))]
sigExon$cmcDirAndP01 = exonStats$dirAndSig[match(sigExon$exon, rownames(exonStats))]
sigJxn$cmcDirAndP01 = jxnStats$dirAndSig[match(sigJxn$jxn, rownames(jxnStats))]
sigTrans$cmcDirAndP01 = transStats$dirAndSig[match(sigTrans$tx, rownames(transStats))]
sigDer$cmcDirAndP01 = derStats$dirAndSig[match(sigDer$der, rownames(derStats))]

### add 
sigGene$cmcDirAndPe5 = geneStats$dirAndSig5[match(sigGene$gene, rownames(geneStats))]
sigExon$cmcDirAndPe5 = exonStats$dirAndSig5[match(sigExon$exon, rownames(exonStats))]
sigJxn$cmcDirAndPe5 = jxnStats$dirAndSig5[match(sigJxn$jxn, rownames(jxnStats))]
sigTrans$cmcDirAndPe5 = transStats$dirAndSig5[match(sigTrans$tx, rownames(transStats))]
sigDer$cmcDirAndPe5 = derStats$dirAndSig5[match(sigDer$der, rownames(derStats))]
	
#### add GTEX #####
xx = load("/users/ajaffe/Lieber/Projects/RNAseq/GTEX/rdas/eqtl_replication_stats_byRegion.rda")

## just to cortex
geneStatsGtex = as.data.frame(geneStatsGtex[,
	grep("Frontal Cortex",colnames(geneStatsGtex))])
exonStatsGtex = as.data.frame(exonStatsGtex[,
	grep("Frontal Cortex",colnames(exonStatsGtex))])
jxnStatsGtex = as.data.frame(jxnStatsGtex[,
	grep("Frontal Cortex",colnames(jxnStatsGtex))])
txStatsGtex = as.data.frame(txStatsGtex[,
	grep("Frontal_Cortex",colnames(txStatsGtex))])
derStatsGtex = as.data.frame(derStatsGtex[,
	grep("Frontal Cortex",colnames(derStatsGtex))])

colnames(geneStatsGtex) = colnames(exonStatsGtex) = colnames(jxnStatsGtex) = c("Dir","DirAndSig05","DirAndSig0001") 
colnames(txStatsGtex) = colnames(derStatsGtex) = c("Dir","DirAndSig05","DirAndSig0001") 

### add 
sigGene$gtexPresent = sigGene$gene %in% rownames(geneStatsGtex)
sigGene$gtexDir = geneStatsGtex$Dir[match(sigGene$gene, rownames(geneStatsGtex))]
sigGene$gtexDirAnd05 = geneStatsGtex$DirAndSig05[match(sigGene$gene, rownames(geneStatsGtex))]

sigExon$gtexPresent = sigExon$exon %in% rownames(exonStatsGtex)
sigExon$gtexDir = exonStatsGtex$Dir[match(sigExon$exon, rownames(exonStatsGtex))]
sigExon$gtexDirAnd05 = exonStatsGtex$DirAndSig05[match(sigExon$exon, rownames(exonStatsGtex))]

sigJxn$gtexPresent = sigJxn$jxn %in% rownames(jxnStatsGtex)
sigJxn$gtexDir = jxnStatsGtex$Dir[match(sigJxn$jxn, rownames(jxnStatsGtex))]
sigJxn$gtexDirAnd05 = jxnStatsGtex$DirAndSig05[match(sigJxn$jxn, rownames(jxnStatsGtex))]

sigTrans$gtexPresent = sigTrans$tx %in%  rownames(txStatsGtex)
sigTrans$gtexDir = txStatsGtex$Dir[match(sigTrans$tx, rownames(txStatsGtex))]
sigTrans$gtexDirAnd05 = txStatsGtex$DirAndSig05[match(sigTrans$tx, rownames(txStatsGtex))]

sigDer$gtexPresent = sigDer$der %in% rownames(derStatsGtex)
sigDer$gtexDir = derStatsGtex$Dir[	match(sigDer$der, rownames(derStatsGtex))]
sigDer$gtexDirAnd05 = derStatsGtex$DirAndSig05[	match(sigDer$der, rownames(derStatsGtex))]

save(sigGene, sigExon, sigJxn, sigTrans, sigDer,
	file="rdas/annotated_and_testedForRep_eQTLs.rda")
	
#############################
## check replication ########

sigList = list(gene = sigGene, exon = sigExon,
	transcript = sigTrans, junction = sigJxn, der=  sigDer)
	
###########################
## numbers for table s8 ###

## CMC
cmcFdrRepTable = t(data.frame(numSigDisc = sapply(sigList,nrow),
	numPresent =sapply(sigList, function(x) sum(!is.na(x$cmcTested))),
	numTested =sapply(sigList, function(x) sum(x$cmcTested, na.rm=TRUE)),
	dirCons =sapply(sigList, function(x) sum(x$cmcDir, na.rm=TRUE)),
	dirAnd01 =sapply(sigList, function(x) sum(x$cmcDirAndP01, na.rm=TRUE)),
	dirAndPe5 =sapply(sigList, function(x) sum(x$cmcDirAndPe5, na.rm=TRUE))))

sigListBonf = lapply(sigList, function(x) x[x$bonf < 0.05,])
cmcBonfRepTable = t(data.frame(numSigDisc = sapply(sigListBonf,nrow),
	numPresent =sapply(sigListBonf, function(x) sum(!is.na(x$cmcTested))),
	numTested =sapply(sigListBonf, function(x) sum(x$cmcTested, na.rm=TRUE)),
	dirCons =sapply(sigListBonf, function(x) sum(x$cmcDir, na.rm=TRUE)),
	dirAnd01 =sapply(sigListBonf, function(x) sum(x$cmcDirAndP01, na.rm=TRUE)),
	dirAndPe5 =sapply(sigListBonf, function(x) sum(x$cmcDirAndPe5, na.rm=TRUE))))

## GTEx
gtexFdrRepTable = t(data.frame(numSigDisc = sapply(sigList,nrow),
	numPresent =sapply(sigList, function(x) sum(x$gtexPresent)),
	numTested =sapply(sigList, function(x) sum(!is.na(x$gtexDir), na.rm=TRUE)),
	dirCons =sapply(sigList, function(x) sum(x$gtexDir, na.rm=TRUE)),
	dirAndP05 =sapply(sigList, function(x) sum(x$gtexDirAnd05, na.rm=TRUE))))

sigListBonf = lapply(sigList, function(x) x[x$bonf < 0.05,])
gtexBonfRepTable = t(data.frame(numSigDisc = sapply(sigListBonf,nrow),
	numPresent =sapply(sigListBonf, function(x) sum(x$gtexPresent)),
	numTested =sapply(sigListBonf, function(x) sum(!is.na(x$gtexDir), na.rm=TRUE)),
	dirCons =sapply(sigListBonf, function(x) sum(x$gtexDir, na.rm=TRUE)),
	dirAndP05 =sapply(sigListBonf, function(x) sum(x$gtexDirAnd05, na.rm=TRUE))))

## both
bothFdrRepTable = t(data.frame(numSigDisc = sapply(sigList,nrow),
	numPresent =sapply(sigList, function(x) sum(!is.na(x$cmcTested & x$gtexPresent))),
	numTested =sapply(sigList, function(x) sum(x$cmcTested &
		!is.na(x$gtexDir), na.rm=TRUE)),
	dirCons =sapply(sigList, function(x) sum(x$cmcDir & x$gtexDir, na.rm=TRUE)),
	dirAndCmc01 =sapply(sigList, function(x) sum(x$cmcDirAndP01 & x$gtexDir, na.rm=TRUE)),
	dirAndCmc01Gtex05 =sapply(sigList, function(x) 
		sum(x$cmcDirAndP01 & x$gtexDirAnd05, na.rm=TRUE))))

bothBonfRepTable = t(data.frame(numSigDisc = sapply(sigListBonf,nrow),
	numPresent =sapply(sigListBonf, function(x) sum(!is.na(x$cmcTested & x$gtexPresent))),
	numTested =sapply(sigListBonf, function(x) sum(x$cmcTested &
		!is.na(x$gtexDir), na.rm=TRUE)),
	dirCons =sapply(sigListBonf, function(x) sum(x$cmcDir & x$gtexDir, na.rm=TRUE)),
	dirAndCmc01 =sapply(sigListBonf, function(x) sum(x$cmcDirAndP01 & x$gtexDir, na.rm=TRUE)),
	dirAndCmc01Gtex05 =sapply(sigListBonf, function(x) 
		sum(x$cmcDirAndP01 & x$gtexDirAnd05, na.rm=TRUE))))
		
repTabOut = rbind(cmcFdrRepTable,cmcBonfRepTable,
	gtexFdrRepTable,gtexBonfRepTable,
	bothFdrRepTable,bothBonfRepTable)
type = rownames(repTabOut)
rownames(repTabOut) = NULL
repTabOut = as.data.frame(repTabOut)
repTabOut$Type = type
repTabOut$Data = rep(c("CMC","GTEx","Both"), times = c(12,10,12))
repTabOut = repTabOut[,c(6:7,1:5)]
repTabOut$Data[c(1,7,13,18,23,29)] = "Discovery"
repTabOut$Type[c(1,7,13,18,23,29)] = ""
write.csv(repTabOut, file="tables/suppTable8_eqtlReplicationMetrics.csv",
	row.names=FALSE)
### other stats	
lapply(sigList, function(x) {
	table(x$cmcDirAndP01, x$gtexDir, 
		useNA="ifany", dnn = c("CMC", "GTEX"))
})

lapply(sigList, function(x) {
	table(x$cmcDirAndP01[x$bonf < 0.05], 
		x$gtexDir[x$bonf < 0.05], 
		useNA="ifany", dnn = c("CMC", "GTEX"))
})	

sigGeneRep = sigGene[which(sigGene$cmcDirAndP01 & 
	sigGene$gtexDir),]
sigExonRep = sigExon[which(sigExon$cmcDirAndP01 & 
	sigExon$gtexDir),]
sigJxnRep = sigJxn[which(sigJxn$cmcDirAndP01 & 
	sigJxn$gtexDir),]
sigTransRep = sigTrans[which(sigTrans$cmcDirAndP01 & 
	sigTrans$gtexDir),]
sigDerRep = sigDer[which(sigDer$cmcDirAndP01 & 
	sigDer$gtexDir),]

## each type	
geneEqtlRep = geneEqtl[geneEqtl$gene %in% sigGeneRep$gene,]
exonEqtlRep = exonEqtl[exonEqtl$gene %in% sigExonRep$exon,]
junctionEqtlRep = junctionEqtl[junctionEqtl$gene %in% sigJxnRep$jxn,]
transcriptEqtlRep = transcriptEqtl[transcriptEqtl$gene %in% sigTransRep$tx,]
derEqtlRep = derEqtl[derEqtl$gene %in% sigDerRep$der,]


### filter to best feature stats
sigEqtlsRep = allEqtlsRep[!duplicated(allEqtlsRep$Feature),]
save(sigEqtlsRep, compress=TRUE,
	file="/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/sigEqtlsRep.rda")


########################	
## best without rep ####

#### add other stats to the discovery for filtering
sigGene$Symbol = geneEqtl$Symbol[match(sigGene$gene, geneEqtl$gene)]
sigGene$EnsemblID = geneEqtl$EnsemblID[match(sigGene$gene, geneEqtl$gene)]
sigGene$Class = geneEqtl$Class[match(sigGene$gene, geneEqtl$gene)]

sigExon$Symbol = exonEqtl$Symbol[match(sigExon$exon, exonEqtl$gene)]
sigExon$EnsemblID = exonEqtl$EnsemblID[match(sigExon$exon, exonEqtl$gene)]
sigExon$Class = exonEqtl$Class[match(sigExon$exon, exonEqtl$gene)]

sigJxn$Symbol = junctionEqtl$Symbol[match(sigJxn$jxn, junctionEqtl$gene)]
sigJxn$EnsemblID = junctionEqtl$EnsemblID[match(sigJxn$jxn, junctionEqtl$gene)]
sigJxn$Class = junctionEqtl$Class[match(sigJxn$jxn, junctionEqtl$gene)]

sigTrans$Symbol = transcriptEqtl$Symbol[match(sigTrans$tx, transcriptEqtl$gene)]
sigTrans$EnsemblID = transcriptEqtl$EnsemblID[match(sigTrans$tx, transcriptEqtl$gene)]
sigTrans$Class = transcriptEqtl$Class[match(sigTrans$tx, transcriptEqtl$gene)]

sigDer$Symbol = derEqtl$Symbol[match(sigDer$der, derEqtl$gene)]
sigDer$EnsemblID = derEqtl$EnsemblID[match(sigDer$der, derEqtl$gene)]
sigDer$Class = derEqtl$Class[match(sigDer$der, derEqtl$gene)]
save(sigExon, sigJxn, sigDer, sigGene, sigTrans, compress=TRUE,
	file = "rdas/discovery_eqtls_forRepStats.rda")


####################################	
### filter smaller for plotting #####

load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/allEqtlsRep.rda")
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/sigEqtlsRep.rda")
