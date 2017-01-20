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

### add number minor alleles
snpMap$numMinorHom = rowSums(snp == 2,na.rm=TRUE)

# load expression annotation
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
exonMap$coord = with(exonMap, paste0(Chr, ":", Start, "-", End, "(", Strand, ")"))
snpMap$numMinorHom_13plus = rowSums(snp[,pd$Age > 13] == 2,na.rm=TRUE)

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
load("rdas/annotated_gene_eqtl_control_13plus_cisOnly.rda")
load("rdas/gene_eqtl_control_13plus_cisOnly.rda")
geneEqtl = meGeneCis$cis$eqtls
geneEqtl = geneEqtl[geneEqtl$FDR < 0.01,]
geneEqtl$gene = as.character(geneEqtl$gene)
geneEqtl$snps = as.character(geneEqtl$snps)
geneEqtl$bonf = p.adjust(geneEqtl$pvalue, "bonf", 
	n = (sigGene$bonf[1]/sigGene$pvalue[1]))

###### EXON
load("rdas/annotated_exon_eqtl_control_13plus_cisOnly.rda")
load("rdas/exon_eqtl_control_13plus_cisOnly.rda")
exonEqtl = meExonCis$cis$eqtls
exonEqtl = exonEqtl[exonEqtl$FDR < 0.01,]
exonEqtl$gene = as.character(exonEqtl$gene)
exonEqtl$snps = as.character(exonEqtl$snps)
exonEqtl$bonf = p.adjust(exonEqtl$pvalue, "bonf", 
	n = (sigExon$bonf[1]/sigExon$pvalue[1]))

######### TX
load("rdas/annotated_transcript_eqtl_control_13plus_cisOnly.rda")
load("rdas/transcript_eqtl_control_13plus_cisOnly.rda")
transcriptEqtl = meTransCis$cis$eqtls
transcriptEqtl = transcriptEqtl[transcriptEqtl$FDR < 0.01,]
transcriptEqtl$gene = as.character(transcriptEqtl$gene)
transcriptEqtl$snps = as.character(transcriptEqtl$snps)
transcriptEqtl$bonf = p.adjust(transcriptEqtl$pvalue, "bonf", 
	n = (sigTrans$bonf[1]/sigTrans$pvalue[1]))
	
###### JUNCTION
load("rdas/annotated_junction_eqtl_control_13plus_cisOnly.rda")
load("rdas/junction_eqtl_control_13plus_cisOnly.rda")
junctionEqtl = meJxnCis$cis$eqtls
junctionEqtl = junctionEqtl[junctionEqtl$FDR < 0.01,]
junctionEqtl$gene = as.character(junctionEqtl$gene)
junctionEqtl$snps = as.character(junctionEqtl$snps)
junctionEqtl$bonf = p.adjust(junctionEqtl$pvalue, "bonf", 
	n = (sigJxn$bonf[1]/sigJxn$pvalue[1]))

###### DER
load("rdas/annotated_der_eqtl_control_13plus_cisOnly.rda")
load("rdas/der_eqtl_control_13plus_cisOnly.rda")
derEqtl = meDerCis$cis$eqtls
derEqtl = derEqtl[derEqtl$FDR < 0.01,]
derEqtl$gene = as.character(derEqtl$gene)
derEqtl$snps = as.character(derEqtl$snps)
derEqtl$bonf = p.adjust(derEqtl$pvalue, "bonf", 
	n = (sigDer$bonf[1]/sigDer$pvalue[1]))

### add number of homozygous minor
geneEqtl$numMinor13 = snpMap$numMinorHom_13plus[
	match(geneEqtl$snps, snpMap$SNP)]
exonEqtl$numMinor13 = snpMap$numMinorHom_13plus[
	match(exonEqtl$snps, snpMap$SNP)]
transcriptEqtl$numMinor13 = snpMap$numMinorHom_13plus[
	match(transcriptEqtl$snps, snpMap$SNP)]
junctionEqtl$numMinor13 = snpMap$numMinorHom_13plus[
	match(junctionEqtl$snps, snpMap$SNP)]
derEqtl$numMinor13 = snpMap$numMinorHom_13plus[
	match(derEqtl$snps, snpMap$SNP)]

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
controlEqtls = rbind(geneEqtl, exonEqtl, 
	transcriptEqtl, junctionEqtl, derEqtl)
controlEqtls$snpRsNum = snpMap$name[match(controlEqtls$snps, snpMap$SNP)]
controlEqtls$Type = rep(c("Gene","Exon", "Transcript", "Junction", "ER"),
	c(nrow(geneEqtl), nrow(exonEqtl), nrow(transcriptEqtl), 
		nrow(junctionEqtl), nrow(derEqtl)))
controlEqtls$statistic = NULL

## make smaller 
controlEqtls = DataFrame(controlEqtls)
controlEqtls$Class = Rle(controlEqtls$Class)
controlEqtls$Type = Rle(controlEqtls$Type)
controlEqtls$Type = Rle(controlEqtls$Type)
controlEqtls = controlEqtls[order(controlEqtls$pvalue),]
names(controlEqtls)[1:2] = c("SNP", "Feature")

#############################
### get transcript info #####
xx=load("/users/ajaffe/Lieber/Projects/RNAseq/ensembl_V75_feature_to_Tx.rda")

## add annotation from tx
theTxs = CharacterList(as.list(tMap$nearest_ref))
names(theTxs) = names(tMap)
allTx = c(allTx, theTxs)

mmTx = match(controlEqtls$Feature, names(allTx))
tx = CharacterList(vector("list", nrow(controlEqtls)))
tx[!is.na(mmTx)] = allTx[mmTx[!is.na(mmTx)]]
controlEqtls$NumTxEqtl = elementLengths(tx)
controlEqtls$WhichTx = tx

mmTxGene = match(controlEqtls$EnsemblID, names(allTx))
txGene = CharacterList(vector("list", nrow(controlEqtls)))
txGene[!is.na(mmTxGene)] = allTx[mmTxGene[!is.na(mmTxGene)]]
controlEqtls$NumTxGene = elementLengths(txGene)

controlEqtls$NumTxEqtl = Rle(controlEqtls$NumTxEqtl)
controlEqtls$NumTxGene = Rle(controlEqtls$NumTxGene)

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
controlEqtls$Coordinates = coords[match(controlEqtls$Feature, names(coords))]

### add distance
starts = c(geneMap$Start, exonMap$Start, start(tMap), 
	start(jMap), start(regions))
ends = c(geneMap$End, exonMap$End, end(tMap), 
	end(jMap), end(regions))
names(starts) = names(ends) = c(rownames(geneMap), rownames(exonMap),
		names(tMap), names(jMap), names(regions))
snpPos = snpMap$POS[match(controlEqtls$SNP, snpMap$SNP)]	
tmp = cbind(snpPos - starts[match(controlEqtls$Feature, names(coords))],
		snpPos - ends[match(controlEqtls$Feature, names(coords))])
mins = matrixStats::rowMins(abs(tmp))
controlEqtls$snpDistToFeature = sign(rowMeans(tmp))*mins

save(controlEqtls, compress=TRUE,
	file="/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/controlEqtls.rda")
	
### filter to best feature stats
sigEqtls = controlEqtls[!duplicated(controlEqtls$Feature),]
save(sigEqtls, compress=TRUE,
	file="/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/sigControlEqtls.rda")


########################	
## best without rep ####

# filter for Age
aIndex= which(pd$Age > 13 & pd$Dx == "Control")
pdSub = pd[aIndex,]

snpSub = as.matrix(snp[rownames(snp) %in% controlEqtls$SNP,aIndex])
snpMapSub = snpMap[rownames(snp) %in% controlEqtls$SNP,]

geneRpkmSub = as.matrix(log2(geneRpkm[
	rownames(geneRpkm) %in% sigEqtls$Feature,aIndex]+1))
geneMapSub = geneMap[rownames(geneMap) %in% sigEqtls$Feature,]

exonRpkmSub = as.matrix(log2(exonRpkm[
	rownames(exonRpkm) %in% sigEqtls$Feature,aIndex]+1))
exonMapSub = exonMap[rownames(exonMap) %in% sigEqtls$Feature,]

jRpkmSub = as.matrix(log2(jRpkm[
	rownames(jRpkm) %in% sigEqtls$Feature,aIndex]+1))
jMapSub = jMap[names(jMap) %in% sigEqtls$Feature]

tFpkmSub = as.matrix(log2(tFpkm[
	rownames(tFpkm) %in% sigEqtls$Feature,aIndex]+1))
tMapSub = tMap[names(tMap) %in% sigEqtls$Feature]

regionMatSub = as.matrix(log2(regionMat[
	rownames(regionMat) %in% sigEqtls$Feature,aIndex]+1))
regionsSub = regions[names(regions) %in% sigEqtls$Feature]

## residualize
mod = model.matrix(~snpPC1 + snpPC2 + snpPC3 + Sex,data=pdSub)

cleanGeneSub = cleaningY(geneRpkmSub, cbind(mod, pcsGene),P=1)
cleanExonSub = cleaningY(exonRpkmSub, cbind(mod, pcsExon),P=1)
cleanJxnSub = cleaningY(jRpkmSub, cbind(mod, pcsJxn),P=1)
cleanTxSub = cleaningY(tFpkmSub, cbind(mod, pcsTrans),P=1)
cleanErSub = cleaningY(regionMatSub, cbind(mod, pcsDer),P=1)

save(pdSub, snpSub, snpMapSub, 
	geneRpkmSub, cleanGeneSub, geneMapSub,
	exonRpkmSub, cleanExonSub, exonMapSub,
	jRpkmSub, cleanJxnSub, jMapSub,
	tFpkmSub, cleanTxSub, tMapSub,
	regionMatSub, cleanErSub, regionsSub, compress=TRUE,
	file = "/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/cleaned_eqtl_data_subset_control_n237.rda")
	
		
#################################
### export database files #######
eqtlStatsOut = controlEqtls[,c(1:6,17)]
eqtlStatsOut = as.data.frame(eqtlStatsOut)
write.table(eqtlStatsOut, row.names=FALSE, sep = "\t",quote=FALSE,
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/controlOnly/eqtl_stats_controlOnly_LIBD.txt.gz"))

## expression anno
exprsAnnoOut = controlEqtls[,c(2,8:16)]
exprsAnnoOut = exprsAnnoOut[!duplicated(exprsAnnoOut$Feature),]
exprsAnnoOut$chr = ss(exprsAnnoOut$Coordinates, ":")
exprsAnnoOut$start = as.integer(ss(ss(exprsAnnoOut$Coordinates, ":",2),"-"))
exprsAnnoOut$end = as.integer(ss(ss(ss(exprsAnnoOut$Coordinates, ":",2),"-",2),"\\("))
exprsAnnoOut = as.data.frame(exprsAnnoOut)
exprsAnnoOut$WhichTx = sapply(exprsAnnoOut$WhichTx, paste, collapse=";")
write.table(exprsAnnoOut, row.names=FALSE, sep = "\t",quote=FALSE,
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/controlOnly/eqtl_exprsAnno_controlOnly_LIBD.txt.gz"))

## SNP anno	
snpAnnoOut = controlEqtls[,c(1,11)]
snpAnnoOut = snpAnnoOut[!duplicated(snpAnnoOut$SNP),]
mm = match(snpAnnoOut$SNP, snpMapSub$SNP)
snpMapTmp = snpMapSub[mm,]
colnames(snpMapTmp)[c(1,3,4:5)] = c("chr","pos","Counted","Ref")
snpMapTmp$chr = paste0("chr", snpMapTmp$chr)

snpTmp = snpSub[mm,]
snpMapTmp$n0 = rowSums(snpTmp == 0,na.rm=TRUE)
snpMapTmp$n1 = rowSums(snpTmp == 1,na.rm=TRUE)
snpMapTmp$n2 = rowSums(snpTmp == 2,na.rm=TRUE)

snpAnnoOut = cbind(snpAnnoOut, 
	snpMapTmp[,c("chr","pos","Counted","Ref","numImp","n0","n1","n2")])
snpAnnoOut = as.data.frame(snpAnnoOut)
write.table(snpAnnoOut, row.names=FALSE, sep = "\t",quote=FALSE,
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/eqtl_snpAnno_controlOnly_LIBD.txt.gz"))

## expression data
cleanExprs = rbind(cleanGeneSub, cleanExonSub, cleanJxnSub,
	cleanTxSub, cleanErSub)
write.table(cleanExprs, sep = "\t", quote=FALSE,
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/eqtl_log2cleanExprsData_controlOnly_LIBD.txt.gz"))

## snp data
write.table(snpSub, sep = "\t",quote=FALSE,
	file=gzfile("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/eqtl_snpData_controlOnly_LIBD.txt.gz"))
