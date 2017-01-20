###
library(GenomicRanges)

source("eqtl_functions.R")

load("phenotype_annotated_szControlEqtl_DLPFC.rda")

junctionFiles = paste0("/dcs01/lieber/ajaffe/Brain/DLPFC_PolyA/Junctions/DLPFC_PolyA_", 
	pd$RNum, "_junctions_primaryOnly.count")
all(file.exists(junctionFiles)) #  TRUE
 names(junctionFiles) = pd$RNum

juncCounts = junctionCount(junctionFiles[1:50],
 	strandSpecific=FALSE, maxCores=12)
juncCounts8 = junctionCount(junctionFiles[1:50],
	minOverhang=8,	maxCores=12)

anno = juncCounts$anno
anno8 = juncCounts8$anno

mm = rowMeans(as.data.frame(juncCounts$countDF))
mm8 = rowMeans(as.data.frame(juncCounts8$countDF))

oo = findOverlaps(anno, anno8, type="equal")

table(mm[queryHits(oo)] - mm8[subjectHits(oo)] == 0)

## annotate junctions
load("/users/ajaffe/Lieber/Projects/RNAseq/ensembl_hg19_v75_junction_annotation.rda")

seqlevels(anno) = seqlevels(anno8) = paste0("chr", c(1:22,"X","Y","M"))

## add additional annotation
anno$inEnsembl = countOverlaps(anno, theJunctions, type="equal") > 0
anno$inEnsemblStart = countOverlaps(anno, theJunctions, type="start") > 0
anno$inEnsemblEnd = countOverlaps(anno, theJunctions, type="end") > 0
anno8$inEnsembl = countOverlaps(anno8, theJunctions, type="equal") > 0
anno8$inEnsemblStart = countOverlaps(anno8, theJunctions, type="start") > 0
anno8$inEnsemblEnd = countOverlaps(anno8, theJunctions, type="end") > 0

anno$code = ifelse(anno$inEnsembl, "InEns", 
	ifelse(anno$inEnsemblStart & anno$inEnsemblEnd, "ExonSkip",
	ifelse(anno$inEnsemblStart | anno$inEnsemblEnd, "AltStartEnd", "Novel")))
anno8$code = ifelse(anno8$inEnsembl, "InEns", 
	ifelse(anno8$inEnsemblStart & anno8$inEnsemblEnd, "ExonSkip",
	ifelse(anno8$inEnsemblStart | anno8$inEnsemblEnd, "AltStartEnd", "Novel")))
