##
####
# R-devel
source("../eqtl_functions.R")

library(GenomicRanges)
library(limma)

# transcript
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/transcript/transcript_data_filtered_n495.rda")
# gene exon junction
load("/dcs01/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
# ERs
load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda")

pd$ageGroup = cut(pd$Age, c(-1,0,1,10,20,50,100))

### add entrez ID
geneMapGR = makeGRangesFromDataFrame(geneMap,keep=TRUE)
jMap$EntrezID = geneMapGR$EntrezID[match(jMap$newGeneID,names(geneMapGR))]
tMap$EntrezID = geneMapGR$EntrezID[match(tMap$EnsemblGeneID,	names(geneMapGR))]

# add annotation to regions
dA = distanceToNearest(regions, geneMapGR)
regions$EnsemblGeneID = names(geneMapGR)[subjectHits(dA)]
regions$EntrezID = geneMapGR$EntrezID[subjectHits(dA)]
regions$EnsemblGeneID[mcols(dA)$distance > 100] = NA
regions$EntrezID[mcols(dA)$distance > 100] = NA

### filter for control
controlIndex = which(pd$Dx == "Control")
pd = pd[controlIndex,]
geneRpkm = geneRpkm[,controlIndex]
exonRpkm = exonRpkm[,controlIndex]
jRpkm = jRpkm[,controlIndex]
tFpkm = tFpkm[,controlIndex]
regionMat = regionMat[,controlIndex]

###### filter gene, exon, junction
exprsGeneIndex = which(rowMeans(geneRpkm) > 0.01)
exprsExonIndex = which(rowMeans(exonRpkm) > 0.1)
exonMap$coord = paste0("chr",exonMap$Chr, ":", exonMap$Start, "-",
	exonMap$End,"(",exonMap$Strand, ")")
exprsJxnIndex = which(rowMeans(jRpkm > 0.2) & jMap$code !="Novel")

##### xform
xformData = list(Gene = log2(geneRpkm+1), Exon = log2(exonRpkm+1), 
	Junction = log2(jRpkm+1), Transcript =  log2(tFpkm+1), 
	DER = log2(regionMat+1))
dataOut = do.call("rbind", xformData)
logDataOut = log2(dataOut+1)
rownames(logDataOut) = ss(rownames(logDataOut), "\\.",2)

### annotation
geneMap$EnsemblGeneID = rownames(geneMap)
geneMap$Class = "InEns"
geneMap$Type = "Gene"

colnames(exonMap)[1] = "EnsemblGeneID"
exonMap$Class = "InEns"
exonMap$Type = "Exon"

jMap$EnsemblGeneID = jMap$newGeneID
jMap$Symbol = jMap$newGeneSymbol
jMap$Type = "Junction"
colnames(mcols(jMap))[9] ="Class"
jMapDf = as.data.frame(jMap)
jMapDf = jMapDf[,c("seqnames","start","end","EnsemblGeneID","Symbol","Class","Type")]
colnames(jMapDf)[1:3] = c("Chr","Start","End")

colnames(mcols(tMap))[3] = "Symbol"
tMap$Type = "Transcript"
tMap$Class = "Novel"
tMap$Class[tMap$class_code == "="] = "InEns"
tMap$Class[tMap$class_code == "j"] = "ExonSkip"
tMap$Class[tMap$class_code %in% c("c", "e", "o", "p", "x") ] = "AltStartEnd"
tMapDf = as.data.frame(tMap)
tMapDf = tMapDf[,c("seqnames","start","end","EnsemblGeneID","Symbol","Class","Type")]
colnames(tMapDf)[1:3] = c("Chr","Start","End")

colnames(mcols(regions))[7] = "Symbol"
regions$Type = "ER"
regions$Class = "Novel"
regions$Class[regions$annoClass == "strictExonic"] = "InEns"
regions$Class[regions$annoClass %in% c("exonIntron", "extendUTR")] = "AltStartEnd"
regionsDf = as.data.frame(regions)
regionsDf = regionsDf[,c("seqnames","start","end","EnsemblGeneID","Symbol","Class","Type")]
colnames(regionsDf)[1:3] = c("Chr","Start","End")

mapOut = rbind(geneMap[,colnames(jMapDf)], exonMap[,colnames(jMapDf)],
	jMapDf, tMapDf, regionsDf)
	
##### write out #####
outPath = "/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/db/szPlusControl/fullData/"

write.csv(logDataOut, file = gzfile(paste0(outPath, "log2_exprs_data_ALL.csv.gz")),
	quote=FALSE)
write.csv(mapOut, file = gzfile(paste0(outPath, "annotation_data_ALL.csv.gz")),
	quote=FALSE)

###################	
#### dev stats ####
###################

load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/devStats_controlSamples.rda")
fullStats = unlist(statList)
fullStatsOut = mcols(fullStats)[,-(1:5)]
rownames(fullStatsOut) = ss(names(fullStats), "\\.",2)
fullStatsOut = as.data.frame(fullStatsOut)
write.csv(fullStatsOut, quote=FALSE,
	file = gzfile(paste0(outPath, "develStats_controls.csv.gz")))
	
###################	
#### Dx stats ####
###################

xx = load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/caseControl/rdas/all_de_features.rda")

dxStats = unlist(outStats)
dxStatsOut = mcols(dxStats)[,-c(1:4,6)]
rownames(dxStatsOut) = ss(names(dxStatsOut), "\\.",2)
dxStatsOut = as.data.frame(dxStatsOut)
write.csv(dxStatsOut, quote=FALSE,
	file = gzfile(paste0(outPath, "dxStats_szControls.csv.gz")))
	
## phenotype
phenoOut = pd[,c(1:8,35:39,45:57)]
write.csv(phenoOut, quote=FALSE,row.names=FALSE,
	file = gzfile(paste0(outPath, "phenotypeData_szControls.csv.gz")))
	