#####
source("eqtl_functions.R")

load("phenotype_annotated_szControlEqtl_DLPFC.rda")

library(GenomicRanges)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(biomaRt)

## gene counts from consortium set

### gene count files
geneFn = paste0("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/Counts/Gene/DLPFC_PolyA_",
	pd$RNum, "_Ensembl_v75_Genes.counts")
names(geneFn) = pd$RNum

### read in annotation ##
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]
rownames(geneMap) = geneMap$Geneid
geneMap$Chr = paste0("chr", ss(geneMap$Chr, ";"))
geneMap$Start = as.numeric(ss(geneMap$Start, ";"))
tmp = strsplit(geneMap$End, ";")
geneMap$End = as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand = ss(geneMap$Strand, ";")
geneMap$Geneid = NULL

### biomart 
ensembl = useMart("ENSEMBL_MART_ENSEMBL", # VERSION 75, hg19
	dataset="hsapiens_gene_ensembl",
	host="feb2014.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
	values=rownames(geneMap), mart=ensembl)
geneMap$Symbol = sym$hgnc_symbol[match(rownames(geneMap), sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene[match(rownames(geneMap), sym$ensembl_gene_id)]
	
## counts
geneCountList = mclapply(geneFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=12)
geneCounts = do.call("cbind", geneCountList)
rownames(geneCounts) = rownames(geneMap)
geneCounts = geneCounts[,pd$RNum] # put in order

# make RPKM
bg = matrix(rep(pd$totalMapped), nc = nrow(pd), 
	nr = nrow(geneCounts),	byrow=TRUE)
widG = matrix(rep(geneMap$Length), nr = nrow(geneCounts), 
	nc = nrow(pd),	byrow=FALSE)
geneRpkm = geneCounts/(widG/1000)/(bg/1e6)

# number of reads assigned
geneStatList = lapply(paste0(geneFn, ".summary"), 
	read.delim,row.names=1)
geneStats = do.call("cbind", geneStatList)
colnames(geneStats) = pd$RNum
pd$totalAssignedGene = as.numeric(geneStats[1,] / colSums(geneStats))

###############
### exon counts

exonFn = paste0("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/Counts/Exon/DLPFC_PolyA_",
	pd$RNum, "_Ensembl_v75_Exons.counts")
names(exonFn) = pd$RNum
all(file.exists(exonFn))

### read in annotation ##
exonMap = read.delim(exonFn[1], skip=1, as.is=TRUE)[,1:6]
exonMap$Symbol = sym$hgnc_symbol[match(exonMap$Geneid, sym$ensembl_gene_id)]
exonMap$EntrezID = sym$entrezgene[match(exonMap$Geneid, sym$ensembl_gene_id)]
rownames(exonMap) = paste0("e", rownames(exonMap))
exonMap$Chr = paste0("chr", exonMap$Chr)

## counts
exonCountList = mclapply(exonFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=12)
exonCounts = do.call("cbind", exonCountList)
rownames(exonCounts) = rownames(exonMap)
exonCounts = exonCounts[,pd$RNum] # put in order

## remove duplicated
eMap = GRanges(exonMap$Chr, IRanges(exonMap$Start, exonMap$End))
keepIndex= which(!duplicated(eMap))
exonCounts = exonCounts[keepIndex,]
exonMap = exonMap[keepIndex,]

## make RPKM
bgE = matrix(rep(pd$totalMapped), nc = nrow(pd), 
	nr = nrow(exonCounts),	byrow=TRUE)
widE = matrix(rep(exonMap$Length), nr = nrow(exonCounts), 
	nc = nrow(pd),	byrow=FALSE)
exonRpkm = exonCounts/(widE/1000)/(bgE/1e6)

#############
##### junctions

## via primary alignments only
junctionFiles = paste0("/dcs01/lieber/ajaffe/Brain/DLPFC_PolyA/Junctions/DLPFC_PolyA_", 
	pd$RNum, "_junctions_primaryOnly_regtool.count")
all(file.exists(junctionFiles)) #  TRUE
 
juncCounts = junctionCount(junctionFiles, pd$RNum,
 	output = "Count", maxCores=12)

## annotate junctions
load("/users/ajaffe/Lieber/Projects/RNAseq/ensembl_hg19_v75_junction_annotation.rda")

anno = juncCounts$anno
seqlevels(anno) = paste0("chr", c(1:22,"X","Y","M"))

## add additional annotation
anno$inEnsembl = countOverlaps(anno, theJunctions, type="equal") > 0
anno$inEnsemblStart = countOverlaps(anno, theJunctions, type="start") > 0
anno$inEnsemblEnd = countOverlaps(anno, theJunctions, type="end") > 0

oo = findOverlaps(anno, theJunctions, type="equal")
anno$ensemblGeneID = NA
anno$ensemblGeneID[queryHits(oo)] = as.character(theJunctions$ensemblID[subjectHits(oo)])
anno$ensemblSymbol = NA
anno$ensemblSymbol[queryHits(oo)] = theJunctions$symbol[subjectHits(oo)]
anno$ensemblStrand = NA
anno$ensemblStrand[queryHits(oo)] = as.character(strand(theJunctions)[subjectHits(oo)])
anno$ensemblTx = CharacterList(vector("list", length(anno)))
anno$ensemblTx[queryHits(oo)] = theJunctions$tx[subjectHits(oo)]
anno$numTx = elementLengths(anno$ensemblTx)

# clean up
anno$ensemblSymbol = geneMap$Symbol[match(anno$ensemblGeneID, rownames(geneMap))]

## junction code
anno$code = ifelse(anno$inEnsembl, "InEns", 
	ifelse(anno$inEnsemblStart & anno$inEnsemblEnd, "ExonSkip",
	ifelse(anno$inEnsemblStart | anno$inEnsemblEnd, "AltStartEnd", "Novel")))

## b/w exons and junctions
exonGR = GRanges( exonMap$Chr,	IRanges(exonMap$Start, exonMap$End))
anno$startExon = match(paste0(seqnames(anno),":",start(anno)-1), 
	paste0(seqnames(exonGR), ":", end(exonGR)))
anno$endExon = match(paste0(seqnames(anno),":",end(anno)+1),
	paste0(seqnames(exonGR), ":", start(exonGR)))
g = data.frame(leftGene = exonMap$Geneid[anno$startExon],
	rightGene = exonMap$Geneid[anno$endExon],
	leftGeneSym = exonMap$Symbol[anno$startExon],
	rightGeneSym = exonMap$Symbol[anno$endExon],
	stringsAsFactors=FALSE)
g$newGene = NA
g$newGeneSym = NA
g$newGene[which(g$leftGene==g$rightGene)] = 
	g$leftGene[which(g$leftGene==g$rightGene)] 
g$newGeneSym[which(g$leftGene==g$rightGene)] = 
	g$leftGeneSym[which(g$leftGene==g$rightGene)] 
g$newGene[which(g$leftGene!=g$rightGene)] = 
	paste0(g$leftGene,"-",g$rightGene)[which(g$leftGene!=g$rightGene)] 
g$newGeneSym[which(g$leftGene!=g$rightGene)] = 
	paste0(g$leftGeneSym,"-",g$rightGeneSym)[which(g$leftGene!=g$rightGene)] 
g$newGene[which(is.na(g$newGene) & is.na(g$leftGene))] = 
	g$rightGene[which(is.na(g$newGene) & is.na(g$leftGene))] 
g$newGene[which(is.na(g$newGene) & is.na(g$rightGene))] = 
	g$leftGene[which(is.na(g$newGene) & is.na(g$rightGene))] 
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] = 
	g$rightGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] 
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] = 
	g$leftGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] 
g$newGeneSym[g$newGeneSym==""] = NA
g$newGeneSym[g$newGeneSym=="-"] = NA
anno$newGeneID = g$newGene
anno$newGeneSymbol = g$newGeneSym
anno$isFusion = grepl("-", anno$newGeneID)

## extract out
jMap = anno
jCounts = juncCounts$countDF
jCounts = jCounts[names(jMap),pd$RNum]

mappedPer80M = pd$totalMapped/80e6
countsM = DataFrame(mapply(function(x,d) x/d, jCounts , mappedPer80M))
rownames(jCounts) = rownames(countsM) = names(jMap)
jRpkm = as.data.frame(countsM)

## sequence of acceptor/donor sites
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
left = right = jMap
end(left) = start(left) +1
start(right) = end(right) -1

jMap$leftSeq  = getSeq(Hsapiens, left)
jMap$rightSeq = getSeq(Hsapiens, right)

### save counts
save(pd, jMap, jCounts, geneCounts, geneMap, exonCounts, exonMap, compress=TRUE,
	file="/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rawCounts_szControlDLPFC.rda")
save(pd, jMap, jRpkm, geneRpkm,	geneMap, exonRpkm, exonMap, compress=TRUE,
	file="/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")

## write out for coverage
write.table(pd[,c("RNum", "bamFile")], 
	"samples_with_bams.txt",
	row.names=FALSE, quote=FALSE, sep="\t")