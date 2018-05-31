##############
## preprocess GTEX Brain Counts
##############

source("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl_functions.R")

### load packages
library(GenomicRanges)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(biomaRt)
library(parallel)
library(stringr)

# load data and filter
load("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/gtexPd.Rdata")
gtexPd$sra_accession = str_trim(gtexPd$sra_accession)
gtexPd = gtexPd[!is.na(gtexPd$sra_accession),]

## drop those w/o SRR
gtexPd = gtexPd[which(gtexPd$SAMPLE_USE == 
	"Seq_RNA_WTSS; Seq_RNA_Expression"),]

## just brain
gtexPd = gtexPd[which(gtexPd$SMTS == "Brain"),]

## add bam and junctions
gtexPd$bamFile = paste0("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/TopHat/",
	gtexPd$sra_accession, "/accepted_hits.bam")
table(file.exists(gtexPd$bamFile))	

## add total mapped
libSize = getTotalMapped(gtexPd$bamFile,mc.cores=16)
gtexPd$totalMapped = libSize$totalMapped
gtexPd$mitoMapped = libSize$mitoMapped

#########################
#### GENE COUNTS ########

### gene count files
geneFn = paste0("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Counts/Gene/",
	gtexPd$sra_accession, "_Ensembl_v75_Genes.counts")
names(geneFn) = gtexPd$sra_accession
all(file.exists(geneFn))

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
geneCounts = geneCounts[,gtexPd$sra_accession] # put in order

# make RPKM
bg = matrix(rep(gtexPd$totalMapped), nc = nrow(gtexPd), 
	nr = nrow(geneCounts),	byrow=TRUE)
widG = matrix(rep(geneMap$Length), nr = nrow(geneCounts), 
	nc = nrow(gtexPd),	byrow=FALSE)
geneRpkm = geneCounts/(widG/1000)/(bg/1e6)

# number of reads assigned
geneStatList = mclapply(paste0(geneFn, ".summary"), 
	read.delim,row.names=1,mc.cores=12)
geneStats = do.call("cbind", geneStatList)
colnames(geneStats) = gtexPd$sra_accession
gtexPd$totalAssignedGene = as.numeric(geneStats[1,] / colSums(geneStats))

#### take eQTLs
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/annotated_gene_eqtl_szControl_cisOnly.rda")
geneRpkmSub = geneRpkm[sigGene$gene,]
geneMapSub = geneMap[sigGene$gene,]

#################################
######### EXON COUNTS ###########

### exon count files
exonFn = paste0("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Counts/Exon/",
	gtexPd$sra_accession, "_Ensembl_v75_Exons.counts")
names(exonFn) = gtexPd$sra_accession
all(file.exists(exonFn))

### read in annotation ##
exonMap = read.delim(exonFn[1], skip=1, as.is=TRUE)[,1:6]
exonMap$Symbol = sym$hgnc_symbol[match(exonMap$Geneid, sym$ensembl_gene_id)]
exonMap$EntrezID = sym$entrezgene[match(exonMap$Geneid, sym$ensembl_gene_id)]
rownames(exonMap) = paste0("e", rownames(exonMap))
mmExon = match(sigExon$exon, rownames(exonMap))
exonMapSub = exonMap[mmExon,]

## counts
exonCountList = lapply(exonFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[mmExon,1]
})
exonCountsSub = do.call("cbind", exonCountList)
rownames(exonCountsSub) = rownames(exonMapSub)
exonCountsSub = exonCountsSub[,gtexPd$sra_accession] # put in order

## make RPKM
bgE = matrix(rep(gtexPd$totalMapped), nc = nrow(gtexPd), 
	nr = nrow(exonCountsSub),	byrow=TRUE)
widE = matrix(rep(exonMap$Length), nr = nrow(exonCountsSub), 
	nc = nrow(gtexPd),	byrow=FALSE)
exonRpkmSub = exonCountsSub/(widE/1000)/(bgE/1e6)

#############################
##### JUNCTIONS #############

junctionFiles = gsub("accepted_hits.bam", "junctions.bed", 
	gtexPd$bamFile)
names(junctionFiles) = gtexPd$sra_accession
# all(file.exists(junctionFiles)) #  TRUE

##### convert using modified bed_to_juncs
junctionCounts = gsub(".bed", ".juncs", junctionFiles)
names(junctionCounts) = names(junctionFiles)
# thecall = paste("/users/ajaffe/Lieber/Projects/RNAseq/bed_to_juncs_withCount <",
	# junctionFiles,">", junctionCounts)[!file.exists(junctionCounts)]
# # parallel::mclapply(thecall, system, mc.cores=12) 

### get junction counts
juncCounts = junctionCount(junctionCounts,maxCores=12)

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
anno$numTx = elementNROWS(anno$ensemblTx)

## junction code
anno$code = ifelse(anno$inEnsembl, "InEns", 
	ifelse(anno$inEnsemblStart & anno$inEnsemblEnd, "NovelTrans",
	ifelse(anno$inEnsemblStart | anno$inEnsemblEnd, "NovelJxn", "NewJxn")))

## only take those junctions that are eqtls
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/eqtl/joint/rdas/annotated_junction_eqtl_szControl_cisOnly.rda")
jxnIndex= which(names(anno) %in% sigJxn$jxn)
jMapSub = anno[jxnIndex,]
jCountsSub = juncCounts$countDF[jxnIndex,gtexPd$sra_accession]

mappedPer80M = gtexPd$totalMapped/80e6
jRpkmSub = as.data.frame(DataFrame(mapply(function(x,d) 	
	x/d, jCountsSub , mappedPer80M)))

### save counts
save(gtexPd, jMapSub, jRpkmSub, geneRpkm, 
	geneMap, exonRpkmSub, exonMapSub, compress=TRUE,
	file="/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/Counts/rpkmCounts_brainGTEX_subsetEqtl.rda")

