###
###
source("../eqtl_functions.R")

library(GenomicRanges)
library(clusterProfiler)

##### summary statistics for all features ####
load("../devel/rdas/devStats_controlSamples.rda")

#### and those with a switch
load("../devel/rdas/isoform_switch_devel_byFeature.rda")

####################################
##### gwas enrichment for switch ###
geneMapGR = statList$Gene
geneFn = "/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/Counts/Gene/DLPFC_PolyA_R4247_Ensembl_v75_Genes.counts"
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]
geneMapGR$codingL = geneMap$Length[match(names(geneMapGR), geneMap$Geneid)]

exonTable = table(statList$Exon$EnsemblGeneID)
geneMapGR$numExons = as.numeric(exonTable[match(names(geneMapGR), names(exonTable))])

## read in GWAS

pgcFinal = read.delim("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/tables/pgc2_loci.txt", as.is=TRUE)
pgcFinal$chr = ss(pgcFinal$Position..hg19.,":")
tmp = ss(pgcFinal$Position..hg19.,":",2)
pgcFinal$start = as.numeric(ss(tmp, "-"))
pgcFinal$end = as.numeric(ss(tmp, "-",2))
gr108 = GRanges(pgcFinal$chr, IRanges(pgcFinal$start,pgcFinal$end))

## genes in 108?
geneMapGR$pgcOverlap = countOverlaps(geneMapGR, gr108) > 0

## compared to expressed background
geneSwitchList = lapply(switchList, rownames)
switchMat = sapply(geneSwitchList, function(x) names(geneMapGR) %in% x)
colnames(switchMat) = paste0("switch_",colnames(switchMat)) 
mcols(geneMapGR) = cbind(mcols(geneMapGR), switchMat)

## just expressed
geneMapGR_exprs = geneMapGR[!is.na(geneMapGR$p_bonf) & 
	!is.na(geneMapGR$EntrezID)]

## gene length in PGC?
summary(glm(pgcOverlap ~ log2(codingL), 
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(pgcOverlap ~ log2(numExons), 
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(pgcOverlap ~ log2(meanRpkm), 
	data= mcols(geneMapGR_exprs),
	family = "binomial"))

## gene length in switch	
summary(glm(switch_Exon ~ log2(codingL), 
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(switch_Exon ~ log2(numExons), 
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(switch_Exon ~ log2(meanRpkm), 
	data= mcols(geneMapGR_exprs),
	family = "binomial"))

summary(glm(switch_Junction ~ log2(codingL), 
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(switch_Junction ~ log2(numExons), 
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(switch_Junction ~ log2(meanRpkm), 
	data= mcols(geneMapGR_exprs),
	family = "binomial"))

## switch on pgc?
summary(glm(pgcOverlap ~ switch_Exon,
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(pgcOverlap ~ switch_Junction,
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
	
## switch on pgc?
summary(glm(pgcOverlap ~ switch_Exon + log2(codingL),
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(pgcOverlap ~ switch_Junction + log2(codingL),
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(pgcOverlap ~ switch_Exon + log2(numExons),
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(pgcOverlap ~ switch_Junction + log2(numExons),
	data= mcols(geneMapGR_exprs),
	family = "binomial"))

####################################
########## other GWAS ##############
xx=load("/users/ajaffe/Lieber/Projects/RNAseq/n36/finalCode/gwasResults_lifted.rda")

geneMapGR_exprs$T2D_GWAS = countOverlaps(geneMapGR_exprs, gwasLift[gwasLift$Dx == "T2D"]) > 0
geneMapGR_exprs$PD_GWAS = countOverlaps(geneMapGR_exprs, gwasLift[gwasLift$Dx == "Park"]) > 0
geneMapGR_exprs$AD_GWAS = countOverlaps(geneMapGR_exprs, gwasLift[gwasLift$Dx == "Alz"]) > 0


## switch on others?
summary(glm(T2D_GWAS ~ switch_Exon+ log2(codingL),
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(T2D_GWAS ~ switch_Junction+ log2(codingL),
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
	
summary(glm(PD_GWAS ~ switch_Exon+ log2(codingL),
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(PD_GWAS ~ switch_Junction+ log2(codingL),
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
		
summary(glm(AD_GWAS ~ switch_Exon+ log2(codingL),
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
summary(glm(AD_GWAS ~ switch_Junction+ log2(codingL),
	data= mcols(geneMapGR_exprs),
	family = "binomial"))
	