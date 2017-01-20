##
library(derfinder)
library(GenomicRanges)

load("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda")
load("/dcs01/ajaffe/Brain/derRuns/szControlEQTL/szControl/regionMatrix/regionMat-cut5.Rdata")

# extract
regions = unlist(GRangesList(lapply(regionMat, '[[', 'regions')))
names(regions) = NULL
regionMat = do.call("rbind", lapply(regionMat, '[[', 'coverageMatrix'))
regionMat = regionMat[,pd$RNum] # put in order

# add rownames
names(regions) = rownames(regionMat) = paste0("er", 1:nrow(regionMat))

# filter
ind = which(width(regions) >= 12)
regions = regions[ind,]
regionMat = regionMat[ind,]

# annotate
load("/users/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
ensemblAnno = annotateRegions(regions,gs)
countTable = ensemblAnno$countTable

## gene annotation?
geneMapGR = makeGRangesFromDataFrame(geneMap, keep=TRUE)

## overlaps
dA = distanceToNearest(regions, geneMapGR)
regions$nearestSymbol = geneMapGR$Symbol[subjectHits(dA)]
regions$distToGene = mcols(dA)$distance
mcols(regions) = cbind(mcols(regions), countTable)

## add additional annotation
regions$annoClass = NA
regions$annoClass[regions$exon > 0 & 
	regions$intron == 0 &
	regions$intergenic == 0] = "strictExonic"
regions$annoClass[regions$exon == 0 & 
	regions$intron > 0 &
	regions$intergenic == 0] = "strictIntronic"
regions$annoClass[regions$exon == 0 & 
	regions$intron == 0 &
	regions$intergenic > 0] = "strictIntergenic"
regions$annoClass[regions$exon > 0 & 
	regions$intron > 0 &
	regions$intergenic == 0] = "exonIntron"
regions$annoClass[regions$exon > 0 & 
	regions$intergenic > 0] = "extendUTR"

save(pd, regions, regionMat, compress=TRUE,
	file = "/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda")
	
rtracklayer::export(regions, "szControl_ERs.bed")