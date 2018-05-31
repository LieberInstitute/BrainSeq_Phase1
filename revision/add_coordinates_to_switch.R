##
library(GenomicRanges)

## load stats as map
load("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/rdas/devStats_controlSamples.rda")

## add coordinates
statList = endoapply(statList, function(x) {
	x$coords = paste0(seqnames(x), ":", start(x), "-", end(x))
	return(x)
})

## read in switch table
table_s5 = read.csv("/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/devel/tables/suppTable5_isoformSwitches.csv",
	as.is=TRUE)

## retain subset of features
uFeatures = c(table_s5$minFeature, table_s5$maxFeature)
statListSub = endoapply(statList, function(x) x[names(x) %in% uFeatures])
coordSub = do.call("rbind", lapply(statListSub, function(x) {
	data.frame(ID = names(x), Coord = x$coords,
		stringsAsFactors=FALSE)
}))

table_s5$minFeature_hg19 = coordSub$Coord[
	match(table_s5$minFeature, coordSub$ID)]
table_s5$maxFeature_hg19 = coordSub$Coord[
	match(table_s5$maxFeature, coordSub$ID)]
write.csv(table_s5, file = "suppTable5_isoformSwitches_withCoord.csv",
	row.names=FALSE)