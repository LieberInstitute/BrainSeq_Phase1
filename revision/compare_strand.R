##

library(derfinder)

## genomic state object
load("/users/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
gsExon = gs[gs$theRegion == "exon"] # just exonic

## total coverage of exonic sequence
totWidth = sum(width(gsExon))
totWidth/1e6

## total width of same-strand overlapping seq
sameStrandOverlap = sum(width(gsExon[lengths(gsExon$gene) > 1]))
sameStrandOverlap/1e6
100*sameStrandOverlap/totWidth

## across strand overlap
gsExonOne = gsExon[lengths(gsExon$gene) == 1]
table(coverage(gsExonOne[strand(gsExonOne) == "+"]))
table(coverage(gsExonOne[strand(gsExonOne) == "-"]))

colSums(table(coverage(gsExonOne)))[3] / 1e6
100*colSums(table(coverage(gsExonOne)))[3] / totWidth

