##

library(GenomicRanges)
library(jaffelab)

## load results
load("rdas/signif_dev_stats_plusBrainSpan.rda")

table(!is.na(statListSig$Junction$meanRpkm_span))
## compare expression levels
corExprsLevel = sapply(statListSig, function(x) 
	cor(log2(x$meanRpkm+1), log2(x$meanRpkm_span+1),use="comp"))
signif(corExprsLevel,3)

## plot
pdf("plots/brainspan_expression_levels.pdf")
par(mar=c(5,6,2,2), cex.lab=2, cex.main=2,cex.axis=2)
for(i in seq(along=statListSig)) {
	x = statListSig[[i]]
	plot(log2(x$meanRpkm+1), log2(x$meanRpkm_span+1),
		xlab="LIBD - Discovery", ylab="BrainSpan - Replication",
		pch = 21, bg="grey", 
		main = paste0(names(statListSig)[i], " - Log2 Exprs Level"))
}
dev.off()

##################
## what proportion are significant?
sapply(statListSig, function(x) mean(x$f_pval_span < 0.05,na.rm=TRUE))
sapply(statListSig, function(x) mean(x$pval_fetal_span < 0.05,na.rm=TRUE))

## and expressed in span?
sapply(statListSig, function(x) 
	mean(x$f_pval_span[x$meanRpkm_span > 0.1] < 0.05,na.rm=TRUE))
sapply(statListSig, function(x) 
	mean(x$pval_fetal_span[x$meanRpkm_span > 0.1]< 0.05,na.rm=TRUE))

########
## age correlations

corAgeCorr = sapply(statListSig, function(x) 
	cor(x$ageCorr, x$ageCorr_span,use="comp"))
kappaAgeCorr = sapply(statListSig, function(x) 
	table(sign(x$ageCorr) == sign(x$ageCorr_span)))
prop.table(kappaAgeCorr, 2)

plot(statListSig$Gene$ageCorr, statListSig$Gene$ageCorr_span)


##########
## iso switch
load("rdas/isoform_switch_devel_byFeature.rda")

switchList = mapply(function(x,y) {
	cat(".")
	x$negCorr_span = y$ageCorr_span[match(x$minFeature, names(y))]
	x$posCorr_span = y$ageCorr_span[match(x$maxFeature, names(y))]
	x$corDiff_span = x$posCorr_span - x$negCorr_span
	return(x)
}, x = switchList, y=statListSig[-1], SIMPLIFY=FALSE)

## check replication rates by cut
corCut= seq(0,1,by=0.01)
checkCor = t(sapply(corCut, function(cc) {
	sapply(switchList, function(x) {
		xx = x[x$corDiff > cc,]
		length(which(xx$negCorr_span < 0 & 
			xx$posCorr_span > 0 ))/nrow(xx)
	})
}))

pdf("plots/correlation_difference_in_brainspan.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2)
palette(brewer.pal(5,"Dark2"))
matplot(corCut, checkCor,type="l", col=1:4, lty=1, lwd=3,
	xlab="Correlation Difference", ylab = "Replication Proportion")
legend("topleft", colnames(checkCor), col = 1:4, 
	pch = 15, cex=1.5, pt.cex=1.5)
dev.off()

num = t(sapply(corCut, function(cc) {
	sapply(switchList, function(x) {
		nrow(x[x$corDiff > cc,])})}))

matplot(corCut, num)

switchTabList = lapply(switchList, function(x) 
	table(x$negCorr_span < 0 & x$posCorr_span > 0,
		useNA="always"))
switchTabList
switchTab = sapply(switchList, function(x) table(x$corDiff_span > 0,useNA="always"))

prop.table(switchTab[-3,],2)
plot(switchList$Exon$corDiff, switchList$Exon$corDiff_span)