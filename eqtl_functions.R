##
getOR = function(x) x[1,1]/x[2,1]/x[1,2]*x[2,2]

# wrapper for string split and sapply
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

# splitting variables into list
splitit = function(x) split(seq(along=x),x) # splits into list
split0 = function(x) splitit(factor(x, levels = unique(x)))

## pca vars
getPcaVars = function(pca)  signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100

## plot transcripts
plotTranscripts = function(gr, txdb,...) {
	require(derfinderPlot)
	require(GenomicFeatures)
	tx <- exonsBy(txdb)
	ov <- findOverlaps(gr, tx)
	txList <- split(tx[subjectHits(ov)], queryHits(ov))
	poly.data <- lapply(txList, derfinderPlot:::.plotData)[[1]]
	yrange <- range(poly.data$exon$y)
	xrange <- c(start(gr),end(gr))
    plot(0,0, type="n", xlim = xrange , ylim = yrange + c(-0.75, 0.75),
		yaxt="n", ylab="", xlab=as.character(seqnames(gr)),
		cex.axis = 1.5, cex.lab =1.5,...)
	derfinderPlot:::.plotPolygon(poly.data$exon, 'blue')
    derfinderPlot:::.plotPolygon(poly.data$intron, 'lightblue')
	yTx <- unique(yrange / abs(yrange))
	if(length(yTx) > 1) {
		axis(2, yTx, c('-', '+')[c(-1, 1) %in% yTx], 
			tick = FALSE, las = 1, cex.axis = 3)
		abline(h = 0, lty = 3)
	}
}

## get total mapped
getTotalMapped = function(bamFile, mc.cores=1, returnM = TRUE) {
	thecall = paste("samtools idxstats",bamFile)
	tmp = parallel::mclapply(thecall, function(x) {
		cat(".")
		xx = system(x,intern=TRUE)
		xx = do.call("rbind", strsplit(xx, "\t"))
		d = data.frame(chr=xx[,1], L=xx[,2], mapped = xx[,3],
			stringsAsFactors=FALSE)
		d
	},mc.cores=mc.cores)
	
	out = list(totalMapped = sapply(tmp, function(x) sum(as.numeric(x$mapped[x$chr %in% paste0("chr", c(1:22,"X","Y"))]))),
		mitoMapped = sapply(tmp, function(x) as.numeric(x$mapped[x$chr=="chrM"])))
	return(out)
}

###########
junctionCount = function(junctionFiles, sampleNames=names(junctionFiles), 
	output = c("Count", "Rail"), minOverhang=0, 
	strandSpecific = FALSE, illuminaStranded=FALSE,
	minCount = 1, maxCores=NULL) {
	
	require(GenomicRanges,quietly=TRUE)
	require(parallel,quietly=TRUE)

	if(is.null(maxCores)) {
		maxCores=1
	}
		
	names(junctionFiles) = sampleNames
	cat("Reading in data.\n")
	if(all(is.character(junctionFiles))) {
		theData = mclapply(junctionFiles, function(x) {
			if(output == "Rail") {
				y = read.delim(x, skip = 1, header=FALSE, 
					colClasses = c("character", "integer", 
					"integer", "integer", "integer", "integer"))
				colnames(y) = c("chr", "start", "end", "leftHang", "rightHang", "count")
				y = y[y$count >= minCount,] # filter based on min number
				y = y[y$leftHang > minOverhang & y$rightHang > minOverhang,]
			} else if(output == "Count") {
				y = read.delim(x, skip = 1, header=FALSE, 
				col.names = c("chr", "start","end", "strand", "count"), 
				colClasses = c("character", "integer", "integer", "character","integer"))
				y = y[y$count >= minCount,] # filter based on min number
			} else stop("Junction formats can only be from Tophat and Rail.\n")
			
			gr = GRanges(y$chr, IRanges(y$start, y$end), 
				strand=y$strand,count = y$count)
			return(gr)
		}, mc.cores=maxCores)
	} else {
		theData = junctionFiles
		stopifnot(all(sapply(theData, class)=="GRanges"))
	}
	cat("Creating master table of junctions.\n")

	## turn into GRangesList
	### THIS STEP IS SLOW...
	grList = GRangesList(theData)

	# each dataset should be checked
	if(illuminaStranded & strandSpecific) {
		grList = GRangesList(mclapply(grList, function(x) {
			strand(x) = ifelse(strand(x)=="+", "-","+")
			return(x)
		},mc.cores=maxCores))
	}
	
	## get into GRanges object of unique junctions
	fullGR = unlist(grList)
	if(!strandSpecific) strand(fullGR) = "*"
	
	fullGR = fullGR[!duplicated(fullGR)] # or unique(fullGR)
	fullGR = sort(fullGR)
	fullGR$count = NULL

	cat(paste("There are", length(fullGR), "total junctions.\n"))
	
	cat("Populating count matrix.\n")

	jNames = paste0(as.character(seqnames(fullGR)),":",
			start(fullGR),"-",end(fullGR),"(",as.character(strand(fullGR)),")")

	## match GRanges
	options(warn=-1)
	mList = mclapply(grList, match, fullGR, 
		ignore.strand = !strandSpecific, mc.cores=maxCores)
	options(warn=0)
	
	countList = mList # initiate 
	M = length(jNames)

	## fill in matrix
	for(i in seq(along=grList)) {
		if(i %% 25 == 0) cat(".")
		cc = rep(0,M)
		cc[mList[[i]]] = theData[[i]]$count
		countList[[i]] = Rle(cc)
	}
	countDF = DataFrame(countList, row.names=jNames)
	
	names(fullGR) = jNames
	## return matrix and GRanges object
	out = list(countDF = countDF, anno = fullGR)
	return(out)
}

### annotate junctions
annotateJunctions = function(juncCounts, build="hg19") {
	if(build == "hg19") {
		load("/users/ajaffe/Lieber/Projects/RNAseq/ensembl_v75_junction_annotation.rda")
	} else stop("Only supports hg19 for now.\n")
	anno = juncCounts$anno

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

	library(bumphunter)
	an = annotateNearest(anno, theJunctions)
	anno$nearestSymbol = as.character(theJunctions$symbol[an$matchIndex])
	anno$nearestDist = an$dist

	# put back in
	juncCounts$anno = anno
	return(juncCounts)
}

## junction stats
junctionStats = function(jRpkm, jMap, 
	cuts = c(0,0.5,1,5,10,20,50,100), output="percent") {

	require(GenomicRanges,quietly = TRUE)
	## junction means
	junctionMeans = rowMeans(jRpkm)
	numAboveCut = sapply(cuts, function(x) sum(junctionMeans >= x))

	## junction stats
	theSeq = paste0(jMap$leftSeq, ":", jMap$rightSeq)
	canJ = ifelse(theSeq %in% c("GT:AG", "CT:AC"), "Canonical", "Not")

	juncList = lapply(cuts, function(x) {
		cat(".")
		expIndex=rep(FALSE, length(junctionMeans))
		expIndex[which(junctionMeans > x)] = TRUE
		
		tabIn = table(jMap$inEnsembl, expIndex, dnn = c("Ensembl", "Exprs"))
		tabNear = table((jMap$inEnsemblEnd | jMap$inEnsemblStart) & 
			!(jMap$inEnsemblEnd & jMap$inEnsemblStart) & !jMap$inEnsembl, 
			expIndex, dnn = c("Ensembl", "Exprs"))
		tabNovel = table(jMap$inEnsemblEnd & jMap$inEnsemblStart & !jMap$inEnsembl, 
			expIndex, dnn = c("Ensembl", "Exprs"))
		tabNew = table(!jMap$inEnsemblEnd & !jMap$inEnsemblStart, 
			expIndex, dnn = c("Ensembl", "Exprs"))
		
		canon = mean(canJ[expIndex] == "Canonical")
		list(tabIn = tabIn, tabNear = tabNear,
			tabNovel = tabNovel,tabNew=tabNew,canon = canon)
	})
	if(output == "percent") {
		ensOut = data.frame(numExprs = numAboveCut, knownTrans = sapply(juncList, 
				function(x) x$tabIn["TRUE","TRUE"]/sum(x$tabIn[,"TRUE"])),
			novelTrans = sapply(juncList, 
				function(x) x$tabNovel["TRUE","TRUE"]/sum(x$tabNovel[,"TRUE"])),
			novelJxn = sapply(juncList, 
				function(x) x$tabNear["TRUE","TRUE"]/sum(x$tabNear[,"TRUE"])),
			newJxn = sapply(juncList, 
				function(x) x$tabNew["TRUE","TRUE"]/sum(x$tabNew[,"TRUE"])),
			canon = sapply(juncList, function(x) x$canon))
		ensOut[,-1] = round(ensOut[,-1]*100,2)
	} else if(output == "number") {
		ensOut = data.frame(numExprs = numAboveCut, 
			knownTrans = sapply(juncList, function(x) x$tabIn["TRUE","TRUE"]),
			novelTrans = sapply(juncList, function(x) x$tabNovel["TRUE","TRUE"]),
			novelJxn = sapply(juncList, function(x) x$tabNear["TRUE","TRUE"]),
			newJxn = sapply(juncList, 	function(x) x$tabNew["TRUE","TRUE"]),
			canon = sapply(juncList, function(x) x$canon))
	} else stop("'output' must be 'percent' or 'number'.\n")		
	rownames(ensOut) = cuts
	return(ensOut)
}

	
## age plotter
agePlotter = function(y, age, mod = matrix(rep(1,length(y)),ncol=1),
	mainText, smoothIt=TRUE, jitter=TRUE, ageLabel = "bottom",
	orderByAge=TRUE,ylim=NULL,ageBreaks = c(-1, 0, 1, 10, 100), 
	ylab="Adjusted Expression",pointColor = 2, lineColor= 1, alreadyFitted=NULL,
		...) {
	
	if(orderByAge) {
		oo = order(age, decreasing = FALSE)
		y = y[oo] ; age = age[oo] ; mod = mod[oo,]
		if(!is.null(alreadyFitted)) alreadyFitted = alreadyFitted[oo]
	}
	
	if(is.null(alreadyFitted)) {
		fit = fitted(lm(y~mod-1))
	} else fit = alreadyFitted
	
	fetal = cut(age, breaks =ageBreaks ,lab=FALSE)
	fIndex = splitit(fetal)	
	
	layout(matrix(c(1,1,1,1,1,2,2,3,3,4,4,4,4),nr = 1,byrow = TRUE))
	palette(brewer.pal(8,"Set1"))
	
	par(mar = c(4,5,3,0.45))
	if(is.null(ylim)) ylims = range(y,na.rm=TRUE) else ylims = ylim

	if(jitter) xx = jitter(age,amount=0.005) else xx=age
	plot(y ~ xx,
		subset=fIndex[[1]],
		main = "",	ylab=ylab,xlab="",
		ylim = ylims,cex.axis = 1.5, cex.lab=1.75,
		pch = 21, cex = 1.4,xaxt="n",bg = pointColor,
		xlim=c(range(age[fIndex[[1]]])+c(-0.01,0.07)),...)
	
	if(smoothIt) lines(age[fIndex[[1]]],fit[fIndex[[1]]],col=lineColor,lwd=6)
	
	axis(1,at=unique(age[fIndex[[1]]]), 
		labels = 40+52*signif(unique(age[fIndex[[1]]]),1), cex.axis=1.5)
	
	if(ageLabel == "bottom") {
		text(x = quantile(age[fIndex[[1]]],0.33), y= min(ylims), "PCW", cex=1.5)
	} else if(ageLabel == "top") {
		text(x = quantile(age[fIndex[[1]]],0.33), y= max(ylims), "PCW", cex=1.5)
	} 
	
	# infant + child
	par(mar = c(4, 0.25,3,0.25))
	for(j in 2:3) {
		plot(y ~ age,	subset=fIndex[[j]],
			main = "",ylab="",xlab="",yaxt = "n", cex=1.4,
			xlim = range(age[fIndex[[j]]])+c(-0.03,0.03),
			ylim = ylims, cex.axis = 1.5,pch = 21,  bg=pointColor)

		if(ageBreaks[2] == 0 & smoothIt) lines(age[fIndex[[j]]],fit[fIndex[[j]]],col=lineColor,lwd=6)
		if(ageBreaks[2] < 0 & smoothIt) lines(age[fIndex[[j]]][-1],fit[fIndex[[j]]][-1],col=lineColor,lwd=6)
	}
	
	# adults
	par(mar = c(4, 0.25,3,1))
	plot(y ~ age,	subset=fIndex[[4]],
			main = "",ylab="",xlab="",yaxt = "n", cex=1.4,
			xlim = range(age[fIndex[[4]]])+c(-0.01,0.01),
			ylim = ylims, cex.axis = 1.5,pch = 21, bg=pointColor)

	if(smoothIt) lines(age[fIndex[[4]]],fit[fIndex[[4]]],col=lineColor,lwd=6)

	mtext(mainText,	outer=T, line=-2.5,cex=1.35)	
	
	mtext("Age", side=1, outer=T, line=-1.5,cex=1.35)
}


# get the f statistic from 2 lmFit objects
getF = function(fit, fit0, theData) {
	
	rss1 = rowSums((fitted(fit)-theData)^2)
	df1 = ncol(fit$coef)
	rss0 = rowSums((fitted(fit0)-theData)^2)
	df0 = ncol(fit0$coef)

	fstat = ((rss0-rss1)/(df1-df0))/(rss1/(ncol(theData)-df1))
	f_pval = pf(fstat, df1-1, ncol(theData)-df1,lower.tail=FALSE)
	fout = cbind(fstat,df1-1,ncol(theData)-df1,f_pval)
	colnames(fout)[2:3] = c("df1","df0")
	fout = data.frame(fout)
	return(fout)
}

# get R2 from matrix via limma
# mod0 must be a nested model, within mod, a la anova
getR2 = function(p,mod,mod0=NULL) {
	require(limma)
	
	fit1 = lmFit(p,mod)
	rss1 = rowSums((p-fitted(fit1))^2)
	n = ncol(p)
	k = ncol(mod)-1
	
	if(is.null(mod0)) {
		rss0 = rowSums((p-rowMeans(p))^2)
	} else {
		fit0 = lmFit(p,mod0)
		rss0 = rowSums((p-fitted(fit0))^2)
	}
	
	r2 = 1 - (rss1/rss0)
	r2adj = 1 - ((1-r2)*(n-1))/(n-k-1)
	out = data.frame(R2 = r2,Adjusted_R2 = r2adj)
	return(out)	
}

## regress out SVs/PCs
cleaningY = function(y, mod, P=ncol(mod)) {
	Hat=solve(t(mod)%*%mod)%*%t(mod)
	beta=(Hat%*%t(y))
	cleany=y-t(as.matrix(mod[,-c(1:P)])%*%beta[-c(1:P),])
	return(cleany)
}

