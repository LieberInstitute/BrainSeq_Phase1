###

source("../eqtl_functions.R")

# load phenotype data
load("../phenotype_annotated_szControlEqtl_DLPFC.rda")

#####################
# make groups #######
pd$metricGroups = "AdultControl"
pd$metricGroups[pd$Dx == "Schizo"] = "SZ"
pd$metricGroups[pd$Age < 0] = "Fetal"
pd$metricGroups[pd$Age > 0 & pd$Age < 17] = "YoungControl"
pd$metricGroups = factor(pd$metricGroups, 
	levels = c("Fetal","YoungControl","AdultControl","SZ"))
	
# make group indices
groupIndexes=splitit(pd$metricGroups)
dxIndex = unlist(groupIndexes[3:4])

### mapping rate ##
mr = tapply(100*pd$mappingRate, pd$metricGroups, mean)
mrsd = tapply(100*pd$mappingRate, pd$metricGroups, sd)
paste0(round(mr,1), "(", round(mrsd,1), ")")
t.test(pd$mappingRate ~ pd$Dx, subset=dxIndex)$p.value

## number of flow cells
pd$Flowcell_1[pd$Flowcell_1 == "NULL"] = NA
pd$Flowcell_2[pd$Flowcell_2 == "NULL"] = NA
pd$Flowcell_3[pd$Flowcell_3 == "NULL"] = NA
pd$numFlow = as.numeric(!is.na(pd$Flowcell_1)) + 
	as.numeric(!is.na(pd$Flowcell_2)) + 
	as.numeric(!is.na(pd$Flowcell_3))
	
nflow = table(pd$numFlow == 1, pd$metricGroups) 
round(100*prop.table(nflow,2),1)[2,]
chisq.test(nflow[,3:4])

## total mapped
tm = tapply(pd$totalAllMapped/1e6, pd$metricGroups, mean)
tmsd = tapply(pd$totalAllMapped/1e6, pd$metricGroups, sd)
paste0(round(tm,1), "(", round(tmsd,1), ")")
t.test(pd$totalAllMapped/1e6 ~ pd$Dx, subset=dxIndex)$p.value

tmstar = tapply(pd$totalMapped/1e6, pd$metricGroups, mean)
tmstarsd = tapply(pd$totalMapped/1e6, pd$metricGroups, sd)
paste0(round(tmstar,1), "(", round(tmstarsd,1), ")")
t.test(pd$totalMapped/1e6 ~ pd$Dx, subset=dxIndex)$p.value

# mito rate
mito = tapply(100*pd$mitoRate, pd$metricGroups, mean)
mitosd = tapply(100*pd$mitoRate, pd$metricGroups, sd)
paste0(round(mito,1), "(", round(mitosd,1), ")")
t.test(pd$mitoRate ~ pd$Dx, subset=dxIndex)$p.value

# gene assignment rate
summ = read.delim("/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/Counts/Ensembl_Exons_n746.counts.summary",
	row.names=1)
colnames(summ) = ss(colnames(summ), "_", 4)
summ = summ[,pd$RNum]
assignRate = as.numeric(summ[1,] / colSums(summ))

ar = tapply(100*assignRate, pd$metricGroups, mean)
arsd = tapply(100*assignRate, pd$metricGroups, sd)
paste0(round(ar,1), "(", round(arsd,1), ")")
t.test(assignRate ~ pd$Dx, subset=dxIndex)$p.value