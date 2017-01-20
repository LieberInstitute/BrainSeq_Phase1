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

# number
N = sapply(groupIndexes,length)
N

# sex
sex = sapply(groupIndexes, function(x) table(pd$Sex[x]))
sexF= signif(prop.table(sex,2),3)[1,]
sexF*100
chisq.test(table(pd$Sex[dxIndex], as.character(pd$Dx[dxIndex])))$p.value

# race
race = sapply(groupIndexes, function(x) table(pd$Race[x]))
raceCauc =  signif(prop.table(race,2),3)["CAUC",]
raceCauc*100
chisq.test(table(pd$Race[dxIndex], as.character(pd$Dx[dxIndex])))$p.value

### age
age = tapply(pd$Age, pd$metricGroups, mean)
agesd = tapply(pd$Age, pd$metricGroups, sd)
paste0(signif(age,3), "(", signif(agesd,3), ")")
t.test(pd$Age ~ pd$Dx, subset=dxIndex)$p.value

### RIN
rin = tapply(pd$RIN, pd$metricGroups, mean)
rinsd = tapply(pd$RIN, pd$metricGroups, sd)
paste0(round(rin,1), "(", round(rinsd,1), ")")
t.test(pd$RIN ~ pd$Dx, subset=dxIndex)$p.value

### PMI
pmi = tapply(pd$PMI, pd$metricGroups, mean,na.rm=TRUE)
pmisd = tapply(pd$PMI, pd$metricGroups,sd,na.rm=TRUE)
paste0(signif(pmi,3), "(", signif(pmisd,3), ")")
t.test(pd$PMI ~ pd$Dx, subset=dxIndex)$p.value

## smoking at death
smoke = sapply(groupIndexes, function(x) table(factor(pd$SmokingEither)[x]))
smokeY= signif(prop.table(smoke,2),3)[2,]
smokeY*100
chisq.test(smoke[,3:4])$p.value

## cause of death
death = table(pd$Manner, pd$metricGroups)
deathNat = signif(prop.table(death,2),3)["Natural",]
deathAcc = signif(prop.table(death,2),3)["Accident",]
deathSuic = signif(prop.table(death,2),3)["Suicide",]
deathNat*100
deathAcc*100
deathSuic*100
chisq.test(death[c(1,3),3:4])$p.value


# age of onset
mean(pd$AgeOnsetSchizo[pd$Dx == "Schizo"],na.rm=TRUE)
sd(pd$AgeOnsetSchizo[pd$Dx == "Schizo"],na.rm=TRUE)
prop.table(table(pd$Antipsychotics[pd$Dx == "Schizo"]))
