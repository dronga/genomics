## PCOS

pheno <- read.table("/well/lindgren/alexd/UKBIOBANK/ICD10/ICD10_main_sec_SELFREP_AGE_DOB.tab", header=T, as.is=T)
geno.samples <- read.table("/well/lindgren/alexd/UKBIOBANK/impv1.ids",as.is=T)

icd <- which(apply(pheno,1,function(x){any(c("E282","L680") %in% x)}))
self <- which(apply(pheno,1,function(x){1350 %in% as.numeric(x)}))

comb.pcos <- unique(c(icd,self))
comb.pcos.rem <- comb.pcos[!apply(pheno[comb.pcos,],1,function(x){any(x %in% c("D352","E22","E24","E25","N643","C56","C74"))})]

summary(as.character(pheno[comb.pcos.rem,1]) %in% as.character(geno.samples[,1]))

## POI

pheno <- read.table("/well/lindgren/alexd/UKBIOBANK/ICD10/ICD10_main_sec_SELFREP_INTAGE_MISC.tab", header=T, as.is=T)

ages <- pheno[,grep("20009",colnames(pheno)),]
diseases <- pheno[,grep("20002",colnames(pheno)),]

age.logic <- ages < 40
diseases.logic <- diseases == 1665

poi <- age.logic & diseases.logic

poi.cases.self <- which(apply(poi,1,any))

poi.in <- read.table("POI_include.txt",as.is=T)$V1
poi.ex <- read.table("POI_exclude.txt",as.is=T)$V1

poi.cases.icd <- which(apply(pheno[,grep("4120",colnames(pheno))],1,function(x){any(as.character(x) %in% poi.in)}))

comb.poi <- unique(c(poi.cases.self,poi.cases.icd))
comb.poi.rem <- comb[!apply(pheno[comb,],1,function(x){any(x %in% poi.ex)})]


summary(as.character(pheno[comb.rem,1]) %in% as.character(geno.samples[,1]))

##Age of Menopuase

pheno <- read.table("/well/lindgren/alexd/UKBIOBANK/ICD10/ICD10_main_sec_SELFREP_INTAGE_MISC_AGEMEN.tab", header=T, as.is=T)

pheno[which(pheno[,60]<0),60] <- NA

cases.aom <- which(pheno[,60]<40)

summary(apply(pheno[cases.aom,grep("4120",colnames(pheno))],1,function(x){any(as.character(x) %in% poi.ex)}))

comb.poi <- unique(c(poi.cases.self,poi.cases.icd,cases.aom))
comb.poi.rem <- comb.poi[!apply(pheno[comb.poi,],1,function(x){any(x %in% poi.ex)})]

summary(as.character(pheno[comb.poi.rem,1]) %in% as.character(geno.samples[,1]))

#miscarriage

misc.col <- grep("3839",colnames(pheno))[1]

sum(pheno[,misc.col]>1 & pheno[,misc.col]<3,na.rm=T)




## using HESIN table

load("ukbtemp.RData")

rownames(phe.cases.icd) <- phe.cases.icd$f.eid

hes <- read.table("ukb.tsv", header=T,sep="\t", as.is=T)

dob <- phe.cases.icd[head(as.character(hes$eid)),2:3]
source("age_calc.R")

age_calc(as.Date(paste(dob[,1],dob[,2],sep="")), as.Date(head(hes$epistart)),units='years')
