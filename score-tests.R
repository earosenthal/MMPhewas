# Mixed Model analysis for one phecode against all other phecodes

library(data.table)
library(tidyverse)
library(Hmisc)
library(GENESIS) 
library(SeqVarTools)
library(Biobase) 

#library(bigparallelr) 
sessionInfo()
# read in the PRS and phenotype data

phenotype.dt <- fread("example-prs.csv")
str(phenotype.dt)
ids <- phenotype.dt$id

load("example-phecodes.RData")

phecodes.dt <- as.data.table(phenotypes)
phecodes.dt <- phecodes.dt[][phecodes.dt$id %in% ids]
setnames(phecodes.dt, make.names(names(phecodes.dt)))
phecodes.dt <- phecodes.dt[,list(id,X153)]

phecode.cols <- names(phecodes.dt)[-1]
# Function to convert T/F,NA to binary
make.ind <- function(vec)
{
  case_when(
    vec==FALSE ~ 0,
    vec==TRUE ~ 1,
    is.na(vec) ~ NA_real_
  )
}

phecodes.dt[,(phecode.cols) := lapply(.SD,make.ind), .SDcols=phecode.cols]
str(phecodes.dt)
setnames(phecodes.dt,c("id","CRC.153"))

setkey(phecodes.dt,"id")
setkey(phenotype.dt, "id")
phenotype.dt <- phenotype.dt[phecodes.dt]

#########
### Need to add in age and sex to phenotype.dt
#########
# add in sex to PRS phenotype file
setkey(phenotype.dt,"id")
sex.dt <- fread("input-files/example-sex.csv")
setkey(sex.dt,"id") 

# add in age to PRS phenotype file
# use max.age for now
phenotype.dt <- phenotype.dt[sex.dt]
setkey(phenotype.dt,"id")

ages.dt <- fread("example-age.csv")
setkey(ages.dt, "id") # created ages.dt in create-phecodes.R
phenotype.dt <- phenotype.dt[ages.dt]

setnames(phenotype.dt,"id","sample.id")
#phenotype.dt[,sample.id:=as.character(sample.id)]
my.dat <- as.data.frame(phenotype.dt)

##################
# Create sparse kinship matrix
#################
relatives.dt <- fread("input-files/example-kin.csv")
relatives.dt <- relatives.dt[][relatives.dt$ID1 %in% ids & 
                                 relatives.dt$ID2 %in% ids]

row.id <- ids
col.id <- ids
kin.value <- rep(1,length(ids))
kin.dat <- cbind(row.id,col.id,kin.value)
kin.dat <- rbind(kin.dat,relatives.dt, use.names=FALSE)
kin.dt <- data.table(kin.dat)
colnames(kin.dt) <- c("ID1","ID2","value")
str(kin.dt)
kin.dt[,ID1:=as.character(ID1)]
kin.dt[,ID2:=as.character(ID2)]
kin.dt[,value:=as.numeric(value)]

kin.mat.gen.sparse <- makeSparseMatrix(kin.dt,thresh=NULL)

all.ids <- dimnames(kin.mat.gen.sparse)[[1]]
length(all.ids)


#read in genotype data
#create gds file for output
gdsfile <- "phecode.gds"
vcffile <- "phecode-example.vcf"

#convert vcf to GDS. Account for imputation
seqVCF2GDS(vcffile, gdsfile, fmt.import="DS", storage.option="ZIP_RA",
           verbose=FALSE)

#open the GDS file
gds <- seqOpen(gdsfile)
seqGetData(gds, "sample.id")
#Convert into SeqData format
seqData <- SeqVarData(gds)
str(seqData)

# get my.dat into same order as in the gds file
ordered.ids <- seqGetData(gds, "sample.id")
my.dat <- my.dat[match(ordered.ids,my.dat$sample.id),]
identical(seqGetData(gds, "sample.id"),my.dat$sample.id)

# add phenotype data to the seqData 
annot <- AnnotatedDataFrame(my.dat)
sampleData(seqData) <- annot


nullmod <- fitNullModel(seqData, outcome = "CRC.153", 
	covars=c("sex","max.age"), cov.mat =kin.mat.gen.sparse, family = "binomial")


# description of the model 
nullmod$model
# fixed effect regression estimates
nullmod$fixef
#head(nullmod)

iterator <- SeqVarBlockIterator(seqData, verbose=TRUE)
assoc <- assocTestSingle(iterator, nullmod, test="Score") #I think score tests are the default, but I want to make it clear
str(assoc)
summary(assoc)
head(assoc)
#assoc
save(assoc, file="CRC-assoc.RData")
write.csv(assoc,"CRC-assoc.csv",quote=FALSE, row.names=FALSE)
seqClose(seqData)
