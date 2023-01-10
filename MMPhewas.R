# Mixed Model Phewas using Phecodes
# step1: Transform phecodes into a VCF file
# step2: Create sparse kinship matrix
# step3: Perform null model
# step4: Perform association
# step5: Convert association results into format for phewas plotting



library(data.table)
library(tidyverse)
library(Hmisc)

#step 1 Transform phecodes into a VCF file
load("example-phecodes.RData")
phecodes.dt = data.table(phenotypes) #created with create-phecodes.R
# save the actual names to a vector
phecode.names <- names(phecodes.dt)[-1]
# get usable names for R
setnames(phecodes.dt,make.names(names(phecodes.dt)))
names(phecodes.dt)
phecode.cols <- names(phecodes.dt)[-1]

#create dummy genetic positions
positions <- copy(phecode.names)#the phecode names can be treated as positions.
chr <- rep(1,length(positions))

#create a phecode map that will be used later in step 5
map.dt <- data.table(data.frame(pos=format(1000*as.numeric(positions),scientific=FALSE),PHECODE=phecode.names))
map.dt[,pos:=as.numeric(pos)]

#####
# tidyverse function case_when()
#####

# Function to convert T/F,NA to genotype
make.gt <- function(vec)
{
    case_when( 
    vec==FALSE ~ "0/0",
    vec==TRUE ~ "1/0",
    is.na(vec) ~ "./."
  )
}

phecodes.dt[,(phecode.cols) := lapply(.SD,make.gt), .SDcols=phecode.cols]
str(phecodes.dt)

#need to add columns so that when I output the data, it looks like a vcf

map <- data.frame(CHROM=chr,
            POS=format(1000*as.numeric(positions),scientific=FALSE),
            ID=positions,REF="A", ALT="G", QUAL=10, FILTER="PASS",
            INFO=".",FORMAT="GT")
str(map)

#create the GT portion of the VCF file
ids <- phecodes.dt$id
ids <- as.character(ids)

phecodes.dt[,id:=NULL]
geno <- data.frame(t(as.matrix(phecodes.dt)))
colnames(geno) <- ids
str(geno)
output <- cbind(map,geno)
setnames(output,"CHROM","#CHROM")

outfile <- "phecode-example.vcf"
write("##fileformat=VCFv4.0",file=outfile)
write.table(output,file=outfile,sep="\t",
            quote=FALSE, row.names=FALSE, append=TRUE)
# ignore the warning. How do I make it so that the warning about 
# column names doesn't appear?
##################
# step 2 Create sparse kinship matrix
#################
library(GENESIS)
relatives.dt <- fread("input-files/example-kin.csv")
relatives.dt <- relatives.dt[][relatives.dt$ID1 %in% ids & 
                              relatives.dt$ID2 %in% ids]
# need the following when IDs are long
#relatives.dt[,ID1:=format(ID1,scientific=FALSE)]
#relatives.dt[,ID2:=format(ID2,scientific=FALSE)]

# create the triplet form that is needed for a sparse matrix
#start with the diagonals, and then add on

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
# when IDs are long use
#kin.dt[,ID1:=as.character(format(ID1,scientific=FALSE))]
#kin.dt[,ID2:=as.character(format(ID2,scientific=FALSE))]

kin.dt[,value:=as.numeric(value)]


kin.mat.gen.sparse <- makeSparseMatrix(kin.dt,thresh=NULL)

##################
# step 3 Perform null model
##################
# requires some sub-steps
# step 3A: read in phenotype data
# step 3B: read in genotype data
# step 3C: Combined phenotype and genotype data into a GDS file
library(SeqVarTools)
library(Biobase) 
# read in PRS and covariate data
phenotype.dt <- fread("example-prs.csv")
str(phenotype.dt)


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
#setdiff(ids,my.dat$sample.id)
#setdiff(seqGetData(gds, "sample.id"),my.dat$sample.id)
#setdiff(my.dat$sample.id,seqGetData(gds, "sample.id"))
#identical(seqGetData(gds, "sample.id"),my.dat$sample.id)

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


#run null model
nullmod <- fitNullModel(seqData, outcome = "prs.total",
                        covars=c("sex","max.age"), 
                        cov.mat =kin.mat.gen.sparse, 
                        family = "gaussian")

# get fixed effects
output <- data.table(cbind(row.names(nullmod$fixef),nullmod$fixef))
setnames(output,"row.names(nullmod$fixef)","variable")

output[,Est:=round(Est,3)]
output[,pval:=formatC(pval,format="e",digits=2)]
(output <- output[,list(variable,Est,pval)])
write.csv(output,"fixed-effects.csv",quote=FALSE, row.names=FALSE)

##################
# step 4 Perform association
##################
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
assoc <- assocTestSingle(iterator, nullmod)
assoc.dt <- data.table(assoc)

seqClose(seqData)

##################
# step 5 Convert assocation results into format for phewas plotting
##################

### If PheWAS package needs to be installed use the following
#install.packages("devtools")
#It may be necessary to install required as not all package dependencies are installed by devtools:
#install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT"))
#devtools::install_github("PheWAS/PheWAS")

library(PheWAS)

# merge the phecode map with the association results
setkey(assoc.dt,"pos")
setkey(map.dt,"pos")
mmphewas.dt <- map.dt[assoc.dt]

# set the names needed for the parts of the data frame 
# used by phewasManhattan
setnames(mmphewas.dt,c("PHECODE","Est","Est.SE","Score.pval","n.obs"),
         c("phenotype","beta","SE","p","n_total"))
results_d=addPhecodeInfo(mmphewas.dt)
write.csv(results_d,"mmphewas-output.csv",quote=FALSE, row.names=FALSE)

# plot the association results
# Note: max.y is set to 0.75 for the example data.
# Remove the max.y, or change it, to accomodate real data.
png("mmphewas.png", width=3300,height=2550)
phewasManhattan(mmphewas.dt, max.y=0.75,size.x.labels=40, size.y.labels=40, 
                annotate.size=15, point.size=10) +
  labs(title="Mixed Model PheWAS Example") + 
  theme(title=element_text(size=40))
dev.off()



### removed code
# this didn't work for me
#kin.mat.gen.sparse <- makeSparseMatrix(relatives.dt, 
#                                       sample.include=ids,
#                                       diag.value=1, thresh=NULL)
