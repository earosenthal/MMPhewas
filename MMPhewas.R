# Mixed Model Phewas using Phecodes
# step1: Transform phecodes into a VCF file
# step2: Create sparse kinship matrix
# step3: Perform null model
# step4: Perform association

install.packages("devtools")
#It may be necessary to install required as not all package dependencies are installed by devtools:
install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT"))
devtools::install_github("PheWAS/PheWAS")
library(PheWAS)


library(data.table)
library(Hmisc)

#step 1 Transform phecodes into a VCF file
phecodes.dt <- fread("phecode-example.txt")
str(phecodes.dt)
# save the actual names to a vector
phecode.names <- names(phecodes.dt)[-1]
# get usable names for R
setnames(phecodes.dt,make.names(names(phecodes.dt)))
names(phecodes.dt)
phecode.cols <- names(phecodes.dt)[-1]

#create dummy genetic positions
positions <- copy(phecode.names)#the phecode names can be treated as positions.
chr <- rep(1,length(positions))

# Function to replace -9 with NA, make it 1's and 0's 
# instead of 2's and 1's, make it logical
fix.input <- function(vec)
{
  (ifelse(vec<0,NA,vec-1))
}

# Fucntion to convert 0,1,NA to genotype
make.gt <- function(vec)
{
  vec <- ifelse(vec=="0","0/0",vec)
  vec <- ifelse(vec=="1","1/0",vec)
  ifelse(is.na(vec),"./.",vec)
}

phecodes.dt[,(phecode.cols) := lapply(.SD,fix.input), .SDcols=phecode.cols]
phecodes.dt[,(phecode.cols) := lapply(.SD,make.gt), .SDcols=phecode.cols]
str(phecodes.dt)

#need to add columns so that when I output the data, it looks like a vcf
map <- data.frame(CHROM=chr,
                  POS=format(1000*as.numeric(positions),scientific=FALSE),
                  ID=positions,REF=rep("A",length(positions)),
                  ALT=rep("G",length(positions)),
                  QUAL=rep(10,length(positions)),
                  FILTER=rep("PASS",length(positions)),
                  INFO=rep(".",length(positions)),
                  FORMAT=rep("GT",length(positions)))
str(map)

#create the GT portion of the VCF file
ids <- phecodes.dt$id.sample
phecodes.dt[,id.sample:=NULL]
geno <- data.frame(t(as.matrix(phecodes.dt)))
colnames(geno) <- ids
str(geno)
output <- cbind(map,geno)
setnames(output,"CHROM","#CHROM")

write("##fileformat=VCFv4.0",file="phecode-example.vcf")
write.table(output,file="phecode-example.vcf",sep="\t",
            quote=FALSE,row.names=FALSE,
            append=TRUE)

##################
# step 2 Create sparse kinship matrix
#################
library(GENESIS)
relatives.dt <- fread("relatives.txt")
relatives.dt <- relatives.dt[][relatives.dt$ID1 %in% ids & 
                              relatives.dt$ID2 %in% ids]
#attach(relatives.dt)

# create the triplet form that is needed for a sparse matrix
#start with the diagonals, and then add on
#related.ids <- unique(sort(c(ID1,ID2)))
#detach(relatives.dt)
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


##################
# step 3 Perform null model
##################
# requires some sub-steps
# step 3A: read in phenotype data
# step 3B: read in genotype data
# step 3C: Combined phenotype and genotype data into a GDS file
library(SeqVarTools)
library(Biobase) 
# read in phenotype data
phenotype.dt <- fread("phenotype-example.csv")
str(phenotype.dt)
setnames(phenotype.dt,"id.sample","sample.id")
phenotype.dt[,sample.id:=as.character(sample.id)]
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

#Convert into SeqData format
seqData <- SeqVarData(gds)
str(seqData)

# add phenotype data to the seqData 
annot <- AnnotatedDataFrame(my.dat)
sampleData(seqData) <- annot

#run null model
nullmod <- fitNullModel(seqData, outcome = "prs.total",
                        covars=c("is.male","age", paste0("PC", 1:4)), 
                        cov.mat =kin.mat.gen.sparse, 
                        family = "gaussian")

##################
# step 4 Perform association
##################
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
assoc <- assocTestSingle(iterator, nullmod)
assoc





#####################
### removed code
for(i in 1:(length(ids)-1)){
  for(j in (i+1):length(ids)){
    min.id <- min(ids[i],ids[j])
    max.id <- max(ids[i],ids[j])
    if(min.id %in% related.ids & max.id %in% related.ids)
    {
      kin.value <- c(kin.value,KIN[ID1==min.id & ID2 == max.id])
      row.id <- c(row.id,min.id)
      col.id <- c(col.id,max.id)
    }
    
  }
}
kin.dt <- as.data.table(cbind(as.character(row.id),as.character(col.id),
                              as.numeric(kin.value)))
colnames(kin.dt) <- c("ID1","ID2","value")
str(kin.dt)
kin.dt[,value:=as.numeric(value)]
