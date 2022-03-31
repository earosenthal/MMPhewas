# calculate PRS
library(data.table)
library(Hmisc)
library(bigsnpr)
library(dplyr)
library(stringr)

# function to get the complement allele
complement.allele <- function(allele)
{
  out <- ifelse(allele=="A","T",NA)
  out <- ifelse(allele=="T","A",out)
  out <- ifelse(allele=="C","G",out)
  out <- ifelse(allele=="G","C",out)
  out
}

infile <- "example"
bedfile <- paste0(infile,".bed")
rds <- snp_readBed(bedfile,backingfile = tempfile())
obj.bigSNP <- snp_attach(rds)


#get the alpha genotype data
ID <- obj.bigSNP$fam$sample.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
G   <- obj.bigSNP$genotypes
A1 <- obj.bigSNP$map$allele1
A2 <- obj.bigSNP$map$allele2
alpha.snps.dt <- data.table(cbind(CHR,POS,A1,A2))
alpha.snps.dt[,CHR:=str_replace(CHR,"chr","")]
alpha.snps.dt[,CHR:=as.integer(CHR)]
alpha.snps.dt[,POS:=as.integer(POS)]


# get the SNP betas
infile <- "example-beta.csv"
prs.snps <- fread(infile)

#a1 is risk allele
#a0 is other allele

#combine SNP info from alpha and beta files
setkeyv(alpha.snps.dt,c("CHR","POS"))
setkeyv(prs.snps,c("CHR","POS"))
betas.dt <- prs.snps[alpha.snps.dt]

# add in the complement alleles
betas.dt[,a0c:=complement.allele(a0)]
betas.dt[,a1c:=complement.allele(a1)]

# zero out the betas for the indels so that they are not included in the analysis
# the example does not include indels at the moment
betas.dt[,num.nucs:=str_length(A1) + str_length(A2)]
betas.dt[,beta.fix:=ifelse(num.nucs>2,0,beta)]

# determine which positions have duplicates. some of these will correspond to the indels
dup.pos <- betas.dt$POS[duplicated(betas.dt$POS)]

# remove (zero out) the duplicated positions with mismatched alleles.
betas.dt[,rm.dup:=ifelse( POS %in% dup.pos & !( (a0==A1 & a1==A2) | (a0==A2 & a1==A1) | 
                                                  (a0c==A1 & a1c==A2) | (a0c==A2 & a1c==A1)), 1, 0)]
betas.dt[,beta.fix:=ifelse(rm.dup==1,0,beta.fix)]

# set up the data tables that go into snp_match(), which corrects the betas automatically
sumstats <- betas.dt[,list(chr=CHR,pos=POS,a0,a1,beta=beta.fix)]
info_snp <- betas.dt[,list(chr=CHR,pos=POS,a0=A2,a1=A1)] 

# run the snp_match()
beta.match.dt <- data.table(snp_match(sumstats,info_snp,return_flip_and_rev=TRUE,remove_dups=FALSE))

#check for SNPs that don't match
if(length(betas.dt$POS) > length(beta.match.dt$pos))
{ #need to zero out the beta for these SNPs and put back into the data table in correct order
  missing.pos <- setdiff(info_snp$pos,beta.match.dt$pos)
  keep.cols <- c("CHR","POS","A1","A2")
  add.rows <- betas.dt[,.SD,.SDcols=keep.cols][betas.dt$POS %in% missing.pos]
  setnames(add.rows,c("chr","pos","a0","a1"))
  add.rows[,beta:=0]
  add.rows[,"_NUM_ID_.ss":=-1]
  add.rows[,"_FLIP_":=FALSE]
  add.rows[,"_REV_":=FALSE]
  add.rows[,"_NUM_ID_":=-1]
  beta.match.dt <- rbind(add.rows,beta.match.dt)
  setkeyv(beta.match.dt,c("pos","_NUM_ID_.ss"))
}

# impute missing genoytpes useing the mode genotype
G2 <- snp_fastImputeSimple(G,method="mode")

#calculate the chr. PRS using the betas from beta.match.dt
prs <- snp_PRS(G=G2, betas.keep = beta.match.dt$beta)[,1]
prs.dt <- data.table(cbind(ID,prs))
setnames(prs.dt,c("id","prs.total"))

outfile <- "example-prs.csv"
write.csv(prs.dt,outfile,quote=FALSE,row.names=FALSE)

rm(rds)
system("rm -r Rtmp*")
