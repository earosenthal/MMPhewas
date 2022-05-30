# create the phecode file that can be used in the MMPhewas

# remove any variables
rm(list=ls())
library(data.table)
library(PheWAS)

icd.dt <- fread("input-files/example-icdcodes.csv")
#icd.dt[,vocabulary_id:=ifelse(ICD_FLAG==9,"ICD9CM",NA)]
#icd.dt[,vocabulary_id:=ifelse(ICD_FLAG==10,"ICD10CM",vocabulary_id)]
#icd.dt[,ICD_FLAG:=NULL]
#setnames(icd.dt,c("id","code","age","vocabulary_id"))

# get counts of each code
setkeyv(icd.dt,c("id","code"))
icd.dt[,count:=length(age),by=c("id","code")]

#pull in the sex
sex.dt <- fread("input-files/example-sex.csv")

# convert the icd codes to phecodes
id.vocab.code.count <- icd.dt[,list(id,vocabulary_id,code,count)]
id.vocab.code.count <- unique(id.vocab.code.count,by=c("id","vocabulary_id","code","count"))
phenotypes = createPhenotypes(id.vocab.code.count, min.code.count=2,
                              id.sex=sex.dt,aggregate.fun=sum)

save(phenotypes, file="example-phecodes.RData")



# get minimum age and maximum age overall
setkeyv(icd.dt, c("id","age"))
icd.dt[,min.age:=min(age),by=id]
icd.dt[,max.age:=max(age),by=id]

ages.dt <- icd.dt[,list(id,min.age,max.age)]
ages.dt <- unique(ages.dt, by="id")
icd.dt[,max.age:=NULL]
icd.dt[,min.age:=NULL]

write.csv(ages.dt,"example-age.csv",quote=FALSE,row.names=TRUE)

# stop here
# Below is code I was playing around with
# create the phenotype table

# get counts of the phecodes for each participant
#setkeyv(phecode.dt,c("id","phecode"))
#phecode.dt[,count.code:=length(age),by=c("id","phecode")]
#phecode.dt <- unique(phecode.dt, by=c("id","phecode"))
#phecode.dt[,age:=NULL]



#phecode.table <- createPhenotypes(icd.dt)

#phecode.dt <- mapCodesToPhecodes((icd.dt))
