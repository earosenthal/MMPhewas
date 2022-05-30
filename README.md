# MMPhewas
This repo contains code and example input files needed to run a mixed model Phewas using a polygenic risk score (PRS) as the response variable and the phecodes (and other covariates) as the predictors. A mixed model is used, rather than a typical linear model, so that relatedness among participants can be accounted for in the data. This improves power, as related participants do not need to be removed from the analysis. 

The phecodes are transformed into a VCF file with genotypes '0/0' for FALSE (phecode not present) and '0/1' for TRUE (phecode is present) in order to use the mixed model provided in the GENESIS package. A mixed model GWAS is then run to get the association statistics, adjusting for relatedness in the data using a kinship covariance matrix. 
The output is then renamed to correspond with the column names used in the PheWAS package. Then plots and other observations can be made with the output. 

R scripts

1. calculate-prs.R calculates the PRS using genotype data in PLINK format and a file with the SNP effects, specifying which allele has the effect.

2. create-phecodes.R calculate the phecodes given ICD 9 and 10 code data, assuming that the code occurs at least twice (different clinic visit). This script also gets the minimum and maximum age for each participant. 

3. MMPhewas.R runs the analysis. 

