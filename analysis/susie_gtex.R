# module load R/4.2.0
# R
# .libPaths()[1]
#
library(data.table)
library(susieR) # 0.14.2

# Read in the covariate data.
cov1 <- read.table("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt.gz",
                   header = TRUE,sep = "\t",stringsAsFactors = FALSE)
cov2 <- read.table("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt.gz",
                   header = TRUE,sep = "\t",quote = "",stringsAsFactors = FALSE)
cov2 <- transform(cov2,SUBJID = substr(SAMPID,1,10))
cov <- merge(cov1,cov2,by = "SUBJID")
cov <- cov[c("SUBJID","SAMPID","SEX","AGE","SMTS","SMGEBTCHT","SMTSD")]
cov <- transform(cov,
                 SEX       = factor(SEX),
                 AGE       = factor(AGE),
                 SMTS      = factor(SMTS),
                 SMGEBTCHT = factor(SMGEBTCHT))

stop()

# Read in the gene expression data.
pheno <- fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
               sep = "\t",skip = 2,header = TRUE,showProgress = TRUE)
class(pheno) <- "data.frame"

# TO DO. Save the fine-mapping data in an .RData file.

# WGS sequencing platform (HiSeq 2000 or HiSeq X), WGS library
# construction protocol (PCR-based or PCR-free) and donor sex were
# included in the set of covariates used in the association
# analyses.
#
# Expression values for each gene were inverse normal transformed
# across samples.
