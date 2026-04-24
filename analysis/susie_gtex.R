# module load R/4.2.0
# R
# .libPaths()[1]
#
library(data.table)
library(susieR) # 0.14.2
set.seed(1)
target_gene   <- "CBX8"
target_tissue <- "Lung"

# Read in the covariate data.
cov1 <- read.table("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt.gz",
                   header = TRUE,sep = "\t",stringsAsFactors = FALSE)
cov2 <- read.table("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt.gz",
                   header = TRUE,sep = "\t",quote = "",stringsAsFactors = FALSE)
cov2 <- transform(cov2,SUBJID = substr(SAMPID,1,10))
cov <- merge(cov1,cov2,by = "SUBJID")
cov <- cov[c("SUBJID","SAMPID","SEX","AGE","SMTS","SMGEBTCHT","SMTSD")]
cov <- transform(cov,
                 SEX       = SEX - 1,
                 AGE       = factor(AGE),
                 SMTS      = factor(SMTS),
                 SMTSD     = factor(SMTSD),
                 SMGEBTCHT = factor(SMGEBTCHT))
rownames(cov) <- cov$SAMPID

# Read in the gene expression data.
pheno <- fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
               sep = "\t",skip = 2,header = TRUE,showProgress = TRUE)
class(pheno) <- "data.frame"
gene_info <- pheno[1:2]
pheno <- pheno[-(1:2)]
pheno <- as.matrix(pheno)
pheno <- t(pheno)
storage.mode(pheno) <- "double"
colnames(pheno) <- gene_info$Name

# Align the gene expression and covariate data.
ids   <- intersect(cov$SAMPID,rownames(pheno))
rows1 <- match(ids,cov$SAMPID)
rows2 <- match(ids,rownames(pheno))
cov   <- cov[rows1,]
pheno <- pheno[rows2,]

# Extract the gene expression data for the target gene.
j <- which(gene_info$Description == target_gene)
pheno <- cbind(cov,data.frame(count = pheno[,j]))

# Extract the data for the target tissue.
pheno <- subset(pheno,SMTS == target_tissue)
pheno <- transform(pheno,SMGEBTCHT = factor(SMGEBTCHT))

# WGS sequencing platform (HiSeq 2000 or HiSeq X), WGS library
# construction protocol (PCR-based or PCR-free) and donor sex were
# included in the set of covariates used in the association
# analyses.
#
# Expression values for each gene were inverse normal transformed
# across samples.
pheno$y <- resid(lm(count ~ SEX,pheno))

# TO DO NEXT: Read in the genotype data for SNPs near the target gene.
# Use PLINK to automate this.

# Save the fine-mapping data in an .RData file.
# TO DO.
