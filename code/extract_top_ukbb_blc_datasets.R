# I ran this script on randi.cri.uchicago.edu.
# module load gcc/12.1.0 R/4.2.1
n <- 100
datasets <- read.csv("../data/plink_ukbb_blood_cell_traits.csv",header = TRUE,
                     skip = 1,stringsAsFactors = FALSE)
datasets <- transform(datasets,
                      chr       = factor(chr),
                      top_trait = factor(top_trait))
datasets <-
  transform(datasets,
            plink_file = sprintf("bloodcells_chr%d.%d.%d.z.rds",chr,start,end),
            ld_file = sprintf("bloodcells_chr%d.%d.%d.matrix.gz",chr,start,end))
rows <- order(datasets$max_log10P,decreasing = TRUE)
rows <- rows[1:n]
rows <- sort(rows)
datasets <- datasets[rows,]

# Save the list of the top data sets.
write.csv(datasets,"plink_ukbb_blood_cell_traits_top100.csv",
          quote = FALSE,row.names = FALSE)

# Copy the PLINK files.
fromdir <- file.path("/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw",
                     "BloodCells/regions_zscores_maf001_info6")
destdir <- "/scratch/pcarbone/ukbb_blood_cell_traits/plink"
for (i in 1:n) {
  system(sprintf("cp -f %s/%s %s/",fromdir,datasets[i,"plink_file"],destdir))
}

# Copy the LD files.
fromdir <- file.path("/gpfs/data/stephens-lab/finemap-uk-biobank/data/raw",
                     "BloodCells/regions_ld_maf001_info6")
destdir <- "/scratch/pcarbone/ukbb_blood_cell_traits/ldstore"
for (i in 1:n) {
  system(sprintf("cp -f %s/%s %s/",fromdir,datasets[i,"ld_file"],destdir))
}
