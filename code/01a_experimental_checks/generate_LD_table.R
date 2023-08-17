library(snpStats)
library(ggplot2)
library(MASS)

base.dir = setwd("~/PROJECTS/PHD/piQTL/code/01_experimental_checks")


genotype_info <- paste0(base.dir,"/../../data/genotype_information/yeast_snps_loc_dec2022.txt")
gsupport <- read.table(genotype_info, header=TRUE, sep=',')


genotype_longfile <- paste0(base.dir,"/../../data/genotype_information/piQTL_genotype_matrix_for_LD_dec2022.txt")
gdata <- read.long(genotype_longfile,
                   fields=c(snp=1, sample=2, genotype=3, confidence=4),
                   gcodes=c("1", "NA", "-1"),
                   threshold=0.95)

summary(gdata)

ld.piQTL <- ld(gdata, stats=c("D.prime", "R.squared"), depth=1000)


#######

pos <- gsupport$position
diags <- vector("list", 1000)
for (i in 1:1000) diags[[i]] <- pos[(i+1):12054] - pos[1:(12054-i)]
dist <- bandSparse(12054, k=1:1000, diagonals=diags)
distance <- dist@x
D.prime <- ld.piQTL$D.prime@x
R.squared <- ld.piQTL$R.squared@x
write_dgCMatrix_csv(ld.piQTL$R.squared, "LD_table.csv", col1_name="SNP")
