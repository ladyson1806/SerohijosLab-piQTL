library(rMVP)
library(data.table)
library(tidyverse)


#### 1. Loading and formatting initial inputs for rMVP

setwd("/home/savvy/PROJECTS/PHD/piQTL_maping/code/pipeline/")

MVP.Data(fileNum="../../data/genotype_information/rMVP/Genotype.txt",
         filePhe="../../results/03_ppi_estimation/logratio/all_PPI_logratio_fitness_before_downsampling_minus_ref_for_rMVP.csv",
         fileMap="../../data/genotype_information/rMVP/Map.txt",
         sep.num="\t",
         sep.map="\t", 
         sep.phe=",",
         fileKin=TRUE,
         filePC=TRUE,
         priority="speed",
         out="../../results/04a_rMVP/formatted_inputs/piQTL_dec_2022/piQTL_dec_2022"
         )


############################################

##### 2. PCA analysis

library(PCAtools)
library(data.table)
library(tidyverse)
library(pheatmap)

setwd("/home/savvy/PROJECTS/PHD/piQTL/results/04a_rMVP/formatted_inputs/piQTL_dec_2022/")


###### Loading inputs
genotype = attach.big.matrix("piQTL_dec_2022.geno.desc")
genotype = genotype[,]
dim(genotype)

phenotypes = fread("piQTL_dec_2022.phe", data.table=FALSE)

kinship_matrix = attach.big.matrix("piQTL_dec_2022.kin.desc")
kinship_matrix = kinship_matrix[,]

colnames(genotype) = phenotypes$Taxa
rownames(phenotypes) = phenotypes$Taxa

###### Run PCA and add results in a dataframe
p_snp = pca(genotype, center = TRUE, scale=TRUE,metadata=phenotypes)

pca.data = data.frame(Genotype=p_snp$metadata$Taxa,
                      PC1 = p_snp$rotated$PC1,
                      PC2 = p_snp$rotated$PC2,
                      PC3 = p_snp$rotated$PC3,
                      stringsAsFactors = FALSE)

###### Print scatterplot of the PCA 

pca.percent = p_snp$variance
head(round(pca.percent, 2))
sum(pca.percent[1:3])

setwd("/home/savvy/PROJECTS/PHD/piQTL/results/04a_rMVP/figures/")
pdf("piQTL_PCA_matrix_dec_2022.pdf", width = 10, height = 9)
lbls = paste("PC", 1:3, "\n", format(pca.percent[1:3], digits=2), "%", sep="")
pairs(pca.data[,2:4], labels=lbls, pch = 19,  cex = 0.9, lower.panel = NULL)
dev.off()


###### Heatmap of the Kinship matrix
pdf("piQTL_Kinship_Matrix_noNAs_dec_2022.pdf", width = 10, height = 9)
pheatmap(kinship_matrix)
dev.off()

############################################

#### 3. Run QTL mapping 

setwd("/home/savvy/PROJECTS/PHD/piQTL/results/04a_rMVP/formatted_inputs/piQTL_dec_2022/")

#### Read in phenotypes
phenotype = read.table("piQTL_dec_2022.phe", head=TRUE)
#### Read in SNPs
genotype = attach.big.matrix("piQTL_dec_2022.geno.desc")
#### Read in SNPs map
map = read.table("piQTL_dec_2022.geno.map", head = TRUE)


setwd("/home/savvy/PROJECTS/PHD/piQTL/results/04a_rMVP/rMVP_outputs/piQTL_dec_2022/")

for(i in 2:ncol(phenotype)){
  imMVP <- MVP(
    phe=phenotype[, c(1, i)],
    geno=genotype,
    map=map,
    nPC.GLM=2,
    priority="speed",
    vc.method="BRENT",
    maxLoop=10,
    method.bin="FaST-LMM",
    bin.size=1e3,
    threshold=0.05,
    method=c("GLM"), 
    file.output=FALSE
  )
  gc()
  title = colnames(phenotype[, c(1, i)])[2]
  all.pvals = cbind(imMVP$map, imMVP$glm.results)
  write.csv(all.pvals, paste0(title,".csv"), row.names=FALSE, quote=FALSE)

  #### QTL Plots
  # MVP.Report(imMVP$map, plot.type='d', col=c("darkgreen", "yellow", "red"), bin.size=1e3, file.type="jpg", dpi=300)

  # MVP.Report(imMVP, plot.type=c("m","q"), LOG10=TRUE, threshold=c(1e-3),col=c("dodgerblue4","deepskyblue"), cex=0.7,
  #              chr.den.col=NULL, file.type="jpg", dpi=300)
 }