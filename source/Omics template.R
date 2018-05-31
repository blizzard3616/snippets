#Omics circular plot template

#Find document online
#"The OmicCircos usages by examples by Ying Hu and Chunhua Yan October 30, 2017"
#Available source online :http://www.cancergenome.nih.gov/

# source("https://bioconductor.org/biocLite.R")
# biocLite("OmicCircos")

options(stringsAsFactors = FALSE)
library(OmicCircos)

setwd("C:/Users/ry147/Desktop/Skin UVB RNA-DNA seq folder 3-9-2018/Paper prepared figures/Methyl seq papers/SFN/")

dt1 <- read.table("combined_dmr_fraction_s32_uvb-skin_renyi_02092018.csv", header = TRUE, sep = "", skip = 0)

head(dt1)

colnames(dt1)

regment1 <- unique(subset(dt1, 
                      select = c("chr",
                                 "start",
                                 "end")))
