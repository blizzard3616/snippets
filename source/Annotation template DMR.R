#Template
#______________________________________________
#Project: Annotitaion by region DMR results

# Workflow: https://www.bioconductor.org/help/workflows/rnaseqGene/
# source("https://bioconductor.org/biocLite.R")
# biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
# biocLite("ChIPseeker")
# biocLite("DO.db")
# biocLite("GenomicRanges")
# biocLite("org.Mm.eg.db")
# biocLite("DSS")
# biocLite("bsseq")

require(data.table)
require(ggplot2)
require(knitr)
require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(DSS)

# Load and view raw counts (no annoation)----

# (version on Ran lab computer)
dt01 <- fread("../Desktop/Skin UVB RNA-DNA seq folder 3-9-2018/Paper prepared figures/Methyl seq papers/SFN/combined_dmr_default__uvb-skin_renyi_02092018.csv")
dt01

# NOTE: there are 14 rows representing mitochondrial DNA
unique(dt01$chr)
dt01[dt01$chr == "chrM",]

# Load and view percent methylated (no annotation)----
dt02 <- fread("../Desktop/Skin UVB RNA-DNA seq folder 3-9-2018/Paper prepared figures/Methyl seq papers/SFN/combined_dmr_fraction_s32_uvb-skin_renyi_02092018.csv")
dt02

# Annotate----
peakAnno1 <- annotatePeak(peak = "../Desktop/Skin UVB RNA-DNA seq folder 3-9-2018/Paper prepared figures/Methyl seq papers/SFN/combined_dmr_default__uvb-skin_renyi_02092018.csv", 
                          tssRegion = c(-3000, 3000), 
                          TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                          annoDb = "org.Mm.eg.db")

head(peakAnno1@detailGenomicAnnotation)
t1 <- peakAnno1@annoStat
t1$Feature <- factor(t1$Feature,
                     levels = as.character(t1$Feature[order(t1$Frequency,
                                                            decreasing = TRUE)]))
t1
p1 <- ggplot(t1,
             aes(x = "",
                 y = Frequency,
                 fill = Feature)) +
  geom_bar(width = 1, 
           stat = "identity",
           color = "black") +
  coord_polar("y",
              start = 0,
              direction = -1) +
  scale_x_discrete("") +
  ggtitle("Methy-seq DMRs annotaed by Region (%)")
p1 

tiff(filename = "tmp/skin_uvb_anno_by_reg.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()