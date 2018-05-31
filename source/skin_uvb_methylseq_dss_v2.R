# |----------------------------------------------------------------------------------|
# | Project: Skin UVB SKH1 mouse model treated with UA/SFN                           |
# | Script: Methyl-seq data analysis and visualization using DSS                     |
# | Coordinator: Ran Yin, Renyi Wu                                                   |
# | Author: Davit Sargsyan                                                           |
# | Created: 03/17/2018                                                              |
# | Modified:04/16/2018, DS: changed hitmaps to donut plots; added more comparisons  |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_skin_uvb_dna_v2.txt")
date()

# NOTE: several packages, e.g. Rcpp, MASS, etc., might be deleted manually and reinstalled
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

# Make data table----
dt1 <- data.table(start = peakAnno1@anno@ranges@start,
                  as.data.frame(peakAnno1@anno@elementMetadata@listData))

# Remove unmapped regions
dt1 <- dt1[!is.na(dt1$SYMBOL == "NA"), ]
# Removed 12 rows
dt1
# Subset data; remove Ursolic Acid samles----
dt1 <- data.table(gene = dt1$SYMBOL,
                  anno = dt1$annotation,
                  geneId = dt1$geneId,
                  chr = dt1$geneChr,
                  pos = dt1$start,
                  reg = NA,
                  dt1[, CpG:X02w_SFN_1.X],
                  dt1[, X02w_UVB_0.N:X15w_SFN_1.X], 
                  dt1[, X15w_UVB_0.N:X25t_SFN_1.X],
                  dt1[, X25t_UVB_0.N:X25w_SFN_1.X],
                  dt1[, X25w_UVB_0.N:X25w_UVB_1.X],
                  geneName = dt1$GENENAME)
dt1

# Regions----
kable(data.table(table(substr(dt1$anno, 1, 9))))
# |V1        |     N|
# |:---------|-----:|
# |3' UTR    |  4869|
# |5' UTR    |   839|
# |Distal In | 69042|
# |Downstrea |  2939|
# |Exon (uc0 | 12726|
# |Intron (u | 60617|
# |Promoter  | 86826|

# Separate Promoter, Body and Downstream----
dt1$reg <- as.character(dt1$anno)

# a. Promoter: up to 3kb upstream
dt1$reg[substr(dt1$anno, 
               1,
               8) == "Promoter"] <- "Promoter"

# b. Body: exons and introns
dt1$reg[substr(dt1$anno, 
               1, 
               4) %in% c("Exon",
                         "Intr")] <- "Body"

# c. Downstream: Distal Intergenic and  Downstream
dt1$reg[substr(dt1$anno, 
               1, 
               4) %in% c("Dist",
                         "Down")] <- "Downstream"

dt1$reg <- factor(dt1$reg,
                  levels = c("Promoter",
                             "5' UTR",
                             "Body",
                             "3' UTR",
                             "Downstream"))
kable(data.table(table(dt1$reg)))
# |V1         |     N|
# |:----------|-----:|
# |Promoter   | 86826|
# |5' UTR     |   839|
# |Body       | 73343|
# |3' UTR     |  4869|
# |Downstream | 71981|

# CpG distribution and coverage----
p2 <- ggplot(dt1,
             aes(x = CpG)) +
  facet_wrap(~ reg,
             scale = "free_y") +
  geom_histogram(color = "black",
                 fill = "grey",
                 binwidth = 5) +
  scale_x_continuous(name = "Number of CpG",
                     breaks = c(3, 
                                seq(from = 5,
                                    to = 60,
                                    by = 5))) +
  coord_cartesian(xlim=c(3, 60)) +
  scale_y_continuous(name = "Counts") +
  ggtitle("Distribution of DMR by Number of CpG and Region")
p2

tiff(filename = "tmp/skin_uvb_CpG_by_reg_hist.tiff",
     height = 8,
     width = 12,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Percent methylation----
tmp <- as.matrix(dt1[, X02w_CON_0.N:X25w_UVB_1.X])
head(tmp)

dtN <- tmp[, seq(1,
                 ncol(tmp) - 1, 
                 by = 2)]
head(dtN)

dtX <- tmp[, seq(2,
                 ncol(tmp), 
                 by = 2)]
head(dtX)

pct <- dtX/dtN
colnames(pct) <- substr(colnames(pct),
                        1,
                        nchar(colnames(pct)) - 2)
head(pct)

# Remove rows with all zeros----
dim(pct[rowSums(pct) == 0, ])
dim(pct[is.na(rowSums(pct)), ])
dim(pct)
# 53,894 (2,548 with repalced zeros) out of 237,858 such rows. Remove.

ndx.keep <- rowSums(pct) != 0 & !is.na(rowSums(pct))
pct <- pct[ndx.keep, ]
dt1 <- dt1[ndx.keep, ]
dtN <- dtN[ndx.keep, ]
dtX <- dtX[ndx.keep, ]
dim(dtX)
# 183,964 remaine

# Hits per CpG average (i.e. vertical coverage)----
t1 <- apply(dtN,
            2,
            function(a) {
              return(round(a/dt1$CpG,
                           1))
            })

mu <- list()
for (i in 1:(ncol(t1)/2)) {
  x1 <- aggregate(x = t1[, 2*i - 1],
                  FUN = mean,
                  by = list(dt1$reg))
  x2 <- aggregate(x = t1[, 2*i],
                  FUN = mean,
                  by = list(dt1$reg))
  x3 <- merge(x1, x2, by = "Group.1")
  mu[[i]] <- data.table(reg = x3[, 1],
                        mu = (x3[, 2] + x3[, 3])/2)
}
names(mu) <- unique(substr(colnames(t1),
                           1,
                           8))
mu

# Average methylation per region per treatment/time
mumth <- list()
for (i in 1:(ncol(pct)/2)) {
  x1 <- aggregate(x = c(pct[, 2*i - 1],
                        pct[, 2*i]),
                  FUN = mean,
                  by = list(rep(dt1$reg, 2)))
  
  x2 <- aggregate(x = c(pct[, 2*i - 1],
                        pct[, 2*i]),
                  FUN = sd,
                  by = list(rep(dt1$reg, 2)))
  mumth[[i]] <- data.table(rep(colnames(pct)[[2*i]],
                               5),
                           merge(x1, 
                                 x2,
                                 by = "Group.1"))
  
  colnames(mumth[[i]]) <- c("trt",
                            "reg",
                            "mu",
                            "std")
}
mumth <- do.call("rbind",
                 mumth)
mumth
mumth$Time <- substr(mumth$trt, 1, 4)
mumth$Time <- mumth$Time <- factor(mumth$Time,
                                   levels = c("X02w",
                                              "X15w",
                                              "X25w",
                                              "X25t"),
                                   labels = c("Week 2",
                                              "Week 15",
                                              "Week 25",
                                              "Week 25 Tumor"))

mumth$Treatment <- substr(mumth$trt, 6, 8)
mumth$Treatment <- factor(mumth$Treatment,
                          levels = c("CON",
                                     "UVB",
                                     "SFN"),
                          labels = c("Control",
                                     "UVB",
                                     "SFN + UVB"))

mumth$`Methylation (%)` <- 100*mumth$mu

p1 <- ggplot(mumth,
             aes(x = reg,
                 y = `Methylation (%)`,
                 group = Treatment,
                 fill = Treatment)) +
  facet_wrap(~ Time) +
  geom_bar(position = position_dodge(),
           stat="identity",
           color = "black") +
  scale_x_discrete("Region") +
  scale_y_continuous(limits = c(0, 60)) +
  ggtitle("Percent of Methylated CpG by Region") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p1
tiff(filename = "tmp/skin_uvb_avg_methyl_by_reg.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

p2 <- ggplot(mumth,
             aes(x = reg,
                 y = mu,
                 group = Treatment,
                 fill = Treatment)) +
  facet_wrap(~ Time) +
  geom_errorbar(aes(ymax = mu + std,
                    ymin = mu),
                width = 0.5,
                position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9),
             size = 1) +
  geom_bar(position = position_dodge(0.9),
           stat="identity",
           color = "black") +
  scale_x_discrete("Region") +
  scale_y_continuous() +
  ggtitle("Proportion of Methylated CpG by Region") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p2
tiff(filename = "tmp/skin_uvb_avg_sd_methyl_by_reg.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Sample names----
tumor <- rep("", ncol(pct))
tumor[substr(colnames(pct),
             4,
             4) == "t"] <- "t"
snames <- paste("Week",
                as.numeric(substr(colnames(pct),
                                  2,
                                  3)),
                tumor,
                "_",
                substr(colnames(pct),
                       6,
                       10),
                sep = "")
snames

# Clustering by Euclidean distances----
dst <- dist(t(pct))
sampleDists <- as.matrix(dst)
rownames(sampleDists) <- colnames(sampleDists) <- snames

# Make zeros less influential
diag(sampleDists) <- min(dst) - 5

tiff(filename = "tmp/skin_uvb_samples_cluster.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
heatmap(sampleDists,
        symm = TRUE,
        col = heat.colors(50),
        margins = c(10, 10))
graphics.off()

# PCA----
m1 <- prcomp(t(pct),
             center = TRUE,
             scale. = TRUE)
summary(m1)
plot(m1)

# Biplot while keep only the most important variables (Javier)----
# Select PC-s to pliot (PC1 & PC2)
choices <- 1:2

# Scores, i.e. points (df.u)
dt.scr <- data.table(m1$x[, choices])
# Add grouping variable
dt.scr$grp <- substr(colnames(pct),
                     6,
                     8)
dt.scr$grp <- factor(dt.scr$grp,
                     levels = c("CON",
                                "UVB",
                                "SFN"),
                     labels = c("Control",
                                "UVB",
                                "SFN + UVB"))

dt.scr$sample <- snames
dt.scr

# Loadings, i.e. arrows (df.v)
dt.rot <- as.data.frame(m1$rotation[, choices])
dt.rot$feat <- rownames(dt.rot)
dt.rot <- data.table(dt.rot)
dt.rot

# Axis labels
u.axis.labs <- paste(colnames(dt.rot)[1:2], 
                     sprintf('(%0.1f%% explained var.)', 
                             100*m1$sdev[choices]^2/sum(m1$sdev^2)))
u.axis.labs

p1 <- ggplot(data = dt.rot,
             aes(x = PC1,
                 y = PC2)) +
  geom_point(data = dt.scr,
             aes(fill = grp),
             shape = 21,
             size = 3,
             alpha = 0.5) +
  geom_text(data = dt.scr,
            aes(x = PC1 + 40,
                y = PC2,
                label = dt.scr$sample),
            size = 2,
            hjust = 0.5) +
  scale_x_continuous(u.axis.labs[1]) +
  scale_y_continuous(u.axis.labs[2]) +
  scale_fill_manual(name = "Treatment",
                    values = c("white", "red", "blue", "green")) +
  ggtitle("PCA of Percent Methylation") +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10))
p1
tiff(filename = "tmp/skin_uvb_pca_plot.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Clustering, no tumor----
ndx.keep <- which(!(colnames(pct) %in% c("X25t_CON_0",
                                         "X25t_CON_1",
                                         "X25t_SFN_0",
                                         "X25t_SFN_1",
                                         "X25t_UVB_0",
                                         "X25t_UVB_1")))

dst <- dist(t(pct[, ndx.keep]))
sampleDists <- as.matrix(dst)
rownames(sampleDists) <- colnames(sampleDists) <- snames[-c(13:18)]
sampleDists

# Make zeros less influential
diag(sampleDists) <- min(dst) - 5

tiff(filename = "tmp/skin_uvb_samples_cluster_no_tumor.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
heatmap(sampleDists,
        symm = TRUE,
        col = heat.colors(50),
        margins = c(10, 10))
graphics.off()

# PCA, no tumor----
xx <- apply(t(pct[, ndx.keep]), 
            2,
            scale)
head(t(xx))

row.keep <- apply(t(xx),
                  1,
                  function(a) {
                    sum(is.na(a)) == 0
                  })

m1 <- prcomp(t(pct[row.keep, ndx.keep]),
             center = TRUE,
             scale. = TRUE)
summary(m1)
plot(m1)

# Biplot while keep only the most important variables (Javier)----
# Select PC-s to pliot (PC1 & PC2)
choices <- 1:2

# Scores, i.e. points (df.u)
dt.scr <- data.table(m1$x[, choices])
# Add grouping variable
dt.scr$grp <- substr(colnames(pct[row.keep, ndx.keep]),
                     6,
                     8)
dt.scr$grp <- factor(dt.scr$grp,
                     levels = c("CON",
                                "UVB",
                                "SFN"),
                     labels = c("Control",
                                "UVB",
                                "SFN + UVB"))

dt.scr$sample <- snames[-c(13:18)]
dt.scr

# Loadings, i.e. arrows (df.v)
dt.rot <- as.data.frame(m1$rotation[, choices])
dt.rot$feat <- rownames(dt.rot)
dt.rot <- data.table(dt.rot)
dt.rot

# Axis labels
u.axis.labs <- paste(colnames(dt.rot)[1:2], 
                     sprintf('(%0.1f%% explained var.)', 
                             100*m1$sdev[choices]^2/sum(m1$sdev^2)))
u.axis.labs

p1 <- ggplot(data = dt.rot,
             aes(x = PC1,
                 y = PC2)) +
  geom_point(data = dt.scr,
             aes(fill = grp),
             shape = 21,
             size = 3,
             alpha = 0.5) +
  geom_text(data = dt.scr,
            aes(x = PC1 + 40,
                y = PC2,
                label = dt.scr$sample),
            size = 2,
            hjust = 0.5) +
  scale_x_continuous(u.axis.labs[1]) +
  scale_y_continuous(u.axis.labs[2]) +
  scale_fill_manual(name = "Treatment",
                    values = c("white", "red", "blue", "green")) +
  ggtitle("PCA of Percent Methylation") +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10))
p1
tiff(filename = "tmp/skin_uvb_pca_plot_no_tumor.tiff",
     height = 6,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Dispersion Shrinkage for Sequencing data (DSS)----
# This is based on Wald test for beta-binomial distribution.
# Source: https://www.bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.pdf
# The DM detection procedure implemented in DSS is based on a rigorous Wald test for betabinomial
# distributions. The test statistics depend on the biological variations (characterized
# by dispersion parameter) as well as the sequencing depth. An important part of the algorithm
# is the estimation of dispersion parameter, which is achieved through a shrinkage estimator
# based on a Bayesian hierarchical model [1]. An advantage of DSS is that the test can be
# performed even when there is no biological replicates. That’s because by smoothing, the
# neighboring CpG sites can be viewed as “pseudo-replicates", and the dispersion can still be
# estimated with reasonable precision.

# DMR Part I: Controls only (aging effect)----
snames <- c("W2Ctrl1",
            "W2Ctrl2",
            "W15Ctrl1",
            "W15Ctrl2",
            "W25Ctrl1",
            "W25Ctrl2")

# Prepare data----
dtl <- list(data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_CON_0.N,
                       X = dt1$X02w_CON_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_CON_1.N,
                       X = dt1$X02w_CON_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_CON_0.N,
                       X = dt1$X15w_CON_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_CON_1.N,
                       X = dt1$X15w_CON_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_CON_0.N,
                       X = dt1$X25w_CON_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_CON_1.N,
                       X = dt1$X25w_CON_1.X))
dtl

BSobj <- makeBSseqData(dat = dtl,
                       sampleNames = snames)
BSobj

design <- data.table(time = rep(rep(c("Week2", 
                                      "Week15", 
                                      "Week25"),
                                    each = 2),
                                1))
design$time <- factor(design$time,
                      levels = c("Week2", 
                                 "Week15", 
                                 "Week25"))
design

DMLfit <- DMLfit.multiFactor(BSobj = BSobj, 
                             design = design,
                             formula = ~ time)
DMLfit$X

# a. Control Week 15 vs. Week 2----
DMLtest.Ctrl.W15.W2 <- DMLtest.multiFactor(DMLfit,
                                           coef = "timeWeek15")

DMLtest.Ctrl.W15.W2$chr <- as.numeric(as.character(DMLtest.Ctrl.W15.W2$chr))
head(pct[, c("X15w_CON_0",
             "X15w_CON_1",
             "X02w_CON_0",
             "X02w_CON_1")])
DMLtest.Ctrl.W15.W2 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                   pct[, c("X15w_CON_0",
                                                           "X15w_CON_1",
                                                           "X02w_CON_0",
                                                           "X02w_CON_1")]),
                                        DMLtest.Ctrl.W15.W2,
                                        by = c("chr",
                                               "pos")))
DMLtest.Ctrl.W15.W2$mu.ctrl15w <- (DMLtest.Ctrl.W15.W2$X15w_CON_0 + DMLtest.Ctrl.W15.W2$X15w_CON_1)/2
DMLtest.Ctrl.W15.W2$mu.ctrl02w <- (DMLtest.Ctrl.W15.W2$X02w_CON_0 + DMLtest.Ctrl.W15.W2$X02w_CON_1)/2

# Mean vs. difference----
DMLtest.Ctrl.W15.W2$mu <- (DMLtest.Ctrl.W15.W2$mu.ctrl15w + DMLtest.Ctrl.W15.W2$mu.ctrl02w)/2
DMLtest.Ctrl.W15.W2$diff <- (DMLtest.Ctrl.W15.W2$mu.ctrl15w - DMLtest.Ctrl.W15.W2$mu.ctrl02w)

# Green: hypomethylated at Week 2; Red: hypermethylated at Week 2
dtp1 <- DMLtest.Ctrl.W15.W2[!is.na(stat), ]
tiff(filename = "tmp/skin_uvb_maplot_ctrl_w15-w02.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dtp1$diff ~ dtp1$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "Proportion of Methylated Reads in Control\nat Week 15 vs. Week 2, FDR < 0.1")
points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff > 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff > 0] ,
       pch = "x",
       col = "green")
points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff < 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff < 0] ,
       pch = "x",
       col = "red")
abline(h = c(-0.2, 0.2),
       lty = 2)
graphics.off()

# b. Control Week 25 vs. Week 2----
DMLtest.Ctrl.W25.W2 <- DMLtest.multiFactor(DMLfit,
                                           coef = "timeWeek25")

DMLtest.Ctrl.W25.W2$chr <- as.numeric(as.character(DMLtest.Ctrl.W25.W2$chr))
head(pct[, c("X25w_CON_0",
             "X25w_CON_1",
             "X02w_CON_0",
             "X02w_CON_1")])
DMLtest.Ctrl.W25.W2 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                   pct[, c("X25w_CON_0",
                                                           "X25w_CON_1",
                                                           "X02w_CON_0",
                                                           "X02w_CON_1")]),
                                        DMLtest.Ctrl.W25.W2,
                                        by = c("chr",
                                               "pos")))
DMLtest.Ctrl.W25.W2$mu.ctrl25w <- (DMLtest.Ctrl.W25.W2$X25w_CON_0 + DMLtest.Ctrl.W25.W2$X25w_CON_1)/2
DMLtest.Ctrl.W25.W2$mu.ctrl02w <- (DMLtest.Ctrl.W25.W2$X02w_CON_0 + DMLtest.Ctrl.W25.W2$X02w_CON_1)/2

# Mean vs. difference----
DMLtest.Ctrl.W25.W2$mu <- (DMLtest.Ctrl.W25.W2$mu.ctrl25w + DMLtest.Ctrl.W25.W2$mu.ctrl02w)/2
DMLtest.Ctrl.W25.W2$diff <- (DMLtest.Ctrl.W25.W2$mu.ctrl25w - DMLtest.Ctrl.W25.W2$mu.ctrl02w)

# Green: hypomethylated at Week 2; Red: hypermethylated at Week 2
dtp1 <- DMLtest.Ctrl.W25.W2[!is.na(stat), ]
tiff(filename = "tmp/skin_uvb_maplot_ctrl_w25-w02.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dtp1$diff ~ dtp1$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "Proportion of Methylated Reads in Control\nat Week 25 vs. Week 2, FDR < 0.1")
points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff > 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff > 0] ,
       pch = "x",
       col = "green")
points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff < 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff < 0] ,
       pch = "x",
       col = "red")
abline(h = c(-0.2, 0.2),
       lty = 2)
graphics.off()

# c. Control Week 25 vs. Week 15----
# NOTE: checked teh contrast by rerunning teh analysis wiht Week 15 as reference
# Identical estimates: proceed with the code below!
DMLtest.Ctrl.W25.W15 <- DMLtest.multiFactor(DMLfit,
                                            Contrast = matrix(c(0, -1, 1), 
                                                              ncol = 1))

DMLtest.Ctrl.W25.W15$chr <- as.numeric(as.character(DMLtest.Ctrl.W25.W15$chr))
head(pct[, c("X25w_CON_0",
             "X25w_CON_1",
             "X15w_CON_0",
             "X15w_CON_1")])
DMLtest.Ctrl.W25.W15 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                    pct[, c("X25w_CON_0",
                                                            "X25w_CON_1",
                                                            "X15w_CON_0",
                                                            "X15w_CON_1")]),
                                         DMLtest.Ctrl.W25.W15,
                                         by = c("chr",
                                                "pos")))
DMLtest.Ctrl.W25.W15$mu.ctrl25w <- (DMLtest.Ctrl.W25.W15$X25w_CON_0 + DMLtest.Ctrl.W25.W15$X25w_CON_1)/2
DMLtest.Ctrl.W25.W15$mu.ctrl15w <- (DMLtest.Ctrl.W25.W15$X15w_CON_0 + DMLtest.Ctrl.W25.W15$X15w_CON_1)/2

# Mean vs. difference----
DMLtest.Ctrl.W25.W15$mu <- (DMLtest.Ctrl.W25.W15$mu.ctrl25w + DMLtest.Ctrl.W25.W15$mu.ctrl15w)/2
DMLtest.Ctrl.W25.W15$diff <- (DMLtest.Ctrl.W25.W15$mu.ctrl25w - DMLtest.Ctrl.W25.W15$mu.ctrl15w)

# Green: hypomethylated at Week 15; Red: hypermethylated at Week 15
dtp1 <- DMLtest.Ctrl.W25.W15[!is.na(stat), ]
tiff(filename = "tmp/skin_uvb_maplot_ctrl_w25-w15.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dtp1$diff ~ dtp1$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "Proportion of Methylated Reads in Control\nat Week 25 vs. Week 15, FDR < 0.1")
points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff > 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff > 0] ,
       pch = "x",
       col = "green")
points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff < 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff < 0] ,
       pch = "x",
       col = "red")
abline(h = c(-0.2, 0.2),
       lty = 2)
graphics.off()

# Heatmaps 1----
sign.15.02 <- c(unique(DMLtest.Ctrl.W15.W2[fdrs < 0.01 & diff > 0.2]$pos),
                unique(DMLtest.Ctrl.W15.W2[fdrs < 0.01 & diff < -0.2]$pos))
l1 <- DMLtest.Ctrl.W15.W2[pos %in% sign.15.02,
                          c("gene", 
                            "pos", 
                            "diff")]

sign.25.02 <- c(unique(DMLtest.Ctrl.W25.W2[fdrs < 0.01 & diff > 0.2]$pos),
                unique(DMLtest.Ctrl.W25.W2[fdrs < 0.01 & diff < -0.2]$pos))
l2 <- DMLtest.Ctrl.W25.W2[pos %in% sign.25.02,
                          c("gene", 
                            "pos", 
                            "diff")]

sign.25.15 <- unique(DMLtest.Ctrl.W25.W15[fdrs < 0.01 & diff > 0]$pos)
l3 <- DMLtest.Ctrl.W25.W15[, c("gene", 
                               "pos", 
                               "diff")]

ll <- l1[l1$pos %in% l2$pos, ]
# ll <- ll[ll$pos %in% l3$pos, ]
ll

t1 <- merge(l1[l1$pos %in% ll$pos],
            l2[l2$pos %in% ll$pos],
            by = c("gene", 
                   "pos"))
t1 <- merge(t1[t1$pos %in% ll$pos],
            l3[l3$pos %in% ll$pos],
            by = c("gene", 
                   "pos"))
write.csv(t1,
          file = "tmp/heatmap1.csv")
ll <- t1

ll <- melt.data.table(data = ll,
                      id.vars = 1:2,
                      measure.vars = 3:5,
                      variable.name = "Comparison",
                      value.name = "Methyl Diff")
ll$reg <- factor(paste(ll$gene, 
                       ll$pos,
                       sep = "_"))
ll$Comparison <- factor(ll$Comparison,
                        levels = rev(c("diff.x",
                                       "diff.y",
                                       "diff")),
                        labels = rev(c("w15w2",
                                       "w25w2",
                                       "w25w15")))
lvls <- ll[ll$Comparison == "w15w2", ]
ll$reg <- factor(ll$reg,
                 levels = lvls$reg[order(lvls$`Methyl Diff`)])

p1 <- ggplot(data = ll) +
  coord_polar("y",
              start = 0,
              direction = -1) +
  geom_tile(aes(x =  as.numeric(Comparison),
                y = reg, 
                fill = `Methyl Diff`),
            color = "white") +
  geom_text(data = ll[Comparison == "w25w15", ],
            aes(x = rep(2.75,
                        nlevels(reg)),
                y = reg,
                label = unique(reg),
                angle = 90 + seq(from = 0,
                            to = 360,
                            length.out = nlevels(reg))[as.numeric(reg)]),
            hjust = 0) +
  scale_fill_gradient2(low = "red", 
                       high = "green", 
                       mid = "grey", 
                       midpoint = 0, 
                       name = "Methyl Diff") +
  scale_x_continuous(limits = c(0, 
                                max(as.numeric(ll$Comparison)) + 0.5),
                     expand = c(0, 0)) + 
  scale_y_discrete("",
                   expand = c(0, 0)) +
  ggtitle("Changes in Controls Over Time") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
# p1

tiff(filename = "tmp/skin_uvb_heatmap_ctrl.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# # Heatmaps 2----
# sign.15.02 <- c(unique(DMLtest.Ctrl.W15.W2[fdrs < 0.01 & diff > 0.2]$pos),
#                 unique(DMLtest.Ctrl.W15.W2[fdrs < 0.01 & diff < -0.2]$pos))
# l1 <- DMLtest.Ctrl.W15.W2[pos %in% sign.15.02,
#                           c("gene", 
#                             "pos", 
#                             "diff")]
# names(l1)[3] <- "w15w2"
# 
# l2 <- DMLtest.Ctrl.W25.W2[pos %in% sign.15.02,
#                           c("gene", 
#                             "pos", 
#                             "diff")]
# names(l2)[3] <- "w25w2"
# 
# l3 <- DMLtest.Ctrl.W25.W15[pos %in% sign.15.02,
#                            c("gene", 
#                              "pos", 
#                              "diff")]
# names(l3)[3] <- "w25w15"
# 
# ll <- merge(l1, l2, by = c("gene", 
#                            "pos"))
# ll <- merge(ll, l3, by = c("gene", 
#                            "pos"))
# ll
# 
# ll <- melt.data.table(data = ll,
#                       id.vars = 1:2,
#                       measure.vars = 3:5,
#                       variable.name = "Comparison",
#                       value.name = "Methyl Diff")
# ll$reg <- factor(paste(ll$gene, 
#                        ll$pos,
#                        sep = "_"))
# ll$Comparison <- factor(ll$Comparison,
#                         levels = unique(ll$Comparison))
# 
# lvls <- ll[ll$Comparison == "w15w2", ]
# ll$reg <- factor(ll$reg,
#                  levels = lvls$reg[order(lvls$`Methyl Diff`)])
# ll
# 
# p1 <- ggplot(data = ll) +
#   geom_tile(aes(x =  Comparison,
#                 y = reg, 
#                 fill = `Methyl Diff`),
#             color = "black") +
#   scale_fill_gradient2(low = "red", 
#                        high = "green", 
#                        mid = "black", 
#                        midpoint = 0, 
#                        # limit = c(-10, 10), 
#                        name = "Methyl Diff") +
#   scale_x_discrete(expand = c(0, 0)) + 
#   scale_y_discrete("Gene Region",
#                    expand = c(0, 0)) +
#   ggtitle("Changes in Controls Over Time") +
#   theme(plot.title = element_text(hjust = 0.5))
# p1
# 
# tiff(filename = "tmp/skin_uvb_heatmap_ctrl.tiff",
#      height = 20,
#      width = 6,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# print(p1)
# graphics.off()

# DMR Part II: Week2 Comparisons----
snames <- c("W2UVB1",
            "W2UVB2",
            "W2Ctrl1",
            "W2Ctrl2",
            "W2SFN1",
            "W2SFN2")

# Multi-factor analysis (treatment + time + interaction)----
dtl <- list(data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_UVB_0.N,
                       X = dt1$X02w_UVB_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_UVB_1.N,
                       X = dt1$X02w_UVB_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_CON_0.N,
                       X = dt1$X02w_CON_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_CON_1.N,
                       X = dt1$X02w_CON_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_SFN_0.N,
                       X = dt1$X02w_SFN_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X02w_SFN_1.N,
                       X = dt1$X02w_SFN_1.X))
dtl

BSobj <- makeBSseqData(dat = dtl,
                       sampleNames = snames)
BSobj

design <- data.table(trt = rep(rep(c("UVB", 
                                     "Ctrl", 
                                     "SFN"),
                                   each = 2),
                               1))
design$trt <- factor(design$trt,
                     levels = unique(design$trt))
design

DMLfit <- DMLfit.multiFactor(BSobj = BSobj, 
                             design = design,
                             formula = ~ trt)

# a. Control vs. UVB at Week 2----
DMLtest.Ctrl.UVB.W2 <- DMLtest.multiFactor(DMLfit,
                                           coef = "trtCtrl")

DMLtest.Ctrl.UVB.W2$chr <- as.numeric(as.character(DMLtest.Ctrl.UVB.W2$chr))
head(pct[, c("X02w_CON_0",
             "X02w_CON_1",
             "X02w_UVB_0",
             "X02w_UVB_1")])
DMLtest.Ctrl.UVB.W2 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                   pct[, c("X02w_CON_0",
                                                           "X02w_CON_1",
                                                           "X02w_UVB_0",
                                                           "X02w_UVB_1")]),
                                        DMLtest.Ctrl.UVB.W2,
                                        by = c("chr",
                                               "pos")))
DMLtest.Ctrl.UVB.W2$mu.ctrl2w <- (DMLtest.Ctrl.UVB.W2$X02w_CON_0 + DMLtest.Ctrl.UVB.W2$X02w_CON_1)/2
DMLtest.Ctrl.UVB.W2$mu.uvb2w <- (DMLtest.Ctrl.UVB.W2$X02w_UVB_0 + DMLtest.Ctrl.UVB.W2$X02w_UVB_1)/2

# Mean vs. difference----
DMLtest.Ctrl.UVB.W2$mu <- (DMLtest.Ctrl.UVB.W2$mu.ctrl2w + DMLtest.Ctrl.UVB.W2$mu.uvb2w)/2
DMLtest.Ctrl.UVB.W2$diff <- (DMLtest.Ctrl.UVB.W2$mu.ctrl2w - DMLtest.Ctrl.UVB.W2$mu.uvb2w)

# Green: hypomethylated in UVB; Red: hypermethylated in UVB
# NOTE: FLIP THE FIGURE!
dtp1 <- DMLtest.Ctrl.UVB.W2[!is.na(stat), ]
tiff(filename = "tmp/skin_uvb_maplot_uvb-ctrl_w2.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot((-1)*dtp1$diff ~ dtp1$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "Proportion of Methylation in UVB vs. Control\nat Week 2, FDR < 0.1")
points((-1)*dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff > 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff > 0] ,
       pch = "x",
       col = "red")
points((-1)*dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff < 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff < 0] ,
       pch = "x",
       col = "green")
abline(h = c(-0.2, 0.2),
       lty = 2)
graphics.off()

# b. SFN vs. UVB at Week 2----
DMLtest.SFN.UVB.W2 <- DMLtest.multiFactor(DMLfit,
                                          coef = "trtSFN")
head(DMLtest.SFN.UVB.W2)
DMLtest.SFN.UVB.W2$chr <- as.numeric(as.character(DMLtest.SFN.UVB.W2$chr))

DMLtest.SFN.UVB.W2 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                  pct[, c("X02w_SFN_0",
                                                          "X02w_SFN_1",
                                                          "X02w_UVB_0",
                                                          "X02w_UVB_1")]),
                                       DMLtest.SFN.UVB.W2,
                                       by = c("chr",
                                              "pos")))
DMLtest.SFN.UVB.W2$mu.sfn2w <- (DMLtest.SFN.UVB.W2$X02w_SFN_0 + DMLtest.SFN.UVB.W2$X02w_SFN_1)/2
DMLtest.SFN.UVB.W2$mu.uvb2w <- (DMLtest.SFN.UVB.W2$X02w_UVB_0 + DMLtest.SFN.UVB.W2$X02w_UVB_1)/2

# Mean vs. difference----
DMLtest.SFN.UVB.W2$mu <- (DMLtest.SFN.UVB.W2$mu.sfn2w + DMLtest.SFN.UVB.W2$mu.uvb2w)/2
DMLtest.SFN.UVB.W2$diff <- (DMLtest.SFN.UVB.W2$mu.sfn2w - DMLtest.SFN.UVB.W2$mu.uvb2w)

# Green: hypomethylated in UVB; Red: hypermethylated in UVB
dtp1 <- DMLtest.SFN.UVB.W2[!is.na(stat), ]
tiff(filename = "tmp/skin_uvb_maplot_sfn_uvb_w2.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dtp1$diff ~ dtp1$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "Proportion of Methylation in SFN vs. UVB\nat Weeks 2, FDR < 0.1")
points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff > 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff > 0] ,
       pch = "x",
       col = "red")
points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff < 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff < 0] ,
       pch = "x",
       col = "green")
abline(h = c(-0.2, 0.2),
       lty = 2)
graphics.off()

# Heatmaps----
uvb.hyper.1 <- unique(DMLtest.Ctrl.UVB.W2[fdrs < 0.1 & diff < 0]$pos)
uvb.hyper.2 <- unique(DMLtest.SFN.UVB.W2[fdrs < 0.1 & dtp1$diff < 0]$pos)
hyper <- uvb.hyper.1[uvb.hyper.1 %in% uvb.hyper.2]
l1 <- DMLtest.Ctrl.UVB.W2[pos %in% hyper, c("gene", 
                                            "pos", 
                                            "diff")]
l1$diff <- (-1)*l1$diff
l2 <- DMLtest.SFN.UVB.W2[pos %in% hyper, c("gene", 
                                           "pos", 
                                           "diff")]
# Venn diagram----
length(uvb.hyper.1)
length(uvb.hyper.2)
length(hyper)

uvb.hypo.1 <- unique(DMLtest.Ctrl.UVB.W2[fdrs < 0.1 & diff > 0]$pos)
uvb.hypo.2 <- unique(DMLtest.SFN.UVB.W2[fdrs < 0.1 & dtp1$diff > 0]$pos)
hypo <- uvb.hypo.1[uvb.hypo.1 %in% uvb.hypo.2]
l3 <- DMLtest.Ctrl.UVB.W2[pos %in% hypo, c("gene", 
                                           "pos", 
                                           "diff")]
l3$diff <- (-1)*l3$diff
l4 <- DMLtest.SFN.UVB.W2[pos %in% hypo, c("gene", 
                                          "pos", 
                                          "diff")]
# Venn diagram----
length(uvb.hypo.1)
length(uvb.hypo.2)
length(hypo)

ll <- merge(l1,
            l2,
            by = c("gene", 
                   "pos"))
ll2 <- merge(l3,
             l4,
             by = c("gene", 
                    "pos"))
ll <- rbind.data.frame(ll,
                       ll2)

ll$reg <- paste(ll$gene,
                ll$pos,
                sep = "_")
write.csv(ll,
          file = "tmp/w2.csv")

ll <- melt.data.table(data = ll,
                      id.vars = c(1, 5),
                      measure.vars = 3:4,
                      variable.name = "Comparison",
                      value.name = "Methyl Diff")
ll$reg <- factor(ll$reg)
ll$Comparison <- factor(ll$Comparison,
                        levels = rev(c("diff.x",
                                       "diff.y")),
                        labels = rev(c("UVB - Ctrl",
                                       "SFN - UVB")))

lvls <- ll[ll$Comparison == "UVB - Ctrl", ]
ll$reg <- factor(ll$reg,
                 levels = lvls$reg[order(lvls$`Methyl Diff`)])
ll

p1 <- ggplot(data = ll) +
  coord_polar("y",
              start = 0,
              direction = -1) +
  geom_tile(aes(x =  as.numeric(Comparison),
                y = reg, 
                fill = `Methyl Diff`),
            color = "white") +
  geom_text(data = ll[Comparison == "UVB - Ctrl", ],
            aes(x = rep(1.75,
                        nlevels(reg)),
                y = reg,
                label = unique(reg),
                angle = 90 + seq(from = 0,
                                 to = 360,
                                 length.out = nlevels(reg))[as.numeric(reg)]),
            hjust = 0,
            color = "black") +
  scale_fill_gradient2(low = "red", 
                       high = "green", 
                       mid = "grey", 
                       midpoint = 0, 
                       name = "Methyl Diff") +
  scale_x_continuous(limits = c(0, 
                                max(as.numeric(ll$Comparison)) + 0.5),
                     expand = c(0, 0)) + 
  scale_y_discrete("",
                   expand = c(0, 0)) +
  ggtitle("Changes in Controls Over Time") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

tiff(filename = "tmp/skin_uvb_heatmap_w2.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# DMR Part III: Week15 Comparisons----
snames <- c("W15UVB1",
            "W15UVB2",
            "W15Ctrl1",
            "W15Ctrl2",
            "W15SFN1",
            "W15SFN2")

# Multi-factor analysis (treatment + time + interaction)----
dtl <- list(data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_UVB_0.N,
                       X = dt1$X15w_UVB_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_UVB_1.N,
                       X = dt1$X15w_UVB_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_CON_0.N,
                       X = dt1$X15w_CON_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_CON_1.N,
                       X = dt1$X15w_CON_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_SFN_0.N,
                       X = dt1$X15w_SFN_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X15w_SFN_1.N,
                       X = dt1$X15w_SFN_1.X))
dtl

BSobj <- makeBSseqData(dat = dtl,
                       sampleNames = snames)
BSobj

design <- data.table(trt = rep(rep(c("UVB", 
                                     "Ctrl", 
                                     "SFN"),
                                   each = 2),
                               1))
design$trt <- factor(design$trt,
                     levels = unique(design$trt))
design

DMLfit <- DMLfit.multiFactor(BSobj = BSobj, 
                             design = design,
                             formula = ~ trt)

# a. Control vs. UVB at Week 15----
DMLtest.Ctrl.UVB.W15 <- DMLtest.multiFactor(DMLfit,
                                            coef = "trtCtrl")

DMLtest.Ctrl.UVB.W15$chr <- as.numeric(as.character(DMLtest.Ctrl.UVB.W15$chr))
head(pct[, c("X15w_CON_0",
             "X15w_CON_1",
             "X15w_UVB_0",
             "X15w_UVB_1")])
DMLtest.Ctrl.UVB.W15 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                    pct[, c("X15w_CON_0",
                                                            "X15w_CON_1",
                                                            "X15w_UVB_0",
                                                            "X15w_UVB_1")]),
                                         DMLtest.Ctrl.UVB.W15,
                                         by = c("chr",
                                                "pos")))
DMLtest.Ctrl.UVB.W15$mu.ctrl15w <- (DMLtest.Ctrl.UVB.W15$X15w_CON_0 + DMLtest.Ctrl.UVB.W15$X15w_CON_1)/2
DMLtest.Ctrl.UVB.W15$mu.uvb15w <- (DMLtest.Ctrl.UVB.W15$X15w_UVB_0 + DMLtest.Ctrl.UVB.W15$X15w_UVB_1)/2

# Mean vs. difference----
DMLtest.Ctrl.UVB.W15$mu <- (DMLtest.Ctrl.UVB.W15$mu.ctrl15w + DMLtest.Ctrl.UVB.W15$mu.uvb15w)/2
DMLtest.Ctrl.UVB.W15$diff <- (DMLtest.Ctrl.UVB.W15$mu.ctrl15w - DMLtest.Ctrl.UVB.W15$mu.uvb15w)

# Green: hypomethylated in UVB; Red: hypermethylated in UVB
# NOTE: FLIP THE FIGURE!
dtp1 <- DMLtest.Ctrl.UVB.W15[!is.na(stat), ]
tiff(filename = "tmp/skin_uvb_maplot_uvb-ctrl_w15.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot((-1)*dtp1$diff ~ dtp1$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "Proportion of Methylation in UVB vs. Control\nat Week 15, FDR < 0.1")
points((-1)*dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff > 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff > 0] ,
       pch = "x",
       col = "red")
points((-1)*dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff < 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff < 0] ,
       pch = "x",
       col = "green")
abline(h = c(-0.2, 0.2),
       lty = 2)
graphics.off()

# b. SFN vs. UVB at Week 15----
DMLtest.SFN.UVB.W15 <- DMLtest.multiFactor(DMLfit,
                                           coef = "trtSFN")
head(DMLtest.SFN.UVB.W15)
DMLtest.SFN.UVB.W15$chr <- as.numeric(as.character(DMLtest.SFN.UVB.W15$chr))

DMLtest.SFN.UVB.W15 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                   pct[, c("X15w_SFN_0",
                                                           "X15w_SFN_1",
                                                           "X15w_UVB_0",
                                                           "X15w_UVB_1")]),
                                        DMLtest.SFN.UVB.W15,
                                        by = c("chr",
                                               "pos")))
DMLtest.SFN.UVB.W15$mu.sfn15w <- (DMLtest.SFN.UVB.W15$X15w_SFN_0 + DMLtest.SFN.UVB.W15$X15w_SFN_1)/2
DMLtest.SFN.UVB.W15$mu.uvb15w <- (DMLtest.SFN.UVB.W15$X15w_UVB_0 + DMLtest.SFN.UVB.W15$X15w_UVB_1)/2

# Mean vs. difference----
DMLtest.SFN.UVB.W15$mu <- (DMLtest.SFN.UVB.W15$mu.sfn15w + DMLtest.SFN.UVB.W15$mu.uvb15w)/2
DMLtest.SFN.UVB.W15$diff <- (DMLtest.SFN.UVB.W15$mu.sfn15w - DMLtest.SFN.UVB.W15$mu.uvb15w)

# Green: hypomethylated in UVB; Red: hypermethylated in UVB
dtp1 <- DMLtest.SFN.UVB.W15[!is.na(stat), ]
tiff(filename = "tmp/skin_uvb_maplot_sfn_uvb_w15.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dtp1$diff ~ dtp1$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "Proportion of Methylation in SFN vs. UVB\nat Weeks 15, FDR < 0.1")
points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff > 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff > 0] ,
       pch = "x",
       col = "red")
points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff < 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff < 0] ,
       pch = "x",
       col = "green")
abline(h = c(-0.2, 0.2),
       lty = 2)
graphics.off()

# Heatmaps----
uvb.hyper.1 <- unique(DMLtest.Ctrl.UVB.W15[fdrs < 0.1 & diff < 0]$pos)
uvb.hyper.2 <- unique(DMLtest.SFN.UVB.W15[fdrs < 0.1 & dtp1$diff < 0]$pos)
hyper <- uvb.hyper.1[uvb.hyper.1 %in% uvb.hyper.2]
l1 <- DMLtest.Ctrl.UVB.W15[pos %in% hyper, c("gene", 
                                             "pos", 
                                             "diff")]
l1$diff <- (-1)*l1$diff
l2 <- DMLtest.SFN.UVB.W15[pos %in% hyper, c("gene", 
                                            "pos", 
                                            "diff")]
# Venn diagram----
length(uvb.hyper.1)
length(uvb.hyper.2)
length(hyper)

uvb.hypo.1 <- unique(DMLtest.Ctrl.UVB.W15[fdrs < 0.1 & diff > 0]$pos)
uvb.hypo.2 <- unique(DMLtest.SFN.UVB.W15[fdrs < 0.1 & dtp1$diff > 0]$pos)
hypo <- uvb.hypo.1[uvb.hypo.1 %in% uvb.hypo.2]
l3 <- DMLtest.Ctrl.UVB.W15[pos %in% hypo, c("gene", 
                                            "pos", 
                                            "diff")]
l3$diff <- (-1)*l3$diff
l4 <- DMLtest.SFN.UVB.W15[pos %in% hypo, c("gene", 
                                           "pos", 
                                           "diff")]
# Venn diagram----
length(uvb.hypo.1)
length(uvb.hypo.2)
length(hypo)

ll <- merge(l1,
            l2,
            by = c("gene", 
                   "pos"))
ll2 <- merge(l3,
             l4,
             by = c("gene", 
                    "pos"))
ll <- rbind.data.frame(ll,
                       ll2)

ll$reg <- paste(ll$gene,
                ll$pos,
                sep = "_")
write.csv(ll,
          file = "tmp/w15.csv")

ll <- melt.data.table(data = ll,
                      id.vars = c(1, 5),
                      measure.vars = 3:4,
                      variable.name = "Comparison",
                      value.name = "Methyl Diff")
ll$reg <- factor(ll$reg)
ll$Comparison <- factor(ll$Comparison,
                        levels = rev(c("diff.x",
                                       "diff.y")),
                        labels = rev(c("UVB - Ctrl",
                                       "SFN - UVB")))

lvls <- ll[ll$Comparison == "UVB - Ctrl", ]
ll$reg <- factor(ll$reg,
                 levels = lvls$reg[order(lvls$`Methyl Diff`)])
ll

p1 <- ggplot(data = ll) +
  coord_polar("y",
              start = 0,
              direction = -1) +
  geom_tile(aes(x =  as.numeric(Comparison),
                y = reg, 
                fill = `Methyl Diff`),
            color = "white") +
  geom_text(data = ll[Comparison == "UVB - Ctrl", ],
            aes(x = rep(1.75,
                        nlevels(reg)),
                y = reg,
                label = unique(reg),
                angle = 90 + seq(from = 0,
                                 to = 360,
                                 length.out = nlevels(reg))[as.numeric(reg)]),
            hjust = 0,
            color = "black") +
  scale_fill_gradient2(low = "red", 
                       high = "green", 
                       mid = "grey", 
                       midpoint = 0, 
                       name = "Methyl Diff") +
  scale_x_continuous(limits = c(0, 
                                max(as.numeric(ll$Comparison)) + 0.5),
                     expand = c(0, 0)) + 
  scale_y_discrete("",
                   expand = c(0, 0)) +
  ggtitle("Changes in Controls Over Time") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

tiff(filename = "tmp/skin_uvb_heatmap_w15.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# DMR Part IV: Week25 Comparisons----
snames <- c("W25UVB1",
            "W25UVB2",
            "W25Ctrl1",
            "W25Ctrl2",
            "W25SFN1",
            "W25SFN2")

# Multi-factor analysis (treatment + time + interaction)----
dtl <- list(data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_UVB_0.N,
                       X = dt1$X25w_UVB_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_UVB_1.N,
                       X = dt1$X25w_UVB_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_CON_0.N,
                       X = dt1$X25w_CON_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_CON_1.N,
                       X = dt1$X25w_CON_1.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_SFN_0.N,
                       X = dt1$X25w_SFN_0.X),
            data.table(dt1[, c("chr", "pos")],
                       N = dt1$X25w_SFN_1.N,
                       X = dt1$X25w_SFN_1.X))
dtl

BSobj <- makeBSseqData(dat = dtl,
                       sampleNames = snames)
BSobj

design <- data.table(trt = rep(rep(c("UVB", 
                                     "Ctrl", 
                                     "SFN"),
                                   each = 2),
                               1))
design$trt <- factor(design$trt,
                     levels = unique(design$trt))
design

DMLfit <- DMLfit.multiFactor(BSobj = BSobj, 
                             design = design,
                             formula = ~ trt)

# a. Control vs. UVB at Week 25----
DMLtest.Ctrl.UVB.W25 <- DMLtest.multiFactor(DMLfit,
                                            coef = "trtCtrl")

DMLtest.Ctrl.UVB.W25$chr <- as.numeric(as.character(DMLtest.Ctrl.UVB.W25$chr))
head(pct[, c("X25w_CON_0",
             "X25w_CON_1",
             "X25w_UVB_0",
             "X25w_UVB_1")])
DMLtest.Ctrl.UVB.W25 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                    pct[, c("X25w_CON_0",
                                                            "X25w_CON_1",
                                                            "X25w_UVB_0",
                                                            "X25w_UVB_1")]),
                                         DMLtest.Ctrl.UVB.W25,
                                         by = c("chr",
                                                "pos")))
DMLtest.Ctrl.UVB.W25$mu.ctrl25w <- (DMLtest.Ctrl.UVB.W25$X25w_CON_0 + DMLtest.Ctrl.UVB.W25$X25w_CON_1)/2
DMLtest.Ctrl.UVB.W25$mu.uvb25w <- (DMLtest.Ctrl.UVB.W25$X25w_UVB_0 + DMLtest.Ctrl.UVB.W25$X25w_UVB_1)/2

# Mean vs. difference----
DMLtest.Ctrl.UVB.W25$mu <- (DMLtest.Ctrl.UVB.W25$mu.ctrl25w + DMLtest.Ctrl.UVB.W25$mu.uvb25w)/2
DMLtest.Ctrl.UVB.W25$diff <- (DMLtest.Ctrl.UVB.W25$mu.ctrl25w - DMLtest.Ctrl.UVB.W25$mu.uvb25w)

# Green: hypomethylated in UVB; Red: hypermethylated in UVB
# NOTE: FLIP THE FIGURE!
dtp1 <- DMLtest.Ctrl.UVB.W25[!is.na(stat), ]
tiff(filename = "tmp/skin_uvb_maplot_uvb-ctrl_w25.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot((-1)*dtp1$diff ~ dtp1$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "Proportion of Methylation in UVB vs. Control\nat Week 25, FDR < 0.1")
points((-1)*dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff > 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff > 0] ,
       pch = "x",
       col = "red")
points((-1)*dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff < 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff < 0] ,
       pch = "x",
       col = "green")
abline(h = c(-0.2, 0.2),
       lty = 2)
graphics.off()

# b. SFN vs. UVB at Week 25----
DMLtest.SFN.UVB.W25 <- DMLtest.multiFactor(DMLfit,
                                           coef = "trtSFN")
head(DMLtest.SFN.UVB.W25)
DMLtest.SFN.UVB.W25$chr <- as.numeric(as.character(DMLtest.SFN.UVB.W25$chr))

DMLtest.SFN.UVB.W25 <- data.table(merge(data.table(dt1[, gene:CpG],
                                                   pct[, c("X25w_SFN_0",
                                                           "X25w_SFN_1",
                                                           "X25w_UVB_0",
                                                           "X25w_UVB_1")]),
                                        DMLtest.SFN.UVB.W25,
                                        by = c("chr",
                                               "pos")))
DMLtest.SFN.UVB.W25$mu.sfn25w <- (DMLtest.SFN.UVB.W25$X25w_SFN_0 + DMLtest.SFN.UVB.W25$X25w_SFN_1)/2
DMLtest.SFN.UVB.W25$mu.uvb25w <- (DMLtest.SFN.UVB.W25$X25w_UVB_0 + DMLtest.SFN.UVB.W25$X25w_UVB_1)/2

# Mean vs. difference----
DMLtest.SFN.UVB.W25$mu <- (DMLtest.SFN.UVB.W25$mu.sfn25w + DMLtest.SFN.UVB.W25$mu.uvb25w)/2
DMLtest.SFN.UVB.W25$diff <- (DMLtest.SFN.UVB.W25$mu.sfn25w - DMLtest.SFN.UVB.W25$mu.uvb25w)

# Green: hypomethylated in UVB; Red: hypermethylated in UVB
dtp1 <- DMLtest.SFN.UVB.W25[!is.na(stat), ]
tiff(filename = "tmp/skin_uvb_maplot_sfn_uvb_w25.tiff",
     height = 6,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
plot(dtp1$diff ~ dtp1$mu,
     pch = ".",
     xlab = "Mean",
     ylab = "Difference",
     main = "Proportion of Methylation in SFN vs. UVB\nat Weeks 25, FDR < 0.1")
points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff > 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff > 0] ,
       pch = "x",
       col = "red")
points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff < 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff < 0] ,
       pch = "x",
       col = "green")
abline(h = c(-0.2, 0.2),
       lty = 2)
graphics.off()

# Heatmaps----
uvb.hyper.1 <- unique(DMLtest.Ctrl.UVB.W25[fdrs < 0.1 & diff < 0]$pos)
uvb.hyper.2 <- unique(DMLtest.SFN.UVB.W25[fdrs < 0.1 & dtp1$diff < 0]$pos)
hyper <- uvb.hyper.1[uvb.hyper.1 %in% uvb.hyper.2]
l1 <- DMLtest.Ctrl.UVB.W25[pos %in% hyper, c("gene", 
                                             "pos", 
                                             "diff")]
l1$diff <- (-1)*l1$diff
l2 <- DMLtest.SFN.UVB.W25[pos %in% hyper, c("gene", 
                                            "pos", 
                                            "diff")]
# Venn diagram----
length(uvb.hyper.1)
length(uvb.hyper.2)
length(hyper)

uvb.hypo.1 <- unique(DMLtest.Ctrl.UVB.W25[fdrs < 0.1 & diff > 0]$pos)
uvb.hypo.2 <- unique(DMLtest.SFN.UVB.W25[fdrs < 0.1 & dtp1$diff > 0]$pos)
hypo <- uvb.hypo.1[uvb.hypo.1 %in% uvb.hypo.2]
l3 <- DMLtest.Ctrl.UVB.W25[pos %in% hypo, c("gene", 
                                            "pos", 
                                            "diff")]
l3$diff <- (-1)*l3$diff
l4 <- DMLtest.SFN.UVB.W25[pos %in% hypo, c("gene", 
                                           "pos", 
                                           "diff")]
# Venn diagram----
length(uvb.hypo.1)
length(uvb.hypo.2)
length(hypo)

ll <- merge(l1,
            l2,
            by = c("gene", 
                   "pos"))
ll2 <- merge(l3,
             l4,
             by = c("gene", 
                    "pos"))
ll <- rbind.data.frame(ll,
                       ll2)

ll$reg <- paste(ll$gene,
                ll$pos,
                sep = "_")
write.csv(ll,
          file = "tmp/w25.csv")

ll <- melt.data.table(data = ll,
                      id.vars = c(1, 5),
                      measure.vars = 3:4,
                      variable.name = "Comparison",
                      value.name = "Methyl Diff")
ll$reg <- factor(ll$reg)
ll$Comparison <- factor(ll$Comparison,
                        levels = rev(c("diff.x",
                                       "diff.y")),
                        labels = rev(c("UVB - Ctrl",
                                       "SFN - UVB")))

lvls <- ll[ll$Comparison == "UVB - Ctrl", ]
ll$reg <- factor(ll$reg,
                 levels = lvls$reg[order(lvls$`Methyl Diff`)])
ll

p1 <- ggplot(data = ll) +
  coord_polar("y",
              start = 0,
              direction = -1) +
  geom_tile(aes(x =  as.numeric(Comparison),
                y = reg, 
                fill = `Methyl Diff`),
            color = "white") +
  geom_text(data = ll[Comparison == "UVB - Ctrl", ],
            aes(x = rep(1.75,
                        nlevels(reg)),
                y = reg,
                label = unique(reg),
                angle = 90 + seq(from = 0,
                                 to = 360,
                                 length.out = nlevels(reg))[as.numeric(reg)]),
            hjust = 0,
            color = "black") +
  scale_fill_gradient2(low = "red", 
                       high = "green", 
                       mid = "grey", 
                       midpoint = 0, 
                       name = "Methyl Diff") +
  scale_x_continuous(limits = c(0, 
                                max(as.numeric(ll$Comparison)) + 0.5),
                     expand = c(0, 0)) + 
  scale_y_discrete("",
                   expand = c(0, 0)) +
  ggtitle("Changes in Controls Over Time") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

tiff(filename = "tmp/skin_uvb_heatmap_w25.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# sessionInfo()
# sink()