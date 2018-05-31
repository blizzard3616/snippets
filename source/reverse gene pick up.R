# |--------------------------------------------------|
# | Script:   Pick up reverse gene from list template|
# | Scientist: Ran                                   |
# | Created:  04/13/2018                             |
# |--------------------------------------------------|

#1 load data table
#2 find out datatset cut line

# Examples, UA+UVb vs UVB group RNA seq reverse pickup
#Step 1: find gene up-and down-regulated in UVB group
#Step 2: find gene up-and down-regulated in UA+UVB group
#Step 3: fine up in (1) and down in (2) or down in (1) and up in (2)



#Control vs UVB


#la means UA downregulated mRNA
#lar means UA upregulated mRNA
la <- dt1$gene[dt1$`log2(fold_change)` >= 1]
lar <- dt1$gene[dt1$`log2(fold_change)` <= -1]

# 2. UA+UVB vs UVB
#lb means UVB upregulated mRNA
#lbr means UVB downregulated mRNA
lb <- dt2$gene[dt2$`log2(fold_change)` <= -1]
lbr <- dt2$gene[dt2$`log2(fold_change)` >= 1]

down.up <- unique(la[la %in% lb])
up.down <- unique(lar[lar %in% lbr])

dt1[dt1$gene %in% down.up, ]
dt1[dt1$gene %in% up.down, ]
dt2[dt2$gene %in% down.up, ]
dt2[dt2$gene %in% up.down, ]

geneList <- unique(c(down.up,
                     up.down))


t1 <- subset(dt1, 
             gene %in% geneList, 
             select = c("gene",
                        "log2(fold_change)"))
t1$comp <- "UVB"

t2 <- subset(dt2, 
             gene %in% geneList, 
             select = c("gene",
                        "log2(fold_change)"))
t2$comp <- "UA+UVB"

tbl1 <- rbindlist(list(t1, t2))
tbl1$comp <- factor(tbl1$comp,
                    levels = unique(tbl1$comp))
summary(tbl1)

tbl1$`log2(fold_change)`[tbl1$`log2(fold_change)` > 5] <- 5
tbl1$`log2(fold_change)`[tbl1$`log2(fold_change)` < -5] <- -5

dcast.data.table(data = tbl1,
                 gene ~ comp,
                 fun.aggregate = mean,
                 value.var = "log2(fold_change)")

lvls <- unique(tbl1$gene)[order(tbl1$`log2(fold_change)`)]
tbl1$gene <- factor(tbl1$gene,
                    levels = lvls)

p1 <- ggplot(data = tbl1) +
  geom_tile(aes(x = comp,
                y = gene,
                fill = `log2(fold_change)`)) +
  scale_fill_gradient2(low = "red", 
                       high = "green", 
                       mid = "black", 
                       midpoint = 0, 
                       limit = c(-5, 5), 
                       space = "Lab", 
                       name="Log2(Fold-Change)") +
  scale_x_discrete("Comparison") + 
  scale_y_discrete("Gene") +
  ggtitle("Week 25 UA gene reverse list")  +
  theme(axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5))
p1

tiff(filename = "",
     height = 25,
     width = 10,
     units = 'in',
     res = 300,
     compression = "")
# png(filename = "tmp/heatmap_ordered.png",
#      height = 15,
#      width = 5,
#      units = 'in',
#      res = 300)
print(p1)
graphics.off()

tmp <- data.frame(gene = unique(tbl1$gene))

write.csv(tmp, file="Week 25 UA reverse list.csv")
