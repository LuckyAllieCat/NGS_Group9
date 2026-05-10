getwd()
setwd("/home/student18/project/")
data <- readRDS("/opt/omicsdata/datasets/SRP041538.rds")

data

assays(data)
# To visualize the raw count matrix we can use assay() function
class(assay(data, "raw_counts"))

# We visualize the first six genes and the first six samples
assay(data, "raw_counts")[1000:1010, ]

colnames(data)

class(colData(data))
colData(data)[1:6, 1:10]
# We use table() function to count living and deceased patients
table(data$sra.sample_attributes)

vars <- colData(data)[, grepl("sra.sample_", colnames(colData(data)))]
summary(data.frame(vars))

# There are 3467 genes whose sum is equal to zero
table(rowSums(assay(data, "raw_counts")) == 0)

# We keep the genes with the sum > 0
data <- data[rowSums(assay(data, "raw_counts")) > 0, ]



raw_counts <- assay(data, "raw_counts")

# For faster plotting, use only the first 20 samples
boxplot(raw_counts[, 1:20])

data$sra.sample_attributes = as.factor(data$sra.sample_attributes)
col <- data$sra.sample_attributes
levels(col)
levels(col) <- c(1:2)

length(data$sra.experiment_attributes)
boxplot(log1p(raw_counts)[, 1:170], col = col, xlab = "Sample", ylab = "log1p expression", ylim = c(0, 20),  las = 2, cex.axis = 0.7)
legend("top", legend = c("COPD", "Normal"), col = 1:2, pch = 20, cex = 0.7, ncol = 2)


limma::plotMA(log1p(raw_counts), main = "", xlab = "A", ylab = "M")
MDPlot(raw_counts, c(10, 12))

plotRLE(raw_counts, outline = FALSE, las = 2,
        ylab = "RLE", main = "RLE of raw counts", col = col)


pca <- prcomp(t(log1p(raw_counts)))
summary(pca)$importance[, 1:20]
autoplot(pca)

screeplot(pca, type = c("lines"))



#GC evaluation
ensid <- as.character(rowRanges(data)$gene_id)
ensid <- substr(ensid, 1, 15)

# we specify genome assembly hg38 and org.db as mode to retrieve the information
gc <- getGeneLengthAndGCContent(id = ensid, org = "hg38", mode = "org.db")

table(is.na(gc[, 2]))

gc <- gc[!is.na(gc[, 2]), ]

rowData(data)$ensid <- ensid
# Remove the genes without gc content information from the SummarizedExperiment object
data <- data[rowData(data)$ensid  %in% rownames(gc), ]

raw_counts <- assay(data, "raw_counts")
row.names(raw_counts) <- row.names(gc)

table(duplicated(row.names(raw_counts)))

set <- newSeqExpressionSet(raw_counts,
                           phenoData = AnnotatedDataFrame(data.frame(conditions = factor(1:ncol(raw_counts)),
                                                                     row.names = colnames(raw_counts))),
                           featureData = AnnotatedDataFrame(data.frame(gc = gc[, 2], l = gc[, 1])))
set

s <- EDASeq::biasPlot(set, "gc", ylim = c(0, 10), log = TRUE)



lrt <- log(raw_counts[, 1] + 1) - log(raw_counts[, 3] + 1)
biasBoxplot(lrt, gc[, 2], outline = FALSE, xlab = "GC-content", ylab = "Gene expression", las = 2, cex.axis = 0.5)

lrt <- log(raw_counts[, 1] + 1)- log(raw_counts[, 122] + 1)
biasBoxplot(lrt, gc[, 2], outline = FALSE, xlab = "GC-content", ylab = "Gene expression", las = 2, cex.axis = 0.5)

#filtering counts<10
data_sel <- data[rowSums(assay(data, "raw_counts")) >= 10, ]
set_sel <- set[rownames(set) %in% rowData(data_sel)$ensid, ]

#within sample normalisation
before <- biasPlot(set_sel, "gc", ylim = c(0, 10), log = TRUE)

wt <- withinLaneNormalization(set_sel, "gc" , which = "upper", offset = TRUE)

after <- biasPlot(wt, "gc", ylim = c(0, 10), log = TRUE)

#between sample
#different types of normalization
# offset = TRUE: normalized value returns as an offset leaving the original counts unchanged
set_full <- betweenLaneNormalization(set, which = "full", offset = TRUE)
set_uq <- betweenLaneNormalization(set, which = "upper", offset = TRUE)

# save the normalized matrix into an object
full <- set_full@assayData$normalizedCounts
uq_edaseq <- set_uq@assayData$normalizedCounts

raw_counts <- counts(set)

lbs <- calcNormFactors(raw_counts, method = "RLE")
rle <- cpm(raw_counts, lib.size = lbs * colSums(raw_counts))

lbs <- calcNormFactors(raw_counts, method = "TMM")
tmm <- cpm(raw_counts, lib.size = lbs * colSums(raw_counts))

lbs <- calcNormFactors(raw_counts, method = "upperquartile")
uq <- cpm(raw_counts, lib.size = lbs * colSums(raw_counts))

tc <- cpm(raw_counts)

#visualising
boxplot(log2(raw_counts + 1), outline = FALSE,  las=2,
        ylab = "log(count+1)", main = "Before normalization", xaxt = 'n')

boxplot(log2(tc + 1), outline = FALSE, las = 2,
        ylab = "log(count+1)", main = "Boxplot Total count", xaxt='n')

boxplot(log2(rle + 1), outline=FALSE,  las=2,
        ylab = "log(count+1)", main = "Boxplot RLE", xaxt = 'n')

boxplot(log2(tmm + 1), outline=FALSE,  las=2,
        ylab = "log(count+1)", main = "Boxplot TMM", xaxt = 'n')

boxplot(log2(uq + 1), outline=FALSE,  las=2,
        ylab = "log(count+1)", main = "Boxplot Upper quartile", xaxt = 'n')

boxplot(log2(full + 1), outline=FALSE,  las=2,
        ylab = "log(count+1)", main = "Boxplot full", xaxt = 'n')

#rle plots
col <- data$sra.sample_attributes
levels(col)
levels(col) <- c(1:2)

EDASeq::plotRLE(raw_counts, outline = FALSE,  las=2, col=col,
                ylab = "RLE", main = "RLE before normalization", xaxt='n')

EDASeq::plotRLE(tc, outline = FALSE,  las=2, col=col,
                ylab = "RLE", main = "RLE total count", xaxt='n')

EDASeq::plotRLE(tmm, outline = FALSE,  las=2, col=col,
                ylab = "RLE", main = "RLE TMM", xaxt='n')

EDASeq::plotRLE(rle, outline = FALSE,  las=2, col=col,
                ylab = "RLE", main = "RLE RLE", xaxt='n')

EDASeq::plotRLE(uq, outline = FALSE,  las=2, col=col,
                ylab = "RLE", main = "RLE upper quartile", xaxt='n')

EDASeq::plotRLE(full, outline = FALSE,  las=2, col=col,
                ylab = "RLE", main = "RLE full", xaxt='n')


#correlation plot
colnames <- colnames(tmm)
colnames(tmm) <- substr(colnames(tmm), 9, 12)
m <- cor(log1p(tmm))

corrplot(m, method = "color", order = 'hclust', tl.cex=0.2)

#pca analysis
data_pca <- prcomp(t(log1p(tmm)))
tmp <- summary(data_pca)
tmp$importance[, 1:20]

autoplot(data_pca)
screeplot(data_pca, type = c("lines"))  

EDASeq::plotPCA(tc, k=3, labels=F, col=as.numeric(col), pch=20)
EDASeq::plotPCA(tmm, k=3, labels=F, col=as.numeric(col), pch=20)
EDASeq::plotPCA(rle, k=3, labels=F, col=as.numeric(col), pch=20)

##################################
# mean-variance relationship
data_prox = data_sel
design <- model.matrix(~ sra.sample_attributes, data = colData(data_prox))
head(design)
# Create a DGEList object
dge <- DGEList(assay(data_prox, "raw_counts"))

dge <- calcNormFactors(dge, method="TMM")
dge <- estimateDisp(dge, design)

plotMeanVar(dge, show.raw.vars = TRUE, show.ave.raw.vars = F, NBline = TRUE, show.tagwise.vars = TRUE)

#identification of degs
fit <- glmFit(dge, design)
fit$offset[,1:20]

res <- glmLRT(fit, coef=2)
topTags(res, sort.by="logFC")

tab <- topTags(res, n=Inf)
tab <- tab$table
head(tab)

de <- tab[tab$FDR <0.05,]
dim(de)

#MD plot
plotMD(res)

#volcano plot
names <- row.names(de)

tab$sign <- NA
tab$sign[tab$logFC < 0] <- "Down"
tab$sign[tab$logFC > 0] <- "Up"
tab$sign[tab$FDR >= 0.05] <- "NoSig"
tab$gene <- rownames(tab)

ggplot(tab, aes(x = logFC, y = -log10(PValue), color = sign)) + geom_point() + theme_classic() +
  scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NoSig" = "black")) +
  geom_text_repel(data = tab[1:20, ], aes(label = gene), size = 2.5, max.overlaps = 20) 

#p-calue distribution
hist(tab$PValue, main="p-values distribution")

###CLUSTER ANALYSIS####
#heatmap#
raw_count <- assay(data_prox, "raw_counts")
lbs <- calcNormFactors(raw_count, method = "TMM")
tmm <- cpm(raw_count, lib.size = lbs * colSums(raw_count))

data_degs <- tmm[names, ]

Col <- data.frame(Subtype=data_prox$sra.sample_attributes)
rownames(Col) <- colnames(data_prox)

ann_colors = list(Subtype = c("disease state;;Normal|source_name;;Normal Lung Tissue|tissue;;lung"="blue","disease state;;COPD|source_name;;COPD lung tissue|tissue;;lung"="lightblue"))

g <- pheatmap(log1p(data_degs), show_colnames = FALSE, scale = "row", show_rownames = FALSE,
              drop_levels = TRUE, annotation_colors = ann_colors, annotation_col = Col)

#Silhouette index

cut <- cutree(g$tree_col, k = 2)

dist_mat <- dist(t(log1p(data_degs)))
plot(silhouette(cut, dist_mat), col=c("#2E86C1", "#EC7063"), border=NA)

########################################
###pathway analysis###
design <- model.matrix(~ sra.sample_attributes, data = colData(data_prox))
dim(data_sel)
dge <- DGEList(raw_counts)
dge <- calcNormFactors(dge, method = "TMM")

dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)

res <- glmLRT(fit, coef = 2)

topTags(res,sort.by = "logFC")

top <- topTags(res, n = Inf)$table
deg <- top[top$FDR <= 0.05,]
dim(deg)

universo <- row.names(top)
sign <- row.names(deg)

universo.entrez <- mapIds(org.Hs.eg.db, # genome wide annotation for Human that contains Entrez Gene identifiers
                          keys = universo, #Column containing Ensembl gene ids
                          column = "ENTREZID",
                          keytype = "ENSEMBL")  

sign.entrez <- mapIds(org.Hs.eg.db,
                      keys = sign, #Column containing Ensembl gene ids
                      column = "ENTREZID",
                      keytype = "ENSEMBL")
universo.entrez <- as.character(universo.entrez)
sign.entrez <- as.character(sign.entrez)

length(sign.entrez)
length(universo.entrez)

#kegg pathways
kk <- enrichKEGG(gene = sign.entrez,
                 universe = universo.entrez,
                 organism = 'hsa',
                 keyType = 'ncbi-geneid',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
dim(kk)
kk@result[1:10,]

barplot(kk) #Quickly rank terms by the number of genes involved
dotplot(kk) #Useful when you want to distinguish between a term that is significant but contains few of your genes vs. one where your DEGs dominate the pathway.
cnetplot(kk) # Useful to see gene overlap between terms and identify multi-pathway genes
heatplot(kk, showCategory = 10) #Reveals which genes are shared across multiple terms and which are unique to one pathway

edo <- pairwise_termsim(kk)
emapplot(edo) #Enrichment Map, terms are nodes; edges connect terms that share genes
obj <- upsetplot(edo, 8)
obj #Use this to understand the unique vs. shared gene membership across enriched categories.

pathview(gene.data = sign.entrez, pathway.id = "hsa03010", species = "hsa")
knitr::include_graphics("hsa03010.pathview.png")

#reactome database
kk <- enrichPathway(gene = sign.entrez,
                    universe = universo.entrez,
                    organism = 'human',
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH")
kk@result[1:10, 1:6]

barplot(kk)
dotplot(kk)
cnetplot(kk)
heatplot(kk, showCategory = 10)

edo <- pairwise_termsim(kk)
emapplot(edo)
upsetplot(edo, 5)

#gene onthology
ego <- enrichGO(gene = sign.entrez,
                universe = universo.entrez,
                OrgDb = 'org.Hs.eg.db',
                keyType = "ENTREZID",
                ont = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                readable = TRUE)

ego@result[1:10,1:6]

barplot(ego)
upsetplot(ego, 5)

#gsea

tt <- top$logFC
names(tt) <- universo.entrez

table(duplicated(names(tt)))
tt <- tt[!duplicated(names(tt))]

sort_tt <- sort(tt, decreasing = TRUE)

kk2 <- gseKEGG(geneList = sort_tt,
               organism = "hsa",
               pvalueCutoff = 0.05)

kk2@result[1:10, 1:7] # The key output is the Normalized Enrichment Score (NES)

gseaplot(kk2, geneSetID = "hsa03010") #Gene set enriched in down-regulated genes (to change)
ridgeplot(kk2, label_format = 50, showCategory = 20) #Shows the distribution of the ranking statistic (logFC) for genes belonging to each enriched gene set. Useful for interpreting the direction of enrichment in GSEA results.
gseaplot(kk2, geneSetID = "hsa00053") #Gene set enriched in top-regulated genes (to change)

#GO categories

kkGO <- gseGO(geneList = sort_tt,
              ont = "BP",
              OrgDb='org.Hs.eg.db',
              keyType = "ENTREZID",
              pAdjustMethod = "BH")

kkGO@result[1:10, 1:6]

gseaplot(kkGO, geneSetID = "GO:0006334")
ridgeplot(kkGO, label_format = 50, showCategory = 20) # all shown gene sets are significant well beyond the resolution of the test --> To rank and interpret results look at the NES instead of the p.adjust
dotplot(kkGO)







































