# =============================================================================
# RNA-seq ANALYSIS PIPELINE — REVISED (INTERACTIVE VERSION)
# COPD vs Normal Lung Tissue — SRP041538
# Reproduction and improvement upon Kim et al. (2015)
#
# Same logic as Cleaner_script_revised.R but all plots go to the active R
# graphics device (e.g. the RStudio Plots pane) instead of being written to
# PDF files. Data outputs (CSV, RDS) are still saved under results/.
# =============================================================================


# =============================================================================
# SECTION 0: LIBRARIES, SEED, OUTPUT DIRECTORY
# =============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(limma)
  library(EDASeq)
  library(ggplot2)
  library(ggfortify)
  library(ggrepel)
  library(corrplot)
  library(pheatmap)
  library(cluster)
  library(RColorBrewer)
  library(dplyr)
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(pathview)
  library(GOSemSim)
  library(simplifyEnrichment)
})

set.seed(42)

out_dir  <- "results"
data_dir <- file.path(out_dir, "intermediate")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# SECTION 1: DATA LOADING AND INITIAL QC
# =============================================================================

# --- Load the dataset ---
data <- readRDS("/opt/omicsdata/datasets/SRP041538.rds")

# --- Clean two-level group factor (Normal = reference) ---
data$group <- factor(
  ifelse(grepl("Normal", data$sra.sample_attributes), "Normal", "COPD"),
  levels = c("Normal", "COPD")
)
cat("Samples per group:\n"); print(table(data$group))

# --- Remove genes with zero counts everywhere ---
data <- data[rowSums(assay(data, "raw_counts")) > 0, ]
raw_counts <- assay(data, "raw_counts")

# --- Named colour palette ---
grp_colors <- c(Normal = "#4DAF4A", COPD = "#E41A1C")
col_vec    <- grp_colors[as.character(data$group)]

# --- QC: raw_data boxplot ---
boxplot(raw_counts[, 1:20])

# --- QC: log1p boxplot ---
boxplot(log1p(raw_counts), col = col_vec, xlab = "Sample", ylab = "log1p expression", ylim = c(0, 20),  las = 2, cex.axis = 0.7)
legend("top", legend = names(grp_colors),
       fill = grp_colors, bty = "o", bg = "white", cex = 0.8)


# --- MA plot of raw counts ---
limma::plotMA(log2(raw_counts + 1), main = "MA plot (raw)",
              xlab = "A", ylab = "M")

# --- MD plot of raw counts ---
MDPlot(raw_counts, c(1, 3))
MDPlot(raw_counts, c(1, 10))

# --- RLE plot of raw counts ---
EDASeq::plotRLE(raw_counts, outline = FALSE, las = 2,
                col = col_vec, ylab = "RLE",
                main = "RLE of raw counts", xaxt = "n")

# --- PCA on raw data, coloured by group ---
pca_raw <- prcomp(t(log2(raw_counts + 1)))
summary(pca_raw)$importance[, 1:20]
autoplot(pca_raw, data = as.data.frame(colData(data)),
         colour = "group")


screeplot(pca_raw, type = "lines", main = "Scree (raw)")

# --- Clean Ensembl IDs ---
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
# =============================================================================
# SECTION 2: Within sample normalisation
# =============================================================================

#filtering counts<10
data_sel <- data[rowSums(assay(data, "raw_counts")) >= 10, ]
set_sel <- set[rownames(set) %in% rowData(data_sel)$ensid, ]

#within sample normalisation
before <- biasPlot(set_sel, "gc", ylim = c(0, 10), log = TRUE)
wt_counts <- withinLaneNormalization(set_sel, "gc" , which = "upper", offset = TRUE)
gc_corrected_counts  <- counts(wt_counts)
after <- biasPlot(wt_counts, "gc", ylim = c(0, 10), log = TRUE)

# =============================================================================
# SECTION 3: FILTERING AND NORMALISATION
# =============================================================================

#TMM normalisation
raw_counts <- counts(wt_counts)
lbs <- calcNormFactors(raw_counts, method = "TMM")
tmm <- cpm(raw_counts, lib.size = lbs * colSums(raw_counts))

#Boxplot
boxplot(log2(raw_counts + 1), outline = FALSE,  las=2,
        ylab = "log(count+1)", main = "Before normalization and after GC correction", xaxt = 'n')
boxplot(log2(tmm + 1), outline=FALSE,  las=2,
        ylab = "log(count+1)", main = "Boxplot TMM", xaxt = 'n')

#MA plot
limma::plotMA(log2(raw_counts + 1), main = "MA plot (raw)",
              xlab = "A", ylab = "M")
limma::plotMA(log2(tmm + 1), main = "MA plot (normalised)",
              xlab = "A", ylab = "M")

#MD plot
MDPlot(raw_counts, c(1, 3))
MDPlot(tmm, c(1, 3))

#RLE plot
EDASeq::plotRLE(raw_counts, outline = FALSE, las = 2,
                col = col_vec, ylab = "RLE",
                main = "RLE of raw counts", xaxt = "n")

EDASeq::plotRLE(tmm, outline = FALSE, las = 2,
                col = col_vec, ylab = "RLE",
                main = "RLE of raw counts", xaxt = "n")

#Correlation plot
tmm_plot <- tmm
colnames(tmm_plot) <- substr(colnames(tmm), 9, 12)
m <- cor(log1p(tmm_plot))
corrplot(m, method = "color", order = 'hclust', tl.cex = 0.2)



# =============================================================================
# SECTION 4: DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

# --- Design matrix with named coefficient ---
design <- model.matrix(~ group, data = colData(data_sel))
colnames(design)[2] <- "COPDvsNormal"

dge  <- DGEList(counts = gc_corrected_counts)
keep <- filterByExpr(dge, design)
dge  <- dge[keep, ]
dge  <- calcNormFactors(dge, method = "TMM")
# --- Dispersion estimation + modern QL test ---
dge <- estimateDisp(dge, design)
plotBCV(dge, main = "Biological coefficient of variation")

fit <- glmQLFit(dge, design, robust = TRUE)
res <- glmQLFTest(fit, coef = "COPDvsNormal")

plotQLDisp(fit, main = "QL dispersion")

# --- Annotated results table ---
tab <- topTags(res, n = Inf)$table
tab$ensembl <- rownames(tab)
tab$symbol  <- mapIds(org.Hs.eg.db, keys = tab$ensembl,
                      column = "SYMBOL", keytype = "ENSEMBL",
                      multiVals = "first")
tab$entrez  <- mapIds(org.Hs.eg.db, keys = tab$ensembl,
                      column = "ENTREZID", keytype = "ENSEMBL",
                      multiVals = "first")

de_05   <- tab[tab$FDR < 0.05, ]
de_01   <- tab[tab$FDR < 0.01, ]
sig_ens <- rownames(de_05)

cat(sprintf("DEG (FDR < 0.05): %d  (up: %d, down: %d)\n",
            nrow(de_05), sum(de_05$logFC > 0), sum(de_05$logFC < 0)))
cat(sprintf("DEG (FDR < 0.01): %d\n", nrow(de_01)))

# --- MD plot with significance highlighting ---
plotMD(res, status = decideTests(res, p.value = 0.05),
       main = "MD plot (FDR < 0.05 highlighted)")

# --- Volcano with HGNC symbols + thresholds ---
tab$sig <- "NoSig"
tab$sig[tab$FDR < 0.05 & tab$logFC > 0] <- "Up"
tab$sig[tab$FDR < 0.05 & tab$logFC < 0] <- "Down"

top_label <- tab %>%
  arrange(PValue) %>%
  filter(sig != "NoSig") %>%
  head(20)

fdr_line <- if (nrow(de_05) > 0) -log10(max(tab$PValue[tab$FDR < 0.05])) else NA

p_volcano <- ggplot(tab, aes(logFC, -log10(PValue), color = sig)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c(Down  = "#377EB8",
                                Up    = "#E41A1C",
                                NoSig = "grey70")) +
  geom_text_repel(data = top_label,
                  aes(label = symbol), size = 3,
                  max.overlaps = 25, show.legend = FALSE) +
  geom_vline(xintercept = c(-1, 1), lty = 2, color = "grey40") +
  theme_classic() +
  ggtitle("Volcano — COPD vs Normal (edgeR QL)")
if (!is.na(fdr_line))
  p_volcano <- p_volcano +
  geom_hline(yintercept = fdr_line, lty = 2, color = "grey40")
print(p_volcano)

# --- p-value histogram ---
hist(tab$PValue, breaks = 50, col = "lightblue",
     main = "p-value distribution (edgeR QL)", xlab = "p-value")

# --- Persist DE results ---
write.csv(tab,   file.path(out_dir, "edgeR_DE_full.csv"),    row.names = FALSE)
write.csv(de_05, file.path(out_dir, "edgeR_DEGs_FDR05.csv"), row.names = FALSE)
saveRDS(tab, file.path(data_dir, "DE_main.rds"))


# =============================================================================
# SECTION 5: CLUSTER ANALYSIS
# =============================================================================

log_tmm <- log1p(tmm)
# --- Heatmap of top-200 DEGs by FDR ---
top200 <- rownames(de_05)[order(de_05$FDR)][1:min(200, nrow(de_05))]
top200 <- top200[top200 %in% rownames(log_tmm)]
mat_h  <- log_tmm[top200, ]
rownames(mat_h) <- tab$symbol[match(top200, tab$ensembl)]

ann_col <- data.frame(
  Group     = data_sel$group,
  row.names = colnames(tmm)
)

ann_col <- ann_col[colnames(mat_h), , drop = FALSE]
ann_colors <- list(Group = grp_colors)

g <- pheatmap(mat_h,
              show_colnames = FALSE,
              show_rownames = FALSE,
              scale = "row",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              annotation_col    = ann_col,
              annotation_colors = ann_colors,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
              main = "Top 200 DEGs (FDR < 0.05)")

# --- Silhouette index of the 2-cluster cut ---
cut_2    <- cutree(g$tree_col, k = 2)
dist_mat <- dist(t(mat_h))
sil_2    <- silhouette(cut_2, dist_mat)
cat("Mean silhouette width (k=2):",
    round(mean(sil_2[, "sil_width"]), 3), "\n")

plot(sil_2, col = c("#2E86C1", "#EC7063"), border = NA,
     main = "Silhouette of 2-cluster partition (top-200 DEGs)")

# --- Within-COPD subtyping ---
copd_samples <- colnames(data_sel)[data_sel$group == "COPD"]
tmm_copd     <- log_tmm[top200, copd_samples]
dist_copd    <- dist(t(tmm_copd), method = "euclidean")
hclust_copd  <- hclust(dist_copd, method = "complete")

# Pick number of subgroups by silhouette
sil_vals <- sapply(2:6, function(k) {
  cl <- cutree(hclust_copd, k = k)
  mean(silhouette(cl, dist_copd)[, "sil_width"])
})
names(sil_vals) <- 2:6
cat("Silhouette by k (COPD subtypes):\n"); print(round(sil_vals, 3))

best_k         <- as.integer(names(which.max(sil_vals)))
copd_subgroups <- cutree(hclust_copd, k = best_k)
cat(sprintf("Chosen k = %d  (best silhouette = %.3f)\n",
            best_k, max(sil_vals)))
cat("Subgroup sizes:\n"); print(table(copd_subgroups))

plot(hclust_copd, labels = FALSE,
     main = sprintf("Hierarchical clustering of COPD samples (k=%d)", best_k))
rect.hclust(hclust_copd, k = best_k, border = "red")

Col_copd <- data.frame(Subgroup  = factor(copd_subgroups),
                       row.names = copd_samples)

pheatmap(tmm_copd,
         show_colnames = FALSE,
         show_rownames = FALSE,
         scale = "row",
         clustering_distance_cols = "euclidean",
         annotation_col = Col_copd,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         main = "Top 200 DEGs across COPD samples (subgroup-coloured)")

# --- DE between COPD subgroups (no double-dipping, joint test) ---
subgroup_factor <- factor(copd_subgroups)
design_copd     <- model.matrix(~ subgroup_factor)

dge_copd <- DGEList(gc_corrected_counts[, copd_samples], group = subgroup_factor)
keep_copd <- filterByExpr(dge_copd, design_copd)
dge_copd  <- dge_copd[keep_copd, , keep.lib.sizes = FALSE]
dge_copd  <- calcNormFactors(dge_copd, method = "TMM")
dge_copd  <- estimateDisp(dge_copd, design_copd)

fit_copd <- glmQLFit(dge_copd, design_copd, robust = TRUE)
res_copd <- glmQLFTest(fit_copd, coef = 2:ncol(design_copd))

tab_copd <- topTags(res_copd, n = Inf)$table
deg_copd <- tab_copd[tab_copd$FDR < 0.05, ]
cat(sprintf("DEG between COPD subgroups (FDR<0.05): %d\n", nrow(deg_copd)))

saveRDS(tab_copd, file.path(data_dir, "DE_copd_subgroups.rds"))


# =============================================================================
# SECTION 6: PATHWAY ANALYSIS
# =============================================================================

# --- ID mapping ---
universe_entrez <- as.character(mapIds(org.Hs.eg.db,
                                       keys    = rownames(tab),
                                       column  = "ENTREZID",
                                       keytype = "ENSEMBL",
                                       multiVals = "first"))

sig_entrez <- as.character(mapIds(org.Hs.eg.db,
                                  keys    = sig_ens,
                                  column  = "ENTREZID",
                                  keytype = "ENSEMBL",
                                  multiVals = "first"))

cat("DEG with Entrez ID:", length(sig_entrez), "\n")
cat("Universe size:",      length(universe_entrez), "\n")

# --- KEGG over-representation ---
kk <- enrichKEGG(gene          = sig_entrez,
                 universe      = universe_entrez,
                 organism      = "hsa",
                 keyType       = "ncbi-geneid",
                 pvalueCutoff  = 0.05,
                 pAdjustMethod = "BH")

kk@result[1:10, ]
barplot(kk)
dotplot(kk)
cnetplot(kk)
heatplot(kk, showCategory = 10)

edo <- pairwise_termsim(kk)
emapplot(edo)
obj <- upsetplot(edo, 8)
obj

pathview(gene.data  = sig_entrez,
         pathway.id = "hsa03010",
         species    = "hsa")

# --- Reactome over-representation ---
kk_reactome <- enrichPathway(gene          = sig_entrez,
                             universe      = universe_entrez,
                             organism      = "human",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH",
                             readable      = TRUE)

kk_reactome@result[1:10, 1:6]
barplot(kk_reactome)
dotplot(kk_reactome)
cnetplot(kk_reactome)
heatplot(kk_reactome, showCategory = 10)

edo_r <- pairwise_termsim(kk_reactome)
emapplot(edo_r)
upsetplot(edo_r, 5)

# --- GO Biological Process over-representation ---
ego <- enrichGO(gene          = sig_entrez,
                universe      = universe_entrez,
                OrgDb         = "org.Hs.eg.db",
                keyType       = "ENTREZID",
                ont           = "BP",
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                readable      = TRUE)

ego@result[1:10, 1:6]
barplot(ego)
upsetplot(ego, 5)

# --- GSEA: ranked by logFC ---
tt <- tab$logFC
names(tt) <- universe_entrez
tt      <- tt[!duplicated(names(tt))]
sort_tt <- sort(tt, decreasing = TRUE)

# GSEA — KEGG
kk2 <- gseKEGG(geneList     = sort_tt,
               organism     = "hsa",
               pvalueCutoff = 0.05)

kk2@result[1:10, 1:7]
gseaplot(kk2, geneSetID = "hsa03010")
ridgeplot(kk2, label_format = 50, showCategory = 20)

# GSEA — GO BP
kkGO <- gseGO(geneList      = sort_tt,
              OrgDb         = "org.Hs.eg.db",
              ont           = "BP",
              keyType       = "ENTREZID",
              pAdjustMethod = "BH")

kkGO@result[1:10, 1:6]
gseaplot(kkGO, geneSetID = "GO:0006334")
ridgeplot(kkGO, label_format = 50, showCategory = 20)
dotplot(kkGO)

saveRDS(list(KEGG = kk, Reactome = kk_reactome, GO_BP = ego,
             GSEA_KEGG = kk2, GSEA_GO = kkGO),
        file.path(data_dir, "enrichments.rds"))
# =============================================================================
# SECTION 7: GO TERM CLUSTERING (BioLattice analogue)
# =============================================================================

sig_go <- ego@result$ID[ego@result$p.adjust < 0.05]
cat("Significant GO BP terms:", length(sig_go), "\n")

if (length(sig_go) > 1) {
  sig_go_use <- head(sig_go, 500)
  d   <- godata("org.Hs.eg.db", ont = "BP")
  mat <- termSim(sig_go_use, sig_go_use,
                 semData = d, method = "Wang")
  
  set.seed(42)
  simplifyGO(mat,
             word_cloud_grob_param = list(max_width = 80),
             fontsize_range = c(6, 12))
}


# =============================================================================
# SECTION 8: ALTERNATIVE SPLICING
# =============================================================================
# Requires exon-level counts not present in SRP041538 RSE.
# Outline kept for completeness:
#
# library(recount3)
# rse_exon <- create_rse_manual(
#   project      = "SRP041538",
#   project_home = "data_sources/sra",
#   organism     = "human",
#   annotation   = "gencode_v26",
#   type         = "exon"
# )
# library(DEXSeq)
# # ... DEXSeq workflow on rse_exon ...


# =============================================================================
# SECTION 9: CONNECTIVITY MAP
# =============================================================================

# --- Top 150 up / 150 down DEGs by logFC, as HGNC symbols ---
up_ids   <- rownames(de_05[de_05$logFC > 0, ])
down_ids <- rownames(de_05[de_05$logFC < 0, ])
up_ids   <- up_ids[order(-de_05[up_ids,   "logFC"])][1:min(150, length(up_ids))]
down_ids <- down_ids[order( de_05[down_ids, "logFC"])][1:min(150, length(down_ids))]

up_symbols   <- tab$symbol[match(up_ids,   tab$ensembl)]
down_symbols <- tab$symbol[match(down_ids, tab$ensembl)]
up_symbols   <- up_symbols[!is.na(up_symbols)]
down_symbols <- down_symbols[!is.na(down_symbols)]

cat("=== UPREGULATED GENES (top 50) ===\n")
cat(paste(head(up_symbols, 50), collapse = "\n"), "\n\n")
cat("=== DOWNREGULATED GENES (top 50) ===\n")
cat(paste(head(down_symbols, 50), collapse = "\n"), "\n")

write.csv(data.frame(gene = up_symbols),
          file.path(out_dir, "cmap_upregulated_genes.csv"),   row.names = FALSE)
write.csv(data.frame(gene = down_symbols),
          file.path(out_dir, "cmap_downregulated_genes.csv"), row.names = FALSE)

# --- After downloading cmap_results.csv from clue.io, visualise top-20 ---
# cmap_results <- read.csv(file.path(out_dir, "cmap_results.csv"))
# print(
#   ggplot(head(cmap_results[order(cmap_results$score), ], 20),
#          aes(x = reorder(pert_iname, score), y = score, fill = score > 0)) +
#     geom_bar(stat = "identity") + coord_flip() + theme_classic() +
#     scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "blue")) +
#     labs(x = "Drug", y = "Connectivity Score",
#          title = "CMap — Top 20 drugs reversing COPD signature") +
#     theme(legend.position = "none")
# )


# =============================================================================
# SESSION INFO
# =============================================================================
sink(file.path(out_dir, "sessionInfo.txt"))
sessionInfo()
sink()