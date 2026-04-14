# NGS_Group9

April 1st 2026

Nikita
Alison

NGS Practical Report

Dataset:  SRP041538
Compare Disease State: COPD vs Control.  The information is stored in the $sra.sample_attributes in colData() object

Publication:  https://pubmed.ncbi.nlm.nih.gov/25834810/
Comprehensive Analysis of Transcriptome Sequencing Data in the Lung Tissues of COPD Subjects

Github:  https://github.com/LuckyAllieCat-web/NGS_Group9?tab=readme-ov-file#ngs_group9
Download data:  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE57148
Biological Problem:  Unclear mapping of molecular pathways in the COPD dysregulation.  

Experimental Problem:
Wet lab protocol:  RNA isolated from lung samples (COPD vs Control)
RNA-seq performed
Processing and filtering of genes
Statistical analysis and tests:  t-test


Analysis:
Exploratory data analysis
Normalization
Identification of differentially expressed genes
(DEGs)
Gene set analysis
Bioconductor Package

Notes (AW 14/4/26)

Computational Pipeline
The pipeline moves through five stages, each building on the outputs of the previous one.
Raw counts (GEO: GSE57148)
        ↓
  Pre-filtering
        ↓
  Exploratory Analysis
        ↓
  Normalization (FPKM + Upper Quartile)
        ↓
  Differential Expression (4 methods)
        ↓
  Gene Set Analysis (ORA + GSEA)
  
Stage 1 — Data Loading & Pre-filtering
Files used: GSE57148_raw_counts_GRCh38_p13_NCBI.tsv, GSE57148_series_matrix.txt, Human_GRCh38_p13_annot.tsv
The raw count matrix contains integer read counts for 39,376 genes across 189 samples. Group labels (COPD/Control) were extracted directly from the series matrix file rather than hardcoded, ensuring correctness. Two filtering steps were applied:
Remove all-zero genes — genes with zero counts across all 189 samples
Remove non-coding genes — using the GeneType == "protein-coding" column from the NCBI annotation file, retaining ~19,400 genes (matching the paper's approach of excluding non-coding genes)
Remove very low-count genes — rowSum < 10, reducing noise before normalization
A SummarizedExperiment object was built to keep counts and metadata coupled throughout the pipeline.

Stage 2 — Exploratory Analysis (before normalization)
Performed entirely on raw counts to establish a baseline and check data quality. Five types of visualization were produced:
Histogram of raw counts and log1p counts — confirms the expected right-skewed distribution
Boxplot of log1p counts across all 189 samples — checks that no sample has drastically different median expression
Library size bar chart — confirms comparable sequencing depth (~38M reads), consistent with the paper
MD/MA plots — checks for mean-expression-dependent bias in individual samples
RLE plot — the most sensitive quality check; medians deviate from zero before normalization as expected
Mean–variance plot — confirms overdispersion (variance > mean), justifying the use of negative binomial models
PCA — checks whether COPD and Control partially separate before any normalization

Stage 3 — Normalization
Within-sample: FPKM
FPKM was computed from scratch using the raw counts and the Length column from the annotation file:
FPKMij=xijLengthi[kb]×Nj[M reads]\text{FPKM}_{ij} = \frac{x_{ij}}{\text{Length}_i[\text{kb}] \times N_j[\text{M reads}]}FPKMij​=Lengthi​[kb]×Nj​[M reads]xij​​
This corrects simultaneously for gene length (longer genes get more reads regardless of expression) and sequencing depth (samples sequenced more deeply get higher counts). This is the within-sample normalization step used by the paper via Cufflinks.
Between-sample: Upper Quartile
The 75th percentile of the non-zero FPKM values in each sample was used as the scaling factor, applied column-wise. This is more robust than dividing by the total count (CPM) because it is less sensitive to the small number of very highly expressed genes that can dominate the sum. This directly reproduces the paper's stated normalization.
Additional methods for comparison
TMM, RLE, CPM, and Full Quantile were also computed (via edgeR and EDASeq) and evaluated using the same boxplot/RLE/PCA battery, allowing direct visual comparison of all methods.


Stage 4 — Differential Expression (4 methods)
Four independent statistical methods were run, ranging from the paper's original approach to modern best-practice alternatives:
Method 1 — t-test on FPKM-UQ (paper method) A gene-wise two-sample Welch t-test on log1p-transformed FPKM-UQ values, with Benjamini–Hochberg FDR correction at q < 0.01. This directly replicates the paper's reported approach and targets the 2,312 DEGs they found.
Method 2 — edgeR GLM (lab method) A negative binomial generalised linear model fitted with glmFit(), tested with a likelihood ratio test (glmLRT()). TMM normalization factors are incorporated into the model offset. This is more statistically appropriate than a t-test on normalized values because it models the count-generating process directly.
Method 3 — DESeq2 (new package) A second NB-based GLM with a different dispersion estimation strategy — maximum a posteriori shrinkage towards a fitted trend, visualized with plotDispEsts(). Log fold-changes are further shrunk with lfcShrink(type="apeglm") to reduce noise for low-count genes. DESeq2 uses its own median-of-ratios normalization internally.
Method 4 — NOISeq (new package)A non-parametric method making no distributional assumption. It builds an empirical noise distribution by comparing all within-condition pairs, then assigns each gene a probability qq q that its fold-change is real signal rather than noise. noiseqbio mode was used because biological replicates are available. Threshold: q≥0.8q \geq 0.8 q≥0.8.
A four-method consensus was computed: genes called DE by ≥ 3 of the 4 methods are the most robust DEG candidates, surviving multiple statistical frameworks.

Stage 5 — Gene Set Analysis (ORA + GSEA)
The gene ID format (Entrez integers) was already compatible with all databases — no conversion needed.
ORA — Over-Representation Analysis (equivalent to paper's DAVID): Fisher exact test asking whether DEGs are over-represented in a given gene set compared to all tested genes. Run with clusterProfiler on three databases: GO (BP, MF, CC), KEGG, and Reactome. Visualized with barplots, dotplots, cnetplots, and enrichment maps.

GSEA — Gene Set Enrichment Analysis (equivalent to paper's GSEA + MSigDB): All tested genes ranked by log fold-change (continuous, no threshold). The algorithm tests whether gene set members cluster at the top (up in COPD) or bottom (down in COPD) of the ranking. Key output is the Normalized Enrichment Score (NES). Run on KEGG and GO-BP. Individual enrichment plots produced for oxidative phosphorylation, protein catabolism, and chromatin modification — the three specific pathways highlighted in the paper.

Wet Lab Pipeline
The study recruited 189 male patients who required surgical lung resection for lung cancer and were registered in the Asan Biobank (Asan Medical Center, Seoul, South Korea) between January 2008 and November 2011.
Patients were divided into two groups based on post-bronchodilator spirometry:


COPD (n=98)

Control (n=91)

Age (years)
67.5 ± 6.4,
60.9 ± 9.5

Smoking (pack-years)
48.0 ± 22.0,
35.2 ± 17.2

FEV1 (%)
71.9 ± 13.4,
91.0 ± 12.4

FEV1/FVC
57.1 ± 7.8,
74.8 ± 4.3

The inclusion criterion for COPD was FEV1/FVC < 0.7; controls had normal spirometry. All subjects were smokers or ex-smokers.

Tissue Collection
RNA was isolated from apparently normal fresh-frozen lung tissue remote from the tumour site — meaning the tissue came from the same resected lung but was sampled far from the cancer, to avoid tumour-associated gene expression changes. This is a practical approach when you cannot ethically obtain lung biopsies from healthy individuals.


The wet lab protocol followed these steps:
RNA extraction — total RNA isolated from fresh-frozen tissue
Quality assessment — RNA integrity checked with an Agilent Bioanalyzer; purity assessed with a NanoDrop spectrophotometer
Library preparation — 1 µg of total RNA per sample processed with the Illumina TruSeq RNA library kit:
Poly-A selection (mRNA enrichment)
RNA fragmentation
Reverse transcription using random hexamer primers
Adapter ligation
Sequencing — 100 bp paired-end sequencing on the Illumina HiSeq 2000
Mean sequencing depth achieved: ~38.7 million reads per sample.


MIAME Compliance
The study satisfies MIAME (Minimum Information About a Microarray Experiment) criteria — though adapted for RNA-seq — because it reports the experimental design, sample characteristics, raw data (deposited in GEO as GSE57148), data processing steps, and analysis methods in sufficient detail to reproduce the results.
