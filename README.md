# sc-rnaseq-shiny-app

A point-and-click Shiny app for single-cell RNA-seq analysis with Seurat v5, covering QC, integration, clustering, subclustering, cell-type annotation, comparative differential expression between samples, and pathway analysis (GSEA and ORA).

---

## Features

### Data Loading

- **10x Genomics MTX**: upload the `barcodes`, `genes`/`features` and `matrix` files (`.tsv`/`.mtx`, gzipped or not). Gene symbols are taken from the second column of the genes file.
- **10x Genomics HDF5**: upload a single `.h5` file (read with `Read10X_h5`).
- **Directory path**: point to a local `filtered_feature_bc_matrix/` folder instead of uploading files (single-sample mode only).
- **Multiple samples**: tick *Multiple Datasets* and add samples one at a time (MTX or H5), each with its own name. Queued samples appear as chips, and at least 2 are required. The samples are merged with sample-prefixed cell barcodes, and RNA layers are split by sample for integration.
- **Automatic object setup**: cells and genes are pre-filtered at load time (genes detected in ≥ 3 cells, cells with ≥ 200 genes). `percent.mt` is calculated from the species-specific mitochondrial prefix. If that finds no mitochondrial genes, the app falls back to a case-insensitive `MT-`/`mt-`/`Mt-` search.
- **Reload a saved `.rds` object**: the app detects how far the object was already processed and unlocks the matching tabs:
  - **QC done**: a `pca` reduction and a `scale.data` layer are present
  - **Integration done**: a `harmony`, `integrated.cca` or `integrated.dr` reduction is present. The integration method (Harmony / CCA / SCTransform+CCA) is inferred from that reduction and set in the Integration tab.
  - **Clustering done**: a `umap` reduction and a `seurat_clusters` column are present
  - **ScType annotation**: a `sctype_classification` column is present
  - **Manual annotation**: a `manual_annotation` column is present. The cluster → label mapping is rebuilt from `seurat_clusters` (or `subcluster_labels`) so the manual-annotation text boxes are pre-filled.
  - A status line summarises what was found, for example `Loaded 12,345 cells | 2 sample(s) | 14 clusters | ScType ✓ | Manual annotation ✓`
- **Uploads of up to 500 MB per file.**

### Data Summary

- **Dataset overview**: total cells and genes, plus median features, UMIs and % mitochondrial reads per cell. After QC filtering, it also shows cells kept, removed and % retained for each sample.
- **Cells per Sample**: bar chart (multi-sample only).
- **Condition groups (optional)**: map each sample to a condition label, stored in the `condition` metadata column (for example, for replicate designs).
- **Cell metadata table**: a paginated table of every cell, downloadable as CSV.

### Counts Matrix

- **Matrix browser**: view raw counts, log-normalized data or scaled data for a chosen number of genes × cells, optionally filtered to one sample.
- **Full-matrix export**: downloads the full selected layer as a `.tsv` file, which avoids Excel's column limit with large cell counts.

### QC Filtering

- **Adjustable thresholds**: minimum features, maximum features and maximum % MT (defaults 200 / 6000 / 5%).
- **Applying filters** runs `NormalizeData` → `FindVariableFeatures` (2,000 genes) → `ScaleData` → `RunPCA` on the cells that pass.
- **Before / After toggle**: shows QC plots for the unfiltered or the filtered object.
- **QC plots**: violin plots of `nFeature_RNA`, `nCount_RNA` and `percent.mt`. Scatters of counts vs. % MT and counts vs. features. Per-sample QC violins (multi-sample only).
- **Filter summary**: cells before and after filtering, with the number and % removed.

### PCA

- **PCA plot**, colored by sample in multi-sample mode
- **Elbow plot** (40 PCs) to choose the number of dimensions for clustering
- **Feature loadings** for a selectable number of PCs
- **PC heatmaps** (`DimHeatmap`, 500 cells) for a selectable number of PCs

### Integration

- **None (PCA)**: no batch correction. This is the only option for a single sample.
- **Harmony**: `RunHarmony` on `sample`, with an **adjustable theta** (default 2, range 0–10).
- **CCA**: Seurat v5 `IntegrateLayers(method = CCAIntegration)`.
- **SCTransform + CCA**: re-runs `SCTransform` (regressing `percent.mt`) from the filtered object, then performs CCA integration on the SCT-normalized data.
- **Pre- vs. post-integration preview**: UMAPs colored by sample, and split to one sample at a time, on the unintegrated PCA. Once integration has run, a toggle switches them to a preview UMAP built on the integrated reduction.
- **Reproducibility**: all stochastic steps (integration, neighbors, clustering, UMAP) use a fixed seed (42).

### Clustering

- **Graph-based clustering**: `FindNeighbors` + `FindClusters` + `RunUMAP` on the chosen reduction, with adjustable **Dims** (default 15) and **Resolution** (default 0.5).
- **Identity selector**: color the UMAP by cluster number, ScType label, manual label or subcluster label. Options appear as each one becomes available, and the selected identity also drives Cluster Markers, Signature Scoring, Manual Annotation and Subclustering.
- **Color by identity or by sample.**
- **Cluster statistics table**: cell count and % per cluster, plus per-sample cell counts in multi-sample mode. Downloadable as CSV.

### Subclustering

The Subclustering panel in the Clustering tab re-analyses a subset of cells in three stages:

1. **Prep subset & run PCA**: pick one or **several parent clusters** from the active identity (cluster numbers, ScType labels, manual labels or a previous subcluster round). The subset gets fresh variable features, scaling and PCA. The app then shows an elbow plot, PC heatmaps and, for multi-sample data, an unintegrated UMAP by sample, with guidance on whether re-integration is needed.
2. **Integrate**: choose *Match global* (reuses the global method and Harmony theta), *None*, *Harmony* (with its own theta), *CCA* or *SCTransform+CCA*. A pre/post-integration toggle previews the effect. If the object was loaded from `.rds` and the original theta is unknown, *Match global* falls back to theta = 2 and shows a warning.
3. **Cluster**: set subcluster Dims and Resolution.

- **Parent-letter naming**: subclusters are named after their parent plus a letter, for example `3_a`, `3_b`, or `cDC1_a`, `cDC1_b` when subclustering a labelled cluster. When several parents are subclustered together, the prefix joins them, for example `cDC1+cDC2_a`.
- **Result outputs**: a labelled subcluster UMAP, a subcluster size table, sizes by sample, and a **parent-of-origin composition table** that shows which parent cluster each subcluster's cells came from (multi-parent runs only). All three tables are downloadable as CSV.
- **Subcluster markers**: `FindAllMarkers` comparing subclusters **against each other** within the subset, not against the whole dataset, with its own Top N / min pct / min log2FC settings, a marker heatmap and a Top Markers table (Top N and all markers as CSV).
- **Merge into main object**: writes the new labels to a `subcluster_labels` metadata column. Cells outside the subset keep their parent label. *Subcluster Label* then becomes available in every identity dropdown, and subclustering can be repeated on the merged labels. Merging clears results that depend on identities (markers, DE cache, GSEA/ORA, pseudobulk, signature scores).
- **Discard and retry**: throws away the subset without touching the main object.

### Cluster QC

- **Per-cluster QC**: violin plots of `nFeature_RNA`, `nCount_RNA` and `percent.mt` by cluster, plus counts vs. % MT and counts vs. features scatters, to catch low-quality or doublet-like clusters.

### Gene Expression

- **Feature, violin or dot plots** for one or more comma-separated genes. Gene names are matched case-insensitively; if there is no exact match, the first partial match is used.
- **Group by any identity** (cluster, ScType, manual, subcluster), with an optional filter to show only selected clusters.
- **Consistent views**: feature plots share the axis range of the full UMAP, and violin colors match the identity colors used elsewhere in the app.

### Cluster Markers

- **Find All Markers**: `FindAllMarkers` (positive markers only) for the active identity, with adjustable **Top N per cluster**, **min pct expressed** and **min log2FC** (Seurat tutorial defaults: 10 / 0.25 / 0.25).
- **Heatmap**: `DoHeatmap` of the top N markers per cluster by average log2FC. Label sizes and plot height scale with the number of genes and clusters.
- **Top Markers table**: filterable, sorted by cluster and log2FC. Download the top N or all markers as CSV.
- **All Genes by Cluster**: every marker for one selected cluster, sorted by log2FC, with an optional top-N cut-off (0 = all). Downloadable per cluster.

### Annotation

- **ScType (automated)**: scores each cluster against the ScType marker database for a chosen tissue: Immune system, Brain, Lung, Liver, Pancreas, Kidney, Eye, Heart, Adrenal, Intestine, Muscle, Placenta, Spleen, Stomach or Thymus. The ScType functions and database are fetched from the [sc-type GitHub repository](https://github.com/IanevskiAleksandr/sc-type) at run time.
  - Clusters whose top score is below a quarter of their cell count are labelled `Unknown`.
  - The results table lists each cluster's assigned type, score, cell count and up to 6 database marker genes that drove the call. Downloadable as CSV.
  - An annotated UMAP is shown, and *ScType Label* is added to the identity dropdowns.
- **Signature scoring**: score clusters against custom gene sets with `AddModuleScore`.
  - Load signatures from a **CSV file** with `Signature` and `Genes` columns (genes comma-separated). Gene names are matched case-insensitively to the dataset.
  - Signatures can also be added one at a time in the app, and removed individually or all at once.
  - Signatures with fewer than 2 genes found in the dataset are skipped, and missing genes are reported.
  - Outputs: a heatmap of average scores per cluster (z-scored per signature), a raw score table, and a **Signature Gene Inspector**. The inspector lists every gene with whether it was found, its mean expression and a High / Medium / Low / Absent contribution tier. The heatmap downloads as PNG; the scores and gene contributions download as CSV.
- **Manual annotation**: type a free-text label for each group of the active identity.
  - Each group shows hints: the ScType suggestion, the top 5 marker genes (after *Find All Markers*) and the top 3 scoring signatures (after Signature Scoring).
  - Text boxes are pre-filled with any existing label or the ScType call. Subclusters inherit their ancestor's manual label: if cluster `3` was labelled `cDC1`, subcluster `3_a` is pre-filled as `cDC1_a`.
  - Buttons: **Apply Labels**, **Reset** (back to ScType suggestions) and **Refresh Hints**.
  - Once applied, *Manual Label* appears in all identity dropdowns. The app also shows a manual-annotation UMAP and per-label cell counts.

### Comparative Analysis (2+ samples)

All sub-tabs share one **Color / group by** identity selector (cluster, ScType, manual or subcluster).

- **UMAP**: colored by sample; split by sample (one sample at a time, colored by identity or sample); colored by identity. Colors stay the same for each sample and identity across every view.
- **Composition**:
  - Stacked bar of the % of each cluster or cell type coming from each sample
  - Grouped bar of the % of each sample's cells in each cluster or cell type
  - Both are downloadable as PNG and CSV
- **Gene Expression**, for one or more genes:
  - **Split feature plot**: one column per sample, on a shared UMAP axis range.
  - **Split violin plot** and **split dot plot**, with a **scope toggle**: *Global* (all identities), *Per-cluster* (one identity, with a layout tuned for a single category) or *Selected* (several chosen identities). Violins drop genes with no variation in the plotted cells. Dot plots use one color gradient per sample, with a matching legend.
- **Conserved Markers**: `FindConservedMarkers` for a selected cluster, grouped by sample, finds genes that mark the cluster consistently across samples. Downloadable as CSV.
- **Pseudobulk scatter**: `AggregateExpression` by identity × sample, then a scatter of sample A vs. sample B for one cluster or cell type, with Pearson *r* and optional highlighted, labelled genes. The app warns if the identity has changed since the pseudobulk was built.
- **Condition DE**: `FindMarkers` (Seurat default test) between sample A and sample B cells **within one cluster or cell type**.
  - With more than 2 samples, pick the numerator and denominator samples. The log2FC direction is spelled out in the app.
  - Results can be viewed as a filterable table, with significant rows highlighted, or as a **volcano plot** with adjustable |log2FC| and adjusted-p thresholds and top-N gene labels per side.
  - **Cached results**: every comparison is kept in memory, so several cluster or sample comparisons can be fed to GSEA and ORA without re-running DE. The cache is cleared when you re-cluster or merge subclusters.
- **Sample display filter**: the sidebar *Active Samples* checkboxes hide samples from the UMAP, PCA-by-sample, violin, feature and comparative dot-plot views. The filter is **display-only**: statistics, composition and exports always use every sample.

### Pathway Analysis

- **GSEA**: `fgseaMultilevel` on any cached DE comparison.
  - Genes are ranked by **sign(log2FC) × −log10(p) × |log2FC|**, with p = 0 replaced by the smallest non-zero p-value.
  - Gene sets come from **MSigDB** via `msigdbr`: Hallmark, GO Biological Process / Cellular Component / Molecular Function, KEGG, or Reactome. Six species are available, chosen separately from the sidebar species.
  - The plot shows the top N significant pathways (padj < 0.05) **per direction** as a dot or bar plot of NES, colored by −log10(padj).
  - The full results, including leading-edge genes, download as CSV.
- **Pathway Explorer**: pick up to 20 pathways from the GSEA results (pre-filled with the top 5 significant).
  - **NES barplot**: colored by direction, with `*` marking padj < 0.05.
  - **Enrichment traces**: `plotEnrichment` running-score curves, one per pathway.
- **ORA**: `clusterProfiler::enricher` on the significant up- and down-regulated genes from a cached comparison, with adjustable adjusted-p cut-off, |log2FC| cut-off and top N terms.
  - Both directions are shown on **one combined plot** (bar or dot), with each term's −log10(adj p) signed by direction and colored by which sample it is higher in.
  - Species and collection are chosen independently of the GSEA tab.

### Figures and Export

- **600 DPI PNG downloads** for every plot.
- **Dynamic sizing**: canvas size scales with the number of labels, genes, identities and pathways, and with label length, so long cell-type or pathway names aren't clipped.
- **Matching plots**: UMAPs use a fixed panel size and colors that stay the same for each sample and identity, so on-screen and downloaded versions match.
- **Session export**: download the full Seurat object (`.rds`) to reload later, and all cell metadata (`.csv`).

---

## Supported Species

The sidebar species setting controls mitochondrial gene detection:

| Option | Mitochondrial prefix |
|---|---|
| Human | `MT-` |
| Mouse | `mt-` |
| Other | `MT-`, `mt-` or `Mt-` |

GSEA and ORA gene sets (MSigDB via `msigdbr`) are available for:

| Species |
|---|
| *Homo sapiens* |
| *Mus musculus* |
| *Rattus norvegicus* |
| *Danio rerio* |
| *Drosophila melanogaster* |
| *Caenorhabditis elegans* |

Genes should be identified by **gene symbol**, which is what Read10X uses by default (second column of the features file).

---

## Input Format

**Raw data:** standard Cell Ranger output for each sample, in one of these forms:
- MTX triplet: `barcodes.tsv(.gz)`, `features.tsv(.gz)` / `genes.tsv(.gz)`, `matrix.mtx(.gz)`
- HDF5: `filtered_feature_bc_matrix.h5`
- A local `filtered_feature_bc_matrix/` directory (single sample only)

**Saved object:** a Seurat v5 `.rds`, such as one previously exported from this app. Sample names are read from the `sample` metadata column, falling back to `orig.ident`.

**Signature CSV (optional, Annotation tab):** two columns, `Signature` and `Genes`, with genes comma-separated in one cell:

```
Signature,Genes
T cell,"CD3D,CD3E,CD2,TRAC"
B cell,"CD79A,MS4A1,CD19"
cDC1,"CLEC9A,XCR1,BATF3"
```

---

## Example Data

`example_data/mouse_spleen_10k_scFFPE_filtered_feature_bc_matrix.h5` is the filtered feature-barcode matrix from the **"10k Mouse Spleen FFPE Tissue Dissociated using gentleMACS Dissociator, Singleplex Sample"** dataset from 10x Genomics. It was produced with Flex and Cell Ranger 7.1.0 and contains 8,268 cells from the spleen of an adult male C57BL/6 mouse.

- **Source:** https://www.10xgenomics.com/datasets/10k-mouse-spleen-ffpe-tissue-dissociated-using-gentlemacs-dissociator-singleplex-sample-1-standard
- **License:** [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/). Dataset credit: 10x Genomics. The dataset is **not** covered by this repository's MIT license.

**How to load it:** in the sidebar, set **Species = Mouse** and **Format = 10X H5**. Leave *Multiple Datasets* unticked, upload the file as the H5 file, and click **Load Data**.

For ScType annotation of this dataset, select **Spleen** or **Immune system** as the tissue.

Because this is a single sample, it demonstrates QC, clustering, markers, annotation and subclustering. The multi-sample features (integration, split-by-sample plots, condition DE, and GSEA/ORA on condition comparisons) require two or more samples.

---

## Installation

### Requirements

- R ≥ 4.1
- Seurat ≥ 5.0
- msigdbr ≥ 10.0
- Internet access when running ScType, which downloads its scripts and marker database from GitHub
- The following R packages:

```r
# CRAN
install.packages(c("shiny", "Seurat", "ggplot2", "dplyr", "patchwork", "DT",
                   "writexl", "hdf5r", "ggrepel", "scales", "msigdbr",
                   "harmony", "HGNChelper", "openxlsx"))

# Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "fgsea", "clusterProfiler",
  "multtest", "metap",   # required by FindConservedMarkers
  "glmGamPoi"            # optional, speeds up SCTransform
))
```

`harmony` and `HGNChelper` are also installed automatically the first time Harmony or ScType runs, if they are missing.

Note: the analysis this app was developed on used harmony 1.2.4. Current CRAN versions (2.x) gave different clustering in testing, though other package versions also differed. To match the tested setup, see [VERSIONS.md](VERSIONS.md).

### Running the App

**1. Clone the repo (in terminal):**
```bash
git clone https://github.com/eabikh/sc-rnaseq-shiny-app.git
cd sc-rnaseq-shiny-app
```

**2. Launch the app (in R or RStudio):**
```r
shiny::runApp("SCGUIvs.R")
```

---

## Usage

**Quick start:** load `example_data/mouse_spleen_10k_scFFPE_filtered_feature_bc_matrix.h5` as a single H5 file with species = Mouse (see [Example Data](#example-data)).

1. **Load data** (sidebar): choose the species and format, then upload one sample, or tick *Multiple Datasets* and queue 2+ named samples. Click **Load Data**. You can also load a saved `.rds` to resume where you left off.
2. **Data Summary**: check cell counts per sample and optionally assign condition groups.
3. **Counts Matrix** (optional): browse or export the raw, normalized or scaled matrix.
4. **QC Metrics**: set the feature and % MT thresholds and click **Apply QC Filters**. Compare the Before and After plots.
5. **PCA Analysis**: use the elbow plot and PC heatmaps to choose the number of dimensions.
6. **Integration**: pick None, Harmony (set theta), CCA or SCTransform+CCA, then click **Run Integration**. Check the pre- and post-integration UMAPs.
7. **Clustering**: set Dims and Resolution and click **Run Clustering**. Optionally subcluster selected clusters (**Prep → Integrate → Cluster**, then **Merge into main object**).
8. **Cluster QC**: check QC metrics per cluster.
9. **Gene Expression**: plot genes of interest as feature, violin or dot plots.
10. **Cluster Markers**: click **Find All Markers** and review the heatmap and tables.
11. **Annotation**: run ScType for your tissue, score custom signatures, and/or apply manual labels.
12. **Comparative Analysis** (2+ samples): explore UMAPs, composition and split gene plots. Then run Conserved Markers, Pseudobulk and **Condition DE** for each cluster and sample pair of interest. Use the cached DE results in **GSEA** (and the Pathway Explorer) and **ORA**.
13. **Export**: download plots (600 DPI PNG) and tables (`.csv`) from each tab, and save the Seurat object (`.rds`) from the sidebar.

---

## Output

| Download | Format | Contents |
|---|---|---|
| Seurat Object | `.rds` | Full processed object (reloadable in the app) |
| Full Metadata | `.csv` | Per-cell metadata (QC metrics, clusters, annotations) |
| Expression Matrix | `.tsv` | Raw, normalized or scaled matrix (all cells, or one sample) |
| Cluster Statistics | `.csv` | Cell count and % per cluster |
| Subcluster Sizes | `.csv` | Overall, by sample, and by parent origin |
| Cluster / Subcluster Markers | `.csv` | Top N, all markers, or one cluster |
| ScType Results | `.csv` | Cluster → cell type, score, cell count, driver genes |
| Signature Scores | `.csv` | Average module score per signature × cluster |
| Signature Gene Contributions | `.csv` | Per-gene found / mean expression / contribution tier |
| Composition | `.csv` | Per-cluster sample composition; per-sample cell-type proportions |
| Conserved Markers | `.csv` | `FindConservedMarkers` results for a cluster |
| DE Results | `.csv` | `FindMarkers` results for the selected comparison |
| GSEA Results | `.csv` | fgsea results with leading-edge genes |
| ORA Results | `.csv` | Enriched terms for both directions |
| All plots | `.png` | Publication-ready figures (600 dpi) |

---

## Citation

If you use this app in your research, please cite:

**This app:**
Abi Khalil, E. & Lee, S. (2026). sc-rnaseq-shiny-app.
DOI to be added

**Underlying tools:**
- **Seurat v5**: Hao et al., *Nature Biotechnology* (2024)
- **CCA integration**: Stuart et al., *Cell* (2019)
- **SCTransform**: Hafemeister & Satija, *Genome Biology* (2019); Choudhary & Satija, *Genome Biology* (2022)
- **Harmony**: Korsunsky et al., *Nature Methods* (2019)
- **ScType**: Ianevski et al., *Nature Communications* (2022)
- **fgsea**: Korotkevich et al., *bioRxiv* (2021)
- **MSigDB**: Liberzon et al., *Cell Systems* (2015)
- **clusterProfiler**: Wu et al., *The Innovation* (2021)

See [VERSIONS.md](VERSIONS.md) for tested package versions and reproducibility notes.
