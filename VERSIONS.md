# Versions and reproducibility notes

## Developed and tested on

- macOS, Apple Silicon (aarch64-apple-darwin20), R 4.4.2
- RNG: Mersenne-Twister / Inversion / Rejection
- shiny 1.14.0, Seurat 5.4.0, SeuratObject 5.3.0, harmony 1.2.4, uwot 0.2.4, irlba 2.3.7, Matrix 1.7.5, hdf5r 1.3.12, fgsea 1.32.4, msigdbr 26.1.1, clusterProfiler 4.14.6, ggplot2 4.0.3, dplyr 1.2.1, patchwork 1.3.2

Also run once on Windows 11 with R 4.3.1, shiny 1.7.5.1, Seurat 5.2.1, SeuratObject 5.0.2, harmony 2.0.5, Matrix 1.6.5, ggplot2 3.5.2. The interface only loaded after the duplicate `de_samp_ui` fix (see below), and the same parameters gave a different cluster count than on the Mac setup.

## Seeds

The app calls `set.seed(42)` before `RunHarmony`/`IntegrateLayers`, `FindNeighbors`, `FindClusters` (`random.seed = 42`) and `RunUMAP` (`seed.use = 42`). Seeds make repeat runs reproducible on one setup; they do not guarantee identical results across different package versions or operating systems.

## Known notes

- Output IDs must be unique. `de_samp_ui` was previously defined in two tabs, which broke the UI on older Shiny versions (conditional panels did not hide and server-rendered panels did not appear). Fixed by removing the duplicate from the Conserved Markers tab.
- harmony 2.x differed from 1.2.4 in testing (cause not verified: other package versions also differed). Use harmony 1.2.4 to match the tested setup.
- Plot sizing and legends were tuned on ggplot2 4.0.3; older ggplot2 may render slightly differently.
- DE, GSEA and ORA results are session-only (not stored in a saved Seurat object); export the CSV tables to keep them.
