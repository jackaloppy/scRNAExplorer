# scRNAExplorer

`scRNAExplorer` is an R package designed to streamline the analysis of single-cell RNA sequencing (scRNA-seq) data. It offers an integrated workflow for loading 10x Genomics Cell Ranger count output into Seurat objects, pre-processing data, and generating various quality control (QC) plots. This package aims to simplify the initial steps of scRNA-seq data analysis, allowing researchers to focus on downstream analyses and discoveries.

## Installation

You can install `scRNAExplorer` from GitHub with:

```r
# Install the devtools package if you haven't already
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install scRNAExplorer from GitHub
devtools::install_github("jackaloppy/scRNAExplorer")
``` 

## Usage

To start using `scRNAExplorer`, load it as you would with any other R package:

```r 
library(scRNAExplorer)
```

### Using `Loadfrom10X` to Load and Preprocess Data

`scRNAExplorer` simplifies the process of loading scRNA-seq data and preprocessing it.. Below is an example of how to use the `Loadfrom10X` function to load data from a 10X Cell Ranger output, optionally use [SoupX](https://github.com/constantAmateur/SoupX) to estimate and remove cell-free mRNA contamination, and prepare it for downstream analysis.

```r
# Assuming you have 10X Genomics data in the specified directory
dataDir <- "path/to/your/cellranger/outs/folder"

# Load the data, preprocess, and calculate QC metrics
seurat_obj <- Loadfrom10X(dataDir = dataDir,
                          project = "ExampleProject",
                          minCells = 3,
                          minFeatures = 200,
                          useSoupX = TRUE)

```
This example demonstrates how to load data, specifying the directory where your Cell Ranger outputs are stored, and how to set project-specific parameters. Adjust the dataDir, project, minCells, minFeatures, and useSoupX parameters as needed for your analysis.

### Merge multiple `seurat_obj` for Comparative Analysis (Optional)

In scenarios where your study involves multiple scRNA-seq samples that you wish to compare, `scRNAExplorer` facilitates merging these samples into a single `Seurat` object. This is particularly useful for comparative analysis across different conditions or time points. After loading each sample with the `Loadfrom10X` function, you can annotate each sample with a group identifier to distinguish them:

```r
# Annotating samples
sample1$condition <- "control"
sample2$condition <- "control"
sample3$condition <- "treatment"
sample4$condition <- "treatment"
```
Then, you can merge these annotated samples into a single Seurat object:
```r
# Merging samples
merged_obj <- merge(x = sample1, y = list(sample2, sample3, sample4), 
                    add.cell.ids = c("s1", "s2", "s3", "s4"))
```

### Plotting QC Plots

It's important to note that if you have merged multiple samples into a single `Seurat` object for comparative analysis, you can specify the `group="group"` argument in each plotting function. This allows you to visualize the plots by the defined groups (e.g., control vs. treatment), giving your ability to compare different conditions or treatments within your dataset.

Here are examples of how to use plotting functions to generate QC plots from a `Seurat` object.

```r
# Plot the number of cells per sample, differentiated by group
PlotNCells(seurat_obj, group = "condition")

# Plot the density of UMI counts
PlotUMIDensity(seurat_obj)

# Plot the density of gene counts
PlotGeneDensity(seurat_obj)

# Plot the density of novelty scores
PlotNovelty(seurat_obj)

# Plot UMI counts against the number of genes
PlotUMIvsGene(seurat_obj)

# Plot the percentage of mitochondrial genes
PlotMitochondrial(seurat_obj)

# Plot the percentage of ribosomal genes
PlotRibosomal(seurat_obj)

# Plot the percentage of hemoglobin genes
PlotHemoglobin(seurat_obj)
```

### Analysis Workflows

`scRNAExplorer` provides  workflows for deeper analysis of scRNA-seq data, including principal component analysis (PCA) and clustering along with UMAP visualization. These workflows streamline the process of identifying and visualizing underlying patterns in your data.

The `RunPCAandElbowPlot` function simplifies the PCA workflow. It normalizes the data, identifies variable features, scales the data, runs PCA, and generates an elbow plot to help determine the number of principal components to use in further analyses (i.e. it provides a reference for the `dims` used in `RunClusteringandUMAP` function.)

```r
# Run PCA and generate an elbow plot
RunPCAandElbowPlot(seurat_obj)
```

Clustering and Visualizing with UMAP

The `RunClustersandUMAP` function takes your Seurat object through a clustering workflow and visualizes the results using UMAP. It allows you to specify the number of dimensions for neighbors and UMAP, as well as the resolution for clustering, offering flexibility in how you interpret the single-cell data.

```r
# Perform clustering and visualize with UMAP
RunClustersandUMAP(seurat_obj, dims = 20, resolution = 0.5)
```

Specify `dims` and `resolution` according to your dataset and analysis needs. This function is key for uncovering and visualizing the cellular heterogeneity within your samples.


```r
# Visualize the UMAP
DimPlot(pbmc, reduction = "umap")

# Visualize the UMAP by condition
DimPlot(pbmc, reduction = "umap", group.by = "condition")
```

