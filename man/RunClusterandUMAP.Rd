% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunClusterandUMAP.R
\name{RunClusterandUMAP}
\alias{RunClusterandUMAP}
\title{Cluster and UMAP Seurat Analysis Workflow}
\usage{
RunClusterandUMAP(seuratObject, dims = 10, resolution = 0.5)
}
\arguments{
\item{seuratObject}{A Seurat object.}

\item{dims}{Vector of dimensions to use.}

\item{resolution}{Resolution parameter for FindClusters.}
}
\value{
A Seurat object with PCA, clustering, and UMAP results.
}
\description{
This function takes a Seurat object, performs clustering and UMAP visualization.
The number of dimensions for neighbors, and UMAP,
and the resolution for clustering can be specified.
}
