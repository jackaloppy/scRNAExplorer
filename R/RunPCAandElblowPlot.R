#' Run PCA Workflow and Generate Elbow Plot for Seurat Object
#'
#' This function takes a Seurat object, normalizes the data, identifies variable features,
#' scales the data, runs PCA, and generates an elbow plot.
#'
#' @param seuratObject A Seurat object.
#' @return Generates an elbow plot showing the variance explained by each principal component.
#' @export
RunPCAandElbowPlot <- function(seuratObject) {
  # Normalize the data
  seuratObject <- NormalizeData(seuratObject)

  # Find variable features
  seuratObject <- FindVariableFeatures(seuratObject)

  # Scale the data
  seuratObject <- ScaleData(seuratObject, features = rownames(seuratObject))

  # Run PCA
  seuratObject <- RunPCA(seuratObject, features = VariableFeatures(seuratObject))

  # Generate Elbow Plot
  ElbowPlot(seuratObject)

  return(seuratObject)
}
