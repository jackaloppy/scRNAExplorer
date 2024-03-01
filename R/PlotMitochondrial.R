#' Plot Mitochondrial Gene Percentage
#'
#' This function takes a Seurat object and plots the density of mitochondrial gene percentages stored in the object's metadata.
#' It requires the `Seurat` and `ggplot2` packages. If the mitochondrial gene percentages have not been pre-calculated and stored
#' in the object using the `PercentageFeatureSet` function with a pattern matching mitochondrial genes (e.g., "^MT-"),
#' this function will first calculate these percentages.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA sequencing data.
#' @param group The name of the metadata column in the Seurat object that identifies the sample or condition.
#'                   This will be used for coloring and filling the density plot.
#' @return A ggplot object representing the density plot of mitochondrial gene percentages.
#' @examples
#' # Assuming `seurat_obj` is your Seurat object and it has a column 'sample' in its metadata:
#' plotMitochondrial(seurat_obj, pattern = "(?i)^mt-", group = "groups")
#' @export
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom ggplot2 ggplot aes geom_density scale_x_log10 theme_classic geom_vline ggtitle
PlotMitochondrial <- function(seurat_obj, group = NULL) {

  # Calculate percent.mt if not already present
  if (!"percent.mt" %in% names(seurat_obj[[]])) {
    seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern="(?i)^mt-")
  }

  # Extract metadata
  metadata <- seurat_obj@meta.data

  # Generate the plot
  p <- ggplot(metadata, aes_string(x = "percent.mt")) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    geom_vline(xintercept = 0.2, linetype="dashed") +
    ggtitle("Mitochondrial Gene Percentage")

  # Check if grouping is specified and exists in metadata
  if (!is.null(group) && group %in% names(metadata)) {
    p <- p + aes_string(color = group, fill = group)
  } else if (!is.null(group)) {
    warning("Specified group column does not exist in metadata. Plotting without grouping.")
  }
  return(p)
}
