#' Plot Number of Cells per Sample from Seurat Object Metadata
#'
#' @param seuratObject A Seurat object with a 'Group' column in the metadata.
#' @param sample Optional; the name of the metadata column to use for grouping.
#' @return A ggplot object showing the number of cells per sample.
#' @examples
#' PlotNCellsPerSample(seuratObject, sample = "Group")
#' @import ggplot2
#' @import dplyr
#' @export
PlotNCells <- function(seuratObject, group = NULL) {
  # Extract metadata
  metadata <- seuratObject@meta.data

  # Define the aesthetic mappings
  if (!is.null(group) && group %in% names(metadata)) {
    # If group column is specified and exists
    aes_mappings <- aes_string(x = group, fill = group)
    x_label <- group
  } else {
    # Default to 'orig.ident' if group is not specified
    aes_mappings <- aes(x = orig.ident, fill = orig.ident)
    x_label <- "orig.ident"
  }

  # Create the plot
  p <- ggplot(metadata, aes_mappings) +
    geom_bar() +
    theme_classic() +
    xlab(x_label) +
    ylab("Count") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("Number of Cells per group")

  return(p)
}
