#' Plot UMI Counts vs. Number of Genes
#'
#' This function takes a Seurat object and generates a scatter plot comparing UMI counts (`nCount_RNA`)
#' to the number of detected genes (`nFeature_RNA`), with points faceted by a specified grouping variable.
#' It leverages the `ggplot2` package for plotting. This visualization is useful for assessing the quality
#' and distribution of single-cell RNA-seq data across different groups.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA sequencing data with pre-calculated
#'                   `nCount_RNA` and `nFeature_RNA` in its metadata.
#' @param group The name of the metadata column in the Seurat object that identifies the grouping
#'                  variable for faceting the plot. This allows visualization of data distribution across different groups.
#' @return A ggplot object representing the scatter plot of UMI counts vs. number of genes, faceted by the specified group.
#' @examples
#' # Assuming `seurat_obj` is your Seurat object and it has a column 'group' in its metadata for faceting:
#' plotUMIvsGenes(seurat_obj, group_col = "group")
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_x_log10 scale_y_log10 theme_classic facet_wrap ggtitle
PlotUMIvsGenes <- function(seurat_obj, group=NULL) {

  # Extract metadata
  metadata <- seurat_obj@meta.data

  # Check for group or use a default column
  facet_column <- if (!is.null(group) && group %in% names(metadata)) {
    group
  } else if (!is.null(group) && 'orig.ident' %in% names(metadata)) {
    warning(paste0(group," is not in metadata. Use orig.ident instead."))
    'orig.ident'
  } else {
    NULL
  }


  # Generate the plot
  p <- ggplot(metadata, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    xlab("UMI Counts") +
    ylab("Number of Genes")+
    theme_classic() +
    ggtitle("UMIs vs. Genes")

  # Add faceting if a valid column name is determined
  if (!is.null(facet_column)) {
    p <- p + facet_wrap(as.formula(paste0("~", facet_column)))
  }

  return(p)
}
