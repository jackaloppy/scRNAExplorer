#' Load Data from 10X Cellranger Output, Optionally Run SoupX, and Create Modified Seurat Object
#'
#' This function reads 10X Cellranger count output from a specified directory
#' and initializes a Seurat object with the provided parameters. Optionally,
#' it can use SoupX to estimate and remove cell-free mRNA contamination. After
#' creating the Seurat object with minimum cells and genes thresholds,
#' it calculates the percentages of mitochondrial genes (percent.mt),
#' ribosomal genes (percent.ribo), and hemoglobin genes
#' (percent.hb) and adds them to the metadata of the Seurat object.
#' You might also initialize the filtering
#'
#'
#' @param dataDir Top level CellRanger count output directory.
#' @param project Name of the project, used to tag the Seurat object.
#' @param minCells Minimum number of cells for a gene to be included in the Seurat object.
#' @param minFeatures Minimum number of features for a cell to be included in the Seurat object.
#' @param useSoupX Logical, whether to use SoupX to adjust for contamination.
#'
#' @return A Seurat object.
#' @examples
#' pbmc <- Loadfrom10X(dataDir = "path/to/your/cellranger/outs/folder",
#'                     project = "ExampleProject",
#'                     min.cells = 3,
#'                     min.features = 200,
#'                     useSoupX = TRUE)
#' @export
#' @importFrom Seurat CreateSeuratObject Read10X PercentageFeatureSet
#' @importFrom SoupX load10X autoEstCont adjustCounts
Loadfrom10X <- function(dataDir, project = "defaultProject", minCells = 3, minFeatures = 200, useSoupX = FALSE) {
  if (useSoupX) {
    if (!requireNamespace("SoupX", quietly = TRUE)) stop("SoupX package is not installed.")
    # Load data using SoupX
    data <- SoupX::load10X(dataDir)
    data <- SoupX::autoEstCont(data)
    data <- SoupX::adjustCounts(data)
  } else {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("SoupX package is not installed.")
    # Adjust path for Seurat's Read10X
    dataDirAdjusted = paste0(dataDir, "/filtered_feature_bc_matrix")
    data <- Seurat::Read10X(data.dir = dataDirAdjusted)
  }

  seurat_obj <- Seurat::CreateSeuratObject(counts = data,
                                           project = project,
                                           min.cells = minCells,
                                           min.features = minFeatures)

  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern="(?i)^mt-")
  seurat_obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seurat_obj,pattern="^(?i)Rp[sl]")
  seurat_obj[["percent.hb"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern="(?i)Hb[^(p)]")

  return(seurat_obj)
}
