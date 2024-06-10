

#' Read in and add some base information to 10x single cell data
#'
#' @param samples A list of sample names
#' @param project.name Name of the project, will become the orig.ident slot
#' @param data.dir Directory the sample folder containing the raw data is located in
#' @param output.dir Directory to export eh seurat object to
#' @param project.dir Root Project directory
#'
#' @return List of all seurat objects
#' @export SW_CPT_n_seurat.rds is a seurat object where n is the relevent capture
#'
#' @examples
#' 
read_10x_add_base_info <- function(samples = samples, 
                                   output.dir = "output/seurat/", 
                                   project.name = project.name,
                                   data.dir = data.dir,
                                   project.dir = project.dir){
  
  if (dir.exists(paths = paste0(project.dir, output.dir)) == F) {
    dir.create(paste0(project.dir, output.dir))
  }
  seurat.list <- list()
  # loop through all samples
  for (i in 1:length(samples)) {
    temp.data <- Read10X(data.dir = paste0(data.dir, samples[i], "/outs/filtered_feature_bc_matrix/"))
    colnames(temp.data) <- stringr::str_remove(colnames(temp.data), "-1")
    temp.seurat <- CreateSeuratObject(temp.data, project = project.name)
    temp.seurat$sample <- paste0(samples[i])
    temp.seurat <- NormalizeData(temp.seurat, verbose = F)
    temp.seurat <- PercentageFeatureSet(temp.seurat, pattern = "^MT-", col.name = "percent.mt")
    temp.seurat <- CellCycleScoring(temp.seurat, s.features = cc.genes.updated.2019$s.genes,
                               g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
    write_rds(temp.seurat, paste0(project.dir, output.dir, samples[i], "_seurat.rds"))
    seurat.list[[paste0(samples[i])]] <- temp.seurat
  }
  return(seurat.list)
}

