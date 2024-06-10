library(tidyverse)
library(Seurat)
library(patchwork)
source(here::here("code/general_r_functions.R"))
source(here::here("code/project_functions.R"))
library(destiny)
seed <- 250395

# function for generating differentiation maps

CalcDiffMap <- function(object, batch.corr = F, batch = NULL ){
  object@assays$SCT@var.features <- c(object@assays$SCT@var.features[object@assays$SCT@var.features %!in% cc.genes.updated.2019$s.genes],
                                      object@assays$SCT@var.features[object@assays$SCT@var.features %!in% cc.genes.updated.2019$g2m.genes])
  object <- object %>% SCTransform(seed.use = 250395, vars.to.regress = c("S.Score", "G2M.Score")) %>% 
    RunPCA(assay = "SCT",npcs=30, verbose=F, seed.use = 250395)
  if (batch.corr==T){
    library(simspec)
    library(Matrix)
    object <- object %>% cluster_sim_spectrum(label_tag = batch, verbose=F)
    dm.out <- destiny::DiffusionMap(object@reductions$css@cell.embeddings)
    return(dm.out)
  }
  else {
    dm.out <- destiny::DiffusionMap(object@reductions$pca@cell.embeddings)
    return(dm.out)
  }
}

# load all data

fs19 <- read_rds(here::here("data/rds/KidOrgTimecourse.rds"))


# by stage

#map(c("Stage1", "Stage2", "Stage3", "Stage4"), ~(fs19[, fs19$stage %in% .x] %>% CalcDiffMap(batch.corr = F))@eigenvectors %>% as.data.frame() %>%
#      write_csv(paste0("data/tables/DiffMap_", .x, ".csv")))

# 2-stages

(fs19[, fs19$stage %in% c("Stage2", "Stage3")] %>% CalcDiffMap(batch.corr = F))@eigenvectors %>% as.data.frame() %>%
  write_csv(paste0("data/tables/DiffMap_PC_Stages23.csv"))

(fs19[, fs19$stage %in% c("Stage3", "Stage4")] %>% CalcDiffMap(batch.corr = F))@eigenvectors %>% as.data.frame() %>%
  write_csv(paste0("data/tables/DiffMap_PC_Stages34.csv"))



