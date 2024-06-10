# write a script to run on the cluster

if (!require("hdf5r", quietly = TRUE))
  install.packages("BiocManager")
library(Seurat)
library(SeuratObject)
library(tidyverse)

# add the iPSC cells to the whole dataset

# generate a umap reduction of the entire datasets



# the data was originally merged as batches and then integrated by batch

# cell populations

ipsc.cells <- read_csv("data/tables/FS19_all_metadata.csv") %>% filter(condition == "iPSC")
ipsc <- Seurat::Read10X_h5(filename = paste0("data/filtered_feature_bc_matrix_cpt1.h5"))
ipsc.cells$rowname <- gsub("_1_1", "-1", ipsc.cells$rowname)

ipsc <- CreateSeuratObject(ipsc[,ipsc.cells$rowname])

s1 <- read_rds("data/rds/Stage1.rds")
s2 <- read_rds("data/rds/Stage2.rds")
s3 <- read_rds("data/rds/Stage3.rds")
s4 <- read_rds("data/rds/Stage4.rds")

ipsc$batch <- "Batch_1"
ipsc$ann <- "iPSC"
ipsc$stage <- "iPSC"
ipsc$ann <- "iPSC"
ipsc$age <- 0
ipsc$chir <- NA
ipsc$gf <- NA
s1@meta.data <- s1@meta.data %>% select(c("ann", "condition", "age", "chir", "gf")) %>% mutate(batch = "Batch_1", stage = "Stage1")
s2@meta.data <- s2@meta.data %>% select(c("ann", "condition", "age", "chir", "gf")) %>% mutate(batch = "Batch_1", stage = "Stage2")
s3@meta.data <- s3@meta.data %>% select(c("ann", "condition", "age", "chir", "gf")) %>% mutate(batch = "Batch_1", stage = "Stage3")
s4@meta.data <- s4@meta.data %>% select(c("ann", "condition", "age", "chir", "gf")) %>% mutate(batch = "Batch_2", stage = "Stage4")



library(simspec)

options(future.globals.maxSize = Inf)
future::nbrOfWorkers() -> workers
print(paste(workers, " workers", sep = ""))
'%!in%' <- function(x,y) !('%in%'(x,y))

fs19 <- merge(ipsc, list(s1, s2, s3, s4))

# add meta data


# perform whole-data integration
fs19 <- NormalizeData(fs19)
fs19 <- CellCycleScoring(fs19, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
fs19 <- fs19 %>% SCTransform(vars.to.regress = c("S.Score", "G2M.Score"), seed.use = 250395)
fs19@assays$SCT@var.features <- c(fs19@assays$SCT@var.features[fs19@assays$SCT@var.features %!in% cc.genes.updated.2019$s.genes],
                                  fs19@assays$SCT@var.features[fs19@assays$SCT@var.features %!in% cc.genes.updated.2019$g2m.genes])
fs19 <- fs19 %>% RunPCA(assay = "SCT",npcs=50, verbose=F, seed.use = 250395)
fs19 <- fs19 %>% cluster_sim_spectrum(label_tag = "batch", verbose=F)
fs19 <- fs19 %>% RunUMAP(reduction = "css", dims = 1:ncol(Embeddings(fs19, "css")), seed.use = 250395, n.components = 2, 
                         n.epochs = 500, verbose = F, reduction.name = "umap", reduction.key = "UMAP")
fs19 <- fs19 %>% RunUMAP(reduction = "css", dims = 1:ncol(Embeddings(fs19, "css")), seed.use = 250395, n.components = 3, 
                         n.epochs = 500, verbose = F, reduction.name = "umap3d", reduction.key = "UMAP3D")

write_rds(fs19, "data/rds/KidOrgTimecourse.rds")