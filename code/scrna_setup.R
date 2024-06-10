if (!require("hdf5r", quietly = TRUE))
  install.packages("BiocManager")
library(Seurat)
library(SeuratObject)
library(tidyverse)

cpt <- paste0("cpt", 1:13)
seurats <- map(cpt, ~Seurat::Read10X_h5(filename = paste0("data/filtered_feature_bc_matrix_", .x, ".h5")))
seurats <- map(seurats, ~CreateSeuratObject(.x, project = "FS19"))

# temporary save

write_rds(seurats, "data/rds/seurat_list.rds")

# read object back in if needed
seurats <- read_rds("data/rds/seurat_list.rds")

# merge all objects

all <- merge(seurats[[1]], seurats[2:13])

# load meta data

stage.md <- read_rds("data/Stages_metadata.rds")

# the data was originally merged as batches and then integrated by batch

batch <- list(batch1 = merge(seurats[[1]], seurats[2:7]),
              batch2 = merge(seurats[[8]], seurats[9:13]))

batch$batch1$batch <- "Batch_1"
batch$batch2$batch <- "Batch_2"

library(simspec)

options(future.globals.maxSize = Inf)
future::nbrOfWorkers() -> workers
print(paste(workers, " workers", sep = ""))
'%!in%' <- function(x,y) !('%in%'(x,y))

fs19 <- merge(batch$batch1, batch$batch2)

# add meta data

fs19$cell <- gsub("-1_", "_", colnames(fs19))
rownames(stage.md$s2) <- paste0(rownames(stage.md$s2), "_1")
rownames(stage.md$s2) <- gsub("_3_1", "_4_1", rownames(stage.md$s2))
rownames(stage.md$s2) <- gsub("_2_1", "_3_1", rownames(stage.md$s2))
rownames(stage.md$s2) <- gsub("_1_1", "_2_1", rownames(stage.md$s2))
map(stage.md, ~(.x %>% rownames()) %in% fs19$cell %>% table())

fs19 <- fs19[, fs19$cell %in% c(rownames(stage.md$s1), rownames(stage.md$s2), rownames(stage.md$s3), rownames(stage.md$s4))]


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

write_rds(fs19, "data/rds/FS19.rds")



stages <- list(s1 = fs19[, fs19$cell %in% rownames(stage.md$s1)],
               s2 = fs19[, fs19$cell %in% rownames(stage.md$s2)],
               s3 = fs19[, fs19$cell %in% rownames(stage.md$s3)],
               s4 = fs19[, fs19$cell %in% rownames(stage.md$s4)])

md <- map(1:4, ~left_join((stages[[.x]]@meta.data), (stage.md[[.x]] %>% mutate(cell = rownames(stage.md[[.x]]))), by = "cell"))

stages$s1@meta.data <- md[[1]] %>% column_to_rownames("cell")
stages$s2@meta.data <- md[[2]] %>% column_to_rownames("cell")
stages$s3@meta.data <- md[[3]] %>% column_to_rownames("cell")
stages$s4@meta.data <- md[[4]] %>% column_to_rownames("cell")


# clean up the metadata


stage.md$s4.fix <- (fs19[, fs19$Stage=="Stage_4"]@meta.data %>% rownames_to_column()) %>% left_join((stage.md$s4) %>% rownames_to_column(), by = "rowname") %>% column_to_rownames()
stage.md$s4.fix <- stage.md$s4.fix[,1:31]
colnames(stage.md$s4.fix) <- gsub(".x", "", colnames(stage.md$s4.fix))
colkeep <- (intersect(colnames(stage.md$s1), colnames(stage.md$s4.fix)) %>% intersect(colnames(stage.md$s2)) %>% intersect(colnames(stage.md$s3)))

stage.md <- map(stage.md, ~.x %>% select(colkeep))

stage.md.all <- bind_rows(stage.md[[1]], stage.md[[2]]) %>% bind_rows(stage.md[[3]]) %>% bind_rows(stage.md[[5]])

rownames(stage.md.all) <- gsub("A_", "A-1_", rownames(stage.md.all))
rownames(stage.md.all) <- gsub("T_", "T-1_", rownames(stage.md.all))
rownames(stage.md.all) <- gsub("C_", "C-1_", rownames(stage.md.all))
rownames(stage.md.all) <- gsub("G_", "G-1_", rownames(stage.md.all))

fs19@meta.data <- stage.md.all

write_rds(fs19, "data/rds/FS19.rds")

#map(1:4, ~fs19[, fs19$Stage == paste0("Stage_", .x)] %>% write_rds(paste0("data/rds/Stage", .x, ".rds")))
# load the DR from other file.


s1 <- read_rds("data/rds/Stage1.rds")
StrDotPlot(s1, group.by = "Ann", features = c("ZFP42", "SOX2", "NANOG", "IFITM3", "PRDM14", "SYCP3", "DAZL", "PRDM1", "STRA8", "NANOS3", "TCL1A", "TFAP2C", "PDPN"))

