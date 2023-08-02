library(Seurat)
library(SeuratObject)
library(tidyverse)

df[df == "Old Value"] <- "New Value"

fs19 <- read_rds("data/rds/FS19.rds")
s1 <- read_rds("data/rds/Stage1.rds")
s2 <- read_rds("data/rds/Stage2.rds")
s3 <- read_rds("data/rds/Stage3.rds")
s4 <- read_rds("data/rds/Stage4.rds")


md.list <- list(fs19@meta.data, s1@meta.data, s2@meta.data, s3@meta.data, s4@meta.data)
# save the old metadata list

write_rds(md.list, "data/rds/OldMetaDataLists.rds")

s1old <- c("S1.Epi-like", "S1.PS-like_d3", "S1.PS-like_d4", "S1.PS-like_d5")
s1new <- c("S1.Epi-like", "S1.Primitive_Streak", "S1.Nascent_Mesoderm_D4", "S1.Nascent_Mesoderm_D5")

for (i in 1:4){
  md.list[[1]][md.list[[1]] == s1old[i]] <- s1new[i]
}

for (i in 1:4){
  md.list[[2]][md.list[[2]] == s1old[i]] <- s1new[i]
}

md.list[[2]]$Ann %>% table()


s2old <- c("S2.IM_1", "S2.IM_2", "S2.IM_3", "S2.IM_4", "S2.NMP_1", "S2.NMP_2", "S2.NMP_IM", "S2.PM_1", "S2.PM_2", "S2.PSM")
s2new <- c("S2.IM_1", "S2.IM_2", "S2.IM_3", "S2.NMP_3", "S2.NMP_1", "S2.NMP_2", "S2.NMP_IM", "S2.PM_1", "S2.PM_2", "S2.PSM")

for (i in 1:4){
  md.list[[1]][md.list[[1]] == s1old[i]] <- s1new[i]
}

for (i in 1:4){
  md.list[[2]][md.list[[2]] == s1old[i]] <- s1new[i]
}

md.list[[2]]$Ann %>% table()


rnav <- map(1:13, ~read_csv(paste0("data/tables/RNAv_SW_CPT_", .x, ".csv")))

rnav[[2]] %>% rownames() %>% head()
rnav[[2]]$SCT_snn_res.0.3 %>% table
# need to remove all cells in clusters 10 and 11 for rnav objects 2-4

s2.rnav.md <- bind_rows(rnav[[2]], list(rnav[[3]], rnav[[4]]))
s2.rnav.md <- s2.rnav.md %>% filter(SCT_snn_res.0.3 < 10)
s2.rnav.md <- as.data.frame(s2.rnav.md)
rownames(s2.rnav.md) <- rownames(s2@meta.data)

s2@reductions$umap.all <- s2@reductions$umap
s2@reductions$umap <- CreateDimReducObject(as.matrix(s2.rnav.md[, c("UMAP_1", "UMAP_2")]), key = "UMAP_", assay = "SCT")

DimPlot(s2, group.by = "Ann", reduction = "umap")

# add a higher resolution clustering to the metadata

s2$SCT_snn_res.0.9 <- s2.rnav.md$SCT_snn_res.0.9
DimPlot(s2, group.by = "SCT_snn_res.0.9", reduction = "umap", label = T)

s2$SCT_snn_res.0.6 <- s2.rnav.md$SCT_snn_res.0.6
DimPlot(s2, group.by = "SCT_snn_res.0.6", reduction = "umap", label = T)

write_rds(s2, "data/rds/Stage2.rds")

