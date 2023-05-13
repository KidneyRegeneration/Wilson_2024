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

write_rds(all, "data/rds/FS19.rds")

all <- read_rds("data/rds/FS19.rds")

# set up metadata

all.md <- map(seurats, ~.x@meta.data)
all.md <- map(1:13, ~all.md[[.x]] %>% mutate(capture = paste0("CPT", .x)))
all.md <- map_dfr(all.md, ~.x)
rownames(all.md) <- colnames(all)
all@meta.data <- all.md

# stage 1
s1 <- all[, all$capture == "CPT1"]

# load meta data and add to seurat
s1.md <- read_csv("data/tables/Stage1_Annotation.csv")
colnames(s1)[1]
s1.md$rowname[1]
s1.md$rowname <- gsub("_1_1", replacement = "-1_1", s1.md$rowname)
s1 <- s1[, colnames(s1)%in%s1.md$rowname]
s1@meta.data <- bind_cols(s1@meta.data, s1.md)
write_rds(s1, file = "data/rds/Stage1.rds")


# stage 2
s2 <- all[, all$capture %in% c("CPT2", "CPT3", "CPT4")]

# load meta data and add to seurat
s2.md <- read_csv("data/tables/Stage2_Annotation_df.csv")
colnames(s2)[c(1, 20001, 40001)]
s2.md$rowname[c(1, 20001, 40001)]
s2.md$rowname <- gsub("_3", replacement = "-1_4", s2.md$rowname)
s2.md$rowname <- gsub("_2", replacement = "-1_3", s2.md$rowname)
s2.md$rowname <- gsub("_1", replacement = "-1_2", s2.md$rowname)
s2 <- s2[, colnames(s2)%in%s2.md$rowname]
s2@meta.data <- bind_cols(s2@meta.data, s2.md)
write_rds(s2, file = "data/rds/Stage2.rds")


# stage 3
s3 <- all[, all$capture %in% c("CPT5", "CPT6", "CPT7")]

# load meta data and add to seurat

s3.md <- read_csv("data/tables/Stage3_Assigned_Annotation.csv")

colnames(s3)[c(1, 20001, 40001)]
s3.md$rowname %>% tail()

s3.md$rowname <- gsub("_5_1", replacement = "-1_5", s3.md$rowname)
s3.md$rowname <- gsub("_6_1", replacement = "-1_6", s3.md$rowname)
s3.md$rowname <- gsub("_7_1", replacement = "-1_7", s3.md$rowname)

s3.md <- left_join(s3@meta.data %>% rownames_to_column(), s3.md, by = "rowname") %>% filter(!is.na(stage))
s3 <- s3[, colnames(s3)%in%s3.md$rowname]
s3@meta.data <- s3.md %>% column_to_rownames()
write_rds(s3, file = "data/rds/Stage3.rds")

# stage 4
s4 <- all[, all$capture %in% paste0("CPT", 8:13)]

# load meta data and add to seurat

s4.md <- read_csv("data/tables/Stage4_Assigned_Annotation.csv")

colnames(s4)[seq(1, 120000, 20000)]
s4.md$rowname[seq(1, nrow(s4.md), nrow(s4.md)/6)]

s4.md$rowname <- gsub("_6_2", replacement = "-1_13", s4.md$rowname)
s4.md$rowname <- gsub("_5_2", replacement = "-1_12", s4.md$rowname)
s4.md$rowname <- gsub("_4_2", replacement = "-1_11", s4.md$rowname)
s4.md$rowname <- gsub("_3_2", replacement = "-1_10", s4.md$rowname)
s4.md$rowname <- gsub("_2_2", replacement = "-1_9", s4.md$rowname)
s4.md$rowname <- gsub("_1_2", replacement = "-1_8", s4.md$rowname)

s4.md <- left_join(s4@meta.data %>% rownames_to_column(), s4.md, by = "rowname") %>% filter(!is.na(stage))
s4 <- s4[, colnames(s4)%in%s4.md$rowname]
s4@meta.data <- s4.md %>% column_to_rownames()
write_rds(s4, file = "data/rds/Stage4.rds")

# combine as stages
s3$ann.broad <- s3$ann.long
s4$ann.broad <- s4$ann.long
s2$stage <- "Stage_2"

all <- merge(s1, list(s2, s3, s4))

write_rds(all, "data/rds/FS19.rds")


