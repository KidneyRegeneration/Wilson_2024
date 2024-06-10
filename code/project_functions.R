
subset_fs19 <- function(x = 0.1, by = "capture", times = 1) {
  fs19 <- readr::read_rds(here::here("output/rds/FS19.rds"))
  ds <- caret::createDataPartition(fs19@meta.data[, by], times = times, p = x)
  if (times==1){
    fs19 <- fs19[, ds$Resample1]
  } else {
    fs19 <- map(ds, ~fs19[, .x])
  }
  return(fs19)
}

subset_seu <- function(seu, x = 0.1, by = "capture", times = 1) {
  
  ds <- caret::createDataPartition(seu@meta.data[, by], times = times, p = x)
  if (times==1){
    seu <- seu[, ds$Resample1]
  } else {
    seu <- map(ds, ~seu[, .x])
  }
  return(seu)
}

plotly3d <- function(seu, col) {
  plot_ly(data.frame(cell = colnames(seu),
                     dim1 = seu@reductions$umap3d@cell.embeddings[,1],
                     dim2 = seu@reductions$umap3d@cell.embeddings[,2],
                     dim3 = seu@reductions$umap3d@cell.embeddings[,3]),
          x = ~dim1, 
          y = ~dim2,
          z = ~dim3, type="scatter3d", mode = 'markers',
          marker = list(opacity = 0.7, size=2),
          
          color = seu@meta.data[,col])
}

StrDotPlot <- function(x, features, group.by = "Identity", assay = "RNA", col.min = 0,
                       col.max = 5, dot.min = 0, dot.scale = 6, scale = F, rainbow = F) {
  if (rainbow){
    c <- scale_color_gradientn(colours = rainbow(5))
  } else {
    c <- scale_colour_gradient(low = "lightgrey", high = "red")
  }
  DotPlot(object = x, features = features, col.min = col.min, col.max = col.max,
          cols = c("lightgrey", "red"), dot.scale = dot.scale,
          group.by = group.by, dot.min = dot.min, scale = scale,
          assay = assay) +
    c +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 1),
          panel.grid.major = element_line(colour = "lightgray")) +
    theme(legend.title=element_text(size=rel(0.5)))
}



'%!in%' <- function(x,y)!('%in%'(x,y))



vcol <- viridis::viridis_pal()(10)

get.psudobulk <- function(seu, column){
  pb <- data.frame(row.names = rownames(seu@assays$RNA@counts))
  for (i in 1:length(unique(seu@meta.data[,column]))){
    temp.seurat <- seu[, seu@meta.data[,column] == unique(seu@meta.data[,column])[i]]
    #temp.seurat <- subset(matorg, ident = unique(matorg$DKCC)[i])
    temp.counts <- as.data.frame(temp.seurat@assays$RNA@counts)
    temp.bulk <- data.frame(rowSums(temp.counts))
    colnames(temp.bulk) <- c(unique(as.character(seu@meta.data[, column]))[i])
    pb <- cbind(pb, temp.bulk)
  }
  pb
}

cc <- viridis::viridis(3)

pseudobulk <- function(seu, ident) {
  counts <- data.frame(row.names = rownames(seu@assays$RNA@counts))
  for (i in 1:length(unique(seu@meta.data[,ident]))){
    temp.seurat <- seu[, seu@meta.data[,ident] == unique(seu@meta.data[,ident])[i]]
    #temp.seurat <- subset(matorg, ident = unique(matorg$DKCC)[i])
    temp.counts <- as.data.frame(temp.seurat@assays$RNA@counts)
    temp.bulk <- data.frame(rowSums(temp.counts) %>% edgeR::cpm())
    colnames(temp.bulk) <- c(unique(as.character(seu@meta.data[,ident]))[i])
    counts <- cbind(counts, temp.bulk)
  }
  return(counts)
}



allcols <- read_rds("../data/rds/allcols.rds")

gcols <- ggplotColors(5)
