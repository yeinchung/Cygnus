## Convert excel output from AIVIA to Cygnus file format
## Create Cygnus object

getCygnus <- function(data.dir){
  read.csv(data.dir)
  # Insert checkpoints to ensure consistent file format
}

setClass("CygnusData",
         slots = list(
           exp_matrix = "data.frame",
           metadata = "list"
         ))

createCygnus <- function(file, markers_col = colnames(file), metadata_col = NA) {
  obj <- new("CygnusData",
             exp_matrix = file[markers_col],
             metadata = file[metadata_col])
  return(obj)
}

setMethod("summary", signature(object = "CygnusData"),
          function(object) {
            cat("Summary of CygnusData:\n")
            cat("Data:\n")
            print(head(object@exp_matrix))
            cat("Metadata:\n")
            print(head(object@metadata))
          })

setMethod("show",
          signature = "CygnusData",
          definition = function(object){
            cat("CygnusData Object: \n")
            cat("Number of EVs: ", dim(object@exp_matrix)[1], "\n")
            cat("Number of markers: ", dim(object@exp_matrix)[2], "\n")
            cat("Present Metadata: ", names(object@metadata))
          })

# Get histograms for marker distributions
plotDistribution <- function(obj) {
  num_cols <- ncol(obj@exp_matrix)
  num_rows <- ceiling(num_cols / 3)

  par(mfrow = c(num_rows, 3))
  par(mar = c(2, 2, 2, 1))

  for (i in 1:num_cols) {
    hist(obj@exp_matrix[, i],  breaks = 1000,
         main = paste(colnames(obj@exp_matrix)[i]),
         xlab = "", ylab = "")  # Clear x and y labels
  }

  par(mfrow = c(1, 1))
}


# Get Heatmap
getHeatmap <- function(obj){ #*require pheatmap
  pheatmap()
}

# Get Heatmap for Average Expressions
# required: dplyr and pheatmap
plotAvgHeatmap <- function(data, group_column,
                          clustering_distance = "euclidean",
                          colors = rev(colorRampPalette(RColorBrewer::brewer.pal(n = 11, "RdYlBu"))(100)),
                          fontsize = 8,
                          scale = 'column') {


  expression_matrix <- data@exp_matrix
  group_metadata <- data@metadata[[group_column]]

  combined_data <- cbind(expression_matrix, group_metadata)

  avg_marker_expressions <- combined_data %>%
    group_by(group_metadata) %>%
    summarise_all(mean, na.rm = TRUE) %>%
    select(-group_metadata) %>%
    as.matrix()

  pheatmap(avg_marker_expressions,
           clustering_distance_cols = clustering_distance,
           cluster_rows = F,
           color = colors,
           main = paste("Average Marker Expressions by", group_column),
           fontsize = fontsize,
           scale = scale,
           labels_row = unique(data@metadata[[group_column]]))
}


addLayer <- function(obj, layer_name){
  ob
}

# required: plotly
runPCA <- function(obj, scale){
  pca_output <- prcomp(obj@exp_matrix, scale=scale)

  obj@metadata[['PCA_coord']] <- pca_output # figure this out
  return(obj)

}

viewPCA <- function(obj, group_column){
  ##
  y_factor <- factor(obj@metadata[[group_column]])
  coord <- as.data.frame(obj@metadata$PCA_coord[['x']])
  fig <- plot_ly(data = coord,
                 x = coord$PC1, y = coord$PC2, #z = coord$PC3,
                 color = ~y_factor,
                 colors = c("lightseagreen",
                            "gray50",
                            "darkgreen",
                            "red4",
                            "red",
                            "turquoise4",
                            "black",
                            "yellow4",
                            "royalblue1",
                            "lightcyan3",
                            "peachpuff3",
                            "khaki3",
                            "gray20",
                            "orange2",
                            "royalblue4",
                            "yellow3",
                            "gray80",
                            "darkorchid1",
                            "lawngreen",
                            "plum2",
                            "darkmagenta"),
                 type = "scatter",
                 mode = "markers",
                 marker = list(size = 2, width=2))
  fig

}

getElbowPlot <- function(obj, scale=F){
  pca_output <- prcomp(obj@exp_matrix, scale=scale)

  variance_explained <- pca_output$sdev^2 / sum(pca_output$sdev^2)
  cumulative_variance_explained <- cumsum(variance_explained)

  elbow_plot <- data.frame(PC = 1:length(cumulative_variance_explained),
                           Variance_Explained = variance_explained,
                           Cumulative_Variance_Explained = cumulative_variance_explained)

  library(ggplot2)
  ggplot(elbow_plot, aes(x = PC, y = Cumulative_Variance_Explained)) +
    geom_point() +
    geom_line() +
    xlab("Principal Component") +
    ylab("Cumulative Proportion of Variance Explained") +
    ggtitle("Elbow Plot")}



# required: rTSNE
runTSNE <- function(obj, usePCs=F, PCs=5){
  if(usePCs){
    tsne_output <- Rtsne(obj@metadata$PCA_coord[['x']][,PCs], dims = 3,  verbose = TRUE, max_iter = 25, pca=F)
    obj@metadata[['TSNE_coord']] <- tsne_output # figure this out
  }else{
    tsne_output <- Rtsne(obj@exp_matrix, dims = 3, verbose=T, pca=F)
    obj@metadata[['TSNE_coord']] <- tsne_output
  }

  return(obj)
} # shape changes everytime you cant have that!

viewtSNE <- function(obj, group_column){
  ##
  y_factor <- factor(obj@metadata[[group_column]])
  coord <- as.data.frame(obj@metadata$TSNE_coord[['Y']])
  fig <- plot_ly(data = coord,
                 x = coord$V1, y = coord$V2,#, z = coord$V3,
                 color = ~y_factor,
                 colors = c("lightseagreen",
                            "gray50",
                            "darkgreen",
                            "red4",
                            "red",
                            "turquoise4",
                            "black",
                            "yellow4",
                            "royalblue1",
                            "lightcyan3",
                            "peachpuff3",
                            "khaki3",
                            "gray20",
                            "orange2",
                            "royalblue4",
                            "yellow3",
                            "gray80",
                            "darkorchid1",
                            "lawngreen",
                            "plum2",
                            "darkmagenta"),
                 type = "scatter",
                 mode = "markers",
                 marker = list(size = 2, width=2))
  fig

}

viewtSNE_3D <- function(obj, group_column){
  ##
  y_factor <- factor(obj@metadata[[group_column]])
  coord <- as.data.frame(obj@metadata$TSNE_coord[['Y']])
  fig <- plot_ly(data = coord,
                 x = coord$V1, y = coord$V2, z = coord$V3,
                 color = ~y_factor,
                 colors = c("lightseagreen",
                            "gray50",
                            "darkgreen",
                            "red4",
                            "red",
                            "turquoise4",
                            "black",
                            "yellow4",
                            "royalblue1",
                            "lightcyan3",
                            "peachpuff3",
                            "khaki3",
                            "gray20",
                            "orange2",
                            "royalblue4",
                            "yellow3",
                            "gray80",
                            "darkorchid1",
                            "lawngreen",
                            "plum2",
                            "darkmagenta"),
                 type = "scatter3d",
                 mode = "markers",
                 marker = list(size = 2, width=2))
  fig

}

# now run UMAP stuff
library(umap) # tried this but it didnt work
library(Seurat)
# its better to have Seurat here


runUMAP_Seurat<- function(obj, usePCs = FALSE, numPCs = 10, n_neighbors = 15, min_dist = 0.3, n_epochs = 200) {
  data_to_embed <- if (usePCs) {
    obj@metadata$PCA_coord[['x']][, 1:numPCs]
  } else {
    obj@exp_matrix
  }

  seu <- CreateSeuratObject(t(data_to_embed))
  seu <- Seurat::RunUMAP(seu, di)
  return(obj)
}

viewUMAP <- function(obj, group_column = "clustering") {
  coord <- as.data.frame(obj@metadata$UMAP_coord)
  group_labels <- factor(obj@metadata[[group_column]])

  fig <- plot_ly(data = coord,
                 x = coord$V1, y = coord$V2,
                 color = ~group_labels,
                 colors = c("lightseagreen", "gray50", "darkgreen", "red4", "red"),
                 type = "scatter", mode = "markers",
                 marker = list(size = 2, width = 2))
  fig
}

runUMAP_3D <- function(obj, usePCs = FALSE, numPCs = 10, n_neighbors = 15, min_dist = 0.3, n_epochs = 200) {
  data_to_embed <- if (usePCs) {
    obj@metadata$PCA_coord[['x']][, 1:numPCs]
  } else {
    obj@exp_matrix
  }

  umap_output <- umap::umap(data_to_embed,
                            n_neighbors = n_neighbors,
                            min_dist = min_dist,
                            n_epochs = n_epochs,
                            n_components = 3,   # Request 3D embedding
                            metric = "euclidean")

  obj@metadata[['UMAP_coord_3D']] <- umap_output$layout
  return(obj)
}

viewUMAP_3D <- function(obj, group_column = "clustering") {
  coord <- as.data.frame(obj@metadata$UMAP_coord_3D)
  group_labels <- factor(obj@metadata[[group_column]])

  fig <- plot_ly(data = coord,
                 x = coord$V1, y = coord$V2, z = coord$V3,
                 color = ~group_labels,
                 colors = c("lightseagreen", "gray50", "darkgreen", "red4", "red"),
                 type = "scatter3d", mode = "markers",
                 marker = list(size = 2, width = 2))
  fig
}

runClustering <- function(obj, k=10){ # please modify this
  pbmc <- CreateSeuratObject(t(obj@exp_matrix))
  pbmc <- NormalizeData(pbmc)
  pbmc <- ScaleData(object = pbmc)
  pbmc <- RunPCA(object = pbmc, features = rownames(pbmc))
  pbmc <- FindNeighbors(object = pbmc, dims = 1:10)
  pbmc <- FindClusters(object = pbmc)
  pbmc <- RunUMAP(object = pbmc, dims = 1:10)
}

