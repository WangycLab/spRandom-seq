library(tidyverse)
library(Seurat)
library(loomR)

visium_coord <- read.table("./Data/visium-v4_coordinates.txt")
colnames(visium_coord) <- c("Barcode", "X", "Y")
visium_coord <- visium_coord %>%
  arrange(X,Y)

cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#8B2DB2", "#0072B2", "#D55E00", "#CC79A7",
          "#CD534C", "#FFBF7F", "#FF5500", "#E51932",
          "#FF99BF", "#654CFF", "#CCBFFF", "#19B2FF",
          "#A5EDFF", "#CC5151", "#260F99", "#85B22C")

############## Brain ####################
dataObj_brain <- readRDS("./Input/10X_ST_Mouse_Brain_add_res.rds")

dataObj_brain@meta.data <- dataObj_brain@meta.data %>%
  rownames_to_column(var = "Barcode") %>%
  inner_join(., visium_coord, by = "Barcode") %>%
  column_to_rownames(var = "Barcode")

ggplot(data = dataObj_brain@meta.data, 
       mapping = aes(x = X, y = Y, color = seurat_clusters)) + 
  geom_point(size = 2.5) +
  scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_brain@meta.data$X)+2) + 
  ylim(0,max(dataObj_brain@meta.data$Y)+2) + 
  theme_bw()
  
ggplot(data = dataObj_brain@meta.data %>%
         filter(nFeature_RNA > 300) %>%
         filter(!seurat_clusters %in% c(12,1,8)) %>%
         filter(X < 125 | Y > 17) %>%
         filter(X < 104 | Y > 15) %>%
         filter(X < 82 | Y > 13) %>%
         filter(Y > 4) %>%
         filter(X > 48 | Y > 7) %>%
         filter(X > 35 | Y < 60) %>%
         filter(X > 20 | Y < 50) %>%
         filter(X > 35 | Y > 10) %>%
         filter(X > 28 | Y > 15) %>%
         filter(X > 25 | Y > 16), 
       mapping = aes(x = X, y = Y, color = seurat_clusters)) + 
  geom_point(size = 2.5) +
  scale_color_manual(values = cbp2) + 
  # geom_hline(yintercept = c(10,15,16)) +
  # geom_vline(xintercept = c(35,28,25)) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  theme_bw()

dataObj_brain_filtered <- subset(dataObj_brain, subset = nFeature_RNA > 300 & seurat_clusters %in% c(0,2:7,9:11,13:18)) %>%
  subset(., subset = X < 125 | Y > 17) %>%
  subset(., subset = X < 104 | Y > 15) %>%
  subset(., subset = X < 82 | Y > 13) %>%
  subset(., subset = Y > 4) %>%
  subset(., subset = X > 48 | Y > 7) %>%
  subset(., subset = X > 35 | Y < 60) %>%
  subset(., subset = X > 20 | Y < 50) %>%
  subset(., subset = X > 35 | Y > 10) %>%
  subset(., subset = X > 28 | Y > 15) %>%
  subset(., subset = X > 25 | Y > 16)


dataObj_brain_filtered <- CreateSeuratObject(counts = dataObj_brain_filtered@assays$RNA@counts) %>%
  SCTransform(., ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(verbose = FALSE)

dataObj_brain_filtered@meta.data <- dataObj_brain_filtered@meta.data %>%
  rownames_to_column(var = "Barcode") %>%
  inner_join(., visium_coord, by = "Barcode") %>%
  column_to_rownames(var = "Barcode")

saveRDS(dataObj_brain_filtered, "./Filtered/10X_ST_Mouse_brain_filtered.rds")

ggplot(data = dataObj_brain_filtered@meta.data %>%
         filter(!seurat_clusters %in% c(13)) %>%
         filter(X < 75 | Y > 10), 
       mapping = aes(x = X, y = Y, color = seurat_clusters)) + 
  geom_point(size = 2.5) +
  scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_brain@meta.data$X)+2) + 
  ylim(0,max(dataObj_brain@meta.data$Y)+2) + 
  theme_bw()

dataObj_brain_filtered <- subset(dataObj_brain_filtered, subset = seurat_clusters %in% c(0:12,14:16)) %>%
  subset(., subset = X < 75 | Y > 10) %>%
  SCTransform(., ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(verbose = FALSE)

dataObj_brain_filtered@meta.data <- dataObj_brain_filtered@meta.data %>%
  dplyr::select(-c("X", "Y")) %>%
  rownames_to_column(var = "Barcode") %>%
  inner_join(., visium_coord, by = "Barcode") %>%
  column_to_rownames(var = "Barcode")

saveRDS(dataObj_brain_filtered, "./Filtered/10X_ST_Mouse_brain_filtered.rds")

ggplot(data = dataObj_brain_filtered@meta.data, 
       mapping = aes(x = X, y = Y, color = seurat_clusters)) + 
  geom_point(size = 2.5) +
  scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_brain@meta.data$X)+2) + 
  ylim(0,max(dataObj_brain@meta.data$Y)+2) + 
  theme_bw()

### Total Gene Number : 27981
### Total spots under tissue : 3190
### Mean Reads per spot : 239202.7
### Median Gene per spot : 2387
### Median UMI per spot : 6743.5
median(dataObj_brain_filtered@meta.data$nFeature_RNA)
median(dataObj_brain_filtered@meta.data$nCount_RNA)
print(778073298*0.9807/3190)


############## Heart ####################
dataObj_heart <- readRDS("./Input/10X_ST_Mouse_Heart_res.rds")

dataObj_heart@meta.data <- dataObj_heart@meta.data %>%
  rownames_to_column(var = "Barcode") %>%
  inner_join(., visium_coord, by = "Barcode") %>%
  column_to_rownames(var = "Barcode")

ggplot(data = dataObj_heart@meta.data, 
       mapping = aes(x = X, y = Y, color = seurat_clusters)) + 
  geom_point(size = 2.5) +
  scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_heart@meta.data$X)+2) + 
  ylim(0,max(dataObj_heart@meta.data$Y)+2) + 
  theme_bw()

ggplot(data = dataObj_heart@meta.data %>%
         filter(nFeature_RNA > 700) %>%
         filter(!seurat_clusters %in% c(7)), 
       mapping = aes(x = X, y = Y, color = seurat_clusters)) + 
  geom_point(size = 2.5) +
  scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_heart@meta.data$X)+2) + 
  ylim(0,max(dataObj_heart@meta.data$Y)+2) + 
  theme_bw()

dataObj_heart_filtered <- subset(dataObj_heart, subset = nFeature_RNA > 700 & seurat_clusters %in% c(0:6,8:10)) 

dataObj_heart_filtered <- CreateSeuratObject(counts = dataObj_heart_filtered@assays$RNA@counts) %>%
  SCTransform(., ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(verbose = FALSE)

dataObj_heart_filtered@meta.data <- dataObj_heart_filtered@meta.data %>%
  rownames_to_column(var = "Barcode") %>%
  inner_join(., visium_coord, by = "Barcode") %>%
  column_to_rownames(var = "Barcode")

saveRDS(dataObj_heart_filtered, "./Filtered/10X_ST_Mouse_Heart_filtered.rds")

ggplot(data = dataObj_heart_filtered@meta.data %>%
         filter(!seurat_clusters %in% c()), 
       mapping = aes(x = X, y = Y, color = seurat_clusters)) + 
  geom_point(size = 2.5) +
  scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_heart@meta.data$X)+2) + 
  ylim(0,max(dataObj_heart@meta.data$Y)+2) + 
  theme_bw()

### Total Gene Number : 29877
### Total spots under tissue : 3048
### Mean Reads per spot : 193663.4
### Median Gene per spot : 2599
### Median UMI per spot : 9924
median(dataObj_heart_filtered@meta.data$nFeature_RNA)
median(dataObj_heart_filtered@meta.data$nCount_RNA)
print(598109406*0.98692/3048)


median(dataObj_heart@meta.data$nCount_RNA)
median(dataObj_brain@meta.data$nCount_RNA)
median(dataObj_kidney@meta.data$nCount_RNA)








