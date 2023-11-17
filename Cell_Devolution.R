library(tidyverse)
library(Seurat)
library(loomR)
library(spacexr)


dataObj <- readRDS("./10X_ST_Mouse_brain_filtered.rds")

# Try reference from Mouse Brain Atlas
MouseBrainAtlas <- connect(filename = "/public/home/wangycgroup/public/02_Data/Internal/ST/10X_test/Integration_with_MBA/l5_all.loom", mode = "r", 
                           skip.validate = T)
MBA_meta <- data.frame(
  CellID = paste0(MouseBrainAtlas$col.attrs$CellID[],seq(160796)),
  CellType = MouseBrainAtlas$col.attrs$TaxonomyRank4[],
  orig.ident = "MBA"
) %>%
  column_to_rownames(var = "CellID")

MBA_matrix <- MouseBrainAtlas[["matrix"]][,]
rownames(MBA_matrix) <- paste0(MouseBrainAtlas$col.attrs$CellID[],seq(160796))
colnames(MBA_matrix) <- MouseBrainAtlas$row.attrs$Gene[]

MBA_seurat <- CreateSeuratObject(
  counts = MBA_matrix %>% t(),
  meta.data = MBA_meta
)

# extract information to pass to the RCTD Reference function
counts <- MBA_seurat@assays$RNA@counts
cluster <- as.factor(MBA_seurat@meta.data$CellType)
names(cluster) <- colnames(MBA_seurat)
nUMI <- MBA_seurat$nCount_RNA
names(nUMI) <- colnames(MBA_seurat)
reference <- Reference(counts, cluster, nUMI)

# set up query with the RCTD function SpatialRNA
st_counts <- dataObj@assays$RNA@counts
coords <- dataObj@meta.data %>%
  dplyr::select("X", "Y")
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, st_counts, colSums(st_counts))

RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD)
dataObj <- AddMetaData(dataObj, metadata = RCTD@results$results_df)
saveRDS(object = dataObj, file = "./10X_ST_Mouse_brain_filtered_devoluted_MBA.rds")










