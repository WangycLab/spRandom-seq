library(tidyverse)
library(Seurat)
library(loomR)
library(ggpubr)

cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#8B2DB2", "#0072B2", "#D55E00", "#CC79A7",
          "#CD534C", "#FFBF7F", "#FF5500", "#E51932",
          "#FF99BF", "#654CFF", "#CCBFFF", "#19B2FF",
          "#A5EDFF", "#CC5151", "#260F99", "#85B22C",
          "#5773CC", "#FFB900", "#350E20", "#3E3E23",
          "#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D",
          "#59A14F", "#8CD17D", "#B6992D", "#499894",
          "#86BCB6", "#E15759", "#FF9D9A", "#79706E",
          "#D37295", "#FABFD2", "#B07AA1", "#9D7660") 

dataObj_ST_brain <- readRDS("./Filtered/10X_ST_Mouse_brain_filtered.rds")
dataObj_10X_Visium <- readRDS("./10X_Data/10X_Visium_Mouse_Brain_v2_res.rds")
dataObj_ST_heart <- readRDS("./Filtered/10X_ST_Mouse_Heart_filtered.rds")

###### Classic Markers ########

### Astrocyte: Slc1a2
p1 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate(Slc1a2 = dataObj_ST_brain@assays$RNA@counts["Slc1a2", ]), 
       mapping = aes(x = X, y = Y, color = log2(Slc1a2+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Slc1a2") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

### Oligodendrocyte: Cldn11
p2 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate(Cldn11 = dataObj_ST_brain@assays$RNA@counts["Cldn11", ]), 
       mapping = aes(x = X, y = Y, color = log2(Cldn11+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Cldn11") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

### Oligodendrocyte: Mobp
p3 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate(Mobp = dataObj_ST_brain@assays$RNA@counts["Mobp", ]), 
       mapping = aes(x = X, y = Y, color = log2(Mobp+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Mobp") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

### Neuroendocrine: Syt1
p4 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate(Syt1 = dataObj_ST_brain@assays$RNA@counts["Syt1", ]), 
       mapping = aes(x = X, y = Y, color = log2(Syt1+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Syt1") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

### stem cell: Zeb1
p5 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate(Zeb1 = dataObj_ST_brain@assays$RNA@counts["Zeb1", ]), 
       mapping = aes(x = X, y = Y, color = log2(Zeb1+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Zeb1") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

### GABAergic: Gad2
p6 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate(Gad2 = dataObj_ST_brain@assays$RNA@counts["Gad2", ]), 
       mapping = aes(x = X, y = Y, color = log2(Gad2+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) +
  labs(color = "Gad2") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

### Glutamatergic: Slc17a6
p7 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate(Slc17a6 = dataObj_ST_brain@assays$RNA@counts["Slc17a6", ]), 
       mapping = aes(x = X, y = Y, color = log2(Slc17a6+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Slc17a6") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

### Glutamatergic: Snap25
p8 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate(Snap25 = dataObj_ST_brain@assays$RNA@counts["Snap25", ]), 
       mapping = aes(x = X, y = Y, color = log2(Snap25+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Snap25") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

### Prefrontal Cortex layer 6 cells: Syt6
p9 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate(Syt6 = dataObj_ST_brain@assays$RNA@counts["Syt6", ]), 
       mapping = aes(x = X, y = Y, color = log2(Syt6+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Syt6") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

### Lamina layer 4 neurons: Rorb
p10 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate(Rorb = dataObj_ST_brain@assays$RNA@counts["Rorb", ]), 
       mapping = aes(x = X, y = Y, color = log2(Rorb+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Rorb") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )


### CA2 pyramidal neurons: Sv2b
p11 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate(Sv2b = dataObj_ST_brain@assays$RNA@counts["Sv2b", ]), 
       mapping = aes(x = X, y = Y, color = log2(Sv2b+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Sv2b") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

### CA1 pyramidal neurons: Satb2
p12 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate(Satb2 = dataObj_ST_brain@assays$RNA@counts["Satb2", ]), 
       mapping = aes(x = X, y = Y, color = log2(Satb2+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Satb2") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

(p1 | p2 | p3 | p4 | p5 | p6)/(p7 | p8 | p9 | p10 | p11 | p12)

ggsave("./Figures_Draft/Classic_Markers.pdf", height = 12, width = 20, device = "pdf")

###### lncRNA Markers #########
p13 <- ggplot(data = dataObj_ST_brain@meta.data %>%
                 mutate("4921539H07Rik" = dataObj_ST_brain@assays$RNA@counts["4921539H07Rik", ]), 
               mapping = aes(x = X, y = Y, color = log2(`4921539H07Rik`+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "4921539H07Rik") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )


p14 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate("Gm6999" = dataObj_ST_brain@assays$RNA@counts["Gm6999", ]), 
       mapping = aes(x = X, y = Y, color = log2(`Gm6999`+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Gm6999") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

p15 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate("B130024G19Rik" = dataObj_ST_brain@assays$RNA@counts["B130024G19Rik", ]), 
       mapping = aes(x = X, y = Y, color = log2(`B130024G19Rik`+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "B130024G19Rik") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )


p16 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate("Gm30382" = dataObj_ST_brain@assays$RNA@counts["Gm30382", ]), 
       mapping = aes(x = X, y = Y, color = log2(`Gm30382`+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Gm30382") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

p17 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate("Gm4876" = dataObj_ST_brain@assays$RNA@counts["Gm4876", ]), 
       mapping = aes(x = X, y = Y, color = log2(`Gm4876`+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Gm4876") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

p18 <- ggplot(data = dataObj_ST_brain@meta.data %>%
         mutate("Rmst" = dataObj_ST_brain@assays$RNA@counts["Rmst", ]), 
       mapping = aes(x = X, y = Y, color = log2(`Rmst`+1))) +
  geom_point(size = 1) +
  scale_color_viridis_c(direction = -1) + 
  #scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "Rmst") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

p13 | p14 | p15 | p16 | p17 | p18

ggsave("./Figures_Draft/lncRNA_Markers.pdf", height = 6, width = 24, device = "pdf")

###### Brain section / Cell Deconvolution from RCTD #########
dataObj_brain_RCTD <- readRDS("./Cell_Deconvolution/10X_ST_Mouse_brain_filtered_devoluted_MBA.rds")
ggplot(data = dataObj_brain_RCTD@meta.data, 
       mapping = aes(x = X, y = Y, color = first_type)) + 
  geom_point(size = 2.5) +
  scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  theme_bw() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

ggplot(data = dataObj_brain_RCTD@meta.data, 
       mapping = aes(x = X, y = Y, color = second_type)) + 
  geom_point(size = 2) +
  scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  theme_bw() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 23),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

ggplot(data = dataObj_brain_RCTD@meta.data, 
       mapping = aes(x = X, y = Y, color = second_type)) + 
  geom_point(size = 2) +
  scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  theme_bw() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 23),
    legend.position = "None"
  ) + 
  facet_wrap(~second_type, ncol = 5)

dataObj_brain_RCTD@meta.data <- dataObj_brain_RCTD@meta.data %>%
  mutate(
    Spot_type = case_when(
      second_type == "Hindbrain neurons" ~ first_type,
      second_type != "Hindbrain neurons" ~ second_type
    )
  )

dataObj_brain_RCTD@meta.data <- dataObj_brain_RCTD@meta.data %>%
  mutate(
    Spot_type = case_when(
      Spot_type == "Enteric glia" ~ first_type,
      Spot_type == "Enteric neurons" ~ first_type,
      Spot_type == "Schwann cells" ~ first_type,
      Spot_type == "Dentate gyrus radial glia-like cells" ~ first_type,
      Spot_type == "Olfactory ensheathing cells" ~ first_type,
      Spot_type == "Olfactory inhibitory neurons" ~ first_type,
      Spot_type == "Spinal cord excitatory neurons" ~ first_type,
      Spot_type == "Spinal cord inhibitory neurons" ~ first_type,
      Spot_type == "Spinal cord inhibitory neurons" ~ first_type,
      Spot_type == first_type ~ first_type,
      Spot_type == second_type ~ second_type 
    )
  ) %>%
  mutate(
    Spot_type = case_when(
      Spot_type == "Subventricular zone radial glia-like cells" | Spot_type == "Satellite glia" ~ "Other glial cells",
      Spot_type == "Perivascular macrophages" ~ "Microglia",
      Spot_type != "Subventricular zone radial glia-like cells" &
        Spot_type != "Satellite glia" &
        Spot_type != "Perivascular macrophages" ~ Spot_type
    )
  )

ggplot(data = dataObj_brain_RCTD@meta.data, 
       mapping = aes(x = X, y = Y, color = Spot_type)) + 
  geom_point(size = 2) +
  scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  theme_bw() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 23),
    legend.position = "None"
  ) + 
  facet_wrap(~Spot_type, ncol = 5)

# Final Annotation
saveRDS(dataObj_brain_RCTD, "./Cell_Deconvolution/10X_ST_Mouse_brain_filtered_annotated.rds")
ggplot(data = dataObj_brain_RCTD@meta.data, 
       mapping = aes(x = X, y = Y, color = Spot_type)) + 
  geom_point(size = 3.5) +
  scale_color_manual(values = cbp2) + 
  coord_flip() + 
  xlim(0,max(dataObj_kidney@meta.data$X)+2) + 
  ylim(0,max(dataObj_kidney@meta.data$Y)+2) + 
  labs(color = "") + 
  theme_void() + 
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    #legend.position = "bottom",
    legend.text = element_text(size = 18)
  ) + 
  guides(color = guide_legend(ncol = 1, byrow = F))
ggsave("./Figures_Draft/Main_Figure_Brain_Section_Annotated.pdf", height = 12, width = 18, device = "pdf")

###### Biotype per Cell Distribution ########
gene_biotype <- read.table("./Data/gene_biotype.gtf", row.names = NULL) %>%
  dplyr::select(gene_symbol, gene_type) %>%
  filter(gene_type %in% c("protein_coding", "lncRNA","miRNA", 
                          "misc_RNA", "snRNA", "snoRNA", "rRNA", "scaRNA","Mt_tRNA"))

# Boxplot
boxplots <- lapply(unique(gene_biotype$gene_type)[1:6], function(x){
  
  temp_meta_STR <- dataObj_ST_brain@meta.data %>%
    mutate(temp = apply(dataObj_ST_brain@assays$RNA@counts[gene_biotype %>% 
                                                      filter(gene_type == x & gene_symbol %in% rownames(dataObj_ST_brain)) %>% 
                                                      pull(gene_symbol), ]>0, 2, sum),
           total_RNA = apply(dataObj_ST_brain@assays$RNA@counts[gene_biotype %>% 
                                                                  filter(gene_symbol %in% rownames(dataObj_ST_brain)) %>% 
                                                                  pull(gene_symbol), ]>0, 2, sum)) %>%
    mutate(
      Project = "Mouse Brain FFPE\n(STR-seq)",
      temp = (temp/total_RNA)*100
    ) %>%
    dplyr::select(Project, temp)
  
  temp_meta_10X <- dataObj_10X_Visium@meta.data %>%
    mutate(temp = apply(dataObj_10X_Visium@assays$RNA@counts[gene_biotype %>% 
                                                             filter(gene_type == x & gene_symbol %in% rownames(dataObj_10X_Visium)) %>% 
                                                             pull(gene_symbol), ]>0, 2, sum)) %>%
    mutate(
      Project = "Mouse Brain FF\n(10X Visium)",
      temp = (temp/nFeature_RNA)*100
    ) %>%
    dplyr::select(Project, temp)
  
  temp_p <- bind_rows(temp_meta_STR, temp_meta_10X) %>%
    ggplot(data = ., 
           mapping = aes(x = Project, y = temp, 
                         color = Project)) +
    geom_boxplot(size = 1, fill = "white", width = 0.5) + 
    labs(y = paste0(x," (%)")) + 
    # stat_compare_means(aes(label = paste0(..method.., "\n", after_stat(p.format))),
    #                    size = 6,
    #                    label.y = max((bind_rows(temp_meta_STR, temp_meta_10X))$temp)) +
    ylim(0,max((bind_rows(temp_meta_STR, temp_meta_10X))$temp)+1) +
    theme_bw() + 
    theme(
      axis.text = element_text(size = 21, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.ticks = element_line(linewidth = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 23, face = "bold"),
      legend.position = "none",
      plot.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_rect(linewidth = 1, color = "black"),
      axis.line = element_line(linewidth = 1, color = "black")
    )
  
  temp_p
  
})
  
(boxplots[[1]] | boxplots[[2]] | boxplots[[3]]) / 
  (boxplots[[4]] | boxplots[[5]] | boxplots[[6]])

ggsave("./Figures_Draft/Biotype_UMI_prop_boxplot_geneCount_v2.pdf", height = 12, width = 15, device = "pdf")

###### Biotype distribution barplot #########
read.table("./Data/gene_biotype.gtf", row.names = NULL) %>%
  dplyr::select(gene_symbol, gene_type) %>%
  filter(gene_symbol %in% rownames(dataObj_ST_brain)) %>%
  mutate(gene_type = case_when(
    gene_type %in% c("protein_coding", "lncRNA","miRNA", "scaRNA",
                     "misc_RNA", "snRNA", "snoRNA", "rRNA", "Mt_tRNA") ~ gene_type,
    ! gene_type %in% c("protein_coding", "lncRNA","miRNA", "scaRNA",
                       "misc_RNA", "snRNA", "snoRNA", "rRNA", "Mt_tRNA") ~ "Others"
  )) %>%
  group_by(gene_type) %>%
  summarise(counts = n()) %>%
  mutate(gene_type = factor(gene_type, levels = c("protein_coding", "lncRNA","miRNA", "scaRNA",
                                                  "misc_RNA", "snRNA", "snoRNA", "rRNA", 
                                                  "Mt_tRNA", "Others"))) %>%
  ggplot(data = ., mapping = aes(x = gene_type, y = counts, fill = gene_type)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = counts),size = 6, nudge_y = 0.2) + 
  #annotate("text", x = .$gene_type, y = .$counts, label = .$counts, fontface = "bold") + 
  scale_y_continuous (trans = "log10") + 
  scale_fill_manual(values = cbp2) + 
  labs(x = "", y = "Gene Counts") + 
  theme_classic() + 
  theme(
    axis.text = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x = element_blank(),
    title = element_text(size = 30),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    #panel.border = element_rect(linewidth = 1, color = "black"),
    axis.line = element_line(linewidth = 1, color = "black"),
    legend.position = "none",
    text = element_text(size = 25, face = "bold")
  )
ggsave("./Figures_Draft/Biotype_gene_counts_barplot.pdf", 
       height = 5, width = 10, device = "pdf")

###### highly expressed RNAs marker #########
highly_expressed_biotype <- lapply(unique(gene_biotype$gene_type), function(x){
  
  temp_data <- dataObj_ST_brain@assays$RNA@counts[gene_biotype %>% 
                                                    filter(gene_type == x & gene_symbol %in% rownames(dataObj_ST_brain)) %>% 
                                                    pull(gene_symbol), ]
  temp_res <- apply(temp_data>0, 1, sum) %>%
    enframe(., name = "Gene_name", value = "Spots") %>%
    arrange(-Spots) %>%
    head(5) %>%
    mutate(Biotype = x)
  
  temp_res
}) %>%
  bind_rows()

highly_expressed_biotype_genes <- c("Actb","Mir6236", "Mir6240", "Rnu12", "Gm24265",
                                    "Malat1", "Neat1", "Gm17344", "Snord118",
                                    "Rny1", "Rny3", "Rn7sk", "mt-Tv", "mt-Th")

heg_STR <- log2(dataObj_ST_brain@assays$RNA@counts[highly_expressed_biotype_genes, ] + 1) %>%
  apply(., 1, mean) %>%
  enframe(., name = "Gene_name", value = "Mean_expr") %>%
  mutate(
    percent = apply(dataObj_ST_brain@assays$RNA@counts[highly_expressed_biotype_genes, ]>0, 1, sum) %>%
      enframe(., name = "Gene_name", value = "Spots") %>%
      mutate(Spots = Spots/ncol(dataObj_ST_brain)*100) %>%
      pull(Spots)
  )

overlapped_heg <- intersect(highly_expressed_biotype_genes, rownames(dataObj_10X_Visium))
heg_10X <- log2(dataObj_10X_Visium@assays$RNA@counts[overlapped_heg, ] + 1) %>%
  apply(., 1, mean) %>%
  enframe(., name = "Gene_name", value = "Mean_expr") %>%
  mutate(
    percent = apply(dataObj_10X_Visium@assays$RNA@counts[overlapped_heg, ]>0, 1, sum) %>%
      enframe(., name = "Gene_name", value = "Spots") %>%
      mutate(Spots = Spots/ncol(dataObj_10X_Visium)*100) %>%
      pull(Spots)
  )

heg_10X <- heg_STR %>%
  dplyr::select(Gene_name) %>%
  left_join(., heg_10X, by = "Gene_name")
heg_10X[is.na(heg_10X)] <- 0

bind_rows(heg_STR %>% mutate(Method = "Mouse Brain FFPE\n(STR-seq)"), 
          heg_10X %>% mutate(Method = "Mouse Brain FF\n(10X Visium)")) %>%
  ggplot(data = .) + 
  geom_point(
    mapping = aes(x = Method, y = factor(Gene_name, levels = rev(highly_expressed_biotype_genes)), 
                  color = Mean_expr, size = percent)
  ) + 
  scale_color_viridis_c() +
  labs(
    color = "Average Log-transformed Expression",
    size = "Percentage Expressed",
    y = ""
  ) + 
  theme_bw() + 
  theme(
    axis.text = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x = element_blank(),
    title = element_text(size = 30),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(linewidth = 1, color = "black"),
    axis.line = element_line(linewidth = 1, color = "black"),
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 25, face = "bold")
  )

ggsave("./Figures_Draft/Biotype_marker_gene_DotPlot_v2.pdf", 
       height = 10, width = 10, device = "pdf")


###### Biotype Marker Featureplot ################
visium_coord_v1 <- read.table("./Data/visium-v1_coordinates.txt")
colnames(visium_coord_v1) <- c("Barcode", "X", "Y")
visium_coord_v1 <- visium_coord_v1 %>%
  arrange(X,Y)

dataObj_10X_Visium@meta.data <- dataObj_10X_Visium@meta.data %>%
  rownames_to_column(var = "Barcode") %>%
  mutate(Barcode = gsub("-1", "", Barcode)) %>%
  inner_join(., visium_coord_v1, by = "Barcode") %>%
  column_to_rownames(var = "Barcode")

ggplot(data = dataObj_10X_Visium@meta.data, 
       mapping = aes(x = X, y = Y, color = seurat_clusters)) + 
  geom_point(size = 2.5) +
  scale_color_manual(values = cbp2) + 
  xlim(0,max(dataObj_10X_Visium@meta.data$X)+2) + 
  ylim(0,max(dataObj_10X_Visium@meta.data$Y)+2) + 
  theme_bw()

biotype_marker_featureplots <- lapply(highly_expressed_biotype_genes[c(1,2,4,6,7,9,10,14)], function(x){
  
  temp_10X <- ""
  if (x %in% rownames(dataObj_10X_Visium)){
    temp_10X <- dataObj_10X_Visium@meta.data %>%
      dplyr::select(X, Y) %>%
      mutate(
        Y = Y + 100,
        temp = dataObj_10X_Visium@assays$RNA@counts[x, ])
  } else{
    temp_10X <- dataObj_10X_Visium@meta.data %>%
      dplyr::select(X, Y) %>%
      mutate(
        Y = Y + 100,
        temp = 0)
  }
  temp_STR <- dataObj_ST_brain@meta.data %>%
    dplyr::select(X, Y) %>%
    mutate(temp = dataObj_ST_brain@assays$RNA@counts[x, ])
  
  
  temp_p <- ggplot(data = bind_rows(temp_10X, temp_STR), 
                   mapping = aes(x = X, y = Y, color = log2(temp+1))) +
    geom_point(size = 0.5) +
    scale_color_viridis_c(direction = -1) + 
    #scale_color_manual(values = cbp2) + 
    xlim(0,max(bind_rows(temp_10X, temp_STR)$X)+2) + 
    ylim(0,max(bind_rows(temp_10X, temp_STR)$Y)+2) + 
    labs(color = x) + 
    theme_void() + 
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",
      title = element_text(size = 23),
      legend.text = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 22)
    )
  temp_p
}) 

(biotype_marker_featureplots[[1]] | biotype_marker_featureplots[[2]] | 
    biotype_marker_featureplots[[3]] | biotype_marker_featureplots[[4]] + 
    theme(legend.text = element_text(size = 10, face = "bold")))/
  ((biotype_marker_featureplots[[5]] | biotype_marker_featureplots[[6]] + 
      theme(legend.text = element_text(size = 10, face = "bold")) | 
      biotype_marker_featureplots[[7]] | biotype_marker_featureplots[[8]]))

ggsave("./Figures_Draft/Biotype_Markers_v2.pdf", height = 12, width = 16, device = "pdf")




##### Heart Sample #######

# load heart markers
heart_markers <- readxl::read_xlsx("./Data/Heart marker genes.xlsx", 
                                   col_names = F)
heart_markers <- heart_markers %>%
  unite(., col = "markers", colnames(.)[2:5], sep = ",", 
        remove = T, na.rm = T) %>%
  separate_rows(., markers, sep = ",") %>%
  dplyr::rename("CellType" = "...1")

# classic markers
heart_DEGs <- FindAllMarkers(dataObj_ST_heart, assay = "RNA", only.pos = T)

selected_markers_heart <- heart_DEGs %>%
  filter(gene %in% heart_markers$markers) %>%
  filter(p_val_adj <= 0.05) %>%
  group_by(cluster) %>%
  arrange(-(pct.1-pct.2)) %>%
  arrange(cluster) %>%
  top_n(5)

heart_marker_featureplots <- lapply(c(unique(selected_markers_heart$gene)[1:11], "Myl4"), function(x){
  
  temp_STR <- dataObj_ST_heart@meta.data %>%
    dplyr::select(X, Y) %>%
    mutate(temp = dataObj_ST_heart@assays$RNA@counts[x, ])
  
  
  temp_p <- ggplot(data = temp_STR, 
                   mapping = aes(x = X, y = Y, color = log2(temp+1))) +
    geom_point(size = 1.5) +
    scale_color_viridis_c(direction = -1) + 
    #scale_color_manual(values = cbp2) + 
    scale_x_reverse(limits = c(max(temp_STR$X)+2,0)) + 
    #scale_y_reverse(limits = c(max(temp_STR$Y)+2,0)) + 
    #xlim(max(temp_STR$X)+2,0) + 
    #ylim(max(temp_STR$Y)+2,0) + 
    coord_flip() + 
    labs(color = x) + 
    theme_void() + 
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",
      title = element_text(size = 23),
      legend.text = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 22)
    )
  temp_p
}) 

(heart_marker_featureplots[[1]] | heart_marker_featureplots[[2]] | heart_marker_featureplots[[3]])/
  (heart_marker_featureplots[[4]] | heart_marker_featureplots[[5]] | heart_marker_featureplots[[6]])/
  (heart_marker_featureplots[[7]] | heart_marker_featureplots[[8]] | heart_marker_featureplots[[9]])/
  (heart_marker_featureplots[[10]] +  
     theme(legend.text = element_text(size = 10, face = "bold")) | heart_marker_featureplots[[11]] | heart_marker_featureplots[[12]])

ggsave("./Figures_Draft/Heart_Classic_Markers.pdf", height = 24, width = 14, device = "pdf")

###### Gene Body #######
STR_brain_genebody <- read.table("./Genebody/10X_Mouse_Brain_add.geneBodyCoverage.txt",
                                 header = T) %>%
  dplyr::select(2:101) %>%
  pivot_longer(cols = colnames(.), names_to = "Percentile", values_to = "Coverage") %>%
  mutate(
    Percentile = as.numeric(gsub("X","", Percentile)),
    Project = "Mouse Brain FFPE\n(STR-seq)",
    Coverage = Coverage/sum(Coverage)*100
  )

STR_heart_genebody <- read.table("./Genebody/10X_Mouse_Heart.geneBodyCoverage.txt",
                                 header = T) %>%
  dplyr::select(2:101) %>%
  pivot_longer(cols = colnames(.), names_to = "Percentile", values_to = "Coverage") %>%
  mutate(
    Percentile = as.numeric(gsub("X","", Percentile)),
    Project = "Mouse Heart FFPE\n(STR-seq)",
    Coverage = Coverage/sum(Coverage)*100
  )

bind_rows(STR_brain_genebody, STR_heart_genebody) %>%
  ggplot(data = ., 
         mapping = aes(x = Percentile, y = Coverage,
                       color = Project, group = Project)) +
  geom_line(size = 3) +
  ylim(0,2) + 
  scale_color_manual(values = cbp2[2:3]) + 
  labs(x = "Gene body Percentile (5' to 3')", y = "Coverage (%)") + 
  theme_classic() + 
  theme(
    axis.text = element_text(size = 25, face = "bold"),
    #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks = element_line(linewidth = 1),
    title = element_text(size = 26),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    #panel.border = element_rect(linewidth = 1, color = "black"),
    axis.line = element_line(linewidth = 1, color = "black"),
    legend.position = c(0.8,0.9),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    text = element_text(size = 25, face = "bold")
  )
ggsave("./Figures_Draft/GeneBody_STR_brain_heart.pdf", 
       height = 5, width = 10, device = "pdf")

###### Correlation Plot ########
DefaultAssay(dataObj_ST_brain) <- "RNA"
Corr_STR_brain <- NormalizeData(dataObj_ST_brain,
                                normalization.method = "RC",
                                scale.factor = 1e6) 
DefaultAssay(dataObj_10X_Visium) <- "RNA"
Corr_10X_brain <- NormalizeData(dataObj_10X_Visium,
                                normalization.method = "RC",
                                scale.factor = 1e6
                                )

# Average and Log-transformation
mean_STR_cpm <- rowMeans(Corr_STR_brain@assays$RNA@data)
#mean_STR_cpm <- mean_STR_cpm[mean_STR_cpm>1]
mean_STR_cpm <- log2(mean_STR_cpm+1) %>%
  enframe(., name = "Gene_name", value = "Expr_val_STR")
mean_10X_cpm <- rowMeans(Corr_10X_brain@assays$RNA@data)
#mean_10X_cpm <- mean_10X_cpm[mean_10X_cpm>1]
mean_10X_cpm <- log2(mean_10X_cpm+1) %>%
  enframe(., name = "Gene_name", value = "Expr_val_10X")

# Merge two dataset
Corr_STR_10X_df <- inner_join(
  mean_STR_cpm, mean_10X_cpm,
  by = "Gene_name"
)

ggplot(Corr_STR_10X_df,
       mapping = aes(x = Expr_val_STR, y = Expr_val_10X)) + 
  geom_point(size = 2, alpha = 0.3) + 
  stat_smooth(
    formula = y ~ x, 
    color = "red", linetype = 6, size = 1.5 
  ) + 
  stat_cor(
    size = 8,
    method = "pearson"
    # label.x = 0.2, label.y = 13,
  ) + 
  labs(
    x = "Log2( CPM + 1) \nMouse Brain FFPE (STR-seq)",
    y = "Log2( CPM + 1) \nMouse Brain FF (10X Visium)"
  ) + 
  theme_classic() + 
  theme(
    axis.text = element_text(size = 25, face = "bold"),
    #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks = element_line(linewidth = 1),
    title = element_text(size = 26),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(linewidth = 1, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    text = element_text(size = 25, face = "bold")
  )

ggsave("./Figures_Draft/Correlation_Plot_STR_vs_10X_Mouse_Brain.pdf", 
       height = 12, width = 16, device = "pdf")

###### Single gene genebody coverage ######
single_gene_coverage <- lapply(c("Actb", "Mir6240", "Rnu12", "Malat1", "Neat1", "Snord118", "Rny3", "mt-Th"), function(x){
  
  temp_sense <- read.table(paste0("./Genebody/",x,"_genebody.txt")) %>%
    pivot_longer(., cols = colnames(.), 
                 names_to = "Pos", values_to = "Coverage") %>%
    mutate(
      Pos = gsub("V", "", Pos) %>% as.numeric(),
      Coverage = log10(Coverage+1),
      Strand = "sense"
    )
  temp_antisense <- read.table(paste0("./Genebody/anti_",x,"_genebody.txt")) %>%
    pivot_longer(., cols = colnames(.), 
                 names_to = "Pos", values_to = "Coverage") %>%
    mutate(
      Pos = gsub("V", "", Pos) %>% as.numeric(),
      Coverage = -log10(Coverage+1),
      Strand = "anti-sense"
    )
  
  temp_p <- ggplot(data = bind_rows(temp_sense, temp_antisense) %>%
           mutate(Strand = factor(Strand, levels = c("sense", "anti-sense"))), 
         mapping = aes(x = Pos, y = Coverage,
                       fill = Strand)) + 
    geom_bar(stat = "identity") + 
    geom_hline(yintercept = 0, size = 2, color = "grey60") + 
    scale_fill_manual(values = c(alpha("red", 1), 
                                 alpha("navyblue", 1))) + 
    labs(
      y = "Log10( Coverage )"
    ) + 
    xlim(0, max(bind_rows(temp_sense, temp_antisense)$Pos)) + 
    theme_classic() + 
    theme(
      axis.text = element_text(size = 21, face = "bold"),
      axis.ticks = element_line(linewidth = 2, color = "grey60"),
      axis.ticks.length = unit(0.3, "cm"),
      title = element_text(size = 21),
      plot.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line = element_line(linewidth = 2, color = "grey60"),
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      text = element_text(size = 21, face = "bold"),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    )
  
  temp_p
  # ggsave(paste0("./Figures_Draft/", x, "_genebody_coverage.pdf"),
  #        device = "pdf", height = 5, width = 10)
  
})

single_gene_coverage_10X <- lapply(c("Actb", "Mir6240", "Rnu12", "Malat1", "Neat1", "Snord118", "Rny3", "mt-Th"), function(x){
  
  temp_sense <- read.table(paste0("./Genebody/",x,"_genebody_10X.txt")) %>%
    pivot_longer(., cols = colnames(.), 
                 names_to = "Pos", values_to = "Coverage") %>%
    mutate(
      Pos = gsub("V", "", Pos) %>% as.numeric(),
      Coverage = log10(Coverage+1),
      Strand = "sense"
    )
  temp_antisense <- read.table(paste0("./Genebody/anti_",x,"_genebody_10X.txt")) %>%
    pivot_longer(., cols = colnames(.), 
                 names_to = "Pos", values_to = "Coverage") %>%
    mutate(
      Pos = gsub("V", "", Pos) %>% as.numeric(),
      Coverage = -log10(Coverage+1),
      Strand = "anti-sense"
    )
  
  temp_p <- ggplot(data = bind_rows(temp_sense, temp_antisense) %>%
           mutate(Strand = factor(Strand, levels = c("sense", "anti-sense"))), 
         mapping = aes(x = Pos, y = Coverage,
                       fill = Strand)) + 
    geom_bar(stat = "identity") + 
    geom_hline(yintercept = 0, size = 2, color = "grey60") + 
    scale_fill_manual(values = c(alpha("red", 1), 
                                 alpha("navyblue", 1))) + 
    labs(
      y = "Log10( Coverage )"
    ) + 
    xlim(0, max(bind_rows(temp_sense, temp_antisense)$Pos)) + 
    theme_classic() + 
    theme(
      axis.text = element_text(size = 21, face = "bold"),
      axis.ticks = element_line(linewidth = 2, color = "grey60"),
      axis.ticks.length = unit(0.3, "cm"),
      title = element_text(size = 21),
      plot.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.line = element_line(linewidth = 2, color = "grey60"),
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      text = element_text(size = 21, face = "bold"),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    )
  
  temp_p
  # ggsave(paste0("./Figures_Draft/", x, "_genebody_coverage_10X.pdf"),
  #        device = "pdf", height = 5, width = 10)
  
})


(single_gene_coverage[[1]] / single_gene_coverage_10X[[1]]) |
  (single_gene_coverage[[4]] / single_gene_coverage_10X[[4]]) |
  (single_gene_coverage[[5]] / single_gene_coverage_10X[[5]])


ggsave(paste0("./Figures_Draft/Single_genebody_coverage_combined.pdf"),
       device = "pdf", height = 8, width = 27)















