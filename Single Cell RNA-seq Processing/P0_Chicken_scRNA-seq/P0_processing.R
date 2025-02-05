# For RStudio, set the working directory to this file's path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary libraries
library(tidyverse)
library(Seurat)
library(BUSpaRse)
library(SeuratWrappers)
library(DropletUtils)
library(velocyto.R)

# Configure directories
dpath <- "P0_counts_unfiltered/"
project <- "P0_output"
dir.create(project)
outtmp <- file.path(project, "tmp")
dir.create(outtmp)
outfinal <- file.path(project, "final")
dir.create(outfinal)
writeLines(capture.output(sessionInfo()), file.path(project, "sessionInfo.txt"))

# Load spliced and unspliced matrices
v_output <- read_velocity_output(
  spliced_dir = file.path(dpath), 
  spliced_name='spliced',
  unspliced_dir = file.path(dpath),
  unspliced_name='unspliced')

spliced <- v_output$spliced
unspliced <- v_output$unspliced

# Combine spliced and unspliced matrices
shared <- intersect(colnames(spliced), colnames(unspliced))
combined <- spliced[, shared]+unspliced[, shared]
combined <- cbind(combined, spliced[, setdiff(colnames(spliced), colnames(unspliced))])
combined <- cbind(combined, unspliced[, setdiff(colnames(unspliced), colnames(spliced))])

# Barcode rank QC
bc_rank_spliced <- barcodeRanks(spliced)
bc_rank_unspliced <- barcodeRanks(unspliced)
bc_rank_combined <- barcodeRanks(combined)

tibble(rank = bc_rank_spliced$rank, total = bc_rank_spliced$total, matrix = "spliced") %>%
  bind_rows(tibble(rank = bc_rank_unspliced$rank, total = bc_rank_unspliced$total, matrix = "unspliced")) %>%
  distinct() %>%
  ggplot(aes(total, rank, color = matrix)) +
  geom_line() +
  scale_color_manual(values = c(spliced = "#0074D9", unspliced = "#FF4136")) +
  geom_vline(xintercept = metadata(bc_rank_spliced)$inflection, color = "#0074D9", linetype = 3) +
  geom_vline(xintercept = metadata(bc_rank_unspliced)$inflection, color = "#FF4136", linetype = 3) +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Rank", x = "Total UMI counts") +
  theme_classic()

# Also QC combined matrix
bc_rank <- bc_rank_combined
qplot(bc_rank$total, bc_rank$rank, geom = "line") +
  geom_vline(xintercept = metadata(bc_rank)$knee, color = "blue", linetype = 2) +
  geom_vline(xintercept = metadata(bc_rank)$inflection, color = "green", linetype = 2) +
  annotate("text", y = 1000, x = 1.5 * c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
           label = c("knee", "inflection"), color = c("blue", "green")) +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Barcode rank", x = "Total UMI count")

bcs_filtered_spliced <- colnames(spliced)[colSums(spliced) > metadata(bc_rank_spliced)$inflection]
bcs_filtered_unspliced <- colnames(unspliced)[colSums(unspliced) > metadata(bc_rank_unspliced)$inflection]
bcs_filtered_combined <- colnames(combined)[colSums(combined) > metadata(bc_rank_combined)$inflection]

# Filter matrics using combined inflection point
spliced_filtered <- spliced[, bcs_filtered_combined]
unspliced_filtered <- unspliced[, bcs_filtered_combined]
combined_filtered <- combined[, bcs_filtered_combined]

features_filtered <- rownames(combined_filtered)[Matrix::rowSums(combined_filtered) > 0]

spliced_filtered <- spliced_filtered[features_filtered,]
unspliced_filtered <- unspliced_filtered[features_filtered,]
combined_filtered <- combined_filtered[features_filtered, ]

# Splicing ratio QC
overall_ratio <- sum(unspliced_filtered@x) / (sum(unspliced_filtered@x) + sum(spliced_filtered@x))
calculate_ratio <- function(x, y) {
  ratios <- x / y
  ratios[x == 0] <- 0
  ratios[y == 0] <- Inf
  ratios[x == 0 & y == 0] <- NA
  ratios
}
ratios <- calculate_ratio(Matrix::rowSums(unspliced_filtered), Matrix::rowSums(spliced_filtered))
ratios <- tibble(ratio = ratios, gene = rownames(unspliced_filtered))
ratios %>% 
  filter(is.finite(ratio) & ratio > 0) %>% 
  ggplot(aes(x = ratio)) +
  geom_density() +
  scale_x_log10("Ratio per gene", breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3), labels = identity) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = overall_ratio, color = "red") +
  theme_classic() +
  scale_y_continuous("Gene density", expand = c(0, 0)) +
  annotate(geom = "text", label = "More unspliced", y = Inf, x = Inf, hjust = 1, vjust = 1) +
  annotate(geom = "text", label = "More spliced", y = Inf, x = 0, hjust = 0, vjust = 1)

ratios %>% 
  mutate(ratio = ifelse(ratio > 0 & ratio < Inf, 1, ratio)) %>% 
  pull(ratio) %>% 
  table(useNA = "ifany")

# Check splice ratio per cell
common_cells <- intersect(colnames(unspliced_filtered), colnames(spliced_filtered))
ratios <- calculate_ratio(Matrix::colSums(unspliced_filtered[, common_cells]), Matrix::colSums(spliced_filtered[, common_cells]))
ratios <- tibble(ratio = ratios, cell = common_cells)

ratios %>% 
  filter(is.finite(ratio) & ratio > 0) %>% 
  ggplot(aes(x = ratio)) +
  geom_density() +
  scale_x_log10("Ratio per cell", breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3), labels = identity) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = overall_ratio, color = "red") +
  theme_classic() +
  scale_y_continuous("Cell density", expand = c(0, 0)) +
  annotate(geom = "text", label = "More unspliced", y = Inf, x = Inf, hjust = 1, vjust = 1) +
  annotate(geom = "text", label = "More spliced", y = Inf, x = 0, hjust = 0, vjust = 1)

ratios %>% 
  mutate(ratio = ifelse(ratio > 0 & ratio < Inf, 1, ratio)) %>% 
  pull(ratio) %>% 
  table(useNA = "ifany")

# Rename genes to satisfy Seurat's requirement
rownames(spliced_filtered) <- gsub('_', '-', rownames(spliced_filtered))
rownames(unspliced_filtered) <- gsub('_', '-', rownames(unspliced_filtered))
rownames(spliced) <- gsub('_', '-', rownames(spliced))
rownames(unspliced) <- gsub('_', '-', rownames(unspliced))
rownames(combined_filtered) <- gsub("_", "-", rownames(combined_filtered))

#################
# Create object #
#################
seu <- CreateSeuratObject(combined_filtered, assay = "RNA", min.cells = 3, min.features = 500)
seu[["spliced"]] <- CreateAssayObject(spliced[rownames(seu@assays$RNA$counts), colnames(seu@assays$RNA$counts)]) %>%
  NormalizeData()
seu[["unspliced"]] <- CreateAssayObject(unspliced[rownames(seu@assays$RNA$counts), colnames(seu@assays$RNA$counts)]) %>%
  NormalizeData()

# Cell filtration by multiple parameters
# Splice ratio
ratio <- colSums(seu@assays$unspliced$counts)/colSums(seu@assays$spliced$counts)
seu@meta.data$Log2ratio <- log(ratio, 2)
# Mitochondrial gene percent
seu <- PercentageFeatureSet(seu, pattern  = "^J6367-", col.name = "percent.mt")
#Remove mitochondrial genes in the subsequent analysis
all.genes <- rownames(seu)
features <- all.genes[!grepl("^J6367-", all.genes)]
seu <- subset(seu, features=features)

# Cell quality QC plots first round
p1 <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "Log2ratio"), ncol = 4)
p2 <- FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend()
p3 <- FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "Log2ratio") + NoLegend()
p4 <- FeatureScatter(seu, feature1 = "percent.mt", feature2 = "Log2ratio") + NoLegend()
p5 <- FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "nCount_RNA") + NoLegend()
plot1 <- p1 + p2 + p3 + p4 + p5
ggsave(file = file.path(project,"featurep_plot1st.png"), plot = plot1, dpi = 100, width=15, height=15)

# Subset cells accordingly
seu <- subset(seu, subset = percent.mt <= 20)
seu <- subset(seu, subset = Log2ratio >= -4)
seu <- subset(seu, subset = nCount_RNA <= 20000)

# Cell quality QC plots second round
p1 <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "Log2ratio"), ncol = 4)
p2 <- FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend()
p3 <- FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "Log2ratio") + NoLegend()
p4 <- FeatureScatter(seu, feature1 = "percent.mt", feature2 = "Log2ratio") + NoLegend()
p5 <- FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "nCount_RNA") + NoLegend()
plot1 <- p1 + p2 + p3 + p4 + p5
ggsave(file = file.path(project,"featurep_plot1st_sel.png"), plot = plot1, dpi = 100, width=15, height=15)

# Data Normalization
seu <- seu %>%
  NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData() %>%
  SCTransform(method = "glmGamPoi", verbose = TRUE, assay='RNA') 

# Clustering
seu <- RunPCA(seu, verbose = FALSE, npcs = 100, features=features)
ElbowPlot(seu, ndims = 100, reduction = "pca")
set.seed(1)

PC_number <- 51
seu <- FindNeighbors(seu, dims = 1:PC_number)
seu <- RunUMAP(seu, dims = 1:PC_number)
res = 0.1
seu <- FindClusters(seu, resolution = res, random.seed = 1)


dp <- DimPlot(seu, label = TRUE, reduction = "umap") + 
  NoLegend() + theme(axis.text.x=element_blank(), 
                     axis.ticks.x=element_blank(), 
                     axis.text.y=element_blank(), 
                     axis.ticks.y=element_blank()) 
ggsave(file = file.path(outfinal, paste("P0_Complete_clusters.png", sep="")), plot = dp, dpi = 600, width = 5, height = 5)

fp <- FeaturePlot(seu, features = c("OPN1LW", "OPN1MSW", "OPN1SW", "OPN2SW", "RHO", "RCVRN"), ncol=3, order=T)
ggsave(file=file.path(outfinal, "P0_Complete_markers.png"), fp, width=15, height=10, unit="in", dpi=600)
saveRDS(seu, file.path(project, "P0_Complete_all_clusters.rds"))

# Subset photoreceptors
PR <- subset(seu, idents = c(1, 3, 6, 9, 13))
# Remove cells with very low mitochondrial gene content, likely without inner and outer segment
PR <- subset(PR, subset = percent.mt > 2)

# Normalization and clustering
PR <- PR %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  SCTransform(method = "glmGamPoi", verbose = TRUE, assay='RNA')

PR <- RunPCA(PR, verbose = FALSE, npcs = 100, features=features)
ElbowPlot(PR, ndims = 100, reduction = "pca")
set.seed(1)
PC_number <- 40
PR <- FindNeighbors(PR, dims = 1:PC_number)
PR <- RunUMAP(PR, dims = 1:PC_number)
res = 0.8
PR <- FindClusters(PR, resolution = res, random.seed = 1)
dp <- DimPlot(PR, label = TRUE, reduction = "umap") + NoLegend()
ggsave(file = file.path(outfinal, paste("P0_Complete_PR_clusters.png", sep="")), plot = dp, dpi = 600, width = 5, height = 5)
fp <- FeaturePlot(PR, features = c("OPN1LW", "OPN1MSW", "OPN1SW", "OPN2SW", "RHO", "RCVRN"), ncol=3, order=T)
ggsave(file=file.path(outfinal, "P0_Complete_PR_features.png"), fp, width=15, height=10, unit="in", dpi=600)
saveRDS(PR, file.path(project, "P0_Complete_PR_clusters.rds"))

# Check markers to identify putative doulbets
markers <- FindAllMarkers(PR, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(markers, file.path(outfinal, "P0_Complete_PR_all_markers.csv"))
write.csv(Top20, file.path(outfinal, "P0_Complete_PR_Top20_markers.csv"))

ids <- c("DC-A", "DC-P", "DC3", "Green", "Red", 
         "Rods", "BCdbl", "Cone-Rod", "MGdbl1", 
         "Violet", "Blue", "MGdbl2", "MGdbl3", 
         "ACdbl1", "ACdbl2")
names(ids) <- levels(PR)
PR <- RenameIdents(PR, ids)
saveRDS(PR, file.path(project, "P0_Complete_PR_clusters_labeled.rds"))

dp <- DimPlot(PR, label = TRUE, reduction = "umap", repel = T) + 
  NoLegend() + theme(axis.text.x=element_blank(), 
                     axis.ticks.x=element_blank(), 
                     axis.text.y=element_blank(), 
                     axis.ticks.y=element_blank()) 
ggsave(file = file.path(outfinal, paste("P0_Complete_PR_clusters_labeled.png", sep="")), plot = dp, dpi = 600, width = 5, height = 5)

# Remove doublets and re-clustering
Clean <- subset(PR, idents = c("DC-A", "Rods", "Green", "Red", "DC-P", "Violet", "Blue"))
Clean <- RunPCA(Clean, verbose = FALSE, npcs = 100, features=features)
ElbowPlot(Clean, ndims = 100, reduction = "pca")
set.seed(1)
PC_number <- 62
Clean <- FindNeighbors(Clean, dims = 1:PC_number)
Clean <- RunUMAP(Clean, dims = 1:PC_number)
res = 0.5
Clean <- FindClusters(Clean, resolution = res, random.seed = 1)
ids <- c("DC-A", "DC-P", "Green", "Red", 'Rods', "Violet", "Blue")
names(ids) <- levels(Clean)
Clean <- RenameIdents(Clean, ids)
cols <- c('Red'='red2', 'Green'=rgb(0, 176, 80, maxColorValue = 255), 'Blue'='blue', 'Violet'='purple', 
          'DC-P'='goldenrod1', 'DC-A'='goldenrod3', 'Rods'='dimgrey', 
          'DevB/V'='deeppink', 'DevCones'='yellow2', 'DevDC-P'='rosybrown1',
          'DevDC-A'='rosybrown3', 'DevRods'='azure3')
dp <- DimPlot(Clean, label = T, 
              reduction = "umap", 
              cols=cols[levels(Clean)],
              repel = T) + NoLegend() + theme(axis.text.x=element_blank(), 
                                              axis.ticks.x=element_blank(), 
                                              axis.text.y=element_blank(), 
                                              axis.ticks.y=element_blank())
ggsave(file = file.path(outfinal, paste("P0_Complete_Clean_clusters_labeled.png", sep="")), plot = dp, dpi = 600, width = 5, height = 5)
fp <- FeaturePlot(Clean, features = c("OPN1LW", "OPN1MSW", "OPN1SW", "OPN2SW", "RHO", "RCVRN"), ncol=3, order=T)
ggsave(file=file.path(outfinal, "P0_Complete_Clean_features.png"), fp, width=15, height=10, unit="in", dpi=600)
saveRDS(Clean, file.path(project, "P0_Complete_Clean_clusters.rds"))