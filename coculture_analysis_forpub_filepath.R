library(dplyr)
library(ggplot2)
library(spatstat.explore)
library(Seurat)
library(Matrix)
library(patchwork)
library(cowplot)
library(openxlsx)


##Read sc27s media control and create Seurat Object
sc27s.data <- Read10X(data.dir = "/filepath/sc27s_filtered_feature_bc_matrix")
sc27s <- CreateSeuratObject(counts = sc27s.data, project = "sc27s2", min.cells = 3, min.features = 200)
sc27s

##Read SC27s CM and create Seurat Object
sc27scm.data <- Read10X(data.dir = "/filepath/sc27scm_filtered_feature_bc_matrix")
sc27scm <- CreateSeuratObject(counts = sc27scm.data, project = "sc27scm2", min.cells = 3, min.features = 200)
sc27scm

##Read coculture and create Seurat Object
coculture.data <- Read10X(data.dir = "/filepath/coculture_filtered_feature_bc_matrix")
coculture <- CreateSeuratObject(counts = coculture.data, project = "coculture2", min.cells = 3, min.features = 200)
coculture

## Merge Seurat Objects (Coculture, SC27s, SC27s in CM)
cosc27s.merge <- merge(coculture, y = c(sc27s, sc27scm), add.cell.ids = c("coculture", "hNSPCs", "hNSPCsCM"), project = "covssc27s")
cosc27s.merge
##QC of merged data
cosc27s.merge[["percent.mt"]] <- PercentageFeatureSet(cosc27s.merge, pattern = "^MT-")
VlnPlot(cosc27s.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cosc27sclean <- subset(cosc27s.merge, subset = nFeature_RNA > 700 & nFeature_RNA < 7500 & percent.mt < 12)
cosc27sclean
VlnPlot(cosc27sclean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
##lognorm approach
cosc27snorm <- NormalizeData(cosc27sclean, normalization.method = "LogNormalize", scale.factor = 10000)
cosc27snorm <- FindVariableFeatures(cosc27snorm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(cosc27snorm), 10)
plot1 <- VariableFeaturePlot(cosc27snorm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes <- rownames(cosc27snorm)
cosc27scale <- ScaleData(cosc27snorm, features = all.genes)
cosc27spca <- RunPCA(cosc27scale, features = VariableFeatures(object = cosc27scale))
VizDimLoadings(cosc27spca, dims = 1:2, reduction = "pca")
DimPlot(cosc27spca, reduction = "pca")
DimHeatmap(cosc27spca, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(cosc27spca) #ID number of PCs to use
cosc27sumap <- FindNeighbors(cosc27spca, dims = 1:6)
cosc27sumap <- FindClusters(cosc27sumap, resolution = 0.1)
cosc27sumap <- RunUMAP(cosc27sumap, dims = 1:6)

#Visualization
DimPlot(cosc27sumap, reduction = "umap", split.by="orig.ident")
FeaturePlot(cosc27sumap, c("PECAM1"), order = TRUE)
FeaturePlot(cosc27sumap, c("MKI67"), order = TRUE)
FeaturePlot(cosc27sumap, c("GFAP"), order = TRUE)
FeaturePlot(cosc27sumap, c("SOX2"), order = TRUE)
FeaturePlot(cosc27sumap, c("DLL4"), order = TRUE)
FeaturePlot(cosc27sumap, c("JAG1"), order = TRUE)
FeaturePlot(cosc27sumap, c("JAG2"), order = TRUE)
FeaturePlot(cosc27sumap, c("NOTCH1"), order = TRUE)
FeaturePlot(cosc27sumap, c("NOTCH2"), order = TRUE)
FeaturePlot(cosc27sumap, c("NOTCH3"), order = TRUE)

# find markers for every cluster compared to all remaining cells, report only the positive ones
cosc27s.markers <- FindAllMarkers(object = cosc27sumap, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
cosc27s.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(cosc27sumap, features = top10$gene) + NoLegend()

saveRDS(cosc27sumap, file = "/filepath/cosc27sumap6pcsres0.1.rds")
saveRDS(cosc27sclean, file = "/filepath/cosc27sclean.rds")


#####HNSPC subset analysis####
#metadata annotation 
head(cosc27sumap@meta.data)
cosc27sumap@meta.data$cell_type <- as.character(cosc27sumap@meta.data$RNA_snn_res.0.1)
cosc27sumap@meta.data$cell_type[grep("0", cosc27sumap@meta.data$RNA_snn_res.0.1)] <- "hNSPCs"
cosc27sumap@meta.data$cell_type[grep("1", cosc27sumap@meta.data$RNA_snn_res.0.1)] <- "hNSPCs"
cosc27sumap@meta.data$cell_type[grep("3", cosc27sumap@meta.data$RNA_snn_res.0.1)] <- "hNSPCs"
cosc27sumap@meta.data$cell_type[grep("4", cosc27sumap@meta.data$RNA_snn_res.0.1)] <- "hNSPCs"
cosc27sumap@meta.data$cell_type[grep("2", cosc27sumap@meta.data$RNA_snn_res.0.1)] <- "hECs"
cosc27sumap@meta.data$cell_type[grep("5", cosc27sumap@meta.data$RNA_snn_res.0.1)] <- "hECs"

#neural stem cell subset analysis
Idents(cosc27sumap) <- "cell_type"
sc27s.subset <- subset(cosc27sumap, idents = "hNSPCs")
sc27s.subset
DimPlot(sc27s.subset)
saveRDS(sc27s.subset, file = "/filepath/sc27s_subset.rds")

#re-do standard workflow 
sc27s.subset <- NormalizeData(sc27s.subset, normalization.method = "LogNormalize", scale.factor = 10000)
sc27s.subset <- FindVariableFeatures(sc27s.subset, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sc27s.subset)
sc27s.subset <- ScaleData(sc27s.subset, features = all.genes)
sc27s.subset <- RunPCA(sc27s.subset, features = VariableFeatures(object = sc27s.subset))
DimHeatmap(sc27s.subset, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(sc27s.subset) #ID number of PCs to use
sc27sumap.subset <- FindNeighbors(sc27s.subset, dims = 1:8)
sc27sumap.subset <- FindClusters(sc27sumap.subset, resolution = 0.2)
sc27sumap.subset <- RunUMAP(sc27sumap.subset, dims = 1:8)
DimPlot(sc27sumap.subset, reduction = "umap", split.by="orig.ident")
saveRDS(sc27sumap.subset, file = "/filepath/sc27sumap_subset.rds")
sc27sumap.subset <- readRDS(file = "/filepath/sc27sumap_subset.rds")

sc27subset.markers <- FindAllMarkers(object = sc27sumap.subset, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
sc27subset.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(sc27sumap.subset, features = top10$gene) + NoLegend()
write.xlsx(sc27subset.markers, row.names=TRUE, "/filepath/all_cluster_markers.xlsx")

VlnPlot(sc27sumap.subset, features = c("GFAP", "SOX2"))


#get number of cells per cluster for each condition
library(data.table)
library(magrittr)

md3 <- sc27sumap.subset@meta.data %>% as.data.table
md3[, .N, by = c("orig.ident", "seurat_clusters")]
md3[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N")

#differential expression of cluster 0 in coculture2 vs sc27s
Idents(sc27sumap.subset) <- sc27sumap.subset$seurat_clusters
DefaultAssay(sc27sumap.subset) <- "RNA"
cluster0covsmono<- FindMarkers(sc27sumap.subset, ident.1 = "coculture2", ident.2 = "sc27s2", verbose = TRUE, group.by="orig.ident", subset.ident = "0", logfc.threshold = log(2))
write.xlsx(cluster0covsmono, row.names=TRUE, "/filepath/cluster0_covsmono_markers.xlsx")

#annotated sc27s subset data
sc27subset.annotated <- sc27sumap.subset
head(sc27subset.annotated@meta.data)
sc27subset.annotated@meta.data$cell_type <- as.character(sc27subset.annotated@meta.data$RNA_snn_res.0.2)
sc27subset.annotated@meta.data$cell_type[grep("0", sc27subset.annotated@meta.data$RNA_snn_res.0.2)] <- "Neural stem cells (Type B)"
sc27subset.annotated@meta.data$cell_type[grep("1", sc27subset.annotated@meta.data$RNA_snn_res.0.2)] <- "Astrocyte progenitors"
sc27subset.annotated@meta.data$cell_type[grep("2", sc27subset.annotated@meta.data$RNA_snn_res.0.2)] <- "Proliferating progenitors"
sc27subset.annotated@meta.data$cell_type[grep("3", sc27subset.annotated@meta.data$RNA_snn_res.0.2)] <- "Neural stem cells"
sc27subset.annotated@meta.data$cell_type[grep("4", sc27subset.annotated@meta.data$RNA_snn_res.0.2)] <- "Proliferating progenitors"
sc27subset.annotated@meta.data$cell_type[grep("5", sc27subset.annotated@meta.data$RNA_snn_res.0.2)] <- "Transitional cells"
sc27subset.annotated@meta.data$cell_type[grep("6", sc27subset.annotated@meta.data$RNA_snn_res.0.2)] <- "Early neuron progenitors"
saveRDS(sc27subset.annotated, file = "/filepath/sc27subset_annotated.rds")

Idents(sc27subset.annotated) <- "cell_type"
DimPlot(sc27subset.annotated, reduction = "umap", split.by="orig.ident", pt.size = 2)

#Dotplot of hNSPC subset data
dp_genes <- c("GFAP", "PROM1", "HES1", "NGFR", "APOE", "SLC1A3", "S100A6", "SOX2", "MKI67", "TOP2A", "FABP7", "AQP4", "SPARCL1", "ASCL1", "DCLK2", "MAP2")
heshey_genes <- c("HES1", "HES4", "HES5", "HES7", "HEY1", "HEY2", "HEYL")
DotPlot(object = sc27subset.annotated, features = dp_genes, cols = "RdYlBu") + RotatedAxis()
DotPlot(object = sc27subset.annotated, features = heshey_genes, cols = "RdYlBu", split.by = "orig.ident") + RotatedAxis()


Idents(sc27subset.annotated) <- "cell_type"
hNSCs.typeB <- subset(sc27subset.annotated, idents = "Neural stem cells (Type B)")
saveRDS(hNSCs.typeB, file = "/filepath/hNSCs_typeB.rds")

hNSCs.typeB.coculture <- subset(sc27subset.annotated, idents = "Neural stem cells (Type B)", subset = orig.ident == "coculture2")
saveRDS(hNSCs.typeB.coculture, file = "/filepath/hNSCs_typeB_coculture.rds")

hNSCs.typeB.monoculture <- subset(sc27subset.annotated, idents = "Neural stem cells (Type B)", subset = orig.ident == "sc27s2")
saveRDS(hNSCs.typeB.monoculture, file = "/filepath/hNSCs_typeB_monoculture.rds")

FeaturePlot(sc27subset.annotated, c("GFAP"), order = TRUE)

#Find markers for every cluster compared to all remaining cells, report only the positive ones
sc27sannotated.markers <- FindAllMarkers(object = sc27subset.annotated, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
sc27sannotated.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(sc27subset.annotated, features = top10$gene, angle=45, size=3, raster = FALSE) + NoLegend()


#compare gene expression co vs monoculture type b cells
typeb_monovco_merge <- merge(hNSCs.typeB.monoculture, y = c(hNSCs.typeB.coculture), add.cell.ids = c("monoTypeB", "coTypeB"), project = "typeb_monovco")
DefaultAssay(typeb_monovco_merge) <- "RNA"
Idents(typeb_monovco_merge) <- "orig.ident"
avg.exp.typeb <- as.data.frame(log1p(AverageExpression(typeb_monovco_merge, verbose = FALSE)$RNA))
avg.exp.typeb$gene <- rownames(avg.exp.typeb)
write.xlsx(avg.exp.typeb, rowNnames=TRUE, "/filepath/avg_typeB_covmono_markers.xlsx")

genes.to.label = c("HSPG2", "COL1A1", "HES4", "HEY1", "NES", "TNC", "IGFBP7")

p1 <- ggplot(avg.exp.typeb, aes(sc27s2, coculture2)) + geom_point(colour = "blue", size=1) + ggtitle("monoculture_typeB vs coculture_typeB")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1

# Cell cycle analysis
# Read in the expression matrix The first row is a header row, the first
# Column is rownames
exp.mat <- read.table(file = "/filepath/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt", 
                      header = TRUE, as.is = TRUE, row.names = 1)

cc.genes <- readLines(con = "/filepath/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")

# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
#assign cell cycle scores 
sc27sumap.subset <- CellCycleScoring(object = sc27sumap.subset, s.features = s.genes, g2m.features = g2m.genes, 
                                     set.ident = TRUE)
# Visualize the distribution of cell cycle markers across
RidgePlot(object = sc27sumap.subset, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
FeaturePlot(sc27sumap.subset, c("S.Score"), 
            cols=c("light gray","dark blue"), 
            pt.size=1)

FeaturePlot(sc27sumap.subset, c("G2M.Score"), 
            cols=c("light gray","dark blue"), 
            pt.size=1) 
DimPlot(sc27sumap.subset, reduction = "umap", group.by = "Phase")
DimPlot(sc27sumap.subset, reduction = "umap", group.by = "seurat_clusters")



#####Type B subset analysis#####
typebsubset <- hNSCs.typeB

#re-do standard workflow 
typebsubset <- NormalizeData(typebsubset, normalization.method = "LogNormalize", scale.factor = 10000)
typebsubset <- FindVariableFeatures(typebsubset, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(typebsubset)
typebsubset <- ScaleData(typebsubset, features = all.genes)
typebsubset <- RunPCA(typebsubset, features = VariableFeatures(object = typebsubset))
DimHeatmap(typebsubset, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(typebsubset) #ID number of PCs to use
typebumap <- FindNeighbors(typebsubset, dims = 1:8)
typebumap <- FindClusters(typebumap, resolution = 0.1)
typebumap <- RunUMAP(typebumap, dims = 1:8)
DimPlot(typebumap, reduction = "umap")
DimPlot(typebumap, reduction = "umap", pt.size = 1, split.by = "orig.ident")
FeaturePlot(typebumap, features = c("NES"), order = TRUE, keep.scale = "feature", pt.size= 0.5) 

# Visualize co-expression of two features simultaneously
FeaturePlot(typebumap, features = c("S100A6", "GFAP"), cols = c("black", "red", "green"), blend = TRUE, order = TRUE, pt.size= 0.5, blend.threshold = 0.5, min.cutoff = 'q10', max.cutoff = 'q90')
FeaturePlot(typebumap, c("S100A6"), order = TRUE, min.cutoff = 'q10', max.cutoff = 'q90')
FeaturePlot(typebumap, c("STXBP6"), order = TRUE, min.cutoff = 'q10', max.cutoff = 'q90')
FeaturePlot(typebumap, c("NEFL"), order = TRUE, min.cutoff = 'q10', max.cutoff = 'q90')
FeaturePlot(typebumap, c("NGFR"), order = TRUE, min.cutoff = 'q10', max.cutoff = 'q90')
FeaturePlot(typebumap, c("COL1A2"), order = TRUE, min.cutoff = 'q10', max.cutoff = 'q90')
FeaturePlot(typebumap, c("NEAT1"), order = TRUE, min.cutoff = 'q10', max.cutoff = 'q90')
FeaturePlot(typebumap, c("IER3"), order = TRUE, min.cutoff = 'q10', max.cutoff = 'q90')
FeaturePlot(typebumap, c("PCDH9"), order = TRUE, min.cutoff = 'q10', max.cutoff = 'q90')

saveRDS(typebumap, file = "/filepath/typeB_subset_umap.rds")

#Getting percentage of cells expressing a gene of interest
my_genes <- c("GFAP", "S100A6", "APOE", "NGFR", "NTRK2", "MT3", "PROM1", "SOX2", "NES", "SOX3", "TNC", "ANXA2", "COL1A2", "EMP2", "IER3", "PCDH9", "SCRG1", "NEAT1")
S100A6 <- c("S100A6")
exp <- FetchData(typebumap, S100A6)
matrix <- as.matrix(colMeans(exp  > 0))*100

GFAP <- c("GFAP")
exp2 <- FetchData(typebumap, GFAP)
matrix2 <- as.matrix(colMeans(exp2  > 0))*100

#Co-expression data
length(WhichCells(object = typebumap, expression = S100A6 > 0 & GFAP > 0))/nrow(typebumap@meta.data)*100

#Calculating the percent of cells expressing a gene that also expresses GFAP
my_gene <- c("S100A6")

#number of cells expressing a particular gene
sum(GetAssayData(object = typebumap, slot = "data")[my_gene,]>0)

#number of cells co-expressing genes
length(WhichCells(object = typebumap, expression = S100A6 > 0 & GFAP > 0))

