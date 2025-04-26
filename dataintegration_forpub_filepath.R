library(dplyr)
library(ggplot2)
library(spatstat.explore)
library(Seurat)
library(Matrix)
library(patchwork)
library(cowplot)
library(openxlsx) 
library(data.table)


##########Human type B cell data-SC27s only###############
resc27s.data <- Read10X(data.dir = "/filepath/sc27s_filtered_feature_bc_matrix_renamed")
resc27s <- CreateSeuratObject(counts = resc27s.data, project = "sc27s2", min.cells = 3, min.features = 200)
resc27s
##QC
resc27s[["percent.mt"]] <- PercentageFeatureSet(resc27s, pattern = "^MT-")
VlnPlot(resc27s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
resc27sclean <- subset(resc27s, subset = nFeature_RNA > 700 & nFeature_RNA < 7500 & percent.mt < 12)
resc27sclean
VlnPlot(resc27sclean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
## Normalize and cluster each sample to determine appropriate concatenation approach
resc27snorm <- NormalizeData(resc27sclean, normalization.method = "LogNormalize", scale.factor = 10000)
resc27svar <- FindVariableFeatures(resc27snorm, selection.method = "vst", nfeatures = 2000)
reall.genes <- rownames(resc27svar)
resc27sscale <- ScaleData(resc27svar, features = reall.genes)
resc27spca <- RunPCA(resc27sscale, features = VariableFeatures(object = resc27sscale))
VizDimLoadings(resc27spca, dims = 1:2, reduction = "pca")
DimPlot(resc27spca, reduction = "pca")
DimHeatmap(resc27spca, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(resc27spca)
resc27sumap<- FindNeighbors(resc27spca, dims = 1:8)
resc27sumap <- FindClusters(resc27sumap, resolution = 0.2)
resc27sumap <- RunUMAP(resc27sumap, dims = 1:8)
DimPlot(resc27sumap, reduction = "umap")
saveRDS(resc27sumap, file = "/filepath/resc27sumap.rds")

resc27.markers <- FindAllMarkers(object = resc27sumap, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
resc27.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(resc27sumap, features = top10$gene) 
write.xlsx(resc27.markers, rowNames=TRUE, "/filepath/all_cluster_markers.xlsx")
# Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(resc27sumap, features = "GFAP")

head(resc27sumap@meta.data)
resc27sumap@meta.data$cell_type <- as.character(resc27sumap@meta.data$RNA_snn_res.0.2)
resc27sumap@meta.data$cell_type[grep("0", resc27sumap@meta.data$RNA_snn_res.0.2)] <- "Neural stem cells (Type B)"
resc27sumap@meta.data$cell_type[grep("1", resc27sumap@meta.data$RNA_snn_res.0.2)] <- "Astrocyte progenitors"
resc27sumap@meta.data$cell_type[grep("2", resc27sumap@meta.data$RNA_snn_res.0.2)] <- "Neural stem cells"
resc27sumap@meta.data$cell_type[grep("3", resc27sumap@meta.data$RNA_snn_res.0.2)] <- "Proliferating progenitors"
resc27sumap@meta.data$cell_type[grep("4", resc27sumap@meta.data$RNA_snn_res.0.2)] <- "Transitional cells"
resc27sumap@meta.data$cell_type[grep("5", resc27sumap@meta.data$RNA_snn_res.0.2)] <- "Proliferating progenitors"
resc27sumap@meta.data$cell_type[grep("6", resc27sumap@meta.data$RNA_snn_res.0.2)] <- "Early neuron progenitors"
saveRDS(resc27sumap, file = "/filepath/resc27sumapannotated.rds")

Idents(resc27sumap) <- "cell_type"
htypeb.b <- subset(resc27sumap, idents = "Neural stem cells (Type B)")
saveRDS(htypeb.b, file = "/filepath/hNSCs_b.rds")




##### Renamed Human RGCs-Nowakowski data #############
mat <- fread("/filepath/kriegstein_r/k_r.tsv.gz")
meta <- read.table("/filepath/kriegstein_r/meta.kriegstein.tsv", header=T, sep="\t", as.is=T, row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
hRGC <- CreateSeuratObject(counts = mat, project = "hRGCs", meta.data=meta, min.cells = 3, min.features = 200)
hRGC
saveRDS(hRGC, file = "/filepath/hRGCs_kriegstein_renamed.rds")

#Once the object is obtained, you want to process them and apply QC steps and normalization steps
##QC of humanRG data
hRGC[["percent.mt"]] <- PercentageFeatureSet(hRGC, pattern = "^MT.")
VlnPlot(hRGC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hRGk.clean <- subset(hRGC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 8) 
hRGk.clean
VlnPlot(hRGk.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Normalization and scaling
hRGk.norm <- NormalizeData(hRGk.clean, normalization.method = "LogNormalize", scale.factor = 10000)
hRGk.var <- FindVariableFeatures(hRGk.norm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(hRGk.var), 10)
plot1 <- VariableFeaturePlot(hRGk.var)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes <- rownames(hRGk.var)
hRGk.var <- ScaleData(hRGk.var, features = all.genes)
#Dimensional reduction
hRGk.var <- RunPCA(hRGk.var, features = VariableFeatures(object = hRGk.var))
VizDimLoadings(hRGk.var, dims = 1:2, reduction = "pca")
DimPlot(hRGk.var, reduction = "pca")
DimHeatmap(hRGk.var, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(hRGk.var) #ID number of PCs to use
#Creating clusters
hRGk.umap2 <- FindNeighbors(hRGk.var, dims = 1:8)
hRGk.umap2 <- FindClusters(hRGk.umap2, resolution = 0.6)
hRGk.umap2 <- RunUMAP(hRGk.umap2, dims = 1:8)
DimPlot(hRGk.umap2, reduction = "umap")

#Another way is to subset by the meta data. sometimes the original group provides names they added for the 
#clusters they created. create the identity to use as a column stored in their meta data
Idents(hRGk.umap2) <- "WGCNAcluster"
hRGsubset.wgcna <- subset(hRGk.umap2, idents = c("tRG", "vRG", "oRG", "RG-early"))
saveRDS(hRGsubset.wgcna, file = "/filepath/hRGksubset_wgcna_renamed.rds")

head(hRGsubset.wgcna@meta.data) 
hRGsubset.wgcna@meta.data$cell_type <- as.character(hRGsubset.wgcna@meta.data$WGCNAcluster)
hRGsubset.wgcna@meta.data$cell_type[grep("tRG", hRGsubset.wgcna@meta.data$WGCNAcluster)] <- "hRGCs"
hRGsubset.wgcna@meta.data$cell_type[grep("vRG", hRGsubset.wgcna@meta.data$WGCNAcluster)] <- "hRGCs"
hRGsubset.wgcna@meta.data$cell_type[grep("oRG", hRGsubset.wgcna@meta.data$WGCNAcluster)] <- "hRGCs"
hRGsubset.wgcna@meta.data$cell_type[grep("RG-early", hRGsubset.wgcna@meta.data$WGCNAcluster)] <- "hRGCs"
saveRDS(hRGsubset.wgcna, file = "/filepath/hRGksubset_wgcna_annotated.rds")




######Human Adult Type B ########
#Using rawdata
library(janitor)
counts_files<- list.files(path = "/filepath/GSE248995_RAW/", full.names = TRUE, pattern = "*normalized_r.tsv")
for (rr in 1:length(counts_files)){
  
  y<-read.table(counts_files[rr], sep="\t",row.names=1);
  y<- y %>% row_to_names(row_number = 1)
  filename<-sprintf("sobj%d",rr)
  
  sobj<-CreateSeuratObject(counts=y, meta.data=y[,-1], min.cells = 3, min.features = 200)
  assign(filename,sobj)
}

hSVZ <- merge(sobj1, y = c(sobj2, sobj3, sobj4, sobj5, sobj6, sobj7, sobj8, sobj9, sobj10, sobj11, sobj12, sobj13, sobj14, sobj15), add.cell.ids = c("GSM7924388", "GSM7924389", "GSM7924390", "GSM7924391", "GSM7924392", "GSM7924393", "GSM7924394", "GSM7924395", "GSM7924396", "GSM7924397", "GSM7924398", "GSM7924399", "GSM7924400", "GSM7924401", "GSM7924402"), project = "hSVZ_adult_merge")
hSVZ
head(colnames(hSVZ))
tail(colnames(hSVZ))
unique(sapply(X = strsplit(colnames(hSVZ), split = "_"), FUN = "[", 1))
saveRDS(hSVZ, file = "/filepath/hSVZ_Baig.rds")


#QC
hSVZ[["percent.mt"]] <- PercentageFeatureSet(hSVZ, pattern = "^MT-")
VlnPlot(hSVZ, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
hSVZ.clean <- subset(hSVZ, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 8)
VlnPlot(hSVZ.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalization and dimensional reduction
hSVZ.norm <- NormalizeData(hSVZ.clean, normalization.method = "LogNormalize", scale.factor = 10000)
hSVZ.norm <- FindVariableFeatures(hSVZ.norm, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(hSVZ.norm)

hSVZ.scale <- ScaleData(hSVZ.norm, features = all.genes)
hSVZ.pca <- RunPCA(hSVZ.scale, features = VariableFeatures(object = hSVZ.scale))
ElbowPlot(hSVZ.pca)
#Clustering k.param original is 10
hSVZ.umap <- FindNeighbors(hSVZ.pca, dims = 1:35, k.param = 10, reduction = "pca")
hSVZ.umap <- FindClusters(hSVZ.umap, resolution = 0.6)
hSVZ.umap <- RunUMAP(hSVZ.umap, dims = 1:35)
DimPlot(hSVZ.umap, reduction = "umap")

#find markers
hSVZ.umap.joined<- JoinLayers(hSVZ.umap)
hSVZ.markers <- FindAllMarkers(hSVZ.umap.joined, only.pos = TRUE)
hSVZ.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(hSVZ.umap.joined, features = top10$gene) + NoLegend()
write.xlsx(hSVZ.markers, rowNnames=TRUE, "/filepath//Baig_35pcres0.6_markers.xlsx")

saveRDS(hSVZ.umap, file = "/filepath/hsVZ_Baig_umap_35pcsres0.6.rds")

#metadata annotation 
head(hSVZ.umap.joined@meta.data)
hSVZ.umap.joined@meta.data$cell_type <- as.character(hSVZ.umap.joined@meta.data$RNA_snn_res.0.6)
hSVZ.umap.joined@meta.data$cell_type[grep("7", hSVZ.umap.joined@meta.data$RNA_snn_res.0.6)] <- "Neural stem cells (SVZ)"


#neural stem cell subset analysis
Idents(hSVZ.umap.joined) <- "seurat_clusters"
htypeb.a.r <- subset(hSVZ.umap.joined, idents = "7")
htypeb.a.r
saveRDS(htypeb.a.r, file = "/filepath/hSVZ_baig_typeb_subset.rds")


####Merging datasets#########
#merge files
hNSCs.merged <- merge(hRGsubset.wgcna, y = c(htypeb.a.r, htypeb.b), add.cell.ids = c("hRGsubset_wgcna", "hSVZtypeBaig", "hTypeB_monoculture"), project = "NSCs_Kriegstein_Baig")
hNSCs.merged
head(colnames(hNSCs.merged))
tail(colnames(hNSCs.merged))
unique(sapply(X = strsplit(colnames(hNSCs.merged), split = "_"), FUN = "[", 1))
saveRDS(hNSCs.merged, file = "/filepath/hNSCs_kriegstein_baig_merged.rds")

Idents(hNSCs.merged) <- "cell_type"
hNSCs.merged <- PercentageFeatureSet(hNSCs.merged, pattern = "^MT-", col.name = "percent.mt")

VlnPlot(hNSCs.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

hNSCs.merged <- subset(hNSCs.merged, subset = nFeature_RNA > 200 
                    & nFeature_RNA < 6000 & nCount_RNA < 1000000 & percent.mt < 8)


###CCA integration
hNSCs_merged_joined<- JoinLayers(hNSCs.merged)
hNSCs_merged_joined[["RNA"]] <- split(hNSCs_merged_joined[["RNA"]], f = hNSCs_merged_joined$cell_type, layers = c('counts', 'data', 'scale.data'))
# run standard analysis workflow
hNSCs_merged_joined <- NormalizeData(hNSCs_merged_joined)
hNSCs_merged_joined <- FindVariableFeatures(hNSCs_merged_joined)
hNSCs_merged_joined <- ScaleData(hNSCs_merged_joined)
hNSCs_merged_joined <- RunPCA(hNSCs_merged_joined)
hNSCs_merged_joined_pca <- FindNeighbors(hNSCs_merged_joined, dims = 1:30, reduction = "pca")
hNSCs_merged_joined_pca <- FindClusters(hNSCs_merged_joined_pca, resolution = 2, cluster.name = "unintegrated_clusters")
hNSCs_merged_joined_pca <- RunUMAP(hNSCs_merged_joined_pca, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(hNSCs_merged_joined_pca, reduction = "umap.unintegrated", group.by = c("cell_type"))

hNSCs_merged_joined_pca <- IntegrateLayers(object = hNSCs_merged_joined_pca, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                           verbose = FALSE)
# re-join layers after integration
hNSCs_merged_joined_pca[["RNA"]] <- JoinLayers(hNSCs_merged_joined_pca[["RNA"]])
hNSCs_merged_joined_pca <- FindNeighbors(hNSCs_merged_joined_pca, reduction = "integrated.cca", dims = 1:20)
hNSCs_merged_joined_pca <- FindClusters(hNSCs_merged_joined_pca, resolution = 2)
hNSCs_merged_joined_pca <- RunUMAP(hNSCs_merged_joined_pca, dims = 1:20, reduction = "integrated.cca")
saveRDS(hNSCs_merged_joined_pca, file = "/filepath/hNSCs_kriegstein_baig_20pcsres2_intg_cca.rds")

# Visualization
DimPlot(hNSCs_merged_joined_pca, reduction = "umap", group.by = c("cell_type"), pt.size=2, shuffle = TRUE, cols = c("#00A9FF", "#0CB702", "#FB7660"))
FeaturePlot(hNSCs_merged_joined_pca, reduction = "umap", c("GFAP"), split.by = "cell_type", order = TRUE, pt.size=2) & theme(legend.position = "right")


hNSCs_merged_joined2<- JoinLayers(hNSCs_merged_joined_pca)
Idents(hNSCs_merged_joined2) <- "cell_type"
baig.kriegstein <- FindMarkers(hNSCs_merged_joined2, ident.1 = "Neural stem cells (SVZ)", ident.2 = "hRGCs", verbose = FALSE)
head(baig.kriegstein, n = 15)
write.xlsx(baig.kriegstein, rowNames=TRUE, "/filepath/baig_v_kriegstein_markers.xlsx")
saveRDS(hNSCs_merged_joined_pca, file = "/filepath/hNSCs_kriegstein_baig_20pcsres2_intg_cca_joined.rds")


##########Mouse adult SVZ########
mSVZ_cebrian_1 <- Read10X(data.dir = "/filepath/mSVZ_cebrian_1_r")
mSVZ_cebrian_1 <- CreateSeuratObject(counts = mSVZ_cebrian_1$`Gene Expression`, project = "mSVZ", min.cells = 3, min.features = 200)
mSVZ_cebrian_2 <- Read10X(data.dir = "/filepath/mSVZ_cebrian_2_r")
mSVZ_cebrian_2 <- CreateSeuratObject(counts = mSVZ_cebrian_2$`Gene Expression`, project = "mSVZ", min.cells = 3, min.features = 200)


mSVZ <- merge(mSVZ_cebrian_1, y = c(mSVZ_cebrian_2), add.cell.ids = c("mSVZ_cebrian_1", "mSVZ_cebrian_2"), project = "mouseSVZ")
saveRDS(mSVZ, file = "/filepath/mouseSVZ_cebrian_merged.rds")

#Once the object is obtained, you want to process them and apply QC steps and normalization steps
##QC of humanRG data
mSVZ[["percent.mt"]] <- PercentageFeatureSet(mSVZ, pattern = "^Mt.")
VlnPlot(mSVZ, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mSVZ.clean <- subset(mSVZ, subset = nFeature_RNA > 600 & nFeature_RNA < 7500 & percent.mt < 5) 
mSVZ.clean
VlnPlot(mSVZ.clean, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
##lognorm approach
mSVZ.norm <- NormalizeData(mSVZ.clean, normalization.method = "LogNormalize", scale.factor = 10000)

#To pick the number of features, test of diff #s and see how they plot. If the selected features make up
###the scattered cells above then it is a good # to pick 
mSVZ.var <- FindVariableFeatures(mSVZ.norm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mSVZ.var), 10)
plot1 <- VariableFeaturePlot(mSVZ.var)
plot1

###Now scale the data and select top features to take into account for clustering and plotting
all.genes <- rownames(mSVZ.var)
mSVZ.scale <- ScaleData(mSVZ.var, features = all.genes)
mSVZ.pca <- RunPCA(mSVZ.scale, features = VariableFeatures(object = mSVZ.scale))
VizDimLoadings(mSVZ.var, dims = 1:2, reduction = "pca")
DimPlot(mSVZ.var, reduction = "pca")

#Creating clusters
DimHeatmap(mSVZ.var, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(mSVZ.pca) #ID number of PCs to use
mSVZ.umap <- FindNeighbors(mSVZ.pca, dims = 1:10)
mSVZ.umap <- FindClusters(mSVZ.umap, resolution = 1.2)
mSVZ.umap <- RunUMAP(mSVZ.umap, dims = 1:10)
DimPlot(mSVZ.umap, reduction = "umap", label = TRUE)
#look at nFeature and nCount for each cluster
FeaturePlot(mSVZ.umap, c("nFeature_RNA", "nCount_RNA"))
FeaturePlot(mSVZ.umap, c("GPX3"), order=TRUE)
saveRDS(mSVZ.umap, file = "/filepath/mSVZ_umap.rds")

#find markers
mSVZ.umap.joined<- JoinLayers(mSVZ.umap)
mSVZ.markers <- FindAllMarkers(mSVZ.umap.joined, only.pos = TRUE)
mSVZ.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(mSVZ.umap.joined, features = top10$gene) + NoLegend()
write.xlsx(mSVZ.markers, rowNnames=TRUE, "/filepath/mSVZ_allmarkers.xlsx")

head(mSVZ.umap.joined@meta.data)
mSVZ.umap.joined@meta.data$cell_type <- as.character(mSVZ.umap.joined@meta.data$RNA_snn_res.1.2)
mSVZ.umap.joined@meta.data$cell_type[grep("1", mSVZ.umap.joined@meta.data$RNA_snn_res.1.2)] <- "quiescent_TypeB"
mSVZ.umap.joined@meta.data$cell_type[grep("24", mSVZ.umap.joined@meta.data$RNA_snn_res.1.2)] <- "quiescent_TypeB"
mSVZ.umap.joined@meta.data$cell_type[grep("15", mSVZ.umap.joined@meta.data$RNA_snn_res.1.2)] <- "Astrocytes"

#Saving cell type clusters as seurat object 
Idents(mSVZ.umap.joined) <- "seurat_clusters"
quiescent_TypeB <- subset(mSVZ.umap.joined, idents = c("1", "24"))
quiescent_TypeB
DimPlot(quiescent_TypeB)
saveRDS(quiescent_TypeB, file = "/filepath/mSVZ_quiescent_typeB.rds")



#########Mouse RGCs############
#Analyzing rodent RG data 
mat2 <- fread("/filepath/mRG_li_r/mRG_Li_exprMatrix_r.tsv.gz")
meta2 <- read.table("/filepath/mRG_li_r/mRG_Li_meta.tsv.gz", header=T, sep="\t", as.is=T, row.names=1)
genes = mat2[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat2 = data.frame(mat2[,-1], row.names=genes)
rodentRGli <- CreateSeuratObject(counts = mat2, project = "rodentRGli", meta.data=meta2, min.cells = 3, min.features = 200)
saveRDS(rodentRGli, file = "/filepath/rodentRGli.rds")

#Once the object is obtained, you want to process them and apply QC steps and normalization steps
##QC of rodentRG data
rodentRGli[["percent.mt"]] <- PercentageFeatureSet(rodentRGli, pattern = "^Mt.")
VlnPlot(rodentRGli, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
rRGli.clean <- subset(rodentRGli, subset = nFeature_RNA > 750 & nFeature_RNA < 3750 & percent.mt < 5) 
rRGli.clean
VlnPlot(rRGli.clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
##lognorm approach
rRGli.norm <- NormalizeData(rRGli.clean, normalization.method = "LogNormalize", scale.factor = 10000)
rRGli.var <- FindVariableFeatures(rRGli.norm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(rRGli.var), 10)
plot1 <- VariableFeaturePlot(rRGli.var)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

###Now scale the data and select top features to take into account for clustering and plotting
all.genes <- rownames(rRGli.var)
rRGli.var <- ScaleData(rRGli.var, features = all.genes)
rRGli.var <- RunPCA(rRGli.var, features = VariableFeatures(object = rRGli.var))
VizDimLoadings(rRGli.var, dims = 1:2, reduction = "pca")
DimPlot(rRGli.var, reduction = "pca")

#Creating clusters
DimHeatmap(rRGli.var, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(rRGli.var) #ID number of PCs to use
rRGli.umap <- FindNeighbors(rRGli.var, dims = 1:8)
rRGli.umap <- FindClusters(rRGli.umap, resolution = 0.1)
rRGli.umap <- RunUMAP(rRGli.umap, dims = 1:8)
DimPlot(rRGli.umap, reduction = "umap")

#Another way is to subset by the meta data. sometimes the original group provides names they added for the 
#clusters they created. create the identity to use as a column stored in their meta data
Idents(rRGli.umap) <- "clusters"
rRGsubset.wgcna <- subset(rRGli.umap, idents = c("mRGC1", "mRGC2", "mRGC3"))
saveRDS(rRGsubset.wgcna, file = "/filepath/rRGsubset.rds")
saveRDS(rRGli.umap, file = "/filepath/rRG_li_umap.rds")



#####Merge Mouse Data########
###Once you have all you datasets and specific populations you want to compare, merge them into a single file. List
#the datasets you want to merge and make sure to write the name correctly.
mouse_merge2 <- merge(rRGsubset.wgcna, y = c(quiescent_TypeB), add.cell.ids = c("mRGC", "mTypeB"), project = "mouse_merge2")
mouse_merge2
head(colnames(mouse_merge2))
tail(colnames(mouse_merge2))
unique(sapply(X = strsplit(colnames(mouse_merge2), split = "_"), FUN = "[", 1))
saveRDS(mouse_merge2, file = "/filepath/mouse_merged.rds")


##New SCT
# split datasets and process without integration
mouse_merged_joined<- JoinLayers(mouse_merge2)
mouse_merged_joined[["RNA"]] <- split(mouse_merged_joined[["RNA"]], f = mouse_merged_joined$orig.ident, layers = c('counts', 'data', 'scale.data'))
mouse_merged_joined <- SCTransform(mouse_merged_joined)
mouse_merged_joined <- RunPCA(mouse_merged_joined)
mouse_merged_joined <- RunUMAP(mouse_merged_joined, dims = 1:30)
DimPlot(mouse_merged_joined, reduction = "umap", group.by = "orig.ident", pt.size=1)


# integrate datasets
intg <- IntegrateLayers(object = mouse_merged_joined, method = CCAIntegration, normalization.method = "SCT", verbose = F)
intg <- FindNeighbors(intg, reduction = "integrated.dr", dims = 1:30)
intg <- FindClusters(intg, resolution = 0.6)
intg <- RunUMAP(intg, dims = 1:30, reduction = "integrated.dr")
DimPlot(intg, reduction = "umap", group.by = c("orig.ident"), pt.size=1.5, shuffle = TRUE, cols = c("#00BFC4", "#C77CFF"))
saveRDS(intg, file = "/filepath/mouse_merge_30pcsres0.6_intgcca.rds")


# For performing differential expression after integration, we switch back to the original data
theme_set(theme_cowplot())

DefaultAssay(intg) <- "RNA"
Idents(intg) <- "orig.ident"
avg.mnscs <- as.data.frame(log1p(AverageExpression(intg, verbose = FALSE)$RNA))
avg.mnscs$gene <- rownames(avg.mnscs)
write.xlsx(avg.mnscs, rowNnames=TRUE, "/filepath/avg_mouse_nsc_markers.xlsx")
genes.to.label = c("Gfap", "S100a6", "Malat1", "Apoe", "Nrxn3", "Ntrk2", "Btg1", "Cebpd", "Tiam2", "Mt3")
p1 <- ggplot(avg.mnscs, aes(ddSEQ, mSVZ)) + geom_point(colour = "blue", size=1) + ggtitle("mTypeB vs mRG")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge=0, ynudge=0)
p1


Idents(intg) <- "orig.ident"
mouse_int_joined <- JoinLayers(intg)
mouse_nscs <- FindMarkers(mouse_int_joined, ident.1 = "mSVZ", ident.2 = "ddSEQ", verbose = FALSE)
write.xlsx(mouse_nscs, rowNames=TRUE, "/filepath/typeB_RGC_markers.xlsx")

FeaturePlot(intg, c("Apoe"), split.by = "orig.ident", order = TRUE, pt.size=1.5) 



