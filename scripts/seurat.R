#cell quality control 
A_nt[["percent.mt"]] <- PercentageFeatureSet(A_nt, pattern = "^MT-")
plot = VlnPlot(A_nt, pt.size = 0.10,
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot
SaveFigure(plot, "vln_QC_subset", width = 24, height = 6)

plot1 = FeatureScatter(A_nt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
plot2 = FeatureScatter(A_nt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
plot2
SaveFigure((plot1 + plot2),"scatter_QC", width = 12, height = 6, res = 20)

#Filtering 

A_nt = subset(A_nt, subset = nFeature_RNA > 200 & nFeature_RNA < 4000  & nCount_RNA < 20000 & percent.mt < 10)
plot3 = VlnPlot(A_nt, pt.size = 0.10,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot3
SaveFigure(plot3, "vln_QC_subset_after_filtering", width = 12, height = 3)

#Normalize Data

A_nt
A_nt = NormalizeData(A_nt, normalization.method = "LogNormalize", scale.factor = 10000)
A_nt

#Find variable features (2000 highly variable genes that have very high expression in some cells and very low in some)

A_nt = FindVariableFeatures(A_nt, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes

top5 = head(VariableFeatures(A_nt), 5)
top5

plot4 = VariableFeaturePlot(A_nt, cols = c("black", "orange"))
plot4

plot5 = LabelPoints(plot = plot4, points = top5, repel = TRUE)
plot5

SaveFigure((plot4 +plot5), "variable_features", width = 12, height = 3)

#Scale the data

all.genes = rownames(A_nt)
all.genes

A_nt = ScaleData(A_nt, features = all.genes)
A_nt

#High dimension reduction (PCA, Heatmap)
#principal components 50

A_nt = RunPCA(AQP4_nt, npcs = 50, features = VariableFeatures(object = A_nt))
print(A_nt[["pca"]], dims = 1:5, nfeatures = 5)

plot6 = VizDimLoadings(A_nt, dims = 1:4, reduction = "pca", nfeatures = 35, ncol = 4)
plot6

SaveObject(A_nt, "seurat_obj_after_PCA")
A_nt = ReadObject("seurat_obj_after_PCA")

#determine dimensionality of the data
ElbowPlot(A_nt)

DimPlot(A_nt, reduction = "pca", dims = (c(1,2)))
DimHeatmap(A_nt, dims = 1:2, cells = 500, balanced = TRUE, fast = FALSE)

#Cell clustering based in gene expression data based on nearest-neighbor method implemented in Seurat 
A_nt = FindNeighbors(A_nt, dims = 1:10)

#Understanding resolution
A_nt = FindClusters(A_nt, resolution = 0.25)
View(A_nt@meta.data)
DimPlot(A_nt, group.by = "RNA_snn_res.0.25", label = TRUE, label.size = 4)

A_nt = BuildClusterTree(A_nt, reorder = TRUE, reorder.numeric = TRUE)


#UMAP

A_nt = RunUMAP(A_nt, dims = 1:10)
plot_umap = DimPlot(A_nt, reduction = "umap")
plot_umap

SaveObject(A_nt, "seurat_obj_clustered")
A_nt = ReadObject("seurat_obj_clustered")
A_nt = readRDS("seurat_obj_clustered.rds")

aqp4_markers = FindAllMarkers(A_nt, min.pct = 0.25, logfc.threshold = 0.25)
aqp4_markers
aqp4_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


#Visualizing the top n genes per cluster

top = aqp4_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top
to_plot = unique(top$gene)
to_plot

plot8 = DotPlot(A_nt, features = to_plot, group.by = "tree.ident") + coord_flip()
plot8

DoHeatmap(A_nt, features = top$gene, group.bar.height = 0.01, size = 3) + NoLegend()
View(AQP4_nt@meta.data)

#ref = celldex::HumanPrimaryCellAtlasData()

#Identification of monocyte subsets
#ref2 = celldex::DatabaseImmuneCellExpressionData()
View(as.data.frame(colData(ref)))
View(A_nt@meta.data)
