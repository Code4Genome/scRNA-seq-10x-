# 🔬 scRNA-seq 10x 

![grafik](https://github.com/user-attachments/assets/2d69a3a8-1455-4295-856b-8a0075fed1c6)

Single cell RNA sequencing (scRNA-seq) measures whole transcriptome gene expression in individual cells, offering a detailed view of how cells function and interact within complex tissues. 
Tracing individual transcripts back to their cells of origin enables researchers to pinpoint unique gene expression profiles in highly heterogeneous samples, ultimately revealing rare cell 
types and dynamic cell states that other methods often miss.

This repository contains scripts, and resources for analyzing single-cell RNA sequencing (scRNA-seq) data generated using the 10x Genomics platform. 
Designed for reproducibility and scalability, the pipeline supports key steps from raw FASTQ processing to downstream visualization and interpretation.

## 🧬 Workflow Overview

1. **Preprocessing & Alignment**
   - Input: FASTQ files
   - Tool: `cellranger count`
   - Output: Gene-barcode matrices (filtered/unfiltered)

2. **Quality Control**
   - Tools: Seurat / Scanpy
   - Metrics: nFeature_RNA, nCount_RNA, mitochondrial content

3. **Normalization & Integration**
   - SCTransform, log-normalization
   - Batch effect correction (Harmony / Seurat integration)

4. **Dimensionality Reduction & Clustering**
   - PCA, UMAP

5. **Differential Expression & Marker Identification**
   - FindMarkers (Seurat), rank_genes_groups (Scanpy)

6. **Cell Type Annotation**
   - Manual curation or reference-based annotation (SingleR, CellTypist)
  
## 📁 Project Structure

scRNAseq-10x/
├── scripts/ # Helper scripts (e.g., QC, clustering)
├── licens
└── README.md # Project overview
