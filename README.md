# 🔬 scRNA-seq 10x 

![grafik](https://github.com/user-attachments/assets/bdebc1c8-c60f-4e46-a011-1a7b1402bd2f)


Single cell RNA sequencing (scRNA-seq) measures whole transcriptome gene expression in individual cells, offering a detailed view of how cells function and interact within complex tissues. 
Tracing individual transcripts back to their cells of origin enables researchers to pinpoint unique gene expression profiles in highly heterogeneous samples, ultimately revealing rare cell 
types and dynamic cell states that other methods often miss.

This repository contains scripts, and resources for analyzing single-cell RNA sequencing (scRNA-seq) data generated using the 10x Genomics platform. 
Designed for reproducibility and scalability, the pipeline supports key steps from raw FASTQ processing to downstream visualization and interpretation.

More info: https://www.10xgenomics.com/
## 🧬 Workflow Overview

1. **Preprocessing & Alignment**
   - Input: Paired-end FASTQ files (typically from Illumina sequencing)
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

├── bash/ # Raw and processed data (FASTQ, Cell Ranger)  
├── scripts/ # Data analysis scripts in R (Seurat, Harmony, scTransform)  
├── README.md  
└── LICENSE

