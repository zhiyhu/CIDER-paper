# preprocessing covid data
# Zhiyuan Hu
# 28 Oct 2020; last modified 24 Jan 2021

library(ggplot2)
library(Seurat)
dirsave <- "covid19"
## - - - - - - - - - - - - - - - #
## Cell annotations              #
## - - - - - - - - - - - - - - - #
dat <- data.table::fread("raw_data/covid19/Lee_et_al_dataset.tsv.gz", data.table = FALSE, skip = 1, nThread = 2, stringsAsFactors = FALSE)
table(dat$V4)

# B cell, IgG-          B cell, IgG+          CD4, EM-like      CD4, non-EM-like          CD8, EM-like 
# 59573                  4345                  2048                  3517                  2380                  3242 
# CD8, non-EM-like    classical Monocyte                    DC intermediate Monocyte               NK cell nonclassical Monocyte 
# 6651                 18465                   650                   704                  9369                  1919 
# Platelet                   RBC        Uncategorized1        Uncategorized2 
# 3716                  1193                  1191                   182


sum(dat$V3=="severe COVID-19")
# [1] 10296

sum(dat$V3=="mild COVID-19")
# [1] 16742

cell_info <- dat[1:59572, 1:4]
colnames(cell_info) <- c("barcode", "cluster", "disease", "cell_type")
head(cell_info)
rm(dat)

## - - - - - - - - - - - - - - - #
## t-SNE coordinates             #
## - - - - - - - - - - - - - - - #

dat <- data.table::fread("raw_data/covid19/Lee_et_al_dataset.tsv.gz", data.table = FALSE, nThread = 2, skip = 59574)
df_tsne <- dat[ , 1:3]
colnames(df_tsne) <- c("barcode", "tSNE_1", "tSNE_2")
df_tsne <- df_tsne[match(cell_info$barcode, df_tsne$barcode), ]
cell_info <- cbind(cell_info, df_tsne[, 2:3])
cell_info$tSNE_1 <- as.numeric(cell_info$tSNE_1)
cell_info$tSNE_2 <- as.numeric(cell_info$tSNE_2)
rm(dat)
rm(df_tsne)

## - - - - - - - - - - - - - - - #
## 10x mtx data                  #
## - - - - - - - - - - - - - - - #
seu.data <- Read10X(data.dir = "raw_data/covid19/GSE149689/")
seu <- CreateSeuratObject(counts = seu.data)
seu <- seu[ , colnames(seu) %in% cell_info$barcode]
rm(seu.data)
all(colnames(seu) == cell_info$barcode)
# [1] TRUE

seu$original_cluster <- cell_info$cluster
seu$disease <- cell_info$disease
seu$cell_type <- cell_info$cell_type
seu$original_tSNE1 <- cell_info$tSNE_1
seu$original_tSNE2 <- cell_info$tSNE_2
seu$Batch <- seu$disease

## - - - - - - - - - - - - - - - #
## Seurat preprocessing pipeline #
## - - - - - - - - - - - - - - - #
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
ElbowPlot(seu)
seu <- RunTSNE(seu, reduction = "pca", dims = 1:15)

## Tidy up Batch, Group names
seu$Batch_name <- seu$Batch
seu$Batch <- paste0("Batch", as.numeric(as.factor(seu$Batch_name)))
seu$Group <- paste0("Group", as.numeric(as.factor(seu$cell_type)))

## save rds
saveRDS(seu, file = paste0("rdata/",dirsave,"/seurat_object_preprocessed.rds"), compress = TRUE)

