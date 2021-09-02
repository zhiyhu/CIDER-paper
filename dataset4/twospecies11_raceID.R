# RaceID on cross species data
# 25 Jan 2021
# Zhiyuan Hu
# Last modified 3 Aug  2021

# Set up ----
library(batchelor)
library(scater)
library(scran)
library(SingleCellExperiment)
library(Seurat)

verbose <- FALSE
dirsave <- "pancreas_twospecies"
source("/home/z/zhu/cider/functions_pipelines.R")

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
seu <- Seurat::CreateSeuratObject(counts = seu@assays$RNA@counts, 
                                  meta.data = data.frame(Batch = seu$Batch,
                                                         Group = seu$Group))
gc()

n_groups <- length(unique(seu$Group))

# RaceID-----
runtime <- system.time({
  # refer to https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html
  sc <- RaceID::SCseq(seu@assays$RNA@counts)
  sc <- RaceID::filterdata(sc,mintotal=2000)
  fdata <- RaceID::getfdata(sc)
  sc <- RaceID::compdist(sc,metric="pearson")
  sc <- RaceID::clustexp(sc, samp = 1000) # for large dataset 
  
  # RaceID::plotsaturation(sc,disp=FALSE)
  
  sc <- RaceID::clustexp(sc,cln=n_groups,sat=FALSE)
  seu$raceID_cluster <- sc@cluster$clb$result$result$pamobject$clustering
  
})[3]

# Save results ----
results <- data.frame(method = "raceid",
                      runtime  = runtime,
                      ARI = mclust::adjustedRandIndex(seu$raceID_cluster, seu$Group),
                      batch_ARI = mclust::adjustedRandIndex(seu$raceID_cluster, seu$Batch))
clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$raceID_cluster)

saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/raceid_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/raceid_clustering_res.rds"))

