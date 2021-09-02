# run conos on simulated data
# 15 Aug 2021
# Zhiyuan Hu
# last modified 15 Aug 2021

library(conos)
library(dplyr)
library(Seurat)
library(pagoda2)
library(SeuratWrappers)

verbose <- FALSE
dirsave <- "simulation"

for(itor in 1:20){
  seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed",itor,".rds"))
  runtime <- system.time({
    seu.list <- SplitObject(seu, split.by = "Batch")
    for (i in 1:length(seu.list)) {
      seu.list[[i]] <- NormalizeData(seu.list[[i]]) %>% FindVariableFeatures() %>% ScaleData() %>% 
        RunPCA(verbose = FALSE)
    }
    con <- Conos$new(seu.list)
    con$buildGraph(k = 15, k.self = 5, space = "PCA", ncomps = 30, n.odgenes = 2000, matching.method = "mNN", 
                   metric = "angular", score.component.variance = TRUE, verbose = FALSE)
    con$findCommunities()
    con$embedGraph()
    seu.integrated <- as.Seurat(con)
  })
  print("Runtime:")
  print(runtime)
  
  head(seu.integrated@meta.data)
  
  # Save results ----
  results <- data.frame(method = "conos",
                        runtime  = as.numeric(runtime[3]),
                        ARI = mclust::adjustedRandIndex(seu.integrated$Group, seu.integrated$leiden),
                        batch_ARI = mclust::adjustedRandIndex(seu.integrated$Batch, seu.integrated$leiden))
  
  clustering_res <- data.frame(sample = colnames(seu.integrated),
                               clusters = seu.integrated$leiden)
  
  results
  
  saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/conos_ARI_res",itor,".rds"))
  saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/conos_clustering_res",itor,".rds"))
}

sessionInfo()