# Assisted CIDER on COVID data
# 24 Jan 2021
# Zhiyuan Hu
# Last modified 3 Aug 2021

# Set up ----
library(Seurat)
library(doParallel)
library(foreach)
library(limma)
library(edgeR)
source("/home/z/zhu/cider/functions.R")
verbose <- FALSE
dirsave <- "covid19"
n_cores <- 16

# Initialise----
runtime <- c() # count the run time
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds"))
seu$disease <- seu$Batch
seu$Batch <- seu$donor

# calculate dist mat ----
runtime_tmp <- system.time({
  
  # get hvgs
  hvgs <- mclapply(Seurat::SplitObject(seu, split.by = "Batch"),function(seu){
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = verbose)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = verbose)
    VariableFeatures(seu)
  })
  hvgs <- unique(unlist(hvgs))
  
  metadata <- data.frame(label = gsub(" ","",paste(seu$Group, seu$Batch, sep = "_")),
                         batch = seu$Batch,
                         disease = seu$disease,
                         ground_truth = seu$Group, stringsAsFactors = FALSE)
  
  n_size <- 35
  select <- downsampling2(metadata = metadata, n_size = n_size, include = TRUE, replace = TRUE)
  
  matrix <- as.matrix(seu@assays$RNA@counts[,select])
  colnames(matrix) <- paste0(colnames(matrix), "_", 1:ncol(matrix))
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- edgeR::DGEList(counts = matrix[keep,,drop=F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]
  logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
  
  df <- data.frame(g = metadata$label[select],
                   b = metadata$batch[select], ## batch
                   d = metadata$disease[select],
                   ground_truth = metadata$ground_truth[select],
                   stringsAsFactors = F) ## label
  df$detrate <- scale(colMeans(matrix > 0))[,1]
  rownames(df) <- colnames(matrix)
  
  N <- length(unique(df$g)) # prepare pair information
  combinations <- data.frame(g1 = rep(unique(df$g), each = N), 
                             g2 = rep(unique(df$g), N), 
                             stringsAsFactors = FALSE)
  # remove pairs of the same group
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  # remove pairs from the same batch
  combinations$b1 <- df$b[match(combinations$g1, df$g)]
  combinations$b2 <- df$b[match(combinations$g2, df$g)]
  combinations <- combinations[combinations$b1 != combinations$b2,]
  
  idx <- c() # remove redundant pairs (order switches)
  for(i in 2:nrow(combinations)){
    if(!combinations$g2[i] %in% combinations$g1[1:(i-1)]) {
      idx <- c(idx, i)
    }
  }
  
  combinations <- combinations[c(1,idx),]
  rownames(combinations) <- 1:nrow(combinations)
  
  dist_coef <- matrix(0, nrow = N, ncol = N)
  colnames(dist_coef) <- rownames(dist_coef) <- unique(df$g)
  registerDoParallel(cores = n_cores)
  
  n.iter <- nrow(combinations) # number of iterations
  i <- NULL
  j <- NULL
  
  rm(dge, matrix, keep)
  gc()
  logCPM <- logCPM[rownames(logCPM) %in% hvgs,]
  
  df_dist <- foreach(i = combinations$g1, j = combinations$g2, df = rep(list(df), n.iter), 
                     logCPM = rep(list(logCPM), n.iter),.combine = "rbind") %dopar% 
    {
      library(limma)
      df$tmp <- "bg"
      df$tmp[df$g == i] <- "g1"
      df$tmp[df$g == j] <- "g2"
      
      ## design and contrast
      design <- model.matrix(~  0 + tmp + d + b + detrate, data = df) 
      contrast_m <- limma::makeContrasts(
        contrasts = c("tmpg1-tmpbg", "tmpg2-tmpbg"),
        levels = design
      )
      getGrouFit.fasterer(logCPM, design, contrast_m)
    }
  stopImplicitCluster()
  for (i in 1:nrow(combinations)) {
    idx1 <- rownames(dist_coef) == combinations$g1[i]
    idx2 <- colnames(dist_coef) == combinations$g2[i]
    dist_coef[idx1, idx2] <- df_dist[i, 1]
  
  }
  
})
runtime <- c(runtime, runtime_tmp[3])
print(paste0("calculate dist mat runtime: ",paste(runtime_tmp, sep=",")))

# saveRDS(dist_coef,  paste0("rdata/", dirsave, "/asCIDER_dist_coef.rds"))

# final clustering ----
runtime_tmp <- system.time({
  tmp <- dist_coef + t(dist_coef)
  hc <- hclust(as.dist(1-(tmp)))
  hcluster <- cutree(hc, h = 0.45)
  df_merge <-  data.frame(initial_clusters = names(hcluster),
                          final_clusters = hcluster)
  
  seu$inicluster <- metadata$label
  
  seu$final_cluster <- df_merge$final_cluster[match(seu$inicluster,df_merge$initial_clusters)]
  seu$final_cluster[is.na(seu$final_cluster)] <- seu$inicluster[is.na(seu$final_cluster)]
})
runtime <- c(runtime, runtime_tmp[3])

# save results -----
results <- data.frame(method = "asCIDER",
                      runtime  = sum(runtime),
                      ARI = mclust::adjustedRandIndex(seu$final_cluster, seu$Group),
                      batch_ARI = mclust::adjustedRandIndex(seu$final_cluster, seu$Batch))

clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$final_cluster,
                             initial_cluster = seu$inicluster)

print(paste0("runtime: ", sum(runtime)))

saveRDS(dist_coef, paste0("/home/z/zhu/cider/rdata/", dirsave, "/asCIDER_dist_coef_sr.rds"))
saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/asCIDER_ARI_res_sr.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/asCIDER_clustering_res_sr.rds"))

