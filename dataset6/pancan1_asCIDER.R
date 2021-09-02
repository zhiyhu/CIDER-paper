# Assisted CIDER on cancer data
# 5 Aug 2021
# Zhiyuan Hu
# Last modified 5 Aug 2021

# Set up ----
library(Seurat)
library(doParallel)
library(foreach)
library(limma)
library(edgeR)
source("/home/z/zhu/cider/functions.R")
verbose <- FALSE
dirsave <- "pancan"
n_cores <- 16

# Initialise----
runtime <- c() # count the run time
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds"))

# calculate dist mat ----
runtime_tmp <- system.time({
  # get hvgs
  hvgs <- mclapply(Seurat::SplitObject(seu, split.by = "Batch"),function(seu){
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = verbose)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 500, verbose = verbose)
    VariableFeatures(seu)
  }, mc.cores = n_cores)
  hvgs <- unique(unlist(hvgs))
  
  metadata <- data.frame(label = gsub(" ","",paste(seu$Group, seu$Batch, sep = "_")),
                         batch = seu$Batch,
                         timepoint = seu$timepoint,
                         ground_truth = seu$Group, stringsAsFactors = FALSE)
  
  n_size <- 35
  select <- downsampling2(metadata =  data.frame(label = paste(seu$Group, seu$Batch, seu$timepoint, sep = "_")), n_size = n_size/2, include = TRUE, replace = TRUE)
  
  matrix <- as.matrix(seu@assays$RNA@counts[,select])
  colnames(matrix) <- paste0(colnames(matrix), "_", 1:ncol(matrix))
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- edgeR::DGEList(counts = matrix[keep,,drop=F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]
  logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
  
  df <- data.frame(g = metadata$label[select],
                   b = metadata$batch[select], ## batch
                   t = metadata$timepoint[select],
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
  
  head(df)
  
  df_dist <- foreach(i = combinations$g1, j = combinations$g2, df = rep(list(df), n.iter), 
                     logCPM = rep(list(logCPM), n.iter),.combine = "rbind") %dopar% 
    {
      library(limma)
      df$tmp <- "bg"
      df$tmp[df$g == i] <- "g1"
      df$tmp[df$g == j] <- "g2"
      
      gidx <- which(df$b %in% df$b[c(df$tmp %in% c("g1","g2"))])
      
      # gidx <- which(df$b %in% df$b[c(df$tmp %in% c("g1","g2"))])
      # gidx2 <- downsampling2(metadata = data.frame(label = df$g), n_size = 10)
      # gidx <- unique(c(gidx, gidx2))
      
      ## design and contrast
      design <- model.matrix(~0 + tmp + b + t + detrate, data = df[gidx,])
      
      contrast_m <- limma::makeContrasts(
        contrasts = c("tmpg1-tmpbg", "tmpg2-tmpbg"),
        levels = design
      )
      getGrouFit.fasterer(logCPM[,gidx], design, contrast_m)
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
  hcluster <- cutree(hc, h = 0.8)
  df_merge <-  data.frame(initial_clusters = names(hcluster),
                          final_clusters = hcluster)
  
  seu$inicluster <- metadata$label
  
  seu$final_cluster <- df_merge$final_cluster[match(seu$inicluster,df_merge$initial_clusters)]
  seu$final_cluster[is.na(seu$final_cluster)] <- seu$inicluster[is.na(seu$final_cluster)]
})
runtime <- c(runtime, runtime_tmp[3])

print("Number of HVGs:")
print(length(hvgs))
# save results -----
results <- data.frame(method = "asCIDER",
                      runtime  = sum(runtime),
                      ARI = mclust::adjustedRandIndex(seu$final_cluster, seu$Group),
                      batch_ARI = mclust::adjustedRandIndex(seu$final_cluster, seu$Batch))

clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$final_cluster,
                             initial_cluster = seu$inicluster)

print(paste0("runtime: ", sum(runtime)))

saveRDS(combinations,  paste0("/home/z/zhu/cider/rdata/", dirsave, "/asCIDER_combinations.rds"))
saveRDS(select, paste0("/home/z/zhu/cider/rdata/", dirsave, "/asCIDER_select.rds"))
saveRDS(dist_coef, paste0("/home/z/zhu/cider/rdata/", dirsave, "/asCIDER_dist_coef.rds"))
saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/asCIDER_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/asCIDER_clustering_res.rds"))

