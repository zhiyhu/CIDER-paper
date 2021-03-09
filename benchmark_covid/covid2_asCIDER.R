# Assisted CIDER on COVID data
# 24 Jan 2021
# Zhiyua Hu
# Last modified 26 Jan 2021

# Set up ----
library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(doParallel)
library(limma)
library(edgeR)

source("functions_simulate.R")
source("functions_initial.R")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
verbose <- FALSE
dirsave <- "covid19"

# number of cores
n.cores <- detectCores(logical = FALSE) - 1

# Initialise----
runtime <- 0 # count the run time
seu <- readRDS(paste0("rdata/",dirsave,"/seurat_object_preprocessed.rds"))

# calculate dist mat ----
runtime_tmp <- system.time({
  
  metadata <- data.frame(label = gsub(" ","",paste(seu$Group, seu$Batch, sep = "_")),
                         batch = seu$Batch,
                         ground_truth = seu$Group,
                         V1 = seu@reductions$tsne@cell.embeddings[,1],
                         V2 = seu@reductions$tsne@cell.embeddings[,2], stringsAsFactors = FALSE)
  
  n_size <- 35
  select <- downsampling(metadata = metadata, n_size = n_size, include = TRUE, replace = TRUE)
  
  matrix <- as.matrix(seu@assays$RNA@counts[,select])
  colnames(matrix) <- paste0(colnames(matrix), "_", 1:ncol(matrix))
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- edgeR::DGEList(counts = matrix[keep,,drop=F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]
  
  df <- data.frame(g = metadata$label[select],
                   b = metadata$batch[select], ## batch
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
  
  dist_p <- dist_coef <- matrix(0, nrow = N, ncol = N)
  colnames(dist_p) <- rownames(dist_p) <- unique(df$g)
  colnames(dist_coef) <- rownames(dist_coef) <- unique(df$g)
  rm(matrix, keep)
  
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  # showConnections()
  n.iter <- nrow(combinations) # number of iterations
  i <- NULL
  j <- NULL
  
  df_dist <- foreach(i = combinations$g1, j = combinations$g2, df = rep(list(df), n.iter), 
                     dge = rep(list(dge), n.iter), .combine = "rbind") %dopar% 
    {
      library(limma)
      df$tmp <- "bg"
      df$tmp[df$g == i] <- "g1"
      df$tmp[df$g == j] <- "g2"
      
    # design # contrast
    design <- model.matrix(~  0 + tmp + b + detrate, data = df)  
    contrast_m <- makeContrasts(contrasts = c("tmpg1-tmpbg", "tmpg2-tmpbg"), levels = design)
    group_fit <- getGroupFit(dge, design, contrast_m, method = "trend") # limma trend
    dist_coef_i <- cor(coef(group_fit)[, 1], coef(group_fit)[, 2])
    dist_p_i    <- cor(-log10(group_fit$p.value)[, 1] * sign(coef(group_fit)[, 1]), 
                       -log10(group_fit$p.value)[, 2] * sign(coef(group_fit)[, 2]))
    print(c(dist_coef_i, dist_p_i))
    }
  
  stopCluster(cl)
  for (i in 1:nrow(combinations)) {
    idx1 <- rownames(dist_coef) == combinations$g1[i]
    idx2 <- colnames(dist_coef) == combinations$g2[i]
    dist_coef[idx1, idx2] <- df_dist[i, 1]
    dist_p   [idx1, idx2] <- df_dist[i, 2]
  }
  
})
runtime <- runtime + runtime_tmp[3]  

saveRDS(dist_coef,  paste0("rdata/", dirsave, "/asCIDER_dist_coef.rds"))
saveRDS(dist_p,  paste0("rdata/", dirsave, "/asCIDER_dist_p.rds"))

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
runtime <- runtime + runtime_tmp[3]  

# save results -----
results <- data.frame(method = "asCIDER",
                      runtime  = runtime,
                      ARI = mclust::adjustedRandIndex(seu$final_cluster, seu$Group),
                      batch_ARI = mclust::adjustedRandIndex(seu$final_cluster, seu$Batch))

clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$final_cluster)

saveRDS(results, paste0("rdata/", dirsave, "/asCIDER_ARI_res.rds"))
saveRDS(clustering_res, paste0("rdata/", dirsave, "/asCIDER_clustering_res.rds"))

