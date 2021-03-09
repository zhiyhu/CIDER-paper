# De Novo CIDER on COVID data
# 24 Jan 2021
# Zhiyua Hu
# Last modified 24 Jan 2021

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

# Initial clustering ----
runtime_tmp <- system.time({
  output <- initialClusteringSeurat(seu, res = 0.5, dim = 14, verbose = verbose) 
})
runtime <- runtime + runtime_tmp[3]
seu_list <- output[[1]]
p_list <- output[[2]]
rm(output)
plot_grid(plotlist = p_list)
ggsave(paste0("plots/",dirsave,"/metaclustering/dnCIDER_tSNE_rawInitialClusters.pdf"), width = 12, height = 9)

# Pre-merge ----
runtime_tmp <- system.time({
  dist_p <- list()
  dist_coef <- list()
  n_size <- 35
  
  for(seu_itor in 1:length(seu_list)){
    
    df_info <- data.frame(label = seu_list[[seu_itor]]$seurat_clusters,
                          batch = seu_list[[seu_itor]]$Batch)
    idx <- downsampling(metadata = df_info, n_size = n_size, include = TRUE, replace = TRUE)
    idx <- sort(idx)
    matrix <- as.matrix(seu_list[[seu_itor]]@assays$RNA@counts[, idx])
    colnames(matrix) <- paste0(colnames(matrix), "_", 1:ncol(matrix))
    df2 <- df_info[idx, ]
    rownames(df2) <- colnames(matrix)
    
    if(length(unique(df_info$label[idx])) > 2) { # more than 2 grousp in batch
      keep <- rowSums(matrix > 0.5) > 5 
      dge <- edgeR::DGEList(counts = matrix[keep, , drop = F]) # make a edgeR object
      dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
      dge <- dge[!grepl("MT-", rownames(dge)),]
      
      df <- data.frame(g = paste(df2$label, df2$batch, sep = "_"),
                       b = df2$batch, ## batch
                       c = df2$label, stringsAsFactors = F) ## label
      df$detrate <- scale(colMeans(matrix > 0))[,1]
      rownames(df) <- colnames(matrix)
      
      N <- length(unique(df$g)) # number of groups
      combinations <- data.frame(g1 = rep(unique(df$g), each = N), 
                                 g2 = rep(unique(df$g), N), stringsAsFactors = FALSE)
      combinations <- combinations[combinations$g1 != combinations$g2, ]
      idx <- c()
      for(i in 2:nrow(combinations)){
        if(!combinations$g2[i] %in% combinations$g1[1:(i-1)]) {
          idx <- c(idx, i)
        }
      }
      combinations <- combinations[c(1,idx),]
      rownames(combinations) <- 1:nrow(combinations)
      combinations <- na.omit(combinations)
      
      dist_p[[seu_itor]] <- dist_coef[[seu_itor]] <- matrix(0, nrow = N, ncol = N)
      colnames(dist_p[[seu_itor]]) <- rownames(dist_p[[seu_itor]]) <- sort(unique(df$g))
      colnames(dist_coef[[seu_itor]]) <- rownames(dist_coef[[seu_itor]]) <- sort(unique(df$g))
      
      # design and contrast
      design <- model.matrix(~  0 + g + detrate, data = df)
      groups <- sort(unique(paste0("g", df$g)))
      n_groups <- length(groups)
      df_contrasts <- data.frame(target_group = groups, contrast = NA)
      for(i in 1:n_groups){
        df_contrasts$contrast[i] <- paste0(groups[i], "-(", paste(groups[-i], collapse = "+"), ")/", (n_groups-1))
      }
      contrast_m <- makeContrasts(contrasts = df_contrasts$contrast, levels = design)
      colnames(contrast_m) <- groups
      
      logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
      fit <- lmFit(logCPM, design)
      fit <- contrasts.fit(fit, contrast_m)
      group_fit <- eBayes(fit, trend = TRUE, robust = TRUE)
      
      for (i in 1:ncol(group_fit$p.value)) {
        group_fit$p.value[,i] <- group_fit$p.value[,i] + 0.00000001
      }
      
      # pairwise comparison
      dist_coef[[seu_itor]] <- cor(coef(group_fit))
      dist_p[[seu_itor]] <- cor(-log10(group_fit$p.value) * sign(coef(group_fit)) ) 
      
    }
  }
})
runtime <- runtime + runtime_tmp[3]

# Create initial clusters ----
runtime_tmp <- system.time({
  for(seu_itor in 1:length(dist_coef)){
    tmp <- as.dist(dist_coef[[seu_itor]] )
    hc <- hclust(1 - tmp)
    hres <- cutree(hc, h = 0.3)
    df_hres <- data.frame(hres)
    df_hres$hres <- paste0(df_hres$hres, "_", unique(seu_list[[seu_itor]]$Batch))
    seu_list[[seu_itor]]$inicluster_tmp <- paste0("g", seu_list[[seu_itor]]$seurat_clusters, "_", seu_list[[seu_itor]]$Batch)
    seu_list[[seu_itor]]$inicluster <- df_hres$hres[match(seu_list[[seu_itor]]$inicluster_tmp, rownames(df_hres))]
  }
})
runtime <- runtime + runtime_tmp[3]

# Calculate Dist mat ----
gc()
runtime_tmp <- system.time({
  
  tmp <- unlist(sapply(seu_list, function(x) return(x$inicluster)))
  names(tmp) <- unlist(sapply(seu_list, function(x) return(colnames(x@assays$RNA@counts))))
  metadata <- data.frame(label = tmp[match(colnames(seu), names(tmp))],
                         batch = seu$Batch,
                         ground_truth = seu$Group,
                         V1 = seu@reductions$tsne@cell.embeddings[,1],
                         V2 = seu@reductions$tsne@cell.embeddings[,2], stringsAsFactors = FALSE)
  
  
  n_size <- 40 # number of cells per initial cluster
  select <- downsampling(metadata = metadata, n_size = n_size, include = TRUE, replace = TRUE)
  
  matrix <- as.matrix(seu@assays$RNA@counts[,select])
  colnames(matrix) <- paste0(colnames(matrix),1:ncol(matrix))
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- edgeR::DGEList(counts = matrix[keep, ,drop = F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]
  
  df <- data.frame(g = metadata$label[select],
                   b = metadata$batch[select], ## batch
                   ground_truth = metadata$ground_truth[select],
                   stringsAsFactors = F) ## label
  df$detrate <- colSums(matrix > 0.5)
  rownames(df) <- colnames(matrix)
  
  N <- length(unique(df$g))
  
  combinations <- data.frame(g1 = rep(unique(df$g), each = N), g2 = rep(unique(df$g), N), stringsAsFactors = FALSE)
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  
  combinations$b1 <- df$b[match(combinations$g1, df$g)]
  combinations$b2 <- df$b[match(combinations$g2, df$g)]
  combinations <- combinations[combinations$b1!=combinations$b2,]
  
  idx <- c()
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
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  
  showConnections()
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
      
      ## design and contrast
      design <- model.matrix(~  0 + tmp + b + detrate, data = df) 
      contrast_m <- limma::makeContrasts(
        contrasts = c("tmpg1-tmpbg", "tmpg2-tmpbg"),
        levels = design
      )
      
      group_fit <- getGroupFit(dge, design, contrast_m, method = "trend")
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
gc()

# save.image(file = paste0("rdata/",dirsave,"/after_DE.RData"), compress = TRUE)

# dendrogram
plot(hclust(as.dist(1-(dist_coef + t(dist_coef)))/2))

# Final clustering ----
runtime_tmp <- system.time({
  tmp <- dist_coef + t(dist_coef)
  diag(tmp) <- 1
  tmp <- 1 - tmp
  hc <- hclust(as.dist(tmp))
  hcluster <- cutree(hc, h = 0.4)
  df_merge <-  data.frame(initial_clusters = names(hcluster),
                          final_clusters = hcluster)
  
  seu$inicluster <- metadata$label
  
  seu$final_cluster <- df_merge$final_cluster[match(seu$inicluster,df_merge$initial_clusters)]
  seu$final_cluster[is.na(seu$final_cluster)] <- seu$inicluster[is.na(seu$final_cluster)]
})
runtime <- runtime + runtime_tmp[3] 

# Save results ----
results <- data.frame(method = "ours_denovo",
                      runtime  = runtime,
                      ARI = mclust::adjustedRandIndex(seu$final_cluster, seu$Group),
                      batch_ARI = mclust::adjustedRandIndex(seu$final_cluster, seu$Batch))

clustering_res <- data.frame(sample = colnames(seu),
                             ours_denovo = seu$final_cluster)

saveRDS(results, paste0("rdata/", dirsave, "/dnCIDER_ARI_res.rds"))
saveRDS(clustering_res, paste0("rdata/", dirsave, "/dnCIDER_clustering_res.rds"))

