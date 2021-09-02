# De Novo CIDER on cross-species data
# 24 Jan 2021
# Zhiyuan Hu
# Last modified 3 AUg 2021
# Set up ----
library(Seurat)
library(doParallel)
library(parallel)
library(foreach)
library(limma)
library(edgeR)

source("/home/z/zhu/cider/functions.R")
verbose <- FALSE
dirsave <- "pancreas_twospecies"
n_cores <- 8

# Initialise----
runtime <- c() # count the run time
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds"))

# Initial clustering ----
gc()
runtime_tmp <- system.time({
  seu_list <- initialClusteringSeurat_faster(seu, res = 0.6, dim = 15, verbose = verbose, 
                                             additional.vars = "sample",n.cores = n_cores) 
  
  hvgs <- mclapply(seu_list, function(seu){
    VariableFeatures(seu)
  }, mc.cores = n_cores)
  hvgs <- unique(unlist(hvgs))
  
})
runtime <- c(runtime, runtime_tmp[3])
print(paste0("initial clustering runtime: ",paste(runtime_tmp, collapse = ",")))

knitr::kable(table(seu_list[[2]]$Group, seu_list[[2]]$seurat_clusters))

# Pre-merge ----
runtime_tmp <- system.time({
  dist_coef <- list()
  n_size <- 50
  
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
      
      dist_coef[[seu_itor]] <- matrix(0, nrow = N, ncol = N)
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
      rm(dge, keep, matrix)
      gc()
      fit <- lmFit(logCPM, design)
      group_fit <- contrasts.fit(fit, contrast_m)

      # pairwise comparison
      dist_coef[[seu_itor]] <- cor(coef(group_fit))
    }
  }
})
runtime <- c(runtime, runtime_tmp[3])
print(paste0("premerge runtime: ", paste(runtime_tmp, collapse = ",")))

# Create initial clusters ----
runtime_tmp <- system.time({
  for(seu_itor in 1:length(dist_coef)){
    tmp <- as.dist(dist_coef[[seu_itor]] )
    hc <- hclust(1 - tmp, method = "average")
    hres <- cutree(hc, h = 0.4)
    df_hres <- data.frame(hres)
    df_hres$hres <- paste0(df_hres$hres, "_", unique(seu_list[[seu_itor]]$Batch))
    seu_list[[seu_itor]]$inicluster_tmp <- paste0("g", seu_list[[seu_itor]]$seurat_clusters, "_", seu_list[[seu_itor]]$Batch)
    seu_list[[seu_itor]]$inicluster <- df_hres$hres[match(seu_list[[seu_itor]]$inicluster_tmp, rownames(df_hres))]
  }
})
runtime <- c(runtime, runtime_tmp[3])


# Calculate Dist mat ----
rm(logCPM)
gc()
runtime_tmp <- system.time({
  
  tmp <- unlist(sapply(seu_list, function(x) return(x$inicluster)))
  names(tmp) <- unlist(sapply(seu_list, function(x) return(colnames(x@assays$RNA@counts))))
  metadata <- data.frame(label = tmp[match(colnames(seu), names(tmp))],
                         batch = seu$Batch,
                         ground_truth = seu$Group, stringsAsFactors = FALSE)
  
  rm(seu_list)
  gc()
  
  n_size <- 35 # number of cells per initial cluster
  select <- downsampling2(metadata = metadata, n_size = n_size, include = TRUE, replace = TRUE)
  
  matrix <- as.matrix(seu@assays$RNA@counts[,select])
  colnames(matrix) <- paste0(colnames(matrix),1:ncol(matrix))
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- edgeR::DGEList(counts = matrix[keep, ,drop = F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]
  logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
  
  rm(dge, keep)
  gc()
  
  df <- data.frame(g = metadata$label[select],
                   b = metadata$batch[select], ## batch
                   ground_truth = metadata$ground_truth[select],
                   stringsAsFactors = F) ## label
  df$detrate <- colSums(matrix > 0.5)
  rownames(df) <- colnames(matrix)
  rm(matrix)
  gc()
  
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
  
  dist_coef <- matrix(0, nrow = N, ncol = N)
  colnames(dist_coef) <- rownames(dist_coef) <- unique(df$g)

  registerDoParallel(cores = 8)

  n.iter <- nrow(combinations) # number of iterations
  i <- NULL
  j <- NULL
  logCPM <- logCPM[rownames(logCPM) %in% hvgs,]
  df_dist <- foreach(i = combinations$g1, j = combinations$g2, df = rep(list(df), n.iter), 
                     logCPM = rep(list(logCPM), n.iter),.combine = "rbind") %dopar% 
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

      getGrouFit.fasterer(logCPM, design, contrast_m)
    }

  for (i in 1:nrow(combinations)) {
    idx1 <- rownames(dist_coef) == combinations$g1[i]
    idx2 <- colnames(dist_coef) == combinations$g2[i]
    dist_coef[idx1, idx2] <- df_dist[i, 1]

  }
  
})
runtime <- c(runtime, runtime_tmp[3])
gc()
print(paste0("calculate dist mat runtime: ",paste(runtime_tmp, collapse = ",")))

saveRDS(dist_coef,  paste0("/home/z/zhu/cider/rdata/", dirsave, "/dnCIDER_dist_coef.rds"))

# save.image(file = paste0("rdata/",dirsave,"/after_DE.RData"), compress = TRUE)

# dendrogram
# plot(hclust(as.dist(1-(dist_coef + t(dist_coef)))/2))

# Final clustering ----
runtime_tmp <- system.time({
  tmp <- dist_coef + t(dist_coef)
  diag(tmp) <- 1
  tmp <- 1 - tmp
  hc <- hclust(as.dist(tmp)/2)
  hcluster <- cutree(hc, h = 0.35)
  df_merge <-  data.frame(initial_clusters = names(hcluster),
                          final_clusters = hcluster)
  
  seu$inicluster <- metadata$label
  
  seu$final_cluster <- df_merge$final_cluster[match(seu$inicluster,df_merge$initial_clusters)]
  seu$final_cluster[is.na(seu$final_cluster)] <- seu$inicluster[is.na(seu$final_cluster)]
})
runtime <- c(runtime, runtime_tmp[3])
print(paste0("final clustering runtime: ",paste(runtime_tmp, collapse = ",")))


knitr::kable(table(seu$final_cluster, seu$Group))

mclust::adjustedRandIndex(seu$final_cluster, seu$Group)

# Save results ----
results <- data.frame(method = "dnCIDER",
                      runtime  = sum(runtime),
                      ARI = mclust::adjustedRandIndex(seu$final_cluster, seu$Group),
                      batch_ARI = mclust::adjustedRandIndex(seu$final_cluster, seu$Batch))

clustering_res <- data.frame(sample = colnames(seu),
                             ours_denovo = seu$final_cluster,
                             initial_cluster = seu$inicluster)

rm(seu)
gc()

print(paste0("runtime: ", sum(runtime)))
saveRDS(dist_coef,  paste0("/home/z/zhu/cider/rdata/", dirsave, "/dnCIDER_dist_coef.rds"))
saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/dnCIDER_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/dnCIDER_clustering_res.rds"))

