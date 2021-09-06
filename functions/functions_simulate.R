## include: seurat_preprocessing
## include: downsampling
## include calculate_ider
## include cca_integration
## include mnn_integration
## Zhiyuan Hu

seurat_preprocess <- function(sce) {
  seu <- as.Seurat(sce)
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst",
                              nfeatures = 1000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
  seu <- RunTSNE(seu, reduction = "pca", dims = 1:10)
  
  return(seu)
}

## replace: using replace = TRUE if the cluster is smaller than required size
## include: using include = TRUE if including the cluster smaller than required size
downsampling <- function(metadata, n_size = 50, seed = 12345, include = FALSE, replace = FALSE, lower.cutoff = 3) {
  cluster <- unique(metadata$label)
  tech <- unique(metadata$batch)
  select <- c()
  for(i in cluster){
    for(j in tech) {
      idx <- which(metadata$label %in% i & metadata$batch %in% j)
      if(length(idx) > n_size){
        set.seed(seed)
        select <- c(select, sample(idx, size = n_size, replace = FALSE))
      } else if (length(idx) == n_size) {
        select <- idx
      } else if (length(idx) < n_size & length(idx) >= lower.cutoff) {        
        if(include & !replace) {
          select <- c(select, idx)
        } else if (include & replace){
          set.seed(seed)
          select <- c(select, sample(idx, size = n_size, replace = TRUE))
        }
      }
    }
  }
  return(select)
}

calculate_ider <- function(matrix, metadata, verbose = TRUE){
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- edgeR::DGEList(counts = matrix[keep,,drop=F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]
  # rm(matrix, keep)
  
  df <- data.frame(g = paste(metadata$label, metadata$batch, sep = "_")[select],
                   b = metadata$batch[select], ## batch
                   c = metadata$label[select], stringsAsFactors = F) ## label
  df$detrate <- colSums(matrix > 0.5)
  rownames(df) <- colnames(matrix)
  
  N <- length(unique(df$g))
  
  combinations <- data.frame(g1 = rep(unique(df$g), each = N), g2 = rep(unique(df$g), N), stringsAsFactors = FALSE)
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  
  idx <- c()
  for(i in 2:nrow(combinations)){
    if(!combinations$g2[i] %in% combinations$g1[1:(i-1)]) {
      idx <- c(idx, i)
    }
  }
  
  combinations <- combinations[c(1,idx),]
  rownames(combinations) <- 1:nrow(combinations)
  
  dist_p <- dist_t <- dist_coef <- matrix(0, nrow = N, ncol = N)
  colnames(dist_p) <- rownames(dist_p) <- unique(df$g)
  colnames(dist_t) <- rownames(dist_t) <- unique(df$g)
  colnames(dist_coef) <- rownames(dist_coef) <- unique(df$g)
  
  # create progress bar
  if(verbose == TRUE) {
    
    message("Generating distance matrix...")
    pb <- txtProgressBar(min = 0, max = nrow(combinations), style = 3)
    k <- 1
  }
  
  for(i in 1:nrow(combinations)){
    df$tmp <- "bg"
    df$tmp[df$g == combinations$g1[i]] <- "g1"
    df$tmp[df$g == combinations$g2[i]] <- "g2"
    
    ## by group
    design <- model.matrix(~  0 + tmp + b + detrate, data = df)  # Use 0 because we do not need intercept for this linear model
    v <- voom(dge, design, plot = F) 
    fit_g <- lmFit(v, design) 
    groups <- paste0("tmp", unique(df$tmp))
    groups <-  groups[groups != "tmpbg"]
    perm_groups <- data.frame(g1 = groups,
                              g2 = "tmpbg", stringsAsFactors = F)
    perm_groups <- perm_groups[perm_groups$g1 != perm_groups$g2,]
    perm_groups$pair <- paste0(perm_groups$g1, "-", perm_groups$g2)
    contrast_m <- makeContrasts(contrasts = perm_groups$pair, 
                                levels = design)
    group_fit <- contrasts.fit(fit_g, contrast_m) # Compute Contrasts from Linear Model Fit
    group_fit <- eBayes(group_fit)
    
    dist_coef[rownames(dist_coef) == combinations$g1[i], colnames(dist_coef) == combinations$g2[i]] <- cor(coef(group_fit)[,1], coef(group_fit)[,2])
    dist_t[rownames(dist_t) == combinations$g1[i], colnames(dist_t) == combinations$g2[i]] <- cor(group_fit$t[,1], group_fit$t[,2])
    dist_p[rownames(dist_p) == combinations$g1[i], colnames(dist_p) == combinations$g2[i]] <- cor(-log10(group_fit$p.value)[,1]*sign(coef(group_fit)[,1]), -log10(group_fit$p.value)[,2]*sign(coef(group_fit)[,2]))
    
    if(verbose == TRUE) {
      setTxtProgressBar(pb, k) # progress bar
      k <- k+1
    }
    
  }
  if(verbose == TRUE) { 
    close(pb) # close progress bar
  }
  
  return(list(dist_coef, dist_t, dist_p))
}



cca_integration <- function(seu) {
  seu.list <- SplitObject(seu, split.by = "Batch")
  for (i in 1:length(seu.list)) {
    seu.list[[i]] <- NormalizeData(seu.list[[i]], verbose = FALSE)
    seu.list[[i]] <- FindVariableFeatures(seu.list[[i]], selection.method = "vst", 
                                          nfeatures = 1000, verbose = FALSE)
  }
  seu.anchors <- FindIntegrationAnchors(object.list = seu.list, dims = 1:15)
  seu.integrated <- IntegrateData(anchorset = seu.anchors, dims = 1:15)
  
  DefaultAssay(seu.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
  seu.integrated <- RunPCA(seu.integrated, npcs = 15, verbose = FALSE)
  seu.integrated <- RunTSNE(seu.integrated, reduction = "pca", dims = 1:15)
  
  return(seu.integrated)
}


mnn_integration <- function(seu){
  sce <- SingleCellExperiment(assays = list(counts = seu@assays$RNA@counts, logcounts = seu@assays$RNA@data))
  # logcounts(sce) <- log2(calculateCPM(sce) + 1)
  sce$Batch  <- seu$Batch
  dec <- modelGeneVar(sce)
  chosen.hvgs <- dec$bio > 0
  
  set.seed(1000)
  f.out <- fastMNN(sce, batch = sce$Batch, subset.row=chosen.hvgs)
  # str(reducedDim(f.out, "corrected"))
  f.out$Group <- seu$Group
  f.out$Batch <- seu$Batch
  f.out <- runTSNE(f.out, dimred="corrected")
}

cosine_similarity <- function(x) {
  y <- t(x) %*% x
  res <- y / (sqrt(diag(y)) %*% t(sqrt(diag(y))))
  return(res)
}

library(limma)
library(edgeR)
getGroupFit <- function(dge, design, contrast_m, method) {
  if (method == "voom") {
    v <- voom(dge, design, plot = FALSE)
    fit <- limma::lmFit(v, design)
    group_fit <- contrasts.fit(fit, contrast_m)
    group_fit <- eBayes(group_fit)
  } else if (method == "trend") {
    logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
    fit <- lmFit(logCPM, design)
    group_fit <- contrasts.fit(fit, contrast_m)
    group_fit <- eBayes(group_fit, trend = TRUE, robust = TRUE)
  }
  
  group_fit$p.value[, 1] <- group_fit$p.value[, 1] + 0.00000001
  group_fit$p.value[, 2] <- group_fit$p.value[, 2] + 0.00000001
  
  return(group_fit)
}
