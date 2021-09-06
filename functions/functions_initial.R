## functions for processing initial clusters
## Include: calculateDistMat, calculateDist2 and mergeInitialClusters
## Zhiyuan Hu

## calculate distance matrix
calculateDistMat <- function(matrix, metadata, downsampling_factor = 5, verbose = TRUE, method = "voom"){
  
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- edgeR::DGEList(counts = matrix[keep,,drop=F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]
  # rm(matrix, keep)
  
  df <- data.frame(g = paste(metadata$label, metadata$batch, sep = "_"),
                   b = metadata$batch, ## batch
                   c = metadata$label, stringsAsFactors = F) ## label
  df$detrate <- colSums(matrix > 0.5)
  rownames(df) <- colnames(matrix)
  
  N <- length(unique(df$g))
  
  combinations <- data.frame(g1 = rep(unique(df$g), each = N), g2 = rep(unique(df$g), N), stringsAsFactors = FALSE)
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  combinations <- sort(combinations)
  
  idx <- c()
  for(i in 2:nrow(combinations)){
    if(!combinations$g2[i] %in% combinations$g1[1:(i-1)]) {
      idx <- c(idx, i)
    }
  }
  
  combinations <- combinations[c(1,idx),]
  rownames(combinations) <- 1:nrow(combinations)
  
  dist_p <- dist_t <- dist_coef <- matrix(0, nrow = N, ncol = N)
  colnames(dist_p) <- rownames(dist_p) <- sort(unique(df$g))
  colnames(dist_t) <- rownames(dist_t) <- sort(unique(df$g))
  colnames(dist_coef) <- rownames(dist_coef) <- sort(unique(df$g))
  
  if(verbose == TRUE) {
    
    message("Generating distance matrix...")
    pb <- txtProgressBar(min = 0, max = nrow(combinations), style = 3)
    k <- 1
  }
    
  # pairwise comparison
  for(i in 1:nrow(combinations)){
    
    df$tmp <- "bg"
    df$tmp[df$g == combinations$g1[i]] <- "g1"
    df$tmp[df$g == combinations$g2[i]] <- "g2"
    
    ## downsampling the bg
    n_bg <- sum(df$tmp == "bg")
    
    if(downsampling_factor > 1) {
      set.seed(12345)
      random.idx <- sample(x = which(df$tmp == "bg"), size = max(n_bg/downsampling_factor, 50), replace = FALSE)
      select2 <- c(random.idx, which(df$tmp %in% c("g1", "g2")) )
      select2 <- sort(select2)
    } else {
      select2  <- 1:nrow(df)
    }
    
    
    design <- model.matrix(~  0 + tmp + detrate, data = df[select2,])
    
    if (method == "voom") {
      
      v <- voom(dge[,select2], design, plot = FALSE)
      fit <- lmFit(v, design)
      groups <- paste0("tmp", unique(df$tmp))
      groups <-  groups[groups != "tmpbg"]
      perm_groups <- data.frame(g1 = groups,
                                g2 = "tmpbg", stringsAsFactors = F)
      perm_groups <- perm_groups[perm_groups$g1 != perm_groups$g2,]
      perm_groups$pair <- paste0(perm_groups$g1, "-", perm_groups$g2)
      contrast_m <- makeContrasts(contrasts = perm_groups$pair,
                                  levels = design)
      group_fit <- contrasts.fit(fit, contrast_m) 
      group_fit <- eBayes(group_fit)
      
    } else if (method == "trend") {
      
      logCPM <- cpm(dge[,select2], log=TRUE, prior.count=3)
      fit <- lmFit(logCPM, design)
      groups <- paste0("tmp", unique(df$tmp))
      groups <-  groups[groups != "tmpbg"]
      perm_groups <- data.frame(g1 = groups,
                                g2 = "tmpbg", stringsAsFactors = F)
      perm_groups <- perm_groups[perm_groups$g1 != perm_groups$g2,]
      perm_groups$pair <- paste0(perm_groups$g1, "-", perm_groups$g2)
      contrast_m <- makeContrasts(contrasts = perm_groups$pair,
                                  levels = design)
      fit <- contrasts.fit(fit, contrast_m)
      group_fit <- eBayes(fit, trend=TRUE, robust = TRUE)
      
    }
    
    group_fit$p.value[,1] <- group_fit$p.value[,1] + 0.00000001
    group_fit$p.value[,2] <- group_fit$p.value[,2] + 0.00000001
    
    idx1 <- rownames(dist_coef) == combinations$g1[i]
    idx2 <- colnames(dist_coef) == combinations$g2[i]
    dist_coef[idx1, idx2] <- cor(coef(group_fit)[,1], coef(group_fit)[,2])
    dist_t[idx1, idx2]    <- cor(group_fit$t[,1], group_fit$t[,2])
    dist_p[idx1, idx2]    <- cor(-log10(group_fit$p.value)[,1]*sign(coef(group_fit)[,1]), 
                                 -log10(group_fit$p.value)[,2]*sign(coef(group_fit)[,2]), 
                                 use = "pairwise.complete.obs")
    
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


calculateDistMatSubBatch <- function(matrix, metadata, downsampling_factor = 5, verbose = TRUE, method = "voom"){
  
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- edgeR::DGEList(counts = matrix[keep,,drop=F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]
  # rm(matrix, keep)
  
  df <- data.frame(g = paste(metadata$label, metadata$batch, sep = "_"),
                   b = metadata$batch, ## batch
                   subb = metadata$donor,
                   c = metadata$label, stringsAsFactors = F) ## label
  ## df$detrate <- colSums(matrix > 0.5) ## detrate updated 15 Sep 2020
  df$detrate <- df$detrate <- scale(colMeans(matrix > 0))[,1]
  rownames(df) <- colnames(matrix)
  
  N <- length(unique(df$g))
  
  combinations <- data.frame(g1 = rep(unique(df$g), each = N), g2 = rep(unique(df$g), N), stringsAsFactors = FALSE)
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  combinations <- sort(combinations)
  
  idx <- c()
  for(i in 2:nrow(combinations)){
    if(!combinations$g2[i] %in% combinations$g1[1:(i-1)]) {
      idx <- c(idx, i)
    }
  }
  
  combinations <- combinations[c(1,idx),]
  rownames(combinations) <- 1:nrow(combinations)
  
  dist_p <- dist_t <- dist_coef <- matrix(0, nrow = N, ncol = N)
  colnames(dist_p) <- rownames(dist_p) <- sort(unique(df$g))
  colnames(dist_t) <- rownames(dist_t) <- sort(unique(df$g))
  colnames(dist_coef) <- rownames(dist_coef) <- sort(unique(df$g))
  
  if(verbose == TRUE) {
    
    message("Generating distance matrix...")
    pb <- txtProgressBar(min = 0, max = nrow(combinations), style = 3)
    k <- 1
  }
  
  # pairwise comparison
  for(i in 1:nrow(combinations)){
    
    df$tmp <- "bg"
    df$tmp[df$g == combinations$g1[i]] <- "g1"
    df$tmp[df$g == combinations$g2[i]] <- "g2"
    
    ## downsampling the bg
    n_bg <- sum(df$tmp == "bg")
    
    if(downsampling_factor > 1) {
      set.seed(12345)
      random.idx <- sample(x = which(df$tmp == "bg"), size = max(n_bg/downsampling_factor, 50), replace = FALSE)
      select2 <- c(random.idx, which(df$tmp %in% c("g1", "g2")) )
      select2 <- sort(select2)
    } else {
      select2  <- 1:nrow(df)
    }
    
    
    design <- model.matrix(~  0 + tmp + subb + detrate, data = df[select2,])
    
    if (method == "voom") {
      
      v <- voom(dge[,select2], design, plot = FALSE)
      fit <- lmFit(v, design)
      groups <- paste0("tmp", unique(df$tmp))
      groups <-  groups[groups != "tmpbg"]
      perm_groups <- data.frame(g1 = groups,
                                g2 = "tmpbg", stringsAsFactors = F)
      perm_groups <- perm_groups[perm_groups$g1 != perm_groups$g2,]
      perm_groups$pair <- paste0(perm_groups$g1, "-", perm_groups$g2)
      contrast_m <- makeContrasts(contrasts = perm_groups$pair,
                                  levels = design)
      group_fit <- contrasts.fit(fit, contrast_m) 
      group_fit <- eBayes(group_fit)
      
    } else if (method == "trend") {
      
      logCPM <- cpm(dge[,select2], log=TRUE, prior.count=3)
      fit <- lmFit(logCPM, design)
      groups <- paste0("tmp", unique(df$tmp))
      groups <-  groups[groups != "tmpbg"]
      perm_groups <- data.frame(g1 = groups,
                                g2 = "tmpbg", stringsAsFactors = F)
      perm_groups <- perm_groups[perm_groups$g1 != perm_groups$g2,]
      perm_groups$pair <- paste0(perm_groups$g1, "-", perm_groups$g2)
      contrast_m <- makeContrasts(contrasts = perm_groups$pair,
                                  levels = design)
      fit <- contrasts.fit(fit, contrast_m)
      group_fit <- eBayes(fit, trend=TRUE, robust = TRUE)
      
    }
    
    group_fit$p.value[,1] <- group_fit$p.value[,1] + 0.00000001
    group_fit$p.value[,2] <- group_fit$p.value[,2] + 0.00000001
    
    idx1 <- rownames(dist_coef) == combinations$g1[i]
    idx2 <- colnames(dist_coef) == combinations$g2[i]
    dist_coef[idx1, idx2] <- cor(coef(group_fit)[,1], coef(group_fit)[,2])
    dist_t[idx1, idx2]    <- cor(group_fit$t[,1], group_fit$t[,2])
    dist_p[idx1, idx2]    <- cor(-log10(group_fit$p.value)[,1]*sign(coef(group_fit)[,1]), 
                                 -log10(group_fit$p.value)[,2]*sign(coef(group_fit)[,2]), 
                                 use = "pairwise.complete.obs")
    
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

## use one model
calculateDistMatOneModel <- function(matrix, metadata, downsampling_factor = 5, subset.row = NULL, verbose = TRUE, method = "voom"){
  
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- edgeR::DGEList(counts = matrix[keep,,drop=F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]
  if(!is.null(subset.row)) {
    subset.row <- subset.row[keep]
  }
    
  df <- data.frame(g = paste(metadata$label, metadata$batch, sep = "_"),
                   b = metadata$batch, ## batch
                   c = metadata$label, stringsAsFactors = F) ## label
  ## df$detrate <- colSums(matrix > 0.5) ## detrate updated 15 Sep 2020
  df$detrate <- scale(colMeans(matrix > 0))[,1]
  rownames(df) <- colnames(matrix)
  
  N <- length(unique(df$g))
  
  combinations <- data.frame(g1 = rep(unique(df$g), each = N), g2 = rep(unique(df$g), N), stringsAsFactors = FALSE)
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  combinations <- sort(combinations)
  
  idx <- c()
  for(i in 2:nrow(combinations)){
    if(!combinations$g2[i] %in% combinations$g1[1:(i-1)]) {
      idx <- c(idx, i)
    }
  }
  
  combinations <- combinations[c(1,idx),]
  rownames(combinations) <- 1:nrow(combinations)
  
  dist_p <- dist_t <- dist_coef <- matrix(0, nrow = N, ncol = N)
  colnames(dist_p) <- rownames(dist_p) <- sort(unique(df$g))
  colnames(dist_t) <- rownames(dist_t) <- sort(unique(df$g))
  colnames(dist_coef) <- rownames(dist_coef) <- sort(unique(df$g))
  
  if("donor" %in% colnames(metadata)){
    df$subb = metadata$donor
    design <- model.matrix(~  0 + g + subb + detrate, data = df)
  } else {
    design <- model.matrix(~  0 + g + detrate, data = df)
  }
  

  groups <- sort(unique(paste0("g", df$g)))
  n_groups <- length(groups)
  
  df_contrasts <- data.frame(target_group = groups, contrast = NA)
  
  for(i in 1:n_groups){
    df_contrasts$contrast[i] <- paste0(groups[i], "-(", paste(groups[-i], collapse = "+"), ")/", (n_groups-1))
  }
  # contrast matrix
  contrast_m <- makeContrasts(contrasts = df_contrasts$contrast, levels = design)
  colnames(contrast_m) <- groups
  
  if (method == "voom") {
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
    group_fit <- contrasts.fit(fit, contrast_m) 
    group_fit <- eBayes(group_fit)
    
  } else if (method == "trend") {
    logCPM <- edgeR::cpm(dge, log=TRUE, prior.count=3)
    fit <- lmFit(logCPM, design)
    fit <- contrasts.fit(fit, contrast_m)
    group_fit <- eBayes(fit, trend=TRUE, robust = TRUE)
    
  }
  
  for (i in 1:ncol(group_fit$p.value)) {
    group_fit$p.value[,i] <- group_fit$p.value[,i] + 0.00000001
  }
  
  if (is.null(subset.row)) { #
    subset.row <- rep(TRUE, nrow(coef(group_fit)))
  }
  
  # pairwise comparison
  for(i in 1:nrow(combinations)){

    idx1 <- rownames(dist_coef) == combinations$g1[i]
    idx2 <- colnames(dist_coef) == combinations$g2[i]
    
    pos1 <- df_contrasts$target_group == paste0("g", combinations$g1[i])
    pos2 <- df_contrasts$target_group == paste0("g", combinations$g2[i])
    
    dist_coef[idx1, idx2] <- cor(coef(group_fit)[subset.row, pos1], coef(group_fit)[subset.row, pos2])
    dist_t[idx1, idx2]    <- cor(group_fit$t[subset.row, pos1], group_fit$t[subset.row, pos2])
    dist_p[idx1, idx2]    <- cor(-log10(group_fit$p.value)[subset.row, pos1]*sign(coef(group_fit)[subset.row, pos1]), 
                                 -log10(group_fit$p.value)[subset.row, pos2]*sign(coef(group_fit)[subset.row, pos2]), 
                                 use = "pairwise.complete.obs")

  }
    
  return(list(dist_coef, dist_t, dist_p))
}


calculateDist2 <- function(matrix, metadata, verbose = TRUE){
  
  keep <- rowSums(matrix > 0.5) > 5 
  dge <- edgeR::DGEList(counts = matrix[keep,,drop=F]) # make a edgeR object
  dge <- dge[!grepl("ERCC-", rownames(dge)),] # remove ERCC
  dge <- dge[!grepl("MT-", rownames(dge)),]
  # rm(matrix, keep)
  
  df <- data.frame(g = paste(metadata$label, metadata$batch, sep = "_"),
                   b = metadata$batch, ## batch
                   c = metadata$label, stringsAsFactors = F) ## label
  df$detrate <- colSums(matrix > 0.5)
  rownames(df) <- colnames(matrix)
  
  N <- length(unique(df$g))
  
  combinations <- data.frame(g1 = rep(unique(df$g), each = N), 
                             g2 = rep(unique(df$g), N), 
                             stringsAsFactors = FALSE)
  combinations <- combinations[combinations$g1 != combinations$g2, ]
  combinations <- sort(combinations)
  
  idx <- c()
  for(i in 2:nrow(combinations)){
    if(!combinations$g2[i] %in% combinations$g1[1:(i-1)]) {
      idx <- c(idx, i)
    }
  }
  
  combinations <- combinations[c(1,idx),]
  rownames(combinations) <- 1:nrow(combinations)
  
  dist_p <- dist_t <- dist_coef <- matrix(0, nrow = N, ncol = N)
  colnames(dist_p) <- rownames(dist_p) <- sort(unique(df$g))
  colnames(dist_t) <- rownames(dist_t) <- sort(unique(df$g))
  colnames(dist_coef) <- rownames(dist_coef) <- sort(unique(df$g))
  
  ## by group
  design <- model.matrix(~  0 + g + detrate, data = df)  # Use 0 because we do not need intercept for this linear model
  # print(paste0("glm start: ", Sys.time()))
  v <- voom(dge, design, plot = FALSE)
  fit_g <- lmFit(v, design)
  # print(paste0("glm end: ", Sys.time()))
  
  for(i in 1:nrow(combinations)){
    bg <- unique(df$g)[!unique(df$g) %in% c(combinations$g1[i], combinations$g2[i])]
    bg <- paste0("g", bg)
    bg_collapsed <- paste0("(", paste(bg, collapse = "+"),")/",length(bg))
    combinations$g1_equ[i] <- paste0("g",combinations$g1[i], "-", bg_collapsed)
    combinations$g2_equ[i] <- paste0("g",combinations$g2[i], "-", bg_collapsed)
  }
  
  pair <- unique(c(combinations$g1_equ, combinations$g2_equ))
  contrast_m <- makeContrasts(contrasts = pair,
                              levels = design)
  
  group_fit <- contrasts.fit(fit_g, contrast_m) # Compute Contrasts from Linear Model Fit
  group_fit <- eBayes(group_fit)

  if(verbose == TRUE) {
    message("Generating distance matrix...")
    pb <- txtProgressBar(min = 0, max = nrow(combinations), style = 3)
    k <- 1
  }
  for(i in 1:nrow(combinations)){

    equ1 <- combinations$g1_equ[i]
    equ2 <- combinations$g2_equ[i]
    
    idx1 <- colnames(coef(group_fit)) == equ1
    idx2 <- colnames(coef(group_fit)) == equ2
    
    dist_coef[rownames(dist_coef) == combinations$g1[i], colnames(dist_coef) == combinations$g2[i]] <- cor(coef(group_fit)[,idx1], 
                                                                                                           coef(group_fit)[,idx2])
    dist_t[rownames(dist_t) == combinations$g1[i], colnames(dist_t) == combinations$g2[i]] <- cor(group_fit$t[,idx1], 
                                                                                                  group_fit$t[,idx2])
    
    dist_p[rownames(dist_p) == combinations$g1[i], 
           colnames(dist_p) == combinations$g2[i]] <- cor(-log10(group_fit$p.value)[,idx1]*sign(coef(group_fit)[,idx1]), 
                                                          -log10(group_fit$p.value)[,idx2]*sign(coef(group_fit)[,idx2]))
    
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

# merge initial clusters based on the distance matrix
mergeInitialClusters <- function(dist, cutoff = 0.7) {
  df <- dist 
  df <- data.frame(df)
  df <- cbind(rownames(df), df)
  colnames(df)[1] <- "g1"
  df2 <- df %>% tidyr::pivot_longer(cols = colnames(df)[2:ncol(df)], names_to = "g2", values_to = "value")
  df2 <- as.data.frame(df2, row.names = 1:nrow(df2))
  df2$g2 <- gsub("X", "", df2$g2)
  df2 <- df2[df2$g1 != df2$g2,]
  df2 <- df2[as.numeric(df2$value) > cutoff,]
  
  return(df2)

}
 

# initial clustering by Seurat Louvain
initialClusteringSeurat <- function(seu, dim = 10, resolution = 0.4, verbose = TURE){
  seu_list <- Seurat::SplitObject(seu, split.by = "Batch")
  p_list <- list()
  for(i in 1:length(seu_list)){
    seu_list[[i]] <- NormalizeData(seu_list[[i]], normalization.method = "LogNormalize", scale.factor = 10000, verbose = verbose)
    seu_list[[i]] <- FindVariableFeatures(seu_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = verbose)
    seu_list[[i]] <- ScaleData(seu_list[[i]], verbose = verbose)
    seu_list[[i]] <- RunPCA(seu_list[[i]], features = VariableFeatures(object = seu_list[[i]]), verbose = verbose)  
    seu_list[[i]] <- FindNeighbors(seu_list[[i]], dims = 1:dim, verbose = verbose)
    seu_list[[i]] <- FindClusters(seu_list[[i]], resolution = resolution, verbose = verbose)
    seu_list[[i]] <- RunTSNE(seu_list[[i]], dims = 1:dim)
    p_list[[i]] <- DimPlot(seu_list[[i]], reduction = "tsne")
  }
  return(list(seu_list, p_list))
  
}