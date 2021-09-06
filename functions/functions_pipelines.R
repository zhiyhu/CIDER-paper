## simulation function

simulation_function <- function(m = c(300,300,300), # number of batches
                                g = 4,
                                out.prob = 0.2,
                                seed,
                                dropout.mid = 0.5) {
  
  dropout.type <- "experiment"
  set.seed(seed)
  x <- runif(g, min = 0.2, max = 1)
  x <- x/sum(x)
  if(sum(x) != 1){
    x[g] <- 1-sum(x[1:(g-1)])
  }
  
  params <- newSplatParams()
  params <- setParams(params, seed = seed, 
                      group.prob = x, 
                      batchCells = m,
                      
                      dropout.type = dropout.type,# dropout
                      dropout.mid = dropout.mid,
                      
                      out.prob = out.prob, # outliear
                      
                      batch.facLoc = 0.15, # batch
                      batch.facScale = 0.10 
  )
  
  sim.groups <- splatSimulate(params, method = "groups",
                              verbose = FALSE)
  logcounts(sim.groups) <- log2(calculateCPM(sim.groups) + 1)
  
  return(sim.groups)
}

## Seurat pipeline

seurat_clustering <- function(seu_obj, res, ndims = 14, nfeatures = 2000) {
  
  # seu <- as.Seurat(sim.groups)
  seu <- seu_obj
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = nfeatures)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, features = VariableFeatures(object = seu),verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:ndims)
  seu <- FindClusters(seu, resolution = res)
  
  seu_obj$seurat_cluster <- seu@meta.data$seurat_clusters
  return(seu_obj)
}

## Race ID pipeline

raceid_clustering <- function(seu_obj, g) {
  # vignette("RaceID")
  sc <- RaceID::SCseq(seu_obj@assays$RNA@counts)
  sc <- RaceID::filterdata(sc,mintotal=2000)
  fdata <- RaceID::getfdata(sc)
  sc <- RaceID::compdist(sc,metric="pearson")
  sc <- RaceID::clustexp(sc)
  
  RaceID::plotsaturation(sc,disp=FALSE)
  
  sc <- RaceID::clustexp(sc,cln=g,sat=FALSE)
  seu_obj$raceID_cluster <- sc@cluster$clb$result$result$pamobject$clustering
  
  return(seu_obj)
}


## SC3 function

sc3_clustering <- function(seu_obj, g) {

  sim.groups <- SingleCellExperiment(assays = list(counts = as.matrix(seu_obj@assays$RNA@counts)))
  logcounts(sim.groups) <- log2(calculateCPM(sim.groups) + 1)
  rowData(sim.groups)$feature_symbol <- rownames(sim.groups)
  sce <- SC3::sc3(sim.groups, ks = g, biology = F)
  
  seu_obj$sc3_cluster <- colData(sce)[,ncol(colData(sce))]
  
  return(seu_obj)
  
}


## MNN function

mnn_ppl <- function(sim.groups, g) {
  
  f.out <- batchelor::fastMNN(sim.groups, batch = sim.groups$Batch)
  
  reducedDim(sim.groups, "corrected") <- reducedDim(f.out, "corrected")
  scalecounts <- t(sim.groups@reducedDims$corrected[,1:15])
  data <- t(scalecounts) # transposed matrix
  
  if(nrow(data) > 1500) {
    nystrom.red <- T
    nystrom.sample <- 1500
  } else {
    nystrom.red <- F
    nystrom.sample <- nrow(data)
  }
  sc <- kernlab::specc(data,  # rows are cells
                       centers = g,
                       iterations = 1000,
                       kernel = "rbfdot",
                       # kpar = kpar,
                       nystrom.red = nystrom.red,
                       nystrom.sample = nystrom.sample
                       # mod.sample = mod.sample,
  )
  
  colData(sim.groups)$mnn_specc<- sc@.Data 
  return(sim.groups)
  
}



