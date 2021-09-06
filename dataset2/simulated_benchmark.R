# benchmark on simulated data
# 10 Aug 2021
# Zhiyuan Hu
# last modified 15 Aug 2021

library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

verbose <- FALSE
dirsave <- "simulation"

methods <- c("asCIDER", "dnCIDER", "CCA", "harmony", "MNN", "scanorama", "liger","combat", "monocle","conos","rpca", "seurat", "sc3", "raceid")
method_names <-  c("asCIDER", "dnCIDER", "CCA", "Harmony", "fastMNN", "Scanorama","LIGER","Combat","Monocle3","conos","RPCA","Seurat\n(Clustering)", "SC3", "RaceID")
df_ari <- c()
df_aribatch <- c()
df_runtime <- c()
for(i in 1:length(methods)){
  for(itor in 1:20){
    df_tmp <- readRDS(paste0("~/Dropbox/CIDER/evaluate-integration/cbrg/rdata/", dirsave, "/", methods[i],"_ARI_res",itor, ".rds"))
    df_ari <- rbind(df_ari, c(df_tmp$ARI, methods[i]))
    df_aribatch <- rbind(df_aribatch , c(df_tmp$batch_ARI, methods[i]))
    df_runtime <- rbind(df_runtime, c(df_tmp$runtime, methods[i]))
  }
}
colnames(df_ari) <- colnames(df_aribatch) <- colnames(df_runtime) <- c("var","method")
df_ari <- as.data.frame(df_ari)
df_aribatch <- as.data.frame(df_aribatch)
df_runtime <- as.data.frame(df_runtime)
df_ari$method_name <- df_aribatch$method_name <- df_runtime$method_name <- method_names[match(df_ari$method, methods)]
df_ari$method_name <- df_aribatch$method_name <- df_runtime$method_name <- factor(df_ari$method_name, levels = rev(method_names))

# barplots -----
df_ari$var <- as.numeric(as.character(df_ari$var))
ggplot(df_ari, aes(x = method_name, y = var)) +
  geom_violin(scale = "width") + geom_boxplot(width=0.3)+ ylab("ARIpopulation")  + 
  theme_linedraw() + xlab("") + coord_flip()
ggsave(paste0("~/Dropbox/CIDER/evaluate-integration/plots/",dirsave,"/",dirsave,"_ariGroup.pdf"), width = 2.8, height = 3.2)

df_aribatch$var <- as.numeric(as.character(df_aribatch$var))
ggplot(df_aribatch, aes(x = method_name, y = (1-var)))  + ylab("1-ARIbatch") +
  geom_violin(scale = "width") + geom_boxplot(width=0.3) + 
  theme_linedraw() + xlab("") + coord_flip() 
ggsave(paste0("~/Dropbox/CIDER/evaluate-integration/plots/",dirsave,"/",dirsave,"_ariBatch.pdf"), width = 2.8, height = 3.2)

df_runtime$var <- as.numeric(as.character(df_runtime$var))
ggplot(df_runtime, aes(x = method_name, y = var))  +
  geom_violin(scale = "width") + geom_boxplot(width=0.3)+ 
  theme_linedraw() + xlab("") + ylab("Runtime (s)")   + coord_flip() 
ggsave(paste0("~/Dropbox/CIDER/evaluate-integration/plots/",dirsave,"/",dirsave,"_runtime.pdf"), width = 2.8, height = 3.2)
