# Benchmarking on pbmc data
# 26 Jan 2021
# Zhiyuan Hu
# Last modified 15 Aug 2021

library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

verbose <- FALSE
dirsave <- "pbmc"

methods <- c("asCIDER", "dnCIDER", "CCA", "harmony", "MNN", "scanorama", "liger","combat", "monocle","conos","rpca", "seurat", "sc3", "raceid")
method_names <-  c("asCIDER", "dnCIDER", "CCA", "Harmony", "fastMNN", "Scanorama","LIGER","Combat","Monocle3","conos","RPCA","Seurat\n(Clustering)", "SC3", "RaceID")

df_ari <- c()
for(i in 1:length(methods)){
  ari_tmp <- readRDS(paste0("~/Dropbox/CIDER/evaluate-integration/cbrg/rdata/", dirsave, "/", methods[i],"_ARI_res.rds"))
  df_ari <- rbind(df_ari, ari_tmp)
}
df_ari$method_name <- method_names
df_ari$method_name <- factor(df_ari$method_name, levels = rev(method_names))
df_ari
# ari-group----
ggplot(df_ari, aes(x = method_name, y = ARI)) +
  geom_bar(stat = "identity", fill = "steelblue")+ ylab("ARIpopulation") +
  geom_text(aes(label = round(ARI, 2)), vjust=0.5,hjust = -0.2, color="black", size=3.5) + 
  theme_linedraw() + xlab("") + coord_flip()
ggsave(paste0("~/Dropbox/CIDER/evaluate-integration/plots/",dirsave,"/",dirsave,"_ariGroup.pdf"), 
       width = 2.8, height = 3.2)

# ari-batch-----
ggplot(df_ari, aes(x = method_name, y = 1-batch_ARI)) +
  geom_bar(stat = "identity", fill = "steelblue") + ylab("1-ARIbatch") +
  geom_text(aes(label = round(1-batch_ARI, 2)), vjust=0.5,hjust = 1.1, color="black", size=3.5) + 
  theme_linedraw() + xlab("") + coord_flip() 
ggsave(paste0("~/Dropbox/CIDER/evaluate-integration/plots/",dirsave,"/",dirsave,"_ariBatch.pdf"), 
       width = 2.8, height = 3.2)

## runtime---
ggplot(df_ari, aes(x = method_name, y = runtime))  +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label = round(runtime, 1)), vjust=0.5,hjust = -0.2, color="black", size=3.5) + 
  theme_linedraw() + xlab("") + ylab("Runtime (s)")   + coord_flip(ylim = c(0, 1600)) 
ggsave(paste0("~/Dropbox/CIDER/evaluate-integration/plots/",dirsave,"/",dirsave,"_runtime.pdf"), 
       width = 2.8, height = 3.2)
