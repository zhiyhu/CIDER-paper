# Benchmarking on cancer data"
# author: "Zhiyuan"
# date: "10 Aug 2021
# last modified 15 Aug 2021

library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

verbose <- FALSE
dirsave <- "pancan"
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

methods <- c("asCIDER", "dnCIDER", "CCA", "harmony", "MNN", "scanorama", "liger","combat", "monocle","conos","rpca", "seurat", "sc3","raceID")
method_names <-  c("asCIDER", "dnCIDER", "CCA","Harmony", "fastMNN", "Scanorama","LIGER","Combat\n(Downsampled)","Monocle3","conos","RPCA","Seurat\n(Clustering)", "SC3\n(Downsampled)", "RaceID\n(Downsampled)")
df_ari <- c()
for(i in 1:length(methods)){
  ari_tmp <- readRDS(paste0("~/Dropbox/CIDER/evaluate-integration/cbrg/rdata/", dirsave, "/", methods[i],"_ARI_res.rds"))
  df_ari <- rbind(df_ari, ari_tmp)
}
df_ari$method_name <- method_names
df_ari$method_name <- factor(df_ari$method_name, levels = rev(method_names))
df_ari
# method   runtime       ARI  batch_ARI           method_name
# 1          asCIDER   232.865 0.8449760 0.06594349               asCIDER
# 2          dnCIDER  1483.582 0.6674808 0.08189789               dnCIDER
# 3          harmony  2102.236 0.4422195 0.03576350               Harmony
# 4              MNN  1709.119 0.4251380 0.05952601               fastMNN
# 5        scanorama  1938.900 0.6528566 0.06411574             Scanorama
# 6            liger 16566.499 0.3436232 0.02378216                 LIGER
# 7           combat  1369.131 0.4019788 0.06613228                Combat
# elapsed     raceid  2893.243 0.3561035 0.06777261              Monocle3
# elapsed1    seurat   359.578 0.4962668 0.09793193  Seurat\n(Clustering)
# elapsed2       SC3  1212.044 0.4376284 0.03350662    SC3\n(Downsampled)
# elapsed3    raceid 48174.399 0.3747453 0.03667719 RaceID\n(Downsampled)

# barplots -----
ggplot(df_ari, aes(x = method_name, y = ARI)) +
  geom_bar(stat = "identity", fill = "steelblue") + ylab("ARIpopulation") +
  geom_text(aes(label = round(ARI, 2)), vjust = 0.5, hjust = -0.2, color="black", size=3.5) + 
  theme_linedraw() + xlab("") + coord_flip(ylim = c(0,1))
ggsave(paste0("~/Dropbox/CIDER/evaluate-integration/plots/",dirsave,"/",dirsave,"_ariGroup.pdf"), 
       width = 2.8, height = 3.8)

ggplot(df_ari, aes(x = method_name, y = 1-batch_ARI)) + 
  geom_bar(stat = "identity", fill = "steelblue") + ylab("1-ARIbatch") +
  geom_text(aes(label = round(1-batch_ARI, 2)), vjust=0.5,hjust = 1.1, color="black", size=3.5) + 
  theme_linedraw() + xlab("") + coord_flip() 
ggsave(paste0("~/Dropbox/CIDER/evaluate-integration/plots/",dirsave,"/",dirsave,"_ariBatch.pdf"), 
       width = 2.8, height = 3.8)

ggplot(df_ari, aes(x = method_name, y = runtime/60))  +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label = round(runtime/60, 1)), vjust=0.5,hjust = -0.2, color="black", size=3.5) + 
  theme_linedraw() + xlab("") + ylab("Runtime (min)")   + coord_flip() 
ggsave(paste0("~/Dropbox/CIDER/evaluate-integration/plots/",dirsave,"/",dirsave,"_runtime.pdf"), 
       width = 2.8, height = 3.8)

# heatmap----
df_cell <- readRDS(paste0("~/Dropbox/CIDER/evaluate-integration/cbrg/rdata/",dirsave,"/seurat_object_metadata.rds"))

res_ascider <- readRDS("cbrg/rdata/pancan/asCIDER_clustering_res.rds")
dist <- readRDS("cbrg/rdata/pancan/asCIDER_dist_coef.rds")
tmp <- dist + t(dist)
diag(tmp) <- 1
hc <- hclust(as.dist(1-(tmp)))

groups <- colnames(tmp)
df_groups <- data.frame(groups  = groups)
df_groups$cellType <- sapply(groups, function(x) unlist(strsplit(x, split = "_"))[1])
df_groups$batch <- sapply(groups, function(x) unlist(strsplit(x, split = "_BIOKEY_"))[2])
df_groups$batch <- paste0("BIOKEY_", df_groups$batch)
df_groups$disease <- df_cell$BC_type[match(df_groups$batch, df_cell$Batch)]
df_groups$timepoint <- df_cell$timepoint[match(df_groups$batch, df_cell$Batch)]
rownames(df_groups) <- groups

my_colour <- list()
cols <- col_vector[1:8]
names(cols) <- unique(df_groups$cellType)
my_colour[[1]] <- cols
cols <- col_vector[9:11]
names(cols) <- unique(df_groups$disease)
my_colour[[2]] <- cols
names(my_colour) <- c("cellType","disease")

pdf("~/Dropbox/CIDER/evaluate-integration/plots/pancan/asCIDER_similarity_matrix.pdf", width = 10, height = 10)
pheatmap::pheatmap(
  tmp,
  annotation_row = df_groups[,c("cellType","disease")],
  annotation_col = df_groups[,c("cellType","disease")],
  annotation_colors = my_colour,
  color = viridis::inferno(100),
  border_color = NA,
  display_numbers = FALSE,
  width = 10,
  height = 10,
  cellwidth = 1.5,
  cellheight = 1.5,
  show_rownames = FALSE, 
  show_colnames = FALSE
)
dev.off()

# Network -----
dist_coef <- readRDS("cbrg/rdata/pancan/asCIDER_dist_coef.rds")
metadata <- data.frame(ground_truth = df_cell$cellType,
                       batch = df_cell$Batch,
                       label = res_ascider$initial_cluster)

df <- data.frame(g = metadata$ground_truth,
                 b = metadata$batch, ## batch
                 c = metadata$label, stringsAsFactors = F) ## label

df$combination <- paste0(df$g, "-", df$c)
freq <- table(df$combination)
df$freq <- freq[match(df$combination, names(freq))]
df <- df[order(df$freq, decreasing = TRUE),]
df <- unique(df)

N <- length(unique(df$g))

edges <- data.frame(from = combinations$g1, to = combinations$g2, weight = NA)

tmp <- dist_coef + t(dist_coef)
for(i in 1:nrow(edges)) {
  edges$weight[i] <- tmp[rownames(tmp) == edges$from[i], colnames(tmp) == edges$to[i]]
}

edges$weight[edges$weight < 0] <- 0
edges <- edges[edges$weight > 0, ]
net <- graph_from_data_frame(edges, directed = FALSE)
E(net)$width <- 2 ^ (E(net)$weight * 5)
vg_names <- attr(V(net), "names")
df_cols <- data.frame(group = unique(df$g),
                      col = col_vector[1:length(unique(df$g))], stringsAsFactors = FALSE)
df$col <- df_cols$col[match(df$g, df_cols$group)]

V(net)$color <- df$col[match(vg_names, df$c)]
V(net)$frame.color <- "#777777"
V(net)$size <- 10
V(net)$label.family <- "Helvetica"

plot(net);legend(x=-1.5, y=-1.1, df_cols$group, pch=21,
                 col="#777777", pt.bg=df_cols$col, pt.cex=2, cex=.8, bty="n", ncol=1)
