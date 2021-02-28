
cluster.subset <- function(so, clusters) {
  if (!is.character(clusters))
    stop("cluster names must be characters")
  
  meta <- so@meta.data
  meta$seurat_clusters <- as.character(meta$seurat_clusters)
  cells <- rownames(meta)[meta$seurat_clusters %in% clusters]
  return(so[, cells])
}

count.by.cluster <- function(so, work.dir=NULL) {
  df <- so@meta.data$seurat_clusters %>% table() %>% as.data.frame()
  colnames(df) <- c("cluster", "freq")
  df$frac <- df$freq/sum(df$freq)
  sum <- df
  sum$cluster <- NULL
  sum <- lapply(sum, sum) %>% as.data.frame()
  sum$cluster <- "sum"
  df <- rbind(df, sum)
  
  if (!is.null(work.dir))
    write.table(df, work.dir %>% paste0("/cluster_freq.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
  return(df)
}


get.cell.names.seurat <- function(so, filter.list, style = "Seurat") {
  meta.df <- so@meta.data %>% scFanc::factor2character.fanc()
  for (i in 1:length(filter.list)) {
    meta <- names(filter.list)[i]
    values <- filter.list[[i]]
    meta.df <- meta.df %>% .[.[, meta] %in% values,]
  }
  if (style == "Seurat") {
    return(rownames(meta.df))
  }
  
  if (style == "ArchR") {
    return(paste0(meta.df$sample, "#", sub("_.+$", "", rownames(meta.df))))
  }
}

add.metrics.seurat <- function(so, bc.metrics.file.list, columns.to.add) {
  # bc.metrics.file.list should be a named list like this:
  ##
  meta.df <- so@meta.data %>% scFanc::factor2character.fanc()
  meta.df$barcode <- rownames(meta.df) %>% sub("_.+$", "", .)
  metrics.df <- meta.df %>% split(., f = factor(.$sample, levels = unique(.$sample)) ) %>% 
    lapply(function(x) {
      sample <- x$sample[1]
      metrics <- read.csv(bc.metrics.file.list[[sample]])
      metrics <- metrics[, c("barcode", columns.to.add)]
      x <- left_join(x, metrics)
      return(x[, c("barcode", columns.to.add)])
    }) %>% `names<-`(NULL) %>% Reduce(rbind, .)
  
  for (i in columns.to.add) {
    so[[i]] <- metrics.df[, i]
  }
  return(so)
}

vln.depth <- function(so, bc.metrics.file.list=NULL, plot.dir, metas = c("nCount_RNA", "atac_fragments"), return.so=F) {
  p.list <- lapply(metas, function(meta) {
    if (! meta %in% colnames(so@meta.data) ) {
      so <- add.metrics.seurat(so, bc.metrics.file.list = bc.metrics.file.list, columns.to.add = meta)
    }
    
    p <- VlnPlot(object = so,features = meta)
    return(p)
  }) 
  p <- scFanc::wrap.plots.fanc(plot.list = p.list, n.col = 1, sub.width = 10, sub.height = 3, plot.out = paste0(plot.dir, "/depth.png"))
  
  if (return.so == T) {
    return(so)
  }
  return()
}
