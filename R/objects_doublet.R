dd.add.doublet <- function(so, dd.results = NULL, dd.dir = NULL, dd.root.name = NULL) {
  if (!is.null(dd.results))
    doublet.df <- dd.results$Final_doublets_groups
  else {
    doublet.file <- system(paste0("ls ", dd.dir), intern = T) %>% 
      .[grepl("Final_doublets_groups.*txt$",.)]
    if (length(doublet.file) != 1) {
      if (length(doublet.file) < 1)
        stop("no doublet file found in dir")
      doublet.file <- doublet.file %>% .[grepl(dd.root.name,.)]
      if (length(doublet.file) != 1)
        stop("length(doublet.file) != 1 after greping dd.root.name")
    }
    
    doublet.df <- read.table(doublet.file)
  }
  doublets <- rownames(doublet.df) %>% sub("\\.", "-",.)
  so@meta.data$doublet <- 0
  so@meta.data[doublets, "doublet"] <- 1
  return(so)
}

t.f.parse.doublet <- function(so, ao, sample.name = "WT", useMatrix = c("TileMatrix", "PeakMatrix"),
                              k = c(10, 50, 100),
                              nTrials = c(5, 10, 20),
                              knnMethod = c("UMAP", "LSI"),
                              plot = T, fracs = NULL, panel.list = c("seurat_clusters"),
                              order, raster = T, stat.max = T, do.highlight = F,
                              master.dir, show.removed = F,
                              plot.dir, reduction = "ArchRUMAP",
                              plot.out,
                              pt.size = 1,
                              limits = NULL) {
  so <- seurat.add.archr.embed(so = so, ao = ao,ao.embedding = "UMAP")
  so <- seurat.add.archr.embed(so = so, ao = ao,ao.embedding = "TSNE")
  so <- seurat.add.archr.meta(so = so, ao = ao)
  pl <- lapply(useMatrix, function(matrix) {
    lapply(knnMethod, function(method) {
      lapply(nTrials, function(n) {
        lapply(k, function(ki) {
          outDir <- paste0(master.dir, "/", matrix, "/", method,"/n_",n,"/k_",ki)
          # system(paste0("mkdir -p ", outDir))
          stats <- readRDS(paste0(outDir, "/", sample.name, "/", sample.name, "-Doublet-Summary.rds"))
          so.meta <- so@meta.data 
          so.meta$cells <- paste0(so.meta$sample, "#", sub("_.+$","",colnames(so)))
          scores <- c("doubletEnrichLSI", "doubletScoreLSI", "doubletEnrichUMAP", "doubletScoreUMAP")
          if (is.null(fracs)) {
            px <- lapply(scores, function(x) {
              if (plot == T) {
                df <- stats$doubletResults@listData[[x]] %>% as.data.frame()
                colnames(df) <- x
                df$cells <- rownames(df)
                so.meta <- left_join(so.meta, df)
                rownames(so.meta) <- colnames(so)
                so@meta.data <- so.meta
                p <- plot.panel.list(panel.list = x, obj = so, order = order ,
                                     assay = "SCT", return.list = T, raster = raster, pt.size = pt.size,
                                     reduction = paste0("ArchR", "UMAP"))[[1]][[1]]
                p <- p + ggtitle(paste0(matrix, "..", method, "..","n__", n, "..", "k__", ki,"..",x))
              } else {
                p <- quantile(stats$doubletResults@listData[[x]], 0.05*(10:20))
              }
              
              return(p)
            }) 
            
            if (plot == F)
              names(px) <- paste0(matrix, "..", method, "..","n__", n, "..", "k__", ki,"..",scores)
          } else {
            px <- lapply(scores, function(x) {
              pl <- lapply(fracs, function(i) {
                df <- stats$doubletResults@listData[[x]] %>% as.data.frame()
                stat <- stats$doubletResults@listData[[x]]
                colnames(df) <- x
                df$cells <- rownames(df)
                so.meta <- left_join(so.meta, df)
                rownames(so.meta) <- colnames(so)
                so@meta.data <- so.meta
                
                lim <- quantile(stat, 1-i)[1]
                max <- max(stat)
                keep <- names(stat[stat <= lim])
                
                if (show.removed == T) {
                  cells <- colnames(so)[! so.meta$cells %in% keep]
                } else {
                  cells <- colnames(so)[so.meta$cells %in% keep]
                }
                
                c8 <- colnames(so)[so.meta$seurat_clusters == "8"]
                cells <- cells[!cells %in% c8]
                panel.list <- c(x, panel.list)
                if (stat.max == T) {
                  limits <- c(0, max)
                }
                if (do.highlight == T) {
                  # p <- DimPlot(object = so, cells.highlight = cells, cells = colnames(so) %>% .[!.%in% c8],
                  #              label = T, pt.size = 0.01, reduction = "ArchRUMAP") + ggtitle(paste0("remove fraction: ", i))
                  # ps <- list(p)
                  so[["highlight"]] <- 0
                  so@meta.data[cells, "highlight"] <- 1
                  panel.list <- "highlight"
                  cells <- colnames(so) %>% .[!.%in% c8]
                } 
                p <- plot.panel.list(panel.list = panel.list, cells = cells, obj = so, order = order ,
                                     assay = "SCT", return.list = T, raster = raster, pt.size = pt.size, limits = limits,
                                     reduction = reduction)
                ps <- lapply(p, function(y) return(y[[1]]))     
                ps <- lapply(ps, function(pi) {
                  pi <- pi + ggtitle(paste0("remove fraction: ", i, " ", panel.list))
                  return(pi)
                })
                return(ps)
              }) %>% Reduce(c, .)
              plot.out <- paste0(plot.dir, "/", paste0(matrix, "..", method, "..","n__", n, "..", "k__", ki,"..",x), ".png")
              p <- wrap.plots.fanc(plot.list = pl, plot.out = plot.out, page.limit = 200, n.split = 2)
              return()
            })
            
          }
          
          
          return(px)
        }) %>% Reduce(c, .) %>% return()
      }) %>% Reduce(c, .) %>% return()
    }) %>% Reduce(c, .) %>% return()
  }) %>% Reduce(c, .) 
  
  if (plot == T) {
    p <- wrap.plots.fanc(plot.list = pl, plot.out = plot.out, page.limit = 200, n.col = 4)
  } else {
    
    p <- Reduce(rbind, pl)
    rownames(p) <- names(pl)
  }
  
  return(p)
  
}

seurat.doublet.removal <- function(so, ao=NULL, Rds=NULL, score, frac) {
  if (!is.null(Rds)) {
    stats <- readRDS(Rds)
    stat <- stats$doubletResults@listData[[score]]
  } else {
    # l <- seurat.sync.archr(so = so, ao =ao)
    # so <- l$so
    # ao <- l$ao
    if (!is.null(ao)) {
      stat <- getCellColData(ao)[, score]
      names(stat) <- rownames(ao)
    } else {
      stat <- so@meta.data[, score]
      names(stat) <- get.cell.names.seurat(so = so, style = "ArchR")
    }
    
  }
  cells.avail <- get.cell.names.seurat(so = so, style = "ArchR")
  stat <- stat[cells.avail]
  so.meta <- so@meta.data 
  so.meta$cells <- paste0(so.meta$sample, "#", sub("_.+$","",colnames(so)))
  
  lim <- quantile(stat, 1-frac)[1]
  keep <- names(stat[stat <= lim])
  cells <- colnames(so)[so.meta$cells %in% keep]
  return(so[, cells])
}

t.f.purity.test <- function(so, Rds, frac, genes, out.file = NULL) {
  pct.list <- lapply(c("RAW","LSI", "UMAP"), function(x) {
    if (x != "RAW")
      so <- seurat.doublet.removal(so = so, Rds = Rds, score = paste0("doubletEnrich", x), frac = frac)
    exp.mat <- so@assays$RNA@counts[genes, ] %>% t()
    pos.mat <- exp.mat %>% aggregate.Matrix(groupings = so$seurat_clusters, fun = "count")
    totals <- table(so$seurat_clusters)
    totals <- totals[rownames(pos.mat)]
    pos.mat <- diag(1/totals) %*% pos.mat
    df <- pos.mat %>% t() %>% as.data.frame()
    colnames(df) <- paste0(names(totals))
    df$gene <- rownames(df)
    df.melt <- reshape2::melt(df, id.vars = "gene", variable.name = "cluster", value.name = paste0("pct.", x))
    df.melt <- df.melt %>% mutate(gene = factor(gene, levels = genes)) %>% arrange(gene)
    return(df.melt)
  })
  df <- Reduce(left_join, pct.list) %>% arrange(gene, pct.RAW)
  if (!is.null(out.file)) {
    write.table(df, out.file, quote = F, row.names = F, col.names = T, sep = "\t")
  }
  return(df)  
}