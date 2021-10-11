archr.titrate.doublets.core <- function(so, cells.to.remove, max.quantile = 1, 
                                        label.vars = c("Clusters", "seurat_clusters")) {
  # cells.to.remove must observe archr naming convention
  if (is.null(so@meta.data$cells)) {
    so[["cells"]] <- get.cell.names.seurat(so = so, style = "ArchR")
  }
  vars.not.found <- label.vars %>% .[!.%in% colnames(so@meta.data)]
  if (length(vars.not.found) > 0) {
    stop(paste0(paste0(vars.not.found, collapse = ", "), " not found in so"))
  }
  pl <- utilsFanc::safelapply(label.vars, function(label.var) {
    pl <- list()
    pl$clusters <- DimPlot(so, group.by = label.var, 
                           pt.size = 0.02, label = T, label.size = 3)
    max.color <- quantile(so$DoubletEnrichment, max.quantile)
    pl$DoubletEnrich <- plot.panel.list(panel.list = "DoubletEnrichment", obj = so, order = F,
                                        assay = "RNA", limits = c(0, max.color), ident = label.var)
    so$remove <- 0
    so$remove[so$cells %in% cells.to.remove] <- 1
    
    pl$remove <- plot.panel.list(panel.list = "remove", obj = so, order = F,
                                 assay = "RNA", ident = label.var)
    so.clean <- so[, so$remove == 0]
    pl$cleaned <- plot.panel.list(panel.list = "DoubletEnrichment", obj = so.clean, order = F,
                                  assay = "RNA", limits = c(0, max.color), ident = label.var)
    pl$clustersCleaned <- DimPlot(so.clean, group.by = label.var, 
                           pt.size = 0.02, label = T, label.size = 3)
    names(pl) <- paste0(label.var, "..", names(pl))
    return(pl)
  }, threads = 2) %>% Reduce(c, .)
  return(pl)
}

archr.titrate.doublets  <- function(ao, so = NULL, plot.dir, root.name, 
                                    cells.to.remove.list = NULL, remove.frac = 0.08,
                                       label.vars = c("Clusters", "seurat_clusters"),
                                       threads = 8) {
  # cells.to.remove.list must be named.
  dir.create(plot.dir, showWarnings = F, recursive = T)
  if (!is.null(so)) {
    sao <- seurat.sync.archr(so = so, ao = ao)
    so <- sao$so
    ao <- sao$ao
    rm(sao)
    if (is.null(so@meta.data$DoubletEnrichment)) {
      so <- seurat.add.archr.meta(so = so, ao = ao, metas = "DoubletEnrichment")
    }
  }
  
  fake.so <- fake.so.gen(ao = ao, ao.embedding = "UMAP")
  
  if (!is.null(cells.to.remove.list)) {
    if (is.null(names(cells.to.remove.list)))
      stop("cells.to.remove.list, if specified, must be named")
  } else if (!is.null(remove.frac)) {
    cells.to.remove.list <- lapply(remove.frac, function(frac) {
      max.doub <- quantile(ao$DoubletEnrichment, 1 - frac)
      cells <- ao$cellNames[ao$DoubletEnrichment >= max.doub]
      return(cells)
    })
    names(cells.to.remove.list) <- as.character(remove.frac)
    
  }
  
  utilsFanc::safelapply(names(cells.to.remove.list), function(list.name) {
    cells.to.remove <- cells.to.remove.list[[list.name]]
    so.list <- list(fake.so = fake.so)
    if (!is.null(so)) {
      so.list$so <- so
    }
    
    pl <- lapply(so.list, function(so) {
      pl <- archr.titrate.doublets.core(so = so, label.vars = label.vars, cells.to.remove = cells.to.remove)
      return(pl)
    }) %>% Reduce(c,.)
    saveRDS(pl, paste0(plot.dir, "/", root.name, "_", list.name, "_pl.Rds"))
    p <- wrap.plots.fanc(plot.list = pl, plot.out = paste0(plot.dir, "/", root.name, "_", list.name, ".png"), 
                         n.col = 5)
    
    ao.new <- ao[ ! ao$cellNames %in% cells.to.remove,]
    saveRDS(ao.new, paste0(plot.dir, "/", root.name, "_", list.name, "_ao.Rds"))
    if (!is.null(so)) {
      so.new <- so[, ! so$cells %in% cells.to.remove]
      saveRDS(so.new, paste0(plot.dir, "/", root.name, "_", list.name, "_so.Rds"))
    }
    return()
  }, threads = threads)
  return()
}