plot.along.trajec <- function(trajec.df, ptime, y.vars, type.var = NULL, scale, meta.vars = NULL, no.melt = F,
                              plot.out = NULL, return.p.list = F, pt.size = 0.1,
                              meta.brewer = "Spectral") {
  n.meta <- 0
  df <- trajec.df[,c(ptime, y.vars, meta.vars, type.var)]
  if (scale == T) {
    for (v in y.vars) {
      df[, v] <- scale(df[, v])
    }
  }
  if (no.melt == T) {
    df.melt <- df
    df.melt$y <- df.melt[, y.vars]
    df.melt$type <- df.melt[, type.var]
  } else {
    df.melt <- reshape2::melt(data = df, id.vars = c(ptime, meta.vars),
                              measure.vars = y.vars, variable.name = "type",
                              value.name = "y")
  }

  p.main <- ggplot(df.melt, aes_string(x = ptime, y = "y", color = "type")) +
    geom_point(size = pt.size, alpha = 0.4) +
    geom_smooth() +
    ggtitle(paste0(y.vars, collapse = " "))
  if (!is.null(meta.vars)) {
    n.meta <- length(meta.vars)
    df.melt$y.m <- "meta"
    p.list.meta <- lapply(meta.vars, function(meta) {
      p <- ggplot(df.melt, aes_string(x = ptime, y = "y.m", color = meta)) +
        geom_point(alpha = 0.3) +
        scale_color_brewer(palette = meta.brewer)
      return(p)
    })
  }
  p.list <- c(list(p.main), p.list.meta)
  if (!is.null(meta.vars)) {
    p <- cowplot::plot_grid(plotlist = p.list, align = "v", rel_heights = c(5, rep(1, n.meta)),
                            axis = "lr", 
                            ncol = 1)
  } else {
    p <- p.main
  }
  if (!is.null(plot.out)) {
    system(paste0("mkdir -p ", dirname(plot.out)))
    ggsave(plot.out, p, width = 7, height = 5 + n.meta, dpi = 100, units = "in")
  }
  
  if (return.p.list == T) 
    return(p.list)
  else
    return(p)
  
}

plot.along.trajec.all <- function(trajec.df, ptimes, y.vars.list, plot.dir, meta.vars, ..., 
                                  threads.ptime = 4, threads.yvar = 4) {
  mclapply(ptimes, function(ptime) {
    p.list <- mclapply(names(y.vars.list), function(y.vars.name) {
      y.vars <- y.vars.list[[y.vars.name]]
      # browser()
      p <- plot.along.trajec(trajec.df = trajec.df , ptime = ptime,  
                             meta.vars = meta.vars,
                             y.vars = y.vars, ...) 
      return(p)
    }, mc.cores = threads.yvar, mc.cleanup = T)
    # browser()
    out.file <- paste0(plot.dir, "/", ptime, ".png")
    trash <- wrap.plots.fanc(plot.list = p.list, plot.out = out.file, 
                    sub.width = 7, sub.height = 5 + length(meta.vars))
    return()
  }, mc.cores = threads.ptime, mc.cleanup = T)
}

slingshot.archr <- function(ao, cluster.ident, root = NULL, ends=NULL, approx_points,
                            mat = NULL, embedding, scale.LSI = T, work.dir = NULL, root.name = NULL) {
  if (is.null(mat)) {
    if (embedding == "LSI") {
      mat <- ao@reducedDims$IterativeLSI@listData$matSVD 
    } else if (embedding == "UMAP") {
      mat <- ao@embeddings$UMAP$df
    } else {
      stop("only LSI and UMAP are supported at this point")
    }
  }

  
  cl <- getCellColData(ArchRProj = ao, select = cluster.ident)[, cluster.ident]
  lin <- getLineages(mat, cl, start.clus = root, end.clus = ends)
  crv <- getCurves(lin, approx_points = approx_points)
  
  if (is.null(work.dir)) {
    work.dir <- getOutputDirectory(ArchRProj = ao)
  }
  system(paste0("mkdir -p ", work.dir))
  if (is.null(root.name))
    root.name <- "slingshot"
  saveRDS(crv, paste0(work.dir, "/", root.name, "_", embedding, "_", cluster.ident, "_", approx_points, ".Rds"))
  return(crv)
}

archr.build.trajec.df <- function(ao, crv, values.df, metas = NULL) {
  if (!is.data.frame(crv) && !is.matrix(crv))
    crv <- slingPseudotime(crv) 
  crv.df <- crv %>% as.data.frame() %>% mutate(., cell = rownames(.))
  values.df <- values.df %>% as.data.frame() %>% mutate(., cell = rownames(.))
  trajec.df <- inner_join(crv.df, values.df)
  if (!is.null(metas)) {
    meta.df <- ArchR::getCellColData(ao, select = metas) %>% 
      as.data.frame() %>% mutate(., cell = rownames(.))
    trajec.df <- left_join(trajec.df, meta.df)
  }
  return(trajec.df)
}


slingshot.qc <- function(crv.o = NULL, mat.list = NULL, umap = NULL, ao = NULL, root.name, plot.dir) {
  if (is.null(umap)) {
    umap <- ao@embeddings$UMAP$df
  }
  umap <- as.data.frame(umap)
  colnames(umap) <- c("UMAP1", "UMAP2")
  if (is.null(mat.list)) {
    mat.list <- list(slingPseudotime = slingPseudotime(crv.o),
                     slingCurveWeights = slingCurveWeights(crv.o))
  } 
  cells.avail <- lapply(mat.list, rownames) %>% Reduce(intersect, .)
  umap <- umap[rownames(umap) %in% cells.avail, ]
  
  mat.list <- lapply(mat.list, function(mat)  {
    return(as.matrix(mat)[rownames(umap), ])
  })
  lapply(names(mat.list), function(mat.name) {
    mat <- mat.list[[mat.name]]
    p.list <- lapply(colnames(mat), function(curve) {
      p <- ggplot(umap, aes(x = UMAP1, y = UMAP2, color = mat[, curve, drop = T])) + 
        geom_point(size = 0.05) + 
        # scale_color_gradient(low = "midnightblue", high = "green") + 
        scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(5,"Spectral")),na.value = "grey70") +
        theme(legend.title = element_blank()) +
        ggtitle(curve)
      return(p)
    })
    trash <- wrap.plots.fanc(plot.list = p.list, plot.out = paste0(plot.dir, "/", root.name, "_", mat.name, ".png"))
    return()
  })
}

archr.pileup.and.trajec.df <- function(ao, curve.o, peak.mat.full, peak.cat.df, kmo = NULL, 
                                 metas = c("Sample", "Clusters")) {
  if (is.character(peak.mat.full))
    peak.mat.full <- readRDS(peak.mat.full)
  
  pileup.df <- archr.pileup(ao = ao, peak.mat.full = peak.mat.full,
                            peak.cat.df = peak.cat.df, kmo = kmo,
                            normalize = T, return.df = T)
  trajec.df <- archr.build.trajec.df(ao = ao, crv = curve.o, values.df = pileup.df, 
                                     metas = metas)
  return(trajec.df) 
}

slingshot.seurat <- function(so, cluster.ident, reducedDim = "PCA",
                             start.clus, end.clus, 
                             approx_points, plot.dir, root.name = NULL, saveRDS = T,
                             sling.rds = NULL) {
  if (is.null(sling.rds)) {
    sce <- as.SingleCellExperiment(so)
    sce <- getLineages(sce, clusterLabels = cluster.ident, reducedDim = reducedDim, 
                       start.clus = start.clus, end.clus = end.clus)
    sce <- getCurves(sce, approx_points = approx_points)
    
    sling <- SlingshotDataSet(sce)
    if (saveRDS) {
      sling.rds <- paste0(plot.dir, "/", root.name, "_slingshot.Rds")
      dir.create(plot.dir, showWarnings = F, recursive = T)
      saveRDS(sling, sling.rds)
    }
    # if (!is.null(sling.rds)) {
    #   dir.create(dirname(sling.rds), showWarnings = F, recursive = T)
    #   saveRDS(sling, sling.rds)
    # }
    
    # browser()
  } else {
    sling <- readRDS(sling.rds)
  }
  
  trash <- slingshot.qc(crv.o = sling, umap = so@reductions$umap@cell.embeddings, 
                        plot.dir = plot.dir, root.name = root.name)
  return(sling)
}

sling.cmp <- function(sling.1, sling.1.curve, df.1 = NULL,
                      sling.2, sling.2.curve, df.2 = NULL,
                      sling.1.so = NULL, sling.2.so = NULL,
                      plot.out = NULL, pt.size = 0.1) {
  
  
  if (is.null(df.1)) {
    sling.1.curve <- paste0("sling1.", sling.1.curve)
    df.1 <- slingPseudotime(sling.1) %>% as.data.frame() %>% `colnames<-`(., paste0("sling1.", colnames(.))) %>% 
      mutate(., cell = rownames(.)) %>% select(cell, !!as.name(sling.1.curve))
  }
  if (is.null(df.2)) {
    sling.2.curve <- paste0("sling2.", sling.2.curve)
    df.2 <- slingPseudotime(sling.2) %>% as.data.frame() %>% `colnames<-`(., paste0("sling2.", colnames(.))) %>% 
      mutate(., cell = rownames(.)) %>% select(cell, !!as.name(sling.2.curve))
  }
  
  if (!is.null(sling.1.so)) {
    df.1$cell <- get.cell.names.seurat(so = sling.1.so, cells = df.1$cell, style = "ArchR")
  }
  
  if (!is.null(sling.2.so)) {
    df.2$cell <- get.cell.names.seurat(so = sling.2.so, cells = df.2$cell, style = "ArchR")
  }
  df.j <- inner_join(df.1, df.2)
  p <- ggplot(df.j, aes_string(x = sling.1.curve, y = sling.2.curve)) +
    geom_point(size = pt.size, color = "red")
  if (!is.null(plot.out)) {
    dir.create(dirname(plot.out), recursive = T, showWarnings = F)
    ggsave(plot.out, p, width = 7, height = 5, units = "in", dpi = 100)
  }
  return(p)
}

sling.draw.curve <- function(slingo.list, embed.to.list = NULL,
                             color.df = NULL,
                             lineages.list = NULL,
                             width = 100, height = 120, res = 100,
                             plot.dir, root = NULL) {
  # example usage: /bar/cfan/jun/day9_only/revision3.5.1_slingshot_2023-01-05.R
  # embed.to.list: a list of embed.to
  # embed.to: a matrix of UMAP coordinates. rownames are cell names.
  # color.df: must have a column called "cell.name" that matches the cell.name in embed.to
  # also a column called "color"
  if (is.null(names(slingo.list))) {
    stop("slingo.list must be named")
  }
  
  if (!is.null(lineages.list)) {
    if (!identical(sort(names(slingo.list)), sort(names(lineages.list)))) {
      stop("names of slingo.list and lineages.list must be the same")
    }
  }
  
  if (!is.null(embed.to.list)) {
    if (!identical(sort(names(slingo.list)), sort(names(embed.to.list)))) {
      stop("names of slingo.list and lineages.list must be the same")
    }
  }
  
  if (is.null(root)) {
    root <- basename(plot.dir)
  }
  
  lapply(names(slingo.list), function(sample) {
    slingo <- slingo.list[[sample]]
    if (!is.null(lineages.list))
      lineages <- lineages.list[[sample]]
    else
      lineages <- NULL
    
    if (!is.null(embed.to.list))
      embed.to <- embed.to.list[[sample]]
    else
      embed.to <- NULL
    
    file <- paste0(plot.dir, "/", root, "_", sample, ".pdf")
    
    dir.create(dirname(file), showWarnings = F, recursive = T)
    cairo_pdf(file = file, width = width/res, height = height/res)
    try({
      color <- "gold1"
      if (!is.null(embed.to)) {
        slingo <- embedCurves(slingo, newDimRed = as.matrix(embed.to[rownames(slingo@reducedDim), 1:2]))
        to.plot <- embed.to
        if (!is.null(color.df)) {
          umap <- as.data.frame(embed.to) %>% mutate(., cell.name = rownames(.))
          umap <- umap %>% left_join(color.df, by = "cell.name")
          rownames(umap) <- umap$cell.name
          if (nrow(umap) != nrow(embed.to)) {
            stop("nrow(umap) != nrow(embed.to)")
          }
          color <- umap$color
        }
      } else {
        to.plot <- slingo
      }
      par(mai = c(0, 0, 0.2, 0))
      print(plot(to.plot, col = color, asp = 1, cex = 0.1,
                 pch = 18,  bty = "n",
                 xaxt = 'n', yaxt = 'n', ann = FALSE))
      print(plot(slingo, type = "c", add = T, linInd = lineages, lwd = 0.8, col = "black"))
      print(title(sample, cex.main = 0.5))
    })
    dev.off()
  })
  return()
}