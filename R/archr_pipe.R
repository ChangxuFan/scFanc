add.umap.fanc <- function(ao) {
  ao.Rds <- tempfile()
  saveRDS(ao, ao.Rds)
  cmd <- paste0("/bar/cfan/R_for_bash/umap_archr.R ", ao.Rds, " ", ao.Rds, " mm10")
  print(cmd); system(cmd)
  ao <- readRDS(ao.Rds)
  return(ao)

}
ao.gen <- function(frag.files = NULL, arrow.files = NULL, copy.arrow.files = "miao", force = F,
                   minTSS = 4, minFrags = 1000, addGeneScoreMat = T, addTileMat = T,
                   cells.list = NULL,  work.dir) {
  # cells: should be a list, each element corresponding to the cells from each arrow file
  # frag files need to be named as sample name
  
  system(paste0("mkdir -p ", work.dir))
  
  if (sum(grepl(".arrow", list.files(work.dir))) > 0 && force != T && !is.null(frag.files)) {
    stop("there are already arrow files work.dir ")
  }
  
  if (!is.null(frag.files) && is.null(arrow.files)) {
    wd.bk <- getwd()
    setwd(work.dir)
    try(
      arrow.files <- mclapply(seq_along(frag.files), function(i) {
        if (!is.null(cells.list[[i]])) {
          minTSS <- 0
          minFrags <- 0
        }
        arrow <- createArrowFiles(
          inputFiles = frag.files[i],
          sampleNames = names(frag.files[i]),
          # filterTSS = filterTSS, #Dont set this too high because you can always increase later
          # filterFrags = filterFrags, 
          addTileMat = addTileMat,
          addGeneScoreMat = addGeneScoreMat, 
          validBarcodes = cells.list[[i]], 
          minTSS = minTSS,
          minFrags = minFrags,
          maxFrags = 10000000,
          force = force
        )
      }) %>% unlist()
    )
  
    try(print(arrow.files))
    
    try(doubScores <- addDoubletScores(
      input = arrow.files,
      # useMatrix = "TileMatrix",
      # k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
      # nTrials = 5,
      # knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
      # LSIMethod = 1,
      # dimsToUse = 1:30, 
      # scaleDims = T, # this is F by default. changed to T by FANC
      # outDir = paste0(work.dir, "/", getOutputDirectory(arrow.files))[1], 
      force = T
    ))
    
    setwd(wd.bk)
  } else {
    # stop("the option of starting from arrow.files has not been tested")
  }

  if (copy.arrow.files %in% c(T,F)) {
    copy <- copy.arrow.files
  } else {
    if (!is.null(frag.files))
      copy <- F
    else
      copy <- T
    
  }
    
  arrow.files <- paste0(work.dir, "/", arrow.files)
  
  ao <- ArchRProject(
    ArrowFiles = arrow.files, 
    outputDirectory = work.dir,
    copyArrows = copy, showLogo = F
  )
  ao@projectMetadata$outputDirectory <- work.dir
  saveRDS(ao, paste0(work.dir, "/ao.Rds"))
  return(ao)
}

archr.cluster.pipe <- function(ao, so = NULL, plot.only = F,
                               minTSS = 10, 
                               bc.metrics.file.list = NULL,
                               use.matrix = "TileMatrix",
                               add.lsi = T, lsi.iterations = 2, lsi.res = 0.2, do.harmony = F,
                               add.cluster = T, add.umap = T,
                               dimsLSI = 1:30,
                               dimsToUse = NULL, dimsToUse.umap = NULL, dimsToUse.cluster = NULL, cluster.res = 0.8,
                               work.dir, plot.dir = NULL, 
                               use.reduced.dims = c("IterativeLSI", "Harmony"),
                               so.cluster.ident = "seurat_clusters", so.clusters = NULL, 
                               do.confusion.plot = T, do.binary.plot = T,
                               force = T,
                               sample.order = NULL, ...) {
  if (!is.null(dimsToUse)) {
    dimsToUse.umap <- dimsToUse
    dimsToUse.cluster <- dimsToUse
  }
  if (is.null(plot.dir))
    plot.dir <- paste0(work.dir, "/Plots")
  system(paste0("mkdir -p ", work.dir, " ", plot.dir))
  if (plot.only == F) {
    ao <- ao[ao$TSSEnrichment >= minTSS,]
    if (add.lsi == T) {
      utilsFanc::t.stat("adding LSI")
      ao <- addIterativeLSI(
        ArchRProj = ao,
        useMatrix = use.matrix, 
        name = "IterativeLSI", 
        iterations = lsi.iterations, 
        clusterParams = list( #See Seurat::FindClusters
          resolution = lsi.res, 
          # sampleCells = 10000, 
          n.start = 10
        ), 
        dimsToUse = dimsLSI,
        force = force, outDir = work.dir, ...
      )
    }
    
    if (do.harmony == T && length(ArchR::getSampleNames(ao)) > 1 ) {
      utilsFanc::t.stat("adding Harmony")
      ao <- addHarmony(
        ArchRProj = ao,
        reducedDims = "IterativeLSI",
        name = "Harmony",
        groupBy = "Sample",
        force = force
      )
    }
    
    if (!is.null(so)) {
      utilsFanc::t.stat("adding seurat clusters")
      ao <- archr.add.seurat(ao = ao, so = so, meta = so.cluster.ident, as.factor = T)
    }
  }

  reducedDims <- intersect(names(ao@reducedDims), use.reduced.dims)
  if (length(reducedDims) < 1) {
    stop("no reducedDims to use")
  } 
  
  aol <- lapply(reducedDims, function(reducedDim) {
    if (plot.only == F) {
      if (add.cluster == T) {
        utilsFanc::t.stat("adding clusters")
        ao <- addClusters(
          input = ao,
          reducedDims = reducedDim,
          method = "Seurat",
          name = "Clusters",
          resolution = cluster.res,
          force = force, 
          dimsToUse = dimsToUse.cluster, 
        )
        saveRDS(ao, paste0(work.dir,"/ao_",reducedDim,".Rds"))
      }
      
      if (add.umap == T) {
        utilsFanc::t.stat("adding umap")
        ao <- addUMAP(
          ArchRProj = ao, 
          reducedDims = reducedDim, 
          name = "UMAP", 
          nNeighbors = 30, 
          minDist = 0.5, 
          metric = "cosine",
          force = force,
          saveModel = F, 
          seed = 122,
          dimsToUse = dimsToUse.umap,
        )
        saveRDS(ao, paste0(work.dir,"/ao_",reducedDim,".Rds"))
      }
      
    }
    
    utilsFanc::t.stat("plotting")
    # plotEmbedding(ArchRProj = ao, colorBy = "cellColData", name = "Clusters", embedding = "UMAP") +
    #   ggsave(paste0(plot.dir, "/umap_atac_",reducedDim,".png"))
    
    fake.so <- fake.so.gen(ao = ao, ao.embedding = "UMAP")
    seurat.plot.archr(ao = ao, ao.embedding = "UMAP", ao.name = "Clusters", plot.dir = plot.dir,
                      shuffle = T, root.name = reducedDim)
    
    if (do.binary.plot == T) {
      plot.panel.list(panel.list = "Clusters", fake.so, order = T,assay = "RNA", binarize.panel = T, 
                      split.by = "Sample", split.order = sample.order,
                      raster = T, auto.adjust.raster.pt.size = F, pt.size = 0.5, 
                      plot.out = paste0(plot.dir, "/", reducedDim,"_UMAP_Clusters_binary.png"), 
                      page.limit = 100, ident = "Clusters")
      
    }
    
    if (length(unique(ao$Sample)) > 1) {
      seurat.plot.archr(ao = ao, ao.embedding = "UMAP", ao.name = "Sample", plot.dir = plot.dir,
                        shuffle = T, label = F, root.name = reducedDim)
      if (do.binary.plot == T) {
        plot.panel.list(panel.list = "Sample", fake.so, order = T,assay = "RNA", binarize.panel = T, 
                        raster = T, auto.adjust.raster.pt.size = F, pt.size = 0.5, 
                        plot.out = paste0(plot.dir, "/", reducedDim,"_UMAP_Sample_binary.png"), 
                        page.limit = 100, ident = "Clusters")
      }
    }
      
    try({
      if (!is.null(bc.metrics.file.list)) {
        vln.depth.2(so = NULL, ao = ao, bc.metrics.file.list = bc.metrics.file.list, 
                    metas.plot = NULL,
                    ident = "Clusters",
                    plot.out = paste0(plot.dir, "/metrics_Clusters_",reducedDim, ".pdf"),
                    sub.width = 7.5, violin = T, page.limit = 4,
                    max.quantile = 0.99, n.col = 4, numeric.only = T, sub.height = 3)
      }
      
    })

    try({
      p2 <- seurat.plot.archr(ao = ao, 
                              ao.name = so.cluster.ident, ao.embedding = "UMAP", plot.dir = plot.dir, 
                              root.name = reducedDim)
    })
    if (do.confusion.plot == T) {
      try({
        fake.so <- fake.so.gen(ao = ao, ao.embedding = "UMAP")
        trash <- plot.panel.list(panel.list = so.cluster.ident, fake.so, order = F,assay = "RNA", binarize.panel = T, 
                        split.by = "Sample", raster = T, auto.adjust.raster.pt.size = F, pt.size = 1.2, 
                        split.order = sample.order,
                        plot.out = paste0(plot.dir, "/", reducedDim, "_", so.cluster.ident, "_binary.png"), 
                        page.limit = 100, binarize.items = so.clusters)
      })  
    }
  })
  if (length(aol) == 1)
    aol <- aol[[1]]
  
  return(aol)
  
}

archr.crosscheck.pipe <- function(ao, so, study.dir, samples, bc.metrics.file.list) {
  t.f.ao.so.qc(ao = ao, so = so, plot.dir = study.dir)
  
  # ao <- archr.add.seurat(ao = ao, so = so, as.factor = T)
  # seurat.plot.archr(ao = ao, ao.embedding = "UMAP", ao.name = "seurat_clusters", 
  #                   plot.dir = paste0(study.dir, "/"))
  # 
  # seurat.plot.archr(ao = ao[ ! is.na(ao$seurat_clusters),], ao.embedding = "UMAP", ao.name = "seurat_clusters", 
  #                   plot.dir = paste0(study.dir, "/"), root.name = "ao_both")
  trash <- vln.depth.2(so = NULL, ao = ao, bc.metrics.file.list = NULL, 
              metas.plot = NULL,
              ident = "Clusters",
              plot.out = paste0(study.dir, "/metrics","/metrics.pdf"), 
              sub.width = 7.5, violin = T, page.limit = 4,
              max.quantile = 0.99, n.col = 4, numeric.only = T, sub.height = 3)
  return()
}

archr.macs2.pipe <- function(ao, work.dir, cluster.ident = "seurat_clusters", peak.annot = "homer",
                             add.coverage = T, add.peakset = T, add.peakmat = T, get.peakmat = T, add.motifanno = T) {
  system(paste0("mkdir -p ", work.dir))
  if (add.coverage == T) {
    ao <- addGroupCoverages(ArchRProj = ao, groupBy = cluster.ident)
    saveRDS(ao, paste0(work.dir, "/ao.Rds"))
  }
  
  if (add.peakset == T) {
    ao <- addReproduciblePeakSet(
      ArchRProj = ao, 
      groupBy = cluster.ident, 
      pathToMacs2 = "/opt/apps/python2/bin/macs2"
    )
    
    saveRDS(ao, paste0(work.dir, "/ao.Rds"))
  }

  if (add.peakmat == T) {
    ao <- addPeakMatrix(ao)
    saveRDS(ao, paste0(work.dir, "/ao.Rds"))
  }
  
  if (get.peakmat == T) {
    peakmat <- archr.get.peak.mat(ao, mat.name = "PeakMatrix", return.full = T)
    saveRDS(peakmat, paste0(work.dir, "/peakmat_full.Rds"))
  }
  
  if (add.motifanno == T) {
    ao <- addMotifAnnotations(ArchRProj = ao, motifSet = peak.annot, name = peak.annot)
    saveRDS(ao, paste0(work.dir, "/ao.Rds"))
  }
  
  

  # saveRDS(ao, paste0(work.dir, "/ao.Rds"))
  
  return(ao)
  
}

archr.hm.pipe <- function(ao, ao.markersPeaks = NULL, cluster.ident = "seurat_clusters", work.dir, plot.marker.peaks=T, plot.motif=T, peak.annot = "homer",
                          marker.width = 700, marker.height = 600, motif.width = 600 , motif.height = 700, 
                          motif.n = 7, motif.pmax = 50, 
                          plot.dir = NULL) {
  if (is.null(plot.dir))
    plot.dir <- paste0(work.dir, "/Plots")
  system(paste0("mkdir -p ", work.dir, " ", plot.dir))
  
  if (is.null(ao.markersPeaks)) {
    ao.markersPeaks <- getMarkerFeatures(
      ArchRProj = ao, 
      useMatrix = "PeakMatrix", 
      groupBy = cluster.ident,
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon"
    )
    saveRDS(ao.markersPeaks, paste0(work.dir, "/marker_peaks_", Sys.time() %>% sub(" ", "_",.), ".Rds"))
  } else {
    if (is.character(ao.markersPeaks))
      ao.markersPeaks <- readRDS(ao.markersPeaks)
  }

  
  if (plot.marker.peaks == T) {
    utilsFanc::t.stat("preparing to plot marker peaks")
    ao.heatmapPeaks <- plotMarkerHeatmap(
      seMarker = ao.markersPeaks, 
      cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
      transpose = F, 
      labelRows = F,# , nLabel = 0, nPrint = 1
      labelMarkers = NULL,
      nLabel = 1
    )
    
    png(filename = paste0(plot.dir, "/hm_marker_peaks.png"),
        width = marker.width, height = marker.height, res = 100)
    try(print(draw(ao.heatmapPeaks, heatmap_legend_side = "right", annotation_legend_side = "right")))
    dev.off()
  }
  
  if (plot.motif == T) {
    utilsFanc::t.stat("preparing to plot marker motifs")
    ao.enrichMotifs <- peakAnnoEnrichment(
      seMarker = ao.markersPeaks,
      ArchRProj = ao,
      peakAnnotation = peak.annot,
      cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
    )
    
    ao.heatmapEM <- plotEnrichHeatmap(ao.enrichMotifs, n = motif.n, transpose = F, pMax = motif.pmax)
    
    png(filename = paste0(plot.dir, "/hm_",peak.annot,".png"), width = motif.width, height =motif.height, res = 100)
    try(print(ComplexHeatmap::draw(ao.heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")))
    dev.off()
    

  }
  
  return(ao)
  
}



t.f.archr.lsi.resolution <- function(ao = ao, resolutions, plot.dir) {
  mclapply(resolutions, function(res) {
    out.dir <-  paste0(plot.dir, "/LSI_res_", res, "/")
    dir.create(out.dir, recursive = T, showWarnings = F)
    
    # ao <- addIterativeLSI(
    #   ArchRProj = ao,
    #   useMatrix = "TileMatrix", 
    #   name = "IterativeLSI", 
    #   iterations = 2, 
    #   clusterParams = list( #See Seurat::FindClusters
    #     resolution = c(res),
    #     # sampleCells = 10000,
    #     n.start = 10
    #   ),
    #   varFeatures = 25000, 
    #   dimsToUse = 1:30,
    #   force = T, 
    #   threads = 4
    # )
    # 
    # saveRDS(ao, paste0(out.dir, "/ao.Rds"))
    
    ao <- readRDS(paste0(out.dir, "/ao.Rds"))
    so <- readRDS("lodge/rna_check/per_sample/tumor_rep3/soi.Rds")
    ao <- archr.add.seurat(ao = ao, so = so, as.factor = T)
    saveRDS(ao, paste0(out.dir, "/ao.Rds"))
    t.f.archr.umap.seed(ao = ao, seeds = 10*(1:10), plot.dir = out.dir, 
                        root.name = "LSI_res_" %>% paste0(res))
    return()
  }, mc.cores = 1, mc.cleanup = T)
  return()
}

t.f.archr.umap.seed <- function(ao, seeds, plot.dir, root.name) {
  mclapply(seeds, function(seed) {
    ao <- addUMAP(
      ArchRProj = ao, 
      reducedDims = "IterativeLSI", 
      name = "UMAP", 
      nNeighbors = 30, 
      minDist = 0.5, 
      metric = "cosine",
      force = T,
      saveModel = F, 
      seed = seed, threads = 4
    )
    browser()
    p <- seurat.plot.archr(ao = ao,
                           ao.name = "seurat_clusters", ao.embedding = "UMAP", plot.dir = plot.dir,
                           root.name = paste0(root.name, "_seed_", seed))
  }, mc.cores = 4, mc.cleanup = T)
}

t.f.tss.nfrag.qc <- function(tss.nfrag.df, outfile) {
  p <- ggPoint(x = tss.nfrag.df$nfrag, y = tss.nfrag.df$tss, colorDensity = TRUE,
               continuousSet = "sambaNight",
               xlabel = "Log10 Unique Fragments",
               ylabel = "TSS Enrichment",
               xlim = c(log10(500), quantile(tss.nfrag.df$nfrag, probs = 0.99)),
               ylim = c(0, quantile(tss.nfrag.df$tss, probs = 0.99))) +
    ggsave(outfile, width = 5, height = 7, units = "in")
  return(p)
}

tss.nfrag.qc <- function(ao, tss.nfrag.df = NULL, outfile = NULL, 
                         polygon.df = NULL) {
  if (is.null(tss.nfrag.df)) {
    tss.nfrag.df <- data.frame(tss = ao$TSSEnrichment, nfrag = log10(ao$nFrags), 
                               cells = getCellNames(ao))
    
    rownames(tss.nfrag.df) <- tss.nfrag.df$cells
  }
  p <- ggPoint(x = tss.nfrag.df$nfrag, y = tss.nfrag.df$tss, colorDensity = TRUE,
               continuousSet = "sambaNight",
               xlabel = "Log10 Unique Fragments",
               ylabel = "TSS Enrichment",
               xlim = c(log10(500), quantile(tss.nfrag.df$nfrag, probs = 0.99)),
               ylim = c(0, quantile(tss.nfrag.df$tss, probs = 0.99)))
  if (!is.null(polygon.df)) {
    p <- p + geom_polygon(data = polygon.df, mapping = aes(x = x, y = y), 
                     fill = NA, color = "red")
  }
  if (!is.null(outfile)) {
    dir.create(dirname(outfile),showWarnings = F, recursive = T)
    ggsave(outfile, p, width = 5, height = 7, units = "in")
  }
    
  return(p)
}  

subset.cells.in.polygon <- function(ao, tss.nfrag.df = NULL, polygon.df,
                                    key, reduce.pct, keep.higher = T,
                                    plot.dir, root.name,
                                    umap.embedding = "UMAP") {
  suffix <- paste0(key, "_", reduce.pct)
  root.name <- paste0(root.name, "_", suffix)
  if (is.null(tss.nfrag.df)) {
    tss.nfrag.df <- data.frame(tss = ao$TSSEnrichment, nfrag = log10(ao$nFrags), 
                               cells = getCellNames(ao))
    rownames(tss.nfrag.df) <- tss.nfrag.df$cells
  }
  cells.in.poly <- select.cells.by.polygon(embed.df = tss.nfrag.df, x = "nfrag",
                                           y = "tss", polygon.df = polygon.df)
  ao.sub <- ao[ao$cellNames %in% cells.in.poly,]
  bRemove <- getCellColData(ArchRProj = ao.sub, select = key, drop = T) %>% 
    select.by.quantile(frac.cutoff = reduce.pct, na.method = "toF", 
                       larger.than = !keep.higher)
  # we do not remove NA elements
  cells.remove <- ao.sub$cellNames[bRemove]
  cells.left.in.poly <- cells.in.poly %>% .[!.%in% cells.remove]
  rm(ao.sub)
  ao.2 <- ao[! ao$cellNames %in% cells.remove,]
  
  stats.df <- data.frame(total.cell = ao$cellNames %>% length(), 
                         n.in.poly = length(cells.in.poly),
                         n.removed = length(cells.remove),
                         n.left.poly = length(cells.left.in.poly))
  trash <- utilsFanc::write.zip.fanc(df = stats.df, out.file = paste0(plot.dir, "/", root.name, "_stats.tsv"), 
                            zip = F, col.names = T, row.names = F)
  tss.nfrag.df.new <- tss.nfrag.df %>% filter(!cells %in% cells.remove)
  p1 <- tss.nfrag.qc(tss.nfrag.df = tss.nfrag.df,
                        outfile = paste0(plot.dir, "/", root.name, "_tss_nfrag_pre_rm.png"), 
                        polygon.df = polygon.df)
  
  p2 <- tss.nfrag.qc(tss.nfrag.df = tss.nfrag.df.new,
                        outfile = paste0(plot.dir, "/", root.name, "_tss_nfrag_post_rm.png"), 
                        polygon.df = polygon.df)
  trash <- wrap.plots.fanc(plot.list = list(p1, p2),
                           plot.out = paste0(plot.dir, "/", root.name, "_tss_nfrag.png"))
  
  p1 <- plot.panel.list(obj = ao, order = F, assay = "RNA",
                           highlight.list = list(removed = cells.remove), 
                           # plot.out = paste0(plot.dir, "/", root.name, "_umap_rm.png"),
                           ao.embedding = umap.embedding)
  p2 <- plot.panel.list(obj = ao.2, order = F, assay = "RNA",
                           highlight.list = list(left = cells.left.in.poly), 
                          # plot.out = paste0(plot.dir, "/", root.name, "_umap_left.png"),
                           ao.embedding = umap.embedding)
  trash <- wrap.plots.fanc(plot.list = list(p1, p2),
                           plot.out = paste0(plot.dir, "/", root.name, "_umap.png"))
  res.list <- list(cells.in.poly = cells.in.poly, cells.remove = cells.remove,
                   cells.left.in.poly = cells.left.in.poly)
  saveRDS(res.list, paste0(plot.dir, "/", root.name, ".Rds"))
  return(res.list)
}

t.f.ao.so.qc <- function(ao, so, plot.dir) {
  dir.create(plot.dir, showWarnings = F, recursive = T)
  stats <- data.frame(n.ao = length(getCellNames(ao)),
                      n.so = ncol(so))
  so[["cells"]] <- get.cell.names.seurat(so = so, style = "ArchR")
  both <- intersect(getCellNames(ao), so$cells)
  stats$n.both <- length(both)
  write.table(stats, paste0(plot.dir, "/stats.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
  
  tss.nfrag.df <- data.frame(tss = ao$TSSEnrichment, nfrag = log10(ao$nFrags), 
                             cells = getCellNames(ao))
  t.f.tss.nfrag.qc(tss.nfrag.df = tss.nfrag.df, outfile = paste0(plot.dir, "/tss_nfrag_all.png"))
  tss.nfrag.df.both <- tss.nfrag.df %>% filter(cells %in% both)
  t.f.tss.nfrag.qc(tss.nfrag.df = tss.nfrag.df.both, outfile = paste0(plot.dir, "/tss_nfrag_both.png"))
  ao <- ao[getCellNames(ao) %in% both,]
  # ao <- archr.cluster.pipe(ao = ao, so = so, work.dir = paste0(plot.dir, "/both"), force = T)
  
  return(ao)
}

