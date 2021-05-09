add.umap.fanc <- function(ao) {
  ao.Rds <- tempfile()
  saveRDS(ao, ao.Rds)
  cmd <- paste0("/bar/cfan/R_for_bash/umap_archr.R ", ao.Rds, " ", ao.Rds, " mm10")
  print(cmd); system(cmd)
  ao <- readRDS(ao.Rds)
  return(ao)

}
ao.gen <- function(frag.files = NULL, arrow.files = NULL, copy.arrow.files = "miao", force = F,
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
        arrow <- createArrowFiles(
          inputFiles = frag.files[i],
          sampleNames = names(frag.files[i]),
          filterTSS = 0, #Dont set this too high because you can always increase later
          filterFrags = 0, 
          addTileMat = TRUE,
          addGeneScoreMat = TRUE, 
          validBarcodes = cells.list[[i]], 
          minFrags = 0,
          maxFrags = 10000000,
          force = force
        )
      }) %>% unlist()
    )
    try(print(arrow.files))
    
    try(doubScores <- addDoubletScores(
      input = arrow.files,
      k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
      knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
      LSIMethod = 1
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
    

  ao <- ArchRProject(
    ArrowFiles = arrow.files, 
    outputDirectory = work.dir,
    copyArrows = copy, showLogo = F
  )
  saveRDS(ao, paste0(work.dir, "/ao.Rds"))
  return(ao)
}

archr.cluster.pipe <- function(ao, so = NULL, work.dir, plot.dir = NULL, do.harmony = T, force = F) {
  if (is.null(plot.dir))
    plot.dir <- paste0(work.dir, "/Plots")
  system(paste0("mkdir -p ", work.dir, " ", plot.dir))
  utilsFanc::t.stat("adding LSI")
  ao <- addIterativeLSI(
    ArchRProj = ao,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
      resolution = c(0.2), 
      # sampleCells = 10000, 
      n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = force
  )
  
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
  utilsFanc::t.stat("adding clusters")
  ao <- addClusters(
    input = ao,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force = force
  )
  saveRDS(ao, paste0(work.dir, "/ao.Rds"))
  utilsFanc::t.stat("adding umap")
  ao <- addUMAP(
    ArchRProj = ao, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = force
  )
  # ao <- add.umap.fanc(ao)
  
  
  if (!is.null(so)) {
    utilsFanc::t.stat("adding seurat clusters")
    ao <- archr.add.seurat(ao = ao, so = so, meta = "seurat_clusters")
  }
  
  saveRDS(ao, paste0(work.dir, "/ao.Rds"))
  utilsFanc::t.stat("plotting")
  plotEmbedding(ArchRProj = ao, colorBy = "cellColData", name = "Clusters", embedding = "UMAP") +
    ggsave(paste0(plot.dir, "/umap_atac.png"))
  
  if (!is.null(ao)) {
    # p1 <- seurat.plot.archr(so = so, ao = ao, ao.colorBy = "cellColData",
    #                         ao.name = "Clusters", ao.embedding = "UMAP", plot.dir = plot.dir)
    p2 <- seurat.plot.archr(so = so, ao = ao, ao.colorBy = "cellColData",
                            ao.name = "seurat_clusters", ao.embedding = "UMAP", plot.dir = plot.dir)
    # trash <- scFanc::wrap.plots.fanc(plot.list = list(p1, p2), plot.out = paste0(plot.dir, "/seurat_plot_archr.png"))
  }
  return(ao)

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
    saveRDS(peakmat, paste0(work.dir, "/peakmat.Rds"))
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


# archr.pipe.master <- function(ao = NULL, so = NULL, frag.files = NULL, arrow.files = NULL, copy.arrow.files = "miao", force = F,
#                               cells.list = NULL,  work.dir, plot.dir = NULL, do.harmony = T, )