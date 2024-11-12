marker.peak.scatter <- function(a2b.list, mpl, cluster.ident, mpl.clusters = NULL, comp, 
                                plot.out = NULL, peakwatch.dir = NULL, highlight.only=F,
                                add.label = F, threads.a2b = 1, threads.mpl = 1, ...) {
  sample.y <- parse.comp(comp, 1)
  sample.x <- parse.comp(comp, 2)
  if (is.null(mpl.clusters))
    mpl.clusters <- names(mpl)
  
  mpl <- mpl[mpl.clusters]
  p.list <- mclapply(names(a2b.list), function(cluster.name) {
    cluster <- sub(paste0(cluster.ident, "_"), "", cluster.name)
    a2b <- a2b.list[[cluster.name]]
    p.list.sub <- mclapply(names(mpl), function(mpl.cluster) {
      utilsFanc::t.stat(paste0("plotting: ", cluster.name, " and the marker peaks for ", mpl.cluster))
      # if (cluster == "0" )
      #   browser()
      mp <- mpl[[mpl.cluster]]
      mp <- mp %>% `names<-`(NULL) %>% as.data.frame()
      mp <- mp %>% mutate(gene = paste0(seqnames, ":", start, "-", end)) %>% 
        select(gene, Log2FC, FDR, MeanDiff)
      bulkNorm <- a2b$bulkNorm
      bulkNorm <- bulkNorm %>% mutate(id = 1:nrow(bulkNorm)) %>% 
        mutate(plotly = paste0(id, "..", gene))
      if (highlight.only == T) {
        bulkNorm <- bulkNorm %>% filter(gene %in% mp$gene)
      }
      if (!is.null(peakwatch.dir)) {
        colnames(mp)[2:length(colnames(mp))] <- paste0("mp..", colnames(mp)[2:length(colnames(mp))])
        out.root <- paste0(peakwatch.dir, "/", cluster.name, "..mp.", mpl.cluster)
        out.bed <- paste0(out.root, ".bed")
        out.tsv <- paste0(out.root, ".tsv")
        pw.df <- left_join(mp, bulkNorm) %>% utilsFanc::loci.2.df(loci.col.name = "gene", return.gr = F,
                                                                         out.bed = out.bed)
        trash <- utilsFanc::write.zip.fanc(df = pw.df, out.file = out.tsv, zip = F, col.names = T, row.names = F)
        
      }
      label.var <- NULL
      if (add.label == T)
        label.var <- "id"
      p <- xy.plot(df = bulkNorm, x = paste0("bulkNorm_", sample.x), y = paste0("bulkNorm_", sample.y), 
                   highlight.var = "gene", highlight.values = mp$gene, plotly.var = "plotly", label.var = label.var) +
        ggtitle(paste0(cluster.name, "  markerPeak_", mpl.cluster))

      return(p)
    }, mc.cores = threads.mpl, mc.cleanup = T)
    names(p.list.sub) <- names(mpl)
    return(p.list.sub)
  }, mc.cores = threads.a2b, mc.cleanup = T) 
  names(p.list) <- names(a2b.list)
  if (!is.null(plot.out)) {
    p.list.melt <- Reduce(c, p.list)
    trash <- wrap.plots.fanc(p.list.melt, plot.out = plot.out, ...)
  }
  return(p.list)
}

marker.peak.hm <- function(ao, mp.mat, cluster.ident, anchor.clusters, other.clusters = NULL, 
                           filter = "FDR <= 0.1 & Log2FC >= 0.5", plot.dir,
                           width = 700, height = 600) {
  mp <- getMarkers(mp.mat, returnGR = F, cutOff = filter)
  peak.df <- as.data.frame(mp.mat@elementMetadata) %>% mutate(peak = paste0(seqnames, ":", start, "-", end))
  lapply(anchor.clusters, function(cluster) {
    markers <- mp[[cluster]] %>% as.data.frame() %>% mutate(peak = paste0(seqnames, ":", start, "-", end))
    bMarkers <- peak.df$peak %in% markers$peak
    mp.mat.sub <- mp.mat[bMarkers, ]
    if (!is.null(other.clusters)) {
      mp.mat.sub <- mp.mat.sub[, c(anchor.clusters, other.clusters) %>% unique()]
    }
    utilsFanc::t.stat("preparing to plot marker peaks")
    hm <- plotMarkerHeatmap(
      seMarker = mp.mat.sub, 
      cutOff = filter,
      transpose = F, 
      labelRows = F,# , nLabel = 0, nPrint = 1
      labelMarkers = NULL,
      nLabel = 1
    )
    system(paste0("mkdir -p ", plot.dir))
    plot.out <- paste0(plot.dir, "/marker_", cluster, ".png")
    
    png(filename = plot.out,
        width = width, height = height, res = 100)
    try(print(draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")))
    
    dev.off()
    return(NULL)
  })
  return()
}

plot.mat.grouping <- function(mat, x.grouping, x.order,
                              plot.out, width = 700, height=700) {
  ha = HeatmapAnnotation(grouping = anno_simple(x.grouping %>% as.character()),
    annotation_name_side = "top", which = "row", annotation_name_rot = 0, show_annotation_name = F)
  
  hm <- ComplexHeatmap::Heatmap(matrix = mat, show_column_names = T, show_row_names = F,
                          show_row_dend = F, show_column_dend = F, cluster_columns = T,
                          right_annotation = ha, row_order = x.order)
  png(filename = plot.out,
      width = width, height = height, res = 100)
  try(print(draw(hm)))
  dev.off()
}

# plot.mat.grouping <- function(mat, x.grouping, plot.out, width = 700, height=700, ...) {
#   ha = HeatmapAnnotation(grouping = anno_simple(x.grouping %>% as.character()),
#                          annotation_name_side = "top", which = "row",
#                          annotation_name_rot = 0, show_annotation_name = F)
#   hm <- ComplexHeatmap::Heatmap(matrix = mat, show_column_names = T, show_row_names = F,
#                                 show_row_dend = F, show_column_dend = F, cluster_columns = T,
#                                 right_annotation = ha, )
#   png(filename = plot.out,
#       width = width, height = height, res = 100)
#   try(print(draw(hm)))
#   dev.off()
# }

plot.mat.auto.clustering <- function(mat, k.cut=NULL, k.m=NULL, cluster_rows, cluster_columns,
                                     show_column_names = T, show_row_names = F,
                                     show_row_dend = F, show_column_dend = F,
                                     row.dist = "pearson", row.method = "complete", 
                                     col.dist = "pearson", col.method = "complete", plot.out, 
                                     width = 700, height = 700) {
  # warning("when using the block annotation of complexheatmap, a k means clustering will first be performed according to the km parameter")
  if (!is.null(k.m)) {
    anno.block <- 1:k.m
    split.param <- list(row_km = k.m)
  } else {
    if (!is.null(k.cut)) {
      anno.block <- unique(k.cut)
      if (length(anno.block) == 1)
        anno.block <- 1:anno.block
      split.param <- list(row_split = k.cut)
    } else {
      anno.block <- 1
      split.param <- NULL
    }
  }
  
  ha = HeatmapAnnotation(grouping = anno_block(gp = gpar(fill = anno.block)),
                         annotation_name_side = "top", which = "row", annotation_name_rot = 0, show_annotation_name = F)

  hm.params <- list(matrix = mat, show_column_names = show_column_names, show_row_names = show_row_names,
                    show_row_dend = show_row_dend, show_column_dend = show_column_dend, 
                    cluster_rows = cluster_rows, clustering_distance_rows = row.dist, 
                    clustering_method_rows = row.method,
                    cluster_columns = cluster_columns, clustering_distance_columns = col.dist, 
                    clustering_method_columns = col.method,
                    right_annotation = ha
                    # row_km = k
                    #row_split = k
  )
  hm.params <- c(hm.params, split.param)
  
  hm <- do.call(what = ComplexHeatmap::Heatmap, args = hm.params)
  system(paste0("mkdir -p ", dirname(plot.out)))
  png(filename = plot.out,
      width = width, height = height, res = 100)
  try(print(draw(hm)))
  dev.off()
  invisible(hm)
}

plot.mat.rank.row <- function(mat, symbol.mat = NULL, no.ranking = F, no.col.cluster = F,
                              title = "", hm.colors = NULL,
                              hm.values,
                              show_column_names = T, show_row_names = F, show_column_dend = T,
                              col.dist = "pearson", col.method = "complete", plot.out, 
                              row_name_fontSize = 5,
                              width = 2, height = 2, 
                              # You can actually cluster rows with this function, but by default we don't
                              cluster_rows = F, row.dist = "pearson", row.method = "complete",
                              ...) {
  fake.hm <- ComplexHeatmap::Heatmap(
    matrix = mat, cluster_columns = T, clustering_distance_columns = col.dist,
    clustering_method_columns = col.method)
  
  if (!no.col.cluster) {
    col.order <- suppressWarnings(ComplexHeatmap::column_order(fake.hm))
    mat <- mat[, col.order]
  }
  
  if (!no.ranking) {
    # newly added
    
    mat <- mat[rev(order(rowMax(mat))),]
    # 
    which.row.max <- apply(mat, 1, which.max)
    row.order <- order(which.row.max)
    mat <- mat[row.order, ]
  }
  
  suppressMessages(extrafont::loadfonts())
  ht_opt$HEATMAP_LEGEND_PADDING <- unit(0.1, "in")
  ht_opt$DENDROGRAM_PADDING <- unit(0, "in")
  
  col_fun <- NULL
  if (!is.null(hm.colors)) {
    col_fun = circlize::colorRamp2(hm.values, hm.colors)
  }
  
  hm.params <- list(matrix = mat, col = col_fun,
                    show_column_names = show_column_names, show_row_names = show_row_names,
                    show_column_dend = show_column_dend, 
                    cluster_columns = !no.col.cluster, clustering_distance_columns = col.dist, 
                    clustering_method_columns = col.method,
                    # You can cluster roles... but by default we don't
                    cluster_rows = cluster_rows, clustering_distance_rows = row.dist, 
                    clustering_method_rows = row.method,
                    #
                    row_names_gp = gpar(fontsize = row_name_fontSize, fontfamily = "Arial", lineheight = 0.6),
                    column_names_gp = gpar(fontsize = 6, fontfamily = "Arial"),
                    column_dend_height = unit(0.1, "in"),
                    column_dend_gp = gpar(lwd = 0.5),
                    row_dend_gp = gpar(lwd = 0.5),
                    row_dend_width = unit(0.1, "in"),
                    show_heatmap_legend = T,
                    heatmap_legend_param = list(
                      title = title, labels_gp = gpar(fontsize = 5), title_gp = gpar(fontsize = 6),
                      legend_height = unit(0.3, "in"), grid_width = unit(0.05, "in"), gap = unit(2, "in")
                    #   legend_direction = "horizontal"
                    ),
                    ...
                    
  )
  
  if (!is.null(symbol.mat)) {
    if (!identical(sort(rownames(symbol.mat)), sort(rownames(mat)) )) {
      stop("!identical(sort(rownames(symbol.mat)), sort(rownames(mat)) )")
    }
    if (!identical(sort(colnames(symbol.mat)), sort(colnames(mat)) )) {
      stop("!identical(sort(colnames(symbol.mat)), sort(colnames(mat)) )")
    }
    symbol.mat <- symbol.mat[rownames(mat), colnames(mat)]
    hm.params <- c(hm.params, list(
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(symbol.mat[i, j], x, y,
                  gp = gpar(fontsize = 5))
      }
    ))
  }

  hm <- do.call(what = ComplexHeatmap::Heatmap, args = hm.params)

  if (!grepl(".pdf$", plot.out)) {
    stop("plot.out must be pdf")
  }
  system(paste0("mkdir -p ", dirname(plot.out)))
  
  cairo_pdf(filename = plot.out, width = width, height = height, family = "Arial")
  try({print(draw(hm))})
  dev.off()
  
  ht_opt(RESET = T)
  
  invisible(hm)
}


peak.watch.kmo <- function(kmo, out.file.all, out.file.watch, n, seed = 42) {
  mp.df <- data.frame(cluster = kmo$cluster, peak = names(kmo$cluster)) %>% 
    utilsFanc::loci.2.df(loci.col.name = "peak", remove.loci.col = T)
  trash <- utilsFanc::write.zip.fanc(mp.df, out.file.all)
  
  watch.df <- mp.df %>% split(., f = .$cluster) %>% 
    lapply(function(df) {
      set.seed(seed = 42)
      if (n > nrow(df))
        n <- nrow(df)
      df <- df[sample(1:nrow(df), n, replace = F), ]
    }) %>% Reduce(rbind, .)
  trash <- utilsFanc::write.zip.fanc(watch.df, out.file.watch)
  return()
}


# get.hm.mat.fanc <- function(seMarker, ) {
#   # basically copying and slightly modifying the code from ArchR::plotMarkerHeatmap to get intermediate stage data
# }

archr.marker.filter <- function(hm.ns.mat = NULL, mp.mat = NULL, in.markers, size.threshold = 3, max.cluster,
                                save.Rds = NULL) {
  if (is.null(hm.ns.mat)) {
    hm.ns.mat <- plotMarkerHeatmap(mp.mat, scaleTo = 10^6,
                                   transpose = F, scaleRows = F, log2Norm = T, 
                                   limits = c(0, 100000), binaryClusterRows = F,
                                   clusterCols = F, 
                                   labelRows = F,# , nLabel = 0, nPrint = 1
                                   labelMarkers = NULL, returnMatrix = T)
  }
  
  mat <- hm.ns.mat[in.markers,]
  
  markers.filtered <- in.markers[mat[in.markers, ] %>% apply(1, max) == 
                                   mat[in.markers, max.cluster] & 
                                   mat[in.markers, max.cluster] >= size.threshold]
  if (!is.null(save.Rds))
    saveRDS(markers.filtered, save.Rds)
  return(markers.filtered)
}


