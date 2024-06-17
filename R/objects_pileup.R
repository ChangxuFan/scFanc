archr.pileup <- function(ao, peak.cat.df, kmo = NULL, peak.mat.full, normalize = T, 
                         return.df = T, return.ao = F) {
  # format of peak.cat.df: peak and cat (category of the peak) are the columns
  if (is.character(peak.mat.full))
    peak.mat.full <- readRDS(peak.mat.full)
  if (!is.null(kmo))
    peak.cat.df <- data.frame(cat = kmo$cluster %>% paste0("peakClus", .),
                              peak = names(kmo$cluster))
  ## first filter:
  peaks.in.mat <- rowRanges(peak.mat.full) %>% utilsFanc::gr.get.loci()
  peak.mat <- assay(peak.mat.full) %>% `rownames<-`(peaks.in.mat)
  depth.vec <- Matrix.utils::aggregate.Matrix(x = peak.mat, fun = "sum", 
                                              groupings = rep(1, nrow(peak.mat)))[1, ]
  
  peaks.shared <- intersect(peaks.in.mat, peak.cat.df$peak)
  if (length(peaks.shared) <  length(peak.cat.df$peak)) {
    pct <- length(peaks.shared) / length(peak.cat.df$peak)
    warning(paste0(100 * pct," percent of the peaks in peak.cat.df is shared with peak.mat.full"))
  }
    
  
  rownames(peak.cat.df) <- peak.cat.df$peak
  
  peak.cat.df <- peak.cat.df[peaks.shared, ]
  
  
  peak.mat <- peak.mat[peaks.shared, ]
  cellNames <- colnames(peak.mat)
  
  pileup.mat <- Matrix.utils::aggregate.Matrix(x = peak.mat, groupings = peak.cat.df$cat, fun = "sum")
  
  if (normalize == T) {
    pileup.mat <- pileup.mat %*% diag(1/depth.vec)
    colnames(pileup.mat) <- cellNames
  }
  if (return.ao) {
    for (i in rownames(pileup.mat)) {
      cellNames <- cellNames[cellNames %in% ao$cellNames]
      ao <- addCellColData(ao, data = pileup.mat[i,cellNames], cells = cellNames, name = i, force = T)
    }
    return(ao)
  }
  if (return.df == T) {
    # browser()
    pileup.df <- pileup.mat  %>% Matrix::t() %>% as.data.frame()
    rownames(pileup.df) <- cellNames
    return(pileup.df)
  } 
  return(pileup.mat)
}

plot.pileup.fanc <- function(ao, pileup.df, cat.to.plot = NULL, plot.out, 
                             max.quantile = NULL, pt.size = 0.05, ...) {
  if (is.null(max.quantile))
    max.quantile <- 1
  if (is.null(cat.to.plot))
    cat.to.plot <- colnames(pileup.df)
  umap <- ao@embeddings$UMAP$df %>% `colnames<-`(c("UMAP1", "UMAP2")) %>% 
    mutate(., cell = rownames(.))
  if (is.null(pileup.df$cell))
    pileup.df$cell <- rownames(pileup.df)
  plot.df <- left_join(umap, pileup.df)
  p.list <- lapply(cat.to.plot, function(cat) {
    value <- plot.df[, cat]
    # browser()
    value[value >= quantile(value, max.quantile)] <- quantile(value, max.quantile)
    plot.df[,cat] <- value
    p <- ggplot(plot.df, aes_string(x = "UMAP1", y = "UMAP2", color = cat)) +
      geom_point(size = pt.size) +
      scale_color_gradientn(colors = rev(brewer.pal(5,"Spectral")),na.value = "grey70") +
      ggtitle(cat)
    return(p)
  })
  trash <- wrap.plots.fanc(plot.list = p.list, plot.out = plot.out, ...)
  return()
}

