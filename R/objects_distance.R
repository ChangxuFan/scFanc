cluster.distance <- function(obj, cluster.ident, clusters = NULL, 
                                group.ident = NULL, groups = NULL,
                                reduction, dims = 1:15,
                                dist.method, n.cells.each = 100,
                                plot.out = NULL) {
  
  if ("ArchRProject" %in% class(obj)) {
    is.ao <- T
  } else if ("Seurat" %in% class(obj)) {
    is.ao <- F
  } else {
    stop("obj must be ao or so")
  }
  cell.list <- get.cell.list(obj = obj, is.ao = is.ao, group.by = cluster.ident,
                             groups = clusters, split.by = group.ident, splits = groups,
                             na.rm = T, return.named.vec = T, n.cells.each = n.cells.each)
  
  if (is.ao)
    pos.mat <- obj@reducedDims[[reduction]]$matSVD[names(cell.list), dims]
  else
    pos.mat <- obj@reductions[[reduction]]@cell.embeddings[names(cell.list), dims]
  dist.mat <- dist(x = pos.mat, method = dist.method, diag = T, upper = T) %>% 
    as.matrix()
  if (sum(is.na(dist.mat)) > 0) {
    stop("sum(is.na(dist.mat)) > 0")
  }
  dist.mean <- aggregate.fanc(mat = dist.mat, groupings = cell.list, margin = 1, 
                             binarize = F, na.rm = T, take.mean = T, sort = F)
  
  dist.mean <- aggregate.fanc(mat = dist.mean, groupings = cell.list, margin = 2, 
                             binarize = F, na.rm = T, take.mean = T, sort = F)
  dist.mean <- as.matrix(dist.mean)
  if (!is.null(plot.out)) {
    p <- Heatmap(dist.mean)
    save.base.plot(p = p, width = 500, height = 450, file = plot.out)
  }
  return(dist.mean)
}