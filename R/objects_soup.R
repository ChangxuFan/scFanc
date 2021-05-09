soup.pipe <- function(so, tod, out.Rds = NULL, rho = NULL) {
  if (is.character(so))
    so <- readRDS(so)
  if (is.character(tod))
    tod <- readRDS(tod)
  
  toc <- so@assays$RNA@counts
  sc <- SoupChannel(tod, toc)
  sc = setClusters(sc, so$seurat_clusters)
  
  sc = setDR(sc, so@reductions$umap@cell.embeddings %>% `colnames<-`(c("RD1", "RD2")))
  # browser()
  if (!is.null(rho)) {
    sc = setContaminationFraction(sc, rho)
  } else {
    sc <- autoEstCont(sc)
  }
  sc$soupProfile = sc$soupProfile[rownames(sc$toc),]
  out = adjustCounts(sc, roundToInt = T)
  if (!is.null(out.Rds)) {
    system(paste0("mkdir -p ", dirname(out.Rds)))
    saveRDS(out, out.Rds)
  }
  return(out)
  
}

soup.pipe.m <- function(sol, tod.l, names, work.dir, rho = NULL) {
  system(paste0("mkdir -p ", work.dir))
  sol.nosoup <- lapply(seq_along(sol), function(i) {
    so <- sol[[i]]
    if (is.character(so))
      so <- readRDS(so)
    tod <- tod.l[[i]]
    if (is.character(tod))
      tod <- readRDS(tod)
    name <- names[i]
    if (i == 2) 
      browser()
    mat <- soup.pipe(so = so, tod = tod, out.Rds = paste0(work.dir, "/",name, "_nosoup_mat.Rds"),
                     rho = rho)
    so.nosoup <- CreateSeuratObject(counts = mat, project.name = "nosoup", meta.data = so@meta.data)
    return(so.nosoup)
  })
  saveRDS(sol.nosoup, paste0(work.dir, "/sol_nosoup.Rds"))
  return(sol.nosoup)
}

soup.assess <- function(so, desoup.mat.vec, rho.vec, cluster.ident, clusters = NULL,
                        do.plot = T, out.dir, threads = 10) {
  if (is.character(so))
    so <- readRDS(so)
  markers.general <- readRDS("~/hmtp/kc/panel.list.Rds")[["general"]] %>% unlist()
  purity.list <- mclapply(seq_along(rho.vec), function(i) {
    rho <- rho.vec[i]
    desoup.mat <- desoup.mat.vec[i]
    if (is.character(desoup.mat))
      desoup.mat <- readRDS(desoup.mat)
    if (!identical(colnames(so@assays$RNA@counts), colnames(desoup.mat)) ||
        !identical(rownames(so@assays$RNA@counts), rownames(desoup.mat))) {
      stop("the colnames and rownames of desoup.mat must match that of so!")
    }
    root.name <- paste0("rho_", rho)
    soup.assay <- CreateAssayObject(counts = desoup.mat)
    so[[root.name]] <- soup.assay
    so <- NormalizeData(object = so, assay = root.name, normalization.method = "LogNormalize",
                        scale.factor = 10000)
    if (do.plot == T) {
      browser()
      trash <- plot.panel.list(panel.list = markers.general, obj = so, order = T, assay = root.name, 
                               plot.out = paste0(out.dir, "/", root.name, "_hmtp_general.png"))
    }

    purity <- purity.check.seurat(so = so, marker.genes = markers.general, assay = root.name,
                                  cluster.ident = cluster.ident, clusters = clusters, return.df = T)
    purity <- lapply(purity, function(x) {
      x <- utilsFanc::add.column.fanc(df1 = x, df2 = data.frame(rho = rep(rho, nrow(x))), pos = 2)
      return(x)
    })
    return(purity)
  }, mc.cores = threads, mc.cleanup = T)
  
  purity.cat <- lapply(names(purity.list[[1]]), function(item) {
    cat <- lapply(purity.list, function(x) return(x[[item]])) %>% Reduce(rbind, .) %>% 
      arrange(gene, rho)
    utilsFanc::write.zip.fanc(df = cat, out.file = paste0(out.dir, "/", item, ".tsv"), 
                              col.names = T, row.names = F, zip = F)
    return(cat)
  })
  names(purity.cat) <- names(purity.list[[1]])
  return(purity.cat)
}


soup.assess.2 <- function(so, add.soup.mats = T, desoup.mat.vec, rho.vec, # cluster.ident, clusters = NULL,
                          do.plot = T, markers = NULL, plot.out, save.rds = NULL, threads = 10) {
  if (is.character(so))
    so <- readRDS(so)
  if (is.null(markers))
    markers <- readRDS("~/hmtp/kc/panel.list.Rds")[["general"]] %>% unlist()
  if (add.soup.mats == T) {
    for (i in 1:length(desoup.mat.vec)) {
      rho <- rho.vec[i]
      desoup.mat <- desoup.mat.vec[i]
      if (is.character(desoup.mat))
        desoup.mat <- readRDS(desoup.mat)
      if (!identical(colnames(so@assays$RNA@counts), colnames(desoup.mat)) ||
          !identical(rownames(so@assays$RNA@counts), rownames(desoup.mat))) {
        stop("the colnames and rownames of desoup.mat must match that of so!")
      }
      root.name <- paste0("rho_", rho)
      soup.assay <- CreateAssayObject(counts = desoup.mat)
      so[[root.name]] <- soup.assay
      so <- NormalizeData(object = so, assay = root.name, normalization.method = "LogNormalize",
                          scale.factor = 10000)
    }
  }
  # browser()
  # markers <- markers[1]
  if (do.plot == T) {
    trash <- plot.panel.list.m(panel.list = markers, obj = rep(list(so), length(rho.vec) + 1), order = T,
                               assay = as.list(c("RNA", paste0("rho_", rho.vec))), 
                               plot.out = plot.out, threads.master = threads)
  }

  if (!is.null(save.rds))
    saveRDS(so, save.rds)
  return(so)
}