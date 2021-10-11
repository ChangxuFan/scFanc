soup.pipe <- function(so, tod, soup.range = c(0, 100), out.Rds = NULL, rho = NULL, cluster.ident = "seurat_clusters") {
  if (length(soup.range) == 1) {
    soup.range <- c(0, soup.range)
  }
  
  if (length(soup.range) != 2) {
    stop("fanc error: soup.range mis-specified")
  }
  if (is.character(so))
    so <- readRDS(so)
  so@meta.data <- so@meta.data %>% factor2character.fanc()
  if (is.character(tod))
    tod <- readRDS(tod)
  
  toc <- so@assays$RNA@counts
  sc <- SoupChannel(tod, toc, calcSoupProfile = F)
  sc <- estimateSoup(sc = sc, soupRange = soup.range, keepDroplets = F)
  
  sc = setClusters(sc, so@meta.data[, cluster.ident])
  sc = setDR(sc, so@reductions$umap@cell.embeddings %>% `colnames<-`(c("RD1", "RD2")))
  png(paste0(sub(".[Rr]ds", "", out.Rds), "_estimate.png"))
  try(sc <- autoEstCont(sc))
  dev.off()
  
  if (!is.null(rho)) {
    sc = setContaminationFraction(sc, rho)
  } 
  sc$soupProfile = sc$soupProfile[rownames(sc$toc),]
  out = adjustCounts(sc, roundToInt = T)
  if (!is.null(out.Rds)) {
    system(paste0("mkdir -p ", dirname(out.Rds)))
    saveRDS(out, out.Rds)
  }
  return(out)
  
}

soup.pipe.m <- function(sol, tod.l, soup.range = c(0, 100),
                        names, work.dir, rho = NULL, cluster.ident = "seurat_clusters") {
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
                     rho = rho, cluster.ident = cluster.ident, soup.range = soup.range)
    so.nosoup <- CreateSeuratObject(counts = mat, project.name = "nosoup", meta.data = so@meta.data)
    return(so.nosoup)
  })
  names(sol.nosoup) <- names
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
                          do.plot = T, markers = NULL, plot.out, save.rds = NULL, threads = 10,
                          transformation = NULL) {
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
      if (!is.null(transformation))
        desoup.mat <- transformation(desoup.mat)
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

add.desoup.mat <- function(so, desoup.mat.vec, assay.name, save.rds = NULL) {
  rho.name <- paste0(assay.name, "_rho")
  if (is.null(names(desoup.mat.vec))) {
    stop("desoup.mat.vec must be named with sample name")
  }
  
  rho.df <- lapply(desoup.mat.vec, function(x) {
    df <- data.frame(sample = sub("_nosoup.*$", "", basename(x)),
                     rho = sub("^.+rho_", "", dirname(x)) %>% as.numeric())
    return(df)
  }) %>% Reduce(rbind, .)
  colnames(rho.df) <- c("sample", rho.name)
  rownames(rho.df) <- rho.df$sample
  
  fake.soi <- utilsFanc::safelapply(names(desoup.mat.vec), function(sample) {
    mat <- readRDS(desoup.mat.vec[sample])
    if (!grepl("#", mat@Dimnames[[2]][1]))
      mat@Dimnames[[2]] <- paste0(sample, "#",  mat@Dimnames[[2]])
    fake.so <- CreateSeuratObject(counts = mat)
    return(fake.so)
  }, threads = length(desoup.mat.vec)) %>% 
    Reduce(Seurat::merge.Seurat, .)
  
  mat <- fake.soi@assays$RNA@counts
  cells.missing <- colnames(so) %>% .[!. %in% colnames(mat)]
  if (length(cells.missing) > 0) {
    stop(paste0(length(cells.missing)," cells in so is not found in desoup.mat.vec\nthe first 5:\n",
                paste0(cells.missing[1:5], collapse = "\n")))
  }
  mat <- mat[, colnames(so)]
  desoup.assay <- CreateAssayObject(counts = mat)
  so[[assay.name]] <- desoup.assay
  so@meta.data[, rho.name] <- rho.df[so$sample, rho.name]
  so <- NormalizeData(object = so, assay = assay.name)
  if (!is.null(save.rds)) {
    system(paste0("mkdir -p ", dirname(save.rds)))
    saveRDS(so, save.rds)
  }
  return(so)
}

soup.profile.comp <- function(so, tod, soup.range.max.vec, plot.out, threads = NULL) {
  print("this function is not fully finished. The plotting part was not written. But the summary() should be good enough to show that there is no difference")
  if (is.null(threads))
    threads <- length(soup.range.max.vec)
  if (threads > 12)
    threads <- 12
  
  if (is.character(so)) 
    so <- readRDS(so)
  if (is.character(tod))
    tod <- readRDS(tod)
  so@meta.data <- so@meta.data %>% factor2character.fanc()
  toc <- so@assays$RNA@counts
  sc <- SoupChannel(tod, toc, calcSoupProfile = F)
  pl <- utilsFanc::safelapply(soup.range.max.vec, function(soup.range.max) {
    sc <- estimateSoup(sc = sc, soupRange = c(0,soup.range.max), keepDroplets = F)
    genes <- rownames(tod)
    tod <- tod[,colSums(tod) < soup.range.max]
    fanc <- Matrix::t(tod) %>% Matrix.utils::aggregate.Matrix(., groupings = rep("fanc", nrow(.))) %>% 
      Matrix::t()
    df <- sc$soupProfile
    colnames(df) <- c("SoupX_pct", "SoupX_counts")
    df$fanc_counts <- fanc[rownames(df), 1]
    summary(df$fanc_counts)
    summary(df$SoupX_counts)
    summary(df$fanc_counts-df$SoupX_counts)
    
  }, threads = threads)
  
}