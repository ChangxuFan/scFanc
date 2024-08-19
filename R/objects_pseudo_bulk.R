rowname.shffl <- function(df, direction, which.col, keep, ...) {
  if (direction == "in") {
    df2 <- data.frame(V1 = rownames(df))
    colnames(df2) <- which.col
    df <- utilsFanc::add.column.fanc(df1 = df, df2 = df2, ...)
    if (keep == F)
      rownames(df) <- NULL
  }
  if (direction == "out") {
    rownames(df) <- df[,which.col]
    if (keep == F)
      df[, which.col] <- NULL
  }

  return(df)
}


# pseudo.bulk.core <- function(so, assay, slot, cell.names, collapse.method = mean) {
#   exp.df <- GetAssayData(object = so, slot = slot, assay = assay)
#
# }


so.2.bulk <- function(so, root.name = NULL, assay, slot, cells=NULL, genes = NULL, ident, sub.idents,
                      group.by= NULL, groups = NULL, take.mean = F, binarize = F,
                      coldata.columns = NULL, debug = F) {
  warning("so.2.bulk function is intended to give you bulk for one cluster. if you give more " %>%
            paste0("clusters in sub.idents, it will just pile these clusters together"))

  meta.df <- so@meta.data
  if (is.null(group.by)) {
    group.by <- "miao"
    meta.df$miao <- "wang"
  }

  if (length(group.by) > 1) {
    meta.df$mult <- Reduce(function(x,y) paste(x, y, sep = ".."),
                           meta.df[, group.by] %>% as.list())
    if (!is.null(groups)) {
      if (!is.list(groups) || length(groups) != length(group.by)) {
        stop("when length(group.by) > 1 and !is.null(groups), groups must be a list, each element corresponding to group.by")
      }
      groups <- lapply(seq_along(groups), function(i) {
        x <- groups[[i]]
        if (is.null(x)) {
          x <- meta.df[, group.by[i]] %>% unique() %>% .[!is.na(.)]
        }
        return(x)
      })
      groups <- Reduce(function(x, y) outer(x, y, FUN = function(x, y) paste(x,y, sep = "..")),
                       groups) %>% as.character()
    }
    group.by <- "mult"
  }
  meta.df <- meta.df %>% utilsFanc::change.name.fanc(cols.from = group.by, cols.to = "tgroup")
  if (is.null(coldata.columns)) {
    meta.df[, group.by] <- meta.df$tgroup
    coldata.columns <- group.by
  }

  meta.df$cell.id <- rownames(meta.df)
  if (is.null(cells)) {
    utilsFanc::check.intersect(sub.idents, "sub.idents", meta.df[, ident], "meta.df[, ident]")
    cells <- meta.df %>% .[.[,ident] %in% sub.idents,] %>% {if (is.null(groups)) .
      else filter(., tgroup %in% groups)} %>%
      pull(cell.id)
  }
  meta.df <- meta.df %>% .[.$cell.id %in% cells,]

  exp.df <- GetAssayData(object = so, slot = slot, assay = assay)
  if (binarize == T) {
    exp.df[exp.df > 0] <- 1
    exp.df[exp.df <= 0] <- 0
  }

  if (!is.null(genes))
    exp.df <- exp.df[genes, ]
  exp.df <- exp.df[, colnames(exp.df) %in% cells] %>%
    as.matrix() %>% t() %>% as.data.frame()
  exp.df$cell.id <- exp.df %>% rownames()

  exp.df <- left_join(meta.df[, c("cell.id", "tgroup")], exp.df)
  exp.df$cell.id <- NULL
  #### the data frame summarise based method seems too slow.
  # bulk.df <- exp.df %>% group_by(tgroup) %>% summarise_all(sum) %>% as.data.frame()
  # bulk.mat <- bulk.df
  # rownames(bulk.mat) <- bulk.mat$tgroup
  # bulk.mat$tgroup <- NULL
  # bulk.mat <- as.matrix(bulk.mat)
  ########

  bulk.mat <- exp.df %>% Matrix.utils::aggregate.Matrix(.,groupings = .$tgroup, fun = "sum") %>% t()
  bulk.mat <- bulk.mat[rownames(bulk.mat) != "tgroup",,drop=F]
  if (take.mean == T) {
    cell.number.df <- meta.df %>% group_by(tgroup) %>% summarise(n = n())
    for (i in 1:nrow(cell.number.df)) {
      bulk.mat[, cell.number.df$tgroup[i]] <- bulk.mat[, cell.number.df$tgroup[i]] / cell.number.df$n[i]
    }
  }

  coldata <- NULL
  if (!is.null(coldata.columns)) {
    coldata <- meta.df %>% mutate(cell.id = NULL)  %>% `[`(, c(coldata.columns, "tgroup"), drop = F) %>% group_by(tgroup) %>%
      summarise_all(function(x) paste0(unique(x), collapse = ":")) %>% as.data.frame()
    rownames(coldata) <- coldata$tgroup
    coldata$tgroup <- NULL
  }

  return.list <- list(root.name = root.name, bulk.mat = bulk.mat %>% as.matrix(), coldata=coldata)
  if (debug == T) {
    return.list <- list(root.name = root.name, bulk.mat = bulk.mat %>% as.matrix(), coldata=coldata, mat = exp.df)
  }
  return(return.list)

}

so.2.bulk.m <- function(so, root.name, work.dir = NULL, assay, slot,
                        cell.list = NULL, genes.include = NULL, cells.include = NULL,
                        cluster.ident, clusters = NULL,
                        group.ident, groups = NULL,
                        take.mean = F, binarize = F,
                        coldata.df = NULL, coldata.columns = NULL) {
  if (is.null(cell.list)) {
    cell.list <- get.cell.list(obj = so, is.ao = F, cells.include = cells.include,
                               split.by = group.ident, splits = groups,
                               group.by = cluster.ident, groups = clusters)
  }
  cell.vec <- list.switch.name(x = cell.list)
  mat <- GetAssayData(so, slot = slot, assay = assay)
  if (!is.null(genes.include))
    mat <- mat %>% .[rownames(.) %in% genes.include,]
  bulk.mat <- aggregate.fanc(mat = mat, margin = 2, groupings = cell.vec, na.rm = T,
                             binarize = binarize, take.mean = take.mean, sort = T)
  bulk.mat <- bulk.mat %>% as.matrix()
  if (is.null(coldata.df)) {
    coldata.df <- so@meta.data[, unique(c(group.ident, cluster.ident, coldata.columns)),
                               drop = F] %>%
     unique() %>% na.omit()
    if (!is.null(clusters))
      coldata.df <- coldata.df %>% .[.[, cluster.ident] %in% clusters,]
    if (!is.null(groups))
      coldata.df <- coldata.df %>% .[.[, group.ident] %in% groups,]
    rownames(coldata.df) <- paste0(coldata.df[, cluster.ident], "..", coldata.df[, group.ident])
  }
  coldata.df <- coldata.df[colnames(bulk.mat), ]
  s2b <- list(root.name = root.name, bulk.mat = bulk.mat, coldata = coldata.df, cell.list = cell.list)
  if (!is.null(work.dir)) {
    dir.create(work.dir, showWarnings = F, recursive = T)
    saveRDS(s2b, paste0(work.dir, "/", root.name, "_s2b.Rds"))
  }
  return(s2b)
}
soup.2.bulk <- function(droplet.counts.vec, soup.range.max.vec,
                        sample.info.tsv = NULL, coldata = NULL,
                        threads = NULL, save.rds = NULL) {
  if (is.null(threads))
    threads <- length(droplet.counts.vec)
  samples <- names(droplet.counts.vec)
  if (is.null(samples))
    stop("is.null(samples)")

  droplet.counts.vec <- utilsFanc::safelapply(samples, function(sample) {
    droplet.counts <- readRDS(droplet.counts.vec[sample])
    return(droplet.counts)
  }, threads = threads)

  names(droplet.counts.vec) <- samples

  pbl <- lapply(soup.range.max.vec, function(soup.range.max) {
    soup.list <- utilsFanc::safelapply(samples, function(sample) {
      tod <- droplet.counts.vec[[sample]]
      tod <- tod[, colSums(tod) < soup.range.max]
      soup <- Matrix::t(tod) %>% Matrix.utils::aggregate.Matrix(., groupings = rep(sample, nrow(.))) %>%
        Matrix::t()
      return(soup)
    }, threads = threads)

    genes <- Reduce(union, lapply(soup.list, function(x) return(rownames(x))))

    soup.list <- lapply(soup.list, function(x) {
      add <- genes %>% .[! .%in% rownames(x)]
      if (length(add) > 0) {
        add.mat <- matrix(rep(0, length(add)), ncol = 1, dimnames = list(add, colnames(x)))
        x <- rbind(x, add.mat)
      }
      return(x)
    })
    soup.mat <- Reduce(cbind, soup.list)

    if (is.null(coldata)) {
      if (! is.null(sample.info.tsv)) {
        sample.info <- read.table(sample.info.tsv, header = T, sep = "\t") %>% filter(sample %in% samples)
        samples.not.in <- samples %>% .[!.%in% sample.info$sample]
        if (length(samples.not.in) > 0) {
          stop(paste0("some samples were not found in sample.info.tsv, top 5 are: \n",
                      paste0(samples.not.in, collapse = "\n")))
        }
        rownames(sample.info) <- sample.info$sample
        coldata <- sample.info
      }
    }

    s2b <- list(bulk.mat = soup.mat, coldata = coldata)
    return(s2b)
  })
  names(pbl) <- paste0("soup_max_", soup.range.max.vec)
  if (!is.null(save.rds)) {
    system(paste0("mkdir -p ", dirname(save.rds)))
    saveRDS(pbl, save.rds)
  }

  return(pbl)
}

# so.2.bulk.fast <- function(so, assay, slot, cells= NULL, genes, cluster.ident, clusters,
#                            group.by=NULL, groups=NULL) {
#   # this function will not separate each clusters in "clusters". They will be merged
#   if (is.null(cells)) {
#     cells <-
#   }
#
#   mat <- GetAssayData(so, assay = assay, slot = slot)
#
# }

ao.2.bulk <- function(ao=NULL, root.name, cell.list = NULL, gr = NULL, gr.sample.cols =NULL, peakmat = NULL,
                      coldata.columns = NULL, coldata.df = NULL,
                      cell.prep.dir,
                      pseudo, seed = NULL, nprep = 2,
                      cluster.ident, cluster, group.ident, groups=NULL,
                      filter.by = NULL, filter.limits, bQuantileFilter = c(F, F),
                      match.bg.df = NULL) {
  # if (pseudo != T)
  #   stop("getting coldata from non-pseudo based setting is not yet developed")

  if (is.null(gr)) {
    if (is.null(cell.list)) {

      if (length(cluster) > 1)
        stop("ao.2.bulk currently only handle 1 cluster")
      cell.list <- archr.get.cells.grid(ao = ao, cluster.ident = cluster.ident, clusters = cluster,
                                        group.ident = group.ident, groups = groups, melt = T,
                                        seed = seed, pseudo.gen = pseudo, nprep = nprep,
                                        match.bg.df = match.bg.df, match.bg.outdir = cell.prep.dir,
                                        filter.by = filter.by, filter.limits = filter.limits,
                                        bQuantileFilter = bQuantileFilter, filter.plot.dir = cell.prep.dir, root.name = root.name)
      names(cell.list) <- names(cell.list) %>% sub(paste0(cluster.ident, "_", cluster, ".."), "", .)
      # if (is.null(coldata.columns) && is.null(coldata.df) && pseudo == F) {
      #   getCellColData()
      # }
    }
    dir.create(cell.prep.dir, recursive = T, showWarnings = F)
    saveRDS(cell.list, paste0(cell.prep.dir, "/", root.name, "_cell.list.Rds"))
    gr.sample.cols <- names(cell.list)
    if (is.null(peakmat))
      stop("peakmat must be specified when gr is not")
    gr <- vectorization.core(cell.list = cell.list, ao = ao, mat = peakmat, mat.name = "PeakMatrix",
                             deseq2.norm = F)
  }
  if (is.null(coldata.df)) {
    if (pseudo == T) {
      coldata.df <- data.frame(sample = sub("..prep.*$", "", gr.sample.cols))
      rownames(coldata.df) <- gr.sample.cols
    } else {
      coldata.df <- getCellColData(ao, select = unique(c(group.ident, coldata.columns)), drop = F) %>%
        as.data.frame() %>% na.omit() %>% unique()
      if (length(group.ident) > 1) {
        rownames(coldata.df) <- Reduce(function(x,y) paste(x, y, sep = ".."),
                                       coldata.df[, group.ident] %>% as.list())
      } else {
        rownames(coldata.df) <- coldata.df[, group.ident]
      }

    }
  }

  # from gr to DESeq2 format:
  df <- gr %>% `names<-`(NULL) %>% as.data.frame()
  if (any(grepl("^\\d", gr.sample.cols))) {
    cat(paste0("Some sample names start from a number.\n",
                 "Trying to make this work by removing the 'X' from the beginning the the column names of gr\n"))
    colnames(df) <- sub("^X(\\d)", "\\1", colnames(df))
  }
  
  utilsFanc::check.intersect(gr.sample.cols, "gr.sample.cols", colnames(df), "colnames(df)")
  
  df <- df %>% .[, c("seqnames", "start", "end", gr.sample.cols)] %>%
    mutate(peak = paste0(seqnames, ":", start, "-", end))
  rownames(df) <- df$peak
  mat <- df[, gr.sample.cols] %>% as.matrix()
  coldata.df <- coldata.df[gr.sample.cols, , drop = F]
  return(list(root.name = root.name, bulk.mat = mat, coldata = coldata.df, gr = gr, cell.list = cell.list))
}


ao.2.bulk.m <- function(ao=NULL, root.name, cell.list = NULL, gr = NULL, gr.sample.cols =NULL, peakmat = NULL,
                      coldata.columns = NULL, coldata.df = NULL,
                      cell.prep.dir,
                      pseudo, seed = NULL, nprep = 2,
                      cluster.ident, clusters, group.ident, groups=NULL,
                      filter.by = NULL, filter.limits, bQuantileFilter = c(F, F),
                      match.bg.df = NULL) {
  # if (pseudo != T)
  #   stop("getting coldata from non-pseudo based setting is not yet developed")

  if (any(grepl("SummarizedExperiment", class(peakmat))))
    peakmat <- assay(peakmat)
  if (is.null(gr)) {
    if (is.null(cell.list)) {

      cell.list <- archr.get.cells.grid(ao = ao, cluster.ident = cluster.ident, clusters = clusters,
                                        group.ident = group.ident, groups = groups, melt = T,
                                        seed = seed, pseudo.gen = pseudo, nprep = nprep,
                                        match.bg.df = match.bg.df, match.bg.outdir = cell.prep.dir,
                                        filter.by = filter.by, filter.limits = filter.limits,
                                        bQuantileFilter = bQuantileFilter, filter.plot.dir = cell.prep.dir, root.name = root.name)
      names(cell.list) <- names(cell.list) %>% sub(paste0(cluster.ident, "_"), "", .)
      # if (is.null(coldata.columns) && is.null(coldata.df) && pseudo == F) {
      #   getCellColData()
      # }
    }
    dir.create(cell.prep.dir, recursive = T, showWarnings = F)
    saveRDS(cell.list, paste0(cell.prep.dir, "/", root.name, "_cell.list.Rds"))
    gr.sample.cols <- names(cell.list)
    if (is.null(peakmat))
      stop("peakmat must be specified when gr is not")
    gr <- vectorization.core(cell.list = cell.list, ao = ao, mat = peakmat, mat.name = "PeakMatrix",
                             deseq2.norm = F)
  }
  if (is.null(coldata.df)) {
    if (pseudo == T) {
      coldata.df <- data.frame(sample = sub("..prep.*$", "", gr.sample.cols))
      rownames(coldata.df) <- gr.sample.cols
    } else {
      coldata.df <- getCellColData(ao, select = c(group.ident, cluster.ident, coldata.columns) %>% unique(), drop = F) %>%
        as.data.frame() %>% unique() %>% na.omit()
      if (!is.null(clusters))
        coldata.df <- coldata.df %>% .[.[, cluster.ident] %in% clusters,]
      if (!is.null(groups))
        coldata.df <- coldata.df %>% .[.[, group.ident] %in% groups,]
      rownames(coldata.df) <- paste0(coldata.df[, cluster.ident], "..", coldata.df[, group.ident])
      # %>% `rownames<-`(., .[, group.ident])
    }
  }

  # from gr to DESeq2 format:
  df <- gr %>% `names<-`(NULL) %>% as.data.frame() %>% .[, c("seqnames", "start", "end", gr.sample.cols)] %>%
    mutate(peak = paste0(seqnames, ":", start, "-", end))
  rownames(df) <- df$peak
  mat <- df[, gr.sample.cols] %>% as.matrix()
  coldata.df <- coldata.df[gr.sample.cols, ]
  return(list(root.name = root.name, bulk.mat = mat, coldata = coldata.df, gr = gr, cell.list = cell.list))
}

t.f.a2b.da <- function(a2b.list, use.peakwatch=F, p.adj.cutoff, p.cutoff,
                       comp, lfc.cutoff,
                       plot.out.root) {
  sample.y <- sub(":.+$", "", comp)
  sample.x <- sub("^.+:", "", comp)
  system(paste0("mkdir -p ", dirname(plot.out.root)))
  p.list <- lapply(seq_along(a2b.list), function(i) {
    a2b <- a2b.list[[i]]
    name <- names(a2b.list)[i]
    if (use.peakwatch == T) {
      if (is.null(a2b$peakwatch )) {
        stop("peakwatch object is not existent")
      }
      hl.list <- a2b$peakwatch %>% split(f = a2b$peakwatch$stratify) %>%
        lapply(function(x) return(x$gene))
    } else {
      hl <- a2b$res.exp %>% filter(pvalue < p.cutoff, padj < p.adj.cutoff, abs(log2FoldChange) > lfc.cutoff) %>%
        pull(gene)
      hl.list <- list(hl)
    }
    p.list <- lapply(seq_along(hl.list), function(i) {
      hl <- hl.list[[i]]
      stratify <- names(hl.list)[i]
      p <- xy.plot(df = a2b$bulkNorm, x = paste0("bulkNorm_", sample.x), y = paste0("bulkNorm_", sample.y),
                   highlight.var = "gene", highlight.values = hl) +
        ggtitle(paste0(name, ", ",stratify))
      return(p)
    })
    return(p.list)
  }) %>% Reduce(c,.)
  if (use.peakwatch == T) {
    plot.out <- paste0(plot.out.root, "..peakwatch", ".png")
  } else {
    plot.out <- paste0(plot.out.root, "..", "p_", p.cutoff, "..padj_", p.adj.cutoff,
                       "..lfc_", lfc.cutoff, ".png")
  }
  trash <- wrap.plots.fanc(p.list, plot.out = plot.out)
  return()
}

ao.2.bulk.list <- function(ao, peakmat, pseudo, cell.prep.dir = NULL, sample.order = NULL,
                           filter.nz = F, filter.size = NULL, filter.samples = NULL,
                           filter.fun = "max", sequential.filter, protected.genes = NULL,
                           independentFiltering = F,
                           quantile.norm = F, 
                           deseq2.norm.method = "ratio", deseq2.locfunc = NULL,
                           coldata.columns = NULL, coldata.df = NULL,
                           pca.ntop = 10000, pca.groupings = NULL,
                           cluster.ident, clusters = NULL, group.ident, groups = NULL,
                           match.bg.df = NULL,
                           filter.by = NULL, filter.limits, bQuantileFilter = c(F, F),
                           design.formula, contrast,
                           do.plot = F,
                           work.dir,threads = 6, plot.dir = NULL,
                           single.sample = F) {
  if (is.null(plot.dir))
    plot.dir <- paste0(work.dir, "/plots")
  system(paste0("mkdir -p ", work.dir, " ", plot.dir))
  if (is.null(cell.prep.dir)) {
    cell.prep.dir <- paste0(work.dir, "/cell_prep/")
  }

  if (is.null(pca.groupings))
    pca.groupings <- coldata.columns

  if (is.null(cluster.ident)) {
    cluster.ident <- "pd_cluster"
    ao <- addCellColData(ArchRProj = ao, data = rep("0", length(ao$cellNames)),
                         cells = ao$cellNames, name = cluster.ident, force = T)
  }

  if (is.null(clusters)) {
    clusters <- getCellColData(ao, select = cluster.ident)[, cluster.ident] %>% unique() %>%
      gtools::mixedsort()
  }
  a2b.list <- mclapply(clusters, function(cluster){
    a2b <- ao.2.bulk(ao = ao, root.name = paste0(cluster.ident, "_", cluster),
                     peakmat = peakmat, pseudo = pseudo, cluster.ident = cluster.ident,
                     cluster = cluster, group.ident = group.ident, groups = groups,
                     coldata.columns = coldata.columns, coldata.df = coldata.df,
                     match.bg.df = match.bg.df,
                     filter.by = filter.by, filter.limits = filter.limits,
                     bQuantileFilter = bQuantileFilter, cell.prep.dir = cell.prep.dir)
    a2b <- s2b.deseq(s2b.obj = a2b, quantile.norm = quantile.norm, 
                     norm.method = deseq2.norm.method, locfunc = deseq2.locfunc,
                     filter.nz = filter.nz, filter.size = filter.size, filter.samples = filter.samples,
                     filter.fun = filter.fun, sequential.filter = sequential.filter,
                     protected.genes = protected.genes,
                     independentFiltering = independentFiltering,
                     pca.ntop = pca.ntop, pca.groupings = pca.groupings,
                     design = design.formula, contrast = contrast,
                     sample.order = sample.order, try.hm = do.plot,
                     force.hm = F, force = F, p.adj.cutoff = 0.1, p.cutoff = 0.05,
                     title = paste0("cluster_", cluster),
                     plot.dir = plot.dir,
                     outlist = list(plot.dir = plot.dir, root.name.external = "cluster_" %>% paste0(cluster, "_")),
                     single.sample = single.sample)
    return(a2b)
  }, mc.cores = threads, mc.cleanup = T)
  try({
    a2b.list <- a2bl.add.shrinkage(a2bl = a2b.list, contrast = contrast, sample.order = sample.order,
                                   threads = threads)
  })
  names(a2b.list) <- paste0(cluster.ident, "_",clusters)
  saveRDS(a2b.list, paste0(work.dir, "/ao.bulk.list.Rds"))
  return(a2b.list)
}

ao.2.bulk.soupmode <- function(ao, peakmat, ref.da, cluster.ident,
                               quantile.norm = F, contrast,
                               threads = 1,
                               work.dir) {
  if (is.null(ao@cellColData[[cluster.ident]])) {
    stop("is.null(ao@cellColData[[cluster.ident]])")
  }
  clusters <- names(ref.da) %>% sub(paste0(cluster.ident, "_", "", .))
  not.found <- clusters %>% .[!. %in% getCellColData(ao, cluster.ident, drop = T)]
  if (length(not.found) > 0) {
    stop(paste0("clusters not found in ao: ",
                paste0(not.found, collapse = ", ")))
  }
  da <- utilsFanc::safelapply(ref.da, function(ref.s2b) {
    cluster <- ref.s2b$root.name %>% sub(paste0(cluster.ident, "_", "", .))
    peaks <- ref.s2b$bulkNorm$gene
    work.dir <- paste0(work.dir, "/", ref.s2b$root.name)
    da <- ao.2.bulk.list(ao = ao, peakmat = peakmat[peaks, ], pseudo = F,
                         filter.nz = F, filter.size = NULL, filter.samples = NULL,
                         quantile.norm = quantile.norm,
                         coldata.columns = colnames(ref.s2b$coldata),
                         cluster.ident = cluster.ident, clusters = cluster,
                         group.ident = "Sample", groups = rownames(ref.s2b$coldata),
                         design.formula = ref.s2b$dds@design,
                         contrast = contrast, do.plot = F,
                         work.dir = work.dir, threads = 1)
    return(da[[1]])
  }, threads = threads)
  names(da) <- paste0("soup_", names(ref.da))
  return(da)
}
s2b.deseq <- function(s2b.obj, quantile.norm = F, use.genewise.disperisons = F,
                      norm.method = "ratio", locfunc = NULL, norm.cr.genes = NULL,
                      design, contrast=NULL,
                      use.ori = T, filter.nz = F, filter.size = NULL,
                      filter.fun = "max", filter.samples = NULL, sequential.filter = F,
                      protected.genes = NULL, independentFiltering = F,
                      pca.ntop = NULL, pca.groupings = NULL,
                      try.hm = T, force.hm = F,
                      pca.blind = T, pca.add.hm = F, pca.add.3d = F,
                      force = F, sample.order = NULL,
                      plot.dir = NULL, outlist = NULL, single.sample = F,
                      p.adj.cutoff, p.cutoff, title = character(0), pivots = NULL) {

  # I highly regret trying to use the outlist style. In the end I really couldn't remember
  # the syntax. plot.dir parameter was added later. And yes, it's redundant partially with outlist

  # filter.size: length 2. element 1: lower bound. element 2: higher bound.
  # if sequential filter is T, and fitler.size is set with quantile, then first filter out
  # those below lower bound and recalculate upper bound based on filtered vector
  if (!is.null(plot.dir)) {
    system(paste0("mkdir -p ", plot.dir))
  }
  if (!is.null(s2b.obj$dds) && force == F)
    stop("dds already present in s2b object. use force = T to overide")
  if (is.null(s2b.obj$ori.bulk.mat)) {
    s2b.obj$ori.bulk.mat <- s2b.obj$bulk.mat
  }

  if (use.ori == T) {
    s2b.obj$bulk.mat <- s2b.obj$ori.bulk.mat
  }

  if (filter.nz == T) {
    s2b.obj$bulk.mat <- s2b.obj$bulk.mat %>% .[rowMin(.) > 0,]
  }

  if (!is.null(filter.size) ) {
    filtered <- s2b.obj$bulk.mat %>%
      deseq2.filter(filter.size = filter.size, filter.fun = filter.fun,
                    filter.samples = filter.samples, sequential.filter = sequential.filter,
                    protected.genes = protected.genes,
                    plot.dir = plot.dir, root.name = s2b.obj$root.name)
    s2b.obj$bulk.mat <- filtered$mat
  }

  if (quantile.norm == T) {
    s2b.obj$bulk.mat.qn <- qn.fanc(s2b.obj$bulk.mat, T)
    mat.to.use <- s2b.obj$bulk.mat.qn %>% round()
    print("quantile normalization")
  } else {
    mat.to.use <- s2b.obj$bulk.mat
  }
  # if (Sys.getenv("HOME") == "/bar/cfan") {
  #   stop("to self: validate the numeric column removal step!")
  # }

  # debugging purpose:
  # s2b.obj$coldata$t <- as.numeric(s2b.obj$coldata$t)

  num.cols <- sapply(1:ncol(s2b.obj$coldata), function(i) return(is.numeric(s2b.obj$coldata[, i])))
  num.cols <- colnames(s2b.obj$coldata)[num.cols]

  if (length(num.cols) > 0) {
    stop(paste0("numeric columns detected: ", paste0(num.cols, collapse = "; "),
                ". Contact Changxu Fan if you think this is a false alarm."))
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat.to.use, colData = s2b.obj$coldata, design = design)

  if (single.sample == T) {
    dds <- DESeq2::estimateSizeFactors(dds, type = norm.method)
  } else {
    if (! is.null(locfunc)) {
      print("using non-median locfunc")
    } else {
      print("using median locfunc")
      locfunc <- stats::median
      # dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "local", sfType = norm.method)
    }
    if (!is.null(norm.cr.genes)) {
      stop("norm.cr.genes currently doesn't work")
    }
    dds <- estimateSizeFactors(dds, type = norm.method, locfunc = locfunc)
    if (use.genewise.disperisons) {
      # this might be needed for microarray data
      dds <- estimateDispersionsGeneEst(dds)
      dispersions(dds) <- mcols(dds)$dispGeneEst
    } else {
      dds <- estimateDispersions(dds, fitType = "local")
    }
    dds <- nbinomWaldTest(object = dds)
  }
  s2b.obj[["dds"]]  <- dds
  # counts <- SummarizedExperiment::assay(dds) %>% as.matrix()
  # norm.counts <- counts %*% diag((1/dds@colData$sizeFactor))
  # df <- as.data.frame(norm.counts)
  # colnames(df) <- as.character(dds@colData$sample)
  df <- DESeq2::counts(object = dds, normalized = T) %>% as.data.frame()
  # merge pseudo replicates if there are any:
  # if (single.sample == T) {
  #   for (i in 2:3) {
  #     df[, contrast[i]] <- df[, grepl(paste0(contrast[i], "..prep"), colnames(df))] %>% rowSums()
  #   }
  #
  #   df$fc <- (df[, contrast[2]]/df[, contrast[3]]) %>% format(digits = 3)
  # }


 #  colnames(df) <- colnames(df) %>% paste0("bulkNorm_", .)
  df$gene <- rownames(df)
  rownames(df) <- NULL
  s2b.obj[["bulkNorm"]] <- df

  if (single.sample == T)
    return(s2b.obj)

  if (!is.null(pca.ntop)) {
    trash <- deseq2.assess.core(dds = s2b.obj$dds,
                                plot.out = paste0(plot.dir, "/", s2b.obj$root.name, "_deseq2_assess_",pca.blind,".png"),
                                pca.ntop = pca.ntop,
                                pca.groupings = pca.groupings,
                                return.plot.list = F, blind = pca.blind,
                                add.hm = pca.add.hm, add.3d = pca.add.3d)

  }
  if (!is.null(contrast)) {
    s2b.obj[["res"]] <- DESeq2::results(object = dds, contrast = contrast,
                                        independentFiltering = independentFiltering)
    # res.exp means "result.expanded"
    s2b.obj[["res.exp"]] <- deseq2.add.value(dds=s2b.obj$dds, res=s2b.obj$res, sample.order = sample.order)
    rownames(s2b.obj[["res.exp"]]) <- NULL

    # s2b.obj[["res.shrink"]] <- DESeq2::lfcShrink(dds = s2b.obj$dds,
    #                                      coef = paste0(contrast[1], "_", contrast[2], "_vs_", contrast[3]),
    #                                      type="apeglm")
    # s2b.obj[["res.shrink.exp"]] <- deseq2.add.value(dds=s2b.obj$dds, res=s2b.obj$res.shrink,
    #                                                 sample.order = sample.order)
    # rownames(s2b.obj[["res.shrink.exp"]]) <- NULL

    # generate assess plots:

    if (try.hm == T) {
      plot.expr <- expression(exp.hm(
        res.exp = s2b.obj[["res.exp"]], outlist = outlist,
        sample.info = rowname.shffl(df = s2b.obj$coldata, direction = "in", which.col = "sample",
                                    keep = F, pos = 1),
        sample.order = sample.order,
        p.adj.cutoff = p.adj.cutoff, p.cutoff = p.cutoff,
        title = title,
        pivots = pivots
      ))

      if (force.hm == T)
        eval(plot.expr)
      else
        try(eval(plot.expr))
    }
  }

  return(s2b.obj)
}

deseq2.add.value <- function(dds, res, sample.order= NULL) {
  res <- res %>% as.data.frame()

  res <- utilsFanc::add.column.fanc(df1 = res, df2 = data.frame(gene = rownames(res)), pos = 1)
  rownames(res) <- NULL

  exp.mat <- SummarizedExperiment::assays(dds)$counts
  size.vec <- DESeq2::sizeFactors(dds)
  if ( !identical(colnames(exp.mat), names(size.vec))) {
    stop("colnames of exp.mat doesn't match names of size.vec")
  }
  samples <- names(size.vec)

  norm.df <- exp.mat %*% diag((1/size.vec)) %>% as.data.frame()
  colnames(norm.df) <- samples

  if (!is.null(sample.order))
    norm.df <- norm.df[, sample.order]

  norm.df <- utilsFanc::add.column.fanc(df1= norm.df, df2 = data.frame(gene = rownames(norm.df)), pos =1)
  rownames(norm.df) <- NULL

  res <- left_join(res, norm.df)

  return(res)

}

deseq2.shrink <- function(pbl, slot.name, clusters = NULL, threads = 1) {
  print("note: a2bl.add.shrinkage does the same thing as deseq2.shrink")
  if (is.null(clusters))
    clusters <- names(pbl)
  pbl <- utilsFanc::safelapply(pbl, function(s2b) {
    if (! s2b$root.name %in% clusters) {
      return(s2b)
    }
    s2b[[slot.name]] <- lfcShrink(dds = s2b$dds, res = s2b$res, type = "ashr")
    return(s2b)
  }, threads = threads)
  return(pbl)
}
hm.core <- function(plot.mat, res.exp, genes, meta.include,
                    samples, sample.info, dense.hm = T,
                    add.basemean = F, sample.order = NULL,
                    gene.order = NULL,   show.row.names = T,
                    show.row.dend = T, ...) {

  if (length(genes) < 1)
    return()
  if (dense.hm == T) {
    show.row.names <- F
    show.row.dend <- F
  }

  plot.mat <- plot.mat[genes,]
  res.exp <- res.exp[genes,]

  if (is.character(gene.order )) {
    gene.order <- plot.mat %>% .[rev(order(utilsFanc::pmean(.[, gene.order]))),] %>% rownames()
  }

  sa <- ComplexHeatmap::HeatmapAnnotation(
    df = sample.info[samples, meta.include]
  )
  hm.param <- list(matrix = plot.mat, top_annotation = sa, show_row_names = show.row.names,
                   show_row_dend = show.row.dend, ...)
  if (!is.null(sample.order))
    hm.param[["column_order"]] <- sample.order
  if (!is.null(gene.order))
    hm.param[["row_order"]] <- gene.order
  if (add.basemean == T) {

    ea <- ComplexHeatmap::rowAnnotation(log2mean = ComplexHeatmap::anno_barplot(log2(res.exp[genes, ]$baseMean + 1)),
                                        annotation_name_rot = 90)
    hm.param[["right_annotation"]] <- ea
  }

  hm <- do.call(ComplexHeatmap::Heatmap, hm.param)



  return(hm)
}


hm.core.2 <- function(s2b = NULL, use.bubble = F, bubble.log2 = F,
                      plot.mat, res.exp, summ.slot = "summary",
                      genes = NULL, samples = NULL,
                      col.data, metas = NULL, split.scaling.by = NULL,
                      sample.order = NULL, gene.order = NULL, pivots = NULL,
                      dense.hm = T, dense.cutoff = 30, show.row.names = T, show.row.dend = T,
                      add.basemean = F,
                      plot.out = NULL,
                      width = NULL, height = NULL, res = 100,
                      ...) {

  if (!is.null(s2b)) {
    exp.mat <- s2b$bulkNorm %>% `rownames<-`(., .$gene) %>% .[, !colnames(.) %in% "gene"]
    res.exp <- s2b$res.exp %>% `rownames<-`(., .$gene)
    if (!is.null(genes)) {
      if (genes[1] == "all")
        genes <- rownames(exp.mat)
    } else {
      if (!is.null(s2b[[summ.slot]])) {
        genes <- s2b[[summ.slot]]$de.genes
      } else {
        stop("when genes are not specified, the summ.slot must be present in s2b")
      }
    }
    # importantly, keep plot.mat and res.exp in sync
    genes <- genes %>% .[. %in% rownames(exp.mat)]
    if (length(genes) < 1) {
      stop("length(genes) < 1")
    }
    exp.mat <- exp.mat[genes, ]
    res.exp <- res.exp[genes, ]
    col.data <- s2b$coldata %>% as.data.frame()
    metas <- colnames(col.data)
    if (!is.null(split.scaling.by)) {
      plot.mat <- col.data %>% filter(sample %in% colnames(exp.mat)) %>%
        split(., f = .[, split.scaling.by]) %>%
        lapply(function(df) {
          exp.mat[, df$sample] %>% t() %>% scale(center = T, scale = T) %>% t() %>%
            return()
        }) %>% Reduce(cbind, .)

    } else {
      plot.mat <- exp.mat %>% t() %>% scale(center = T, scale = T) %>% t()
    }

  }

  if (nrow(plot.mat) < 1) {
    print("nrow(plot.mat) < 1")
    return()
  }

  if (is.null(genes)) {
    genes <- rownames(plot.mat)
    # must have the genes variable to sync things
  }
  utilsFanc::check.intersect(genes, "genes", rownames(plot.mat), "rownames(plot.mat)")
  plot.mat <- plot.mat[genes,]

  if (!is.null(sample.order)) {
    samples <- sample.order
  }
  if (is.null(samples)) {
    samples <- colnames(plot.mat)
  }
  utilsFanc::check.intersect(samples, "samples", colnames(plot.mat), "colnames(plot.mat)")
  plot.mat <- plot.mat[, samples]

  if (dense.hm == T && length(genes) > dense.cutoff) {
    show.row.names <- F
    show.row.dend <- F
  }

  hm.param <- list(matrix = plot.mat, show_row_names = show.row.names,
                   show_row_dend = show.row.dend, ...)
  if (use.bubble) {
    exp.mat <- exp.mat[rownames(plot.mat),colnames(plot.mat), drop = F] %>% as.matrix()
    if (bubble.log2) {
      exp.mat <- log2(exp.mat + 1)
    }
    max.exp <- quantile(exp.mat, 0.95)
    exp.mat[exp.mat > max.exp] <- max.exp
    col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    hm.param$col <- col_fun
    hm.param$rect_gp <- grid::gpar(type = "none")
    hm.param$cell_fun <- function(j, i, x, y, width, height, fill) {
      # if (i == 2 && j == 2)
      # print(paste0("width: ", width, "; height: ", height))
      # print(paste0("x: ", x, "; y: ", y))
      # size <- ifelse(as.numeric(width) < as.numeric(height), width, height)
      # r <- ((exp.mat[i, j])/max.exp) * 5 * size
      size <- grid::unit(20/res, units = "inches")
      grid::grid.circle(x = x, y = y,
                        r = ((exp.mat[i, j])/max.exp) * size,
                        gp = grid::gpar(fill = col_fun(plot.mat[i, j])))
    }

  }

  if (!is.null(metas)) {
    if ("sample" %in% metas && "Sample" %in% metas) {
      metas <- metas[metas != "Sample"]
    }
    utilsFanc::check.intersect(x = samples, "samples", rownames(col.data), "rownames(col.data)")
    if (!is.null(sample.order)) {
      col.data <- col.data[samples,]
    }
    colors <- lapply(col.data[samples, metas], function(meta) {
      meta <- unique(meta)
      color.map <- utilsFanc::gg_color_hue(length(meta))
      names(color.map) <- meta %>% as.character()
      return(color.map)
    })

    sa <- ComplexHeatmap::HeatmapAnnotation(
      df = col.data[samples, metas], # this is a trick to add multiple annotations at the same time
      col = colors
    )
    hm.param[["top_annotation"]] <- sa
  }

  if (is.character(pivots)) {
    stop("pivots not tested")
    gene.order <- plot.mat %>% .[rev(order(utilsFanc::pmean(.[, pivots]))),] %>% rownames()
  }

  if (!is.null(sample.order))
    hm.param[["column_order"]] <- sample.order
  if (!is.null(gene.order))
    hm.param[["row_order"]] <- gene.order
  if (add.basemean == T) {
    utilsFanc::check.intersect(genes, "genes", res.exp$gene, "res.exp$gene")
    rownames(res.exp) <- res.exp$gene
    baseMean <- res.exp[genes, "baseMean"]
    baseMean <- utilsFanc::log2pm(baseMean)
    ea <- ComplexHeatmap::rowAnnotation(log2mean = ComplexHeatmap::anno_barplot(baseMean),
                                        annotation_name_rot = 90)
    # manually validated that the pairing between baseMean bars and rows are correct.
    hm.param[["right_annotation"]] <- ea
  }

  hm <- do.call(ComplexHeatmap::Heatmap, hm.param)
  if (!is.null(plot.out)) {
    if (is.null(height)) {
      n.lines <- nrow(plot.mat) + length(metas) + 1
      height <- 20 * n.lines + 100
    }
    if (is.null(width)) {
      width <- ncol(plot.mat) * 20 + 450
    }

    save.base.plot(p = hm, file = plot.out, width = width, height = height, res = res)
  }


  return(hm)
}

exp.hm <- function(res.exp, outlist = NULL, sample.info, sample.order=NULL,
                   meta.include=NULL, p.adj.cutoff = NULL, p.cutoff = NULL, samples = NULL,
                   title = character(0), pivots = NULL,
                   plot.all = F) {
  if (is.character(sample.info))
    sample.info <- read.table(sample.info, as.is = T, header = T)
  rownames(sample.info) <- sample.info$sample

  if (is.null(meta.include))
    meta.include <- colnames(sample.info) %>% .[!.%in% c("dir", "sample")]
  if (is.null(samples))
    samples <- sample.info$sample

  rownames(res.exp) <- res.exp$gene

  exp.df <- res.exp[, samples]
  rownames(exp.df) <- res.exp$gene
  exp.mat <- exp.df %>% as.matrix()
  # plot z score
  plot.mat <- exp.mat %>% t() %>% scale(center = T, scale = T) %>% t()

  hm.core.ez <- function(genes, dense.hm, add.basemean,sample.order = NULL,
                         gene.order = NULL, outlist = NULL, sub.name, root.name.internal, ...) {
    hm <- hm.core(plot.mat = plot.mat, res.exp = res.exp, genes = genes, sample.order = sample.order,
                  sample.info = sample.info, meta.include = meta.include, samples = samples,
                  gene.order = gene.order,
                  # heatmap_width = unit(5, "inches"),
                  # heatmap_height = unit(8, "inches"),
                  dense.hm = dense.hm, add.basemean = add.basemean, ...)
    if (!is.null(outlist)) {
      png(filename = utilsFanc::plot.name.construct(outlist = outlist,
                                                    sub.name = sub.name,
                                                    root.name.internal = root.name.internal),
          width = 8, height = 10, units = "in", res = 314)
      try(ComplexHeatmap::draw(hm))
      dev.off()
    }
    return(hm)

  }
  if (plot.all == T) {
    genes <- plot.mat %>% na.omit() %>% rownames()
    trash <- hm.core.ez(genes, dense.hm = T, add.basemean = F,
                        column_title = title, outlist = outlist,
                        sub.name = "all_genes_hm.png", root.name.internal = "")
  }
  if (!is.null(p.adj.cutoff)) {
    genes.padj <- res.exp %>% filter(padj <= p.adj.cutoff ) %>% pull(gene)

    hm <- hm.core.ez(genes.padj, dense.hm = F, add.basemean = T,
                     column_title = title, outlist = outlist,
                     sub.name = "p_adj_" %>% paste0(p.adj.cutoff,"_clust_hm.png"),
                     root.name.internal = "")

  }

  if (!is.null(p.cutoff)) {
    gene.p <- res.exp %>% filter(pvalue <= p.cutoff ) %>% pull(gene)
    hm <- hm.core.ez(gene.p, dense.hm = T, add.basemean = T,
                     column_title = title, outlist = outlist,
                     sub.name = "p_" %>% paste0(p.cutoff,"_clust_hm.png"),
                     root.name.internal = "")
  }

  t.f <- function() {
    if (!is.null(pivots) && !is.null(sample.order)) {
      for (s in pivots) {
        if (!is.null(p.adj.cutoff)) {
          trash <- hm.core.ez(genes.padj, dense.hm = F, add.basemean = T,
                              sample.order = sample.order, gene.order = s,
                              column_title = title, outlist = outlist,
                              sub.name = paste0("p_adj_",p.adj.cutoff,"_",s,"_hm.png"),
                              root.name.internal = "")

        }
        if (!is.null(p.cutoff)) {
          trash <- hm.core.ez(gene.p, dense.hm = T, add.basemean = T,
                              sample.order = sample.order, gene.order = s,
                              column_title = title, outlist = outlist,
                              sub.name = paste0("p_",p.adj.cutoff,"_",s,"_hm.png"),
                              root.name.internal = "")
        }
      }
    }

  }

  t.f()

  return()
}


s2b.plus.deseq <- function(so, assay, slot, cells=NULL, ident, sub.idents,
                           group.by, groups = NULL, coldata.columns, design, contrast = NULL) {
  s2b <- so.2.bulk(so = so, assay = assay, slot = slot, cells = cells,
                   ident = ident, sub.idents = sub.idents, group.by = group.by, groups = groups,
                   coldata.columns = coldata.columns)
  try(s2b <- s2b.deseq(s2b, design = design, contrast = contrast))
  return(s2b)
}

deseq.res.view <- function(s2b=NULL, s2b.slot = NULL, dds.res=NULL, rank.by = "padj", desc = F, View = T ) {
  if (is.null(dds.res)) {
    if (is.null(s2b.slot))
      s2b.slot <- "res"
    dds.res <- s2b[[s2b.slot]]
  }
  df <- dds.res %>% utilsFanc::change.name.fanc(rank.by, "trank") %>%
    as.data.frame() %>% {if (desc == T) arrange(., desc(trank)) else arrange(., trank)} %>%
    utilsFanc::change.name.fanc("trank", rank.by)
  if (View == T)
    return(View(df))
  else
    return(head(df))
}


bulk.list <- function(so = NULL, # pbl = NULL,
                      deseq2 = T,
                      assay = "RNA", slot = "counts",
                      group.by = "sample", groups = NULL, sample.order = NULL,
                      coldata.columns = NULL, design.formula, contrast,
                      cluster.ident = "seurat_clusters", clusters = NULL, prepend.cluster.ident = T,
                      filter.nz = F, filter.size = NULL, sequential.filter = T,
                      independentFiltering = F,
                      pca.ntop = 10000, pca.groupings = NULL,
                      p.adj.cutoff = 0.05, p.cutoff =0.05,
                      single.sample = F,
                      do.plot = F, work.dir, plot.dir=NULL,
                      threads = NULL, stop.on.error = T, # biocparallel
                      save.rds = T) {
  # you can either start from so or pbl. pbl: pseudobulk list. in the latter case you are just adding
  #the deseq2 component.
  utilsFanc::check.intersect(clusters, "clusters", so@meta.data[, cluster.ident], "so@meta.data[, cluster.ident]")
  system("mkdir -p " %>% paste0(work.dir, " "))
  if (is.null(plot.dir)) {
    plot.dir <- paste0(work.dir, "/plots/")
  }
  if (is.null(cluster.ident)) {
    cluster.ident <- "pd_cluster"
    so[[cluster.ident]] <- "0"
  }
  meta.df <- so@meta.data %>% factor2character.fanc()

  if (is.null(clusters)) {
    clusters <- meta.df[, cluster.ident] %>% unique() %>% gtools::mixedsort()
  }

  if (is.null(pca.groupings))
    pca.groupings <- coldata.columns


  if (is.null(threads)) {
    threads <- length(clusters)
  }
  threads <- min(12, threads)

  par.param <- MulticoreParam(workers = threads, stop.on.error = stop.on.error)
  bulk.all.list <- bptry({
    bplapply(clusters, function(cluster) {
      bulk.all <- so.2.bulk(so, root.name = paste0(ifelse(prepend.cluster.ident, paste0(cluster.ident, "_"), ""), cluster), 
                            assay = assay, slot= slot,
                            ident = cluster.ident, sub.idents = cluster,
                            group.by = group.by, groups = groups, coldata.columns = coldata.columns)
      if (deseq2 == T) {
        bulk.all <- s2b.deseq(s2b.obj = bulk.all, design = design.formula , contrast = contrast,
                              filter.nz = filter.nz, filter.size = filter.size, sequential.filter = sequential.filter,
                              independentFiltering = independentFiltering,
                              pca.ntop = pca.ntop, pca.groupings = pca.groupings,
                              force = T, sample.order = sample.order, try.hm = do.plot,
                              plot.dir = plot.dir,
                              outlist = list(plot.dir = plot.dir, root.name.external = "cluster_" %>% paste0(cluster, "_")),
                              p.adj.cutoff = p.adj.cutoff, p.cutoff = p.cutoff,
                              title = paste0("cluster_", cluster), single.sample = single.sample,
                              pivots = NULL)

      }
      return(bulk.all)
    }, BPPARAM = par.param)
  })
  if ("error" %in% class(bulk.all.list)) {
    stop(bulk.all.list$message)
  }
  # bulk.all.list <- utilsFanc::safelapply(clusters, function(cluster) {
  # }, threads = threads)
  names(bulk.all.list) <- lapply(bulk.all.list, function(x) return(x$root.name)) %>% unlist()
  attr(bulk.all.list, "work.dir") <- work.dir
  if (save.rds == T)
    saveRDS(bulk.all.list, paste0(work.dir, "/bulk.list.Rds"))
  return(bulk.all.list)
}

aggregate.fanc <- function(mat, se.assay.name = NULL, margin, groupings, na.rm = T, binarize = F, take.mean = F, sort = F) {
  ### groupings: must be a named vector. the names are the colnames or rownames. mat will be
  #rearranged to match the sequence of groupings
  se.flag <- F
  if (any(grepl("SummarizedExperiment", class(mat)))) {
    se <- mat
    se.flag <- T
    if (is.null(se.assay.name)) {
      if (length(assays(se)) == 1)
        se.assay.name <- names(assays(se))[1]
      else
        stop("is.null(se.assay.name)")
    }
    mat <- assays(se)[[se.assay.name]]
  }

  if (is.factor(groupings)) {
    names <- names(groupings)
    groupings <- as.character(groupings)
    names(groupings) <- names
  }

  if (na.rm == T) {
    groupings <- groupings %>% .[!is.na(.)]
  }

  if (sort == T) {
    groupings <- gtools::mixedsort(groupings)
  }
  # mat <- GetAssayData(object = soi.clean.WT, slot = "counts", assay = "RNA")
  if (margin == 2) {
    mat <- Matrix::t(mat)
  } else if (margin != 1) {
    stop("margin has to be 1 or 2")
  }
  utilsFanc::check.intersect(names(groupings), "names(groupings)",
                             rownames(mat), paste0(ifelse(margin == 2, "col", "row"), "names(mat)"))
  # groupings <- groupings %>% .[names(.) %in% rownames(mat)]
  mat <- mat[names(groupings), ] # important step that syncs groupings and the mat
  if (binarize == T) {
    mat[mat > 0] <- 1
    mat[mat <= 0] <- 0
  }

  aggr <- Matrix.utils::aggregate.Matrix(mat, groupings = groupings, fun = "sum")
  if (sort == T) {
    aggr <- aggr[gtools::mixedsort(rownames(aggr)),]
  }

  if (take.mean == T) {
    n <- table(groupings)[rownames(aggr)]
    names <- dimnames(aggr)
    aggr <- diag(1/n) %*% aggr
    dimnames(aggr) <- names
  }

  if (margin == 2) {
    aggr <- Matrix::t(aggr)
  }

  if (se.flag == T) {
    if (margin == 1) {
      res <- SummarizedExperiment(assays = SimpleList(aggr) %>% `names<-`(se.assay.name),
                                  colData = colData(se))
    } else {
      res <- SummarizedExperiment(assays = SimpleList(aggr) %>% `names<-`(se.assay.name),
                                  rowRanges = rowRanges(se))
      rowData(res) <- rowData(se)
    }
  } else {
    res <- aggr
  }

  return(res)
}

deseq2.filter <- function(mat, filter.size, filter.fun, filter.samples, sequential.filter = T, 
                          protected.genes = NULL,
                          plot.dir, root.name) {
  cpm.mat <- edgeR::cpm(mat)
  if (!is.null(filter.samples)) {
    cpm.mat <- cpm.mat %>% .[, colnames(.) %in% filter.samples]
    print(paste0("using these samples for filtering: \n",
                 paste0(filter.samples, collapse = "\n")))
  }
  if (filter.fun == "max") {
    filter.vec <- rowMax(cpm.mat)
  } else if (filter.fun == "mean") {
    filter.vec <- rowMeans(cpm.mat)
  } else {
    stop("filter.fun not recoganized. has to be 'max' or 'mean'")
  }

  if (filter.size[1] < 1)
    filter.size[1] <- quantile(filter.vec, filter.size[1])
  if (filter.size[2] <= 1) {
    if (sequential.filter == T) {
      filter.size[2] <- quantile(filter.vec %>% .[. > filter.size[1]], filter.size[2])
    } else {
      filter.size[2] <- quantile(filter.vec, filter.size[2])
    }
  }

  bPass <- (filter.vec > filter.size[1] & filter.vec < filter.size[2])
  if (!is.null(protected.genes)) {
    bPass <- bPass | (rownames(cpm.mat) %in% protected.genes)
  }
  trash <- rank.plot(df = data.frame(value = filter.vec) %>% `colnames<-`(filter.fun), vars = filter.fun,
                     add.h.line = filter.size, title = sum(bPass),
                     outfile = paste0(plot.dir,"/", root.name,"_", filter.fun,"_",
                                      filter.size[1], "_", filter.size[2],".rank.png"),
                     quantile.limit.y = 0.999)

  trash <- rank.plot(df = data.frame(value = filter.vec) %>% `colnames<-`(filter.fun), vars = filter.fun,
                     add.h.line = filter.size, title = sum(bPass),
                     transformation = function(x) log2(x + 1),
                     outfile = paste0(plot.dir, "/", root.name, "_", filter.fun, "_",
                                      filter.size[1], "_", filter.size[2],".rank.log.png"))
  res <- list(mat = mat[bPass, ], mat.cpm = cpm.mat[bPass, ],
              filter.size = filter.size, filter.fun = filter.fun, filter.samples = filter.samples,
              sequential.filter = sequential.filter,
              root.name = root.name)
  return(res)
}

mat.xy <- function(mat, plot.out = NULL, trend = T, ...) {
  if (trend == T ) {
    mat <- mat %>% .[sample(1:nrow(.), 10000, replace = F),]
  }
  p <- mat %>% as.data.frame() %>%
    xy.plot(x = "BM", y = "SP", is.regex = T, color.density = T, ...)
  if (trend == T) {
    p <- p + geom_smooth(method = "loess")
  }
  if (!is.null(plot.out)) {
    wrap.plots.fanc(plot.list = list(p), plot.out = plot.out)
  }
  return(p)
}


deseq2.topn.norm <- function(dds, topn.vec, plot.dir, root.name) {
  counts <- counts(dds)
  counts.cpm <- edgeR::cpm(counts)
  norm.mats <- lapply(topn.vec, function(n) {
    cutoff <- counts.cpm %>% rowMeans() %>% sort() %>% rev() %>% .[n]
    bCg <- counts.cpm %>% rowMeans() %>% `>=`(cutoff)
    cg <- rownames(counts.cpm)[bCg]
    dds <- estimateSizeFactors(dds, controlGenes=bCg)
    norm.mat <- counts(object = dds, normalized = T)
    df <- norm.mat %>% as.data.frame()
    df$gene <- rownames(df)
    lapply(c("hl", "noHl"), function(type) {
      params <- list(mat = df, plot.out = paste0(plot.dir, "/", root.name, "_", n, "_", type, "_xy.png"), trend = F)
      if (type == "hl") {
        params <- c(params, list(highlight.var = "gene", highlight.values = cg))
      }
      p <- do.call(mat.xy, params)
    })
    return(norm.mat)
  })
  return(norm.mats)
}

qn.fanc <- function(mat, ori.size) {
  qn <- preprocessCore::normalize.quantiles(x = mat)
  if (ori.size == T) {
    qn <- qn %*% diag((colSums(mat)/colSums(qn)))
  }
  rownames(qn) <- rownames(mat)
  colnames(qn) <- colnames(mat)
  return(qn)
}

s2b.hm.single.gene <- function(so, assay, slot,
                               group.ident.1, groups.1 = NULL, group.ident.2, groups.2 = NULL,
                               genes,
                               maintain.group.order = T, maintain.gene.order = F,
                               plot.ind = T, 
                               height.add = 400,
                               plot.dir, root.name = NULL) {
  if (is.null(root.name)) {
    root.name <- basename(plot.dir)
  }
  # group.ident.1: first layer, group.idetn.2: second layer.
  # you can imagine: colnames(mat) to be in the format of group.ident.1..group.ident.2

  #############
  ## incorrect: normalization was done individually by DESeq2...
  # mats <- lapply(s2bl, function(s2b) {
  #   if (is.null(samples)) {
  #     samples <- s2b$bulkNorm %>% colnames() %>%  .[.!= "gene"]
  #   }
  #   samples <- samples %>% .[. %in% colnames(s2b$bulkNorm)]
  #   genes <- genes %>% .[.%in% s2b$bulkNorm$gene]
  #
  #   rownames(s2b$bulkNorm) <- s2b$bulkNorm$gene
  #   mat <- s2b$bulkNorm %>% .[, colnames(.) != "gene"] %>% as.matrix()
  #   mat <- mat[genes, samples]
  #   colnames(mat) <- paste0(s2b$root.name %>% sub("seurat_clusters_", "", .), "_", colnames(mat))
  #   return(mat)
  # })
  # genes <- Reduce(intersect, lapply(mats, rownames))
  # if (length(genes) < 1) {
  #   stop("length(genes) < 1")
  # }
  #
  # mat <- lapply(mats, function(x) return(x[genes, ])) %>%
  #   Reduce(cbind, .)

  groupings <- get.cell.list(obj = so, is.ao = F, split.by = group.ident.2, splits = groups.2,
                             group.by = group.ident.1, groups = groups.1, na.rm = T, return.named.vec = T)

  mat.big <- Seurat::GetAssayData(object = so, slot = slot, assay = assay)

  mat <- aggregate.fanc(mat = mat.big, margin = 2, groupings = groupings,
                             na.rm = T, take.mean = T, sort = T)
  genes <- genes %>% .[. %in% rownames(mat)]
  mat <- mat[genes,]
  # we first write out mat:
  mat.df <- cbind(data.frame(gene = rownames(mat)), as.data.frame(mat))
  dir.create(plot.dir, showWarnings = F, recursive = T)
  write.table(mat.df, paste0(plot.dir, "/", root.name, "_norm_mat.tsv"),
              col.names = T, row.names = F, sep = "\t", quote = F)
  mat.scale <- mat %>% t() %>% scale() %>% t()
  dimnames(mat.scale) <- dimnames(mat)
  
  if (maintain.group.order) {
    if (!is.null(groups.1) && !is.null(groups.2)) {
      groups <- outer(groups.1, groups.2, function(x, y) return(paste0(x, "..", y))) %>% 
        t() %>% as.vector()
      groups <- groups[groups %in% colnames(mat.scale)]
      mat.scale <- mat.scale[, groups]
    } else{
      stop("groups.1 and groups.2 must be supplied to use maintain.group.order")
    }
    
  }
  
  hm <- ComplexHeatmap::Heatmap(matrix = mat.scale, cluster_rows = !maintain.gene.order, cluster_columns = F,
                                show_row_dend = F, clustering_distance_rows = "pearson")
  save.base.plot(hm, file = paste0(plot.dir, "/", root.name, "_norm_scaled_hm.pdf"),
                 width = 600, height = 15 * nrow(mat.scale) + height.add)
  row.order <- suppressWarnings(ComplexHeatmap::row_order(hm))  
  mat <- mat[row.order,]
  # now we plot individual ones:
  if (plot.ind) {
    ind.dir <- paste0(plot.dir, "/", root.name, "_ind/")
    pl <- lapply(genes, function(gene) {
      col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
      mat.gene <- mat.scale[gene, , drop = F]
      hm <- ComplexHeatmap::Heatmap(matrix = mat.gene, cluster_rows = F, cluster_columns = F, col = col_fun,
                                    column_names_side = "top")
      save.base.plot(hm, file = paste0(ind.dir, "/", root.name, "_", gene, "_hm.png"),
                     width = 600, height = 150)
      df <- data.frame(sample = colnames(mat.gene), expr = mat[gene, ])
      df$sample <- factor(df$sample, levels = rev(df$sample))
      p <- ggplot(df, aes(x = sample, y = expr)) +
        geom_bar(stat = "identity", color = "gray50", fill = "gray50") +
        coord_flip() +
        # theme(axis.line.x = element_blank(), axis.title.x = element_blank(), axis.text.x =  element_blank(),
        #       axis.ticks.x = element_blank()) +
        theme_void() +
        ggtitle(gene) +
        ggsave(paste0(ind.dir, "/", root.name, "_", gene, "_bar.png"), width = 1, height = 5)
      return(p)
    })
    trash <- wrap.plots.fanc(plot.list = pl, sub.height = 5, sub.width = 1,
                             plot.out = paste0(plot.dir, "/", root.name, "_norm_bar.png"))
  }
  return(mat)
}

sc.peaks.for.bulk <- function(scde, sample.info, single.sample = T,
                              npar.clusters = 1, npar.bamCount = 1,
                              work.dir = NULL) {
  # programmed for: ~/hmtp/scAR/ery2/da/da1.1_scPeaks_2023-07-27.R
  # it was programed to support single sample comparisons only
  # scde: the da object for single cell
  if (!single.sample) {
    stop("only single sample comparisons have been developed")
  }

  required.cols <- c("sample", "bamReads")
  utilsFanc::check.intersect(required.cols, "required columns",
                             colnames(sample.info), "colnames(sample.info)")
  sample.info$celltype <- "a"
  sample.info$tissue <- paste0("t", 1:nrow(sample.info))

  if (is.character(scde)) {
    scde <- readRDS(scde)
  }

  # bde: bulk de
  bde <- utilsFanc::safelapply(scde, function(s2b) {
    bulk.mat <- bam.count(sample.info = sample.info, bam.col = "bamReads", sample.col = "sample",
                         features.gr = utilsFanc::loci.2.gr(rownames(s2b$bulk.mat)),
                         bSingleEnd = T, sort = T, return.mat = T, threads = npar.bamCount)

    da <- diffbind.2.raw.a2bl(sample.info = sample.info, mat = bulk.mat,
                              cluster.ident = "celltype",sample.col = "sample",
                              filter.nz = F, filter.size = NULL, sequential.filter = T,
                              coldata.columns = c("tissue"),
                              single.sample = T,
                              design.formula = ~ tissue,
                              threads = 1,
                              work.dir = tempdir())
    b2b <- da[[1]]
    b2b$root.name <- s2b$root.name
    return(b2b)
  }, threads = npar.clusters)
  names(bde) <- names(scde)
  if (!is.null(work.dir)) {
    dir.create(work.dir)
    saveRDS(bde, paste0(work.dir, "/bulk.list.Rds"))
  }
  return(bde)
}


