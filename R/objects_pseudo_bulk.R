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
  
  meta.df <- meta.df %>% utilsFanc::change.name.fanc(cols.from = group.by, cols.to = "tgroup")
  if (is.null(coldata.columns)) {
    meta.df[, group.by] <- meta.df$tgroup
    coldata.columns <- group.by
  }
    
  meta.df$cell.id <- rownames(meta.df)
  if (is.null(cells)) {
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
  # browser()
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
      coldata.df <- getCellColData(ao, select = c(group.ident, coldata.columns), drop = F) %>% 
        as.data.frame() %>% unique() %>% `rownames<-`(., .[, group.ident])
    }
  } 

  # from gr to DESeq2 format:
  df <- gr %>% `names<-`(NULL) %>% as.data.frame() %>% .[, c("seqnames", "start", "end", gr.sample.cols)] %>% 
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
      # browser()
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

ao.2.bulk.list <- function(ao, peakmat, pseudo, sample.order = NULL, cell.prep.dir,
                           filter.nz = F, filter.size = NULL, sequential.filter,
                           coldata.columns = NULL, coldata.df = NULL,
                           pca.ntop = 10000, pca.groupings = NULL,
                           cluster.ident, clusters = NULL, group.ident, groups = NULL, 
                           match.bg.df = NULL, 
                           filter.by = NULL, filter.limits, bQuantileFilter = c(F, F), 
                           design.formula, contrast,
                           do.plot = T, 
                           work.dir,threads = 6, plot.dir = NULL,
                           single.sample = F) {
  if (is.null(plot.dir))
    plot.dir <- paste0(work.dir, "/plots")
  system(paste0("mkdir -p ", work.dir, " ", plot.dir))
  
  if (is.null(pca.groupings))
    pca.groupings <- coldata.columns
  
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
    a2b <- s2b.deseq(s2b.obj = a2b, 
                     filter.nz = filter.nz, filter.size = filter.size, sequential.filter = sequential.filter,
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
  # browser()
  names(a2b.list) <- paste0(cluster.ident, "_",clusters)
  saveRDS(a2b.list, paste0(work.dir, "/ao.bulk.list.Rds"))
  return(a2b.list)
}
  
  
s2b.deseq <- function(s2b.obj, 
                      design, contrast=NULL, 
                      use.ori = T, filter.nz = F, filter.size = NULL, sequential.filter = F,
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
    max.vec <- rowMax(edgeR::cpm(s2b.obj$bulk.mat))

    if (filter.size[1] <= 1)
      filter.size[1] <- quantile(max.vec, filter.size[1])
    if (filter.size[2] <= 1) {
      if (sequential.filter == T) {
        filter.size[2] <- quantile(max.vec %>% .[. > filter.size[1]], filter.size[2])
      } else {
        filter.size[2] <- quantile(max.vec, filter.size[2])
      }
    }
      
    bPass <- (max.vec > filter.size[1] & max.vec < filter.size[2])
    trash <- rank.plot(df = data.frame(max = max.vec), vars = "max",
                       add.h.line = filter.size, title = sum(bPass),
                       outfile = paste0(plot.dir,"/", s2b.obj$root.name,"_max.rank.png"),
                       quantile.limit.y = 0.999)
    
    trash <- rank.plot(df = data.frame(max = max.vec), vars = "max", 
                       add.h.line = filter.size, title = sum(bPass),
                       transformation = function(x) log2(x + 1), 
                       outfile = paste0(plot.dir, "/", s2b.obj$root.name, "_max.rank.log.png"))

    s2b.obj$bulk.mat <- s2b.obj$bulk.mat %>% .[max.vec > filter.size[1] & max.vec < filter.size[2],] 
  }
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = s2b.obj$bulk.mat, colData = s2b.obj$coldata, design = design)
  
  if (single.sample == T) {
    dds <- DESeq2::estimateSizeFactors(dds, type = "ratio")
  } else {
    dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "local", sfType = "ratio") 
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
    s2b.obj[["res"]] <- DESeq2::results(object = dds, contrast = contrast)
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

# s2b.add.shrink <- function(s2b, ) {
#   
# }

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
    gene.order <- plot.mat %>% .[rev(order(.[, gene.order])),] %>% rownames()
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

exp.hm <- function(res.exp, outlist = NULL, sample.info, sample.order=NULL, 
                   meta.include=NULL, p.adj.cutoff, p.cutoff, samples = NULL,
                   title = character(0), pivots = NULL) {
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
                  gene.order = gene.order, heatmap_width = unit(5, "inches"), 
                  heatmap_height = unit(8, "inches"),
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
  
  genes <- plot.mat %>% na.omit() %>% rownames()
  trash <- hm.core.ez(genes, dense.hm = T, add.basemean = F,
                      column_title = title, outlist = outlist,
                      sub.name = "all_genes_hm.png", root.name.internal = "")
  
  genes.padj <- res.exp %>% filter(padj <= p.adj.cutoff ) %>% pull(gene)
  
  hm <- hm.core.ez(genes.padj, dense.hm = F, add.basemean = T,
                   column_title = title, outlist = outlist,
                   sub.name = "p_adj_" %>% paste0(p.adj.cutoff,"_clust_hm.png"), 
                   root.name.internal = "")
  
  gene.p <- res.exp %>% filter(pvalue <= p.cutoff ) %>% pull(gene)
  hm <- hm.core.ez(gene.p, dense.hm = T, add.basemean = T,
                   column_title = title, outlist = outlist,
                   sub.name = "p_" %>% paste0(p.cutoff,"_clust_hm.png"), 
                   root.name.internal = "")
  
  t.f <- function() {
    if (!is.null(pivots) && !is.null(sample.order)) {
      for (s in pivots) {
        trash <- hm.core.ez(genes.padj, dense.hm = F, add.basemean = T,
                            sample.order = sample.order, gene.order = s,
                            column_title = title, outlist = outlist,
                            sub.name = paste0("p_adj_",p.adj.cutoff,"_",s,"_hm.png"), 
                            root.name.internal = "")
        trash <- hm.core.ez(gene.p, dense.hm = T, add.basemean = T,
                            sample.order = sample.order, gene.order = s,
                            column_title = title, outlist = outlist,
                            sub.name = paste0("p_",p.adj.cutoff,"_",s,"_hm.png"),  
                            root.name.internal = "")
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
                      cluster.ident = "seurat_clusters", clusters = NULL,
                      filter.nz = T, filter.size = NULL, sequential.filter = T,
                      pca.ntop = 10000, pca.groupings = NULL,
                      p.adj.cutoff = 0.05, p.cutoff =0.05, 
                      single.sample = F,
                      do.plot = F, work.dir, plot.dir=NULL, threads = NULL, save.rds = T) {
  # you can either start from so or pbl. pbl: pseudobulk list. in the latter case you are just adding
  #the deseq2 component.
  system("mkdir -p " %>% paste0(work.dir, " "))
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
  
  bulk.all.list <- utilsFanc::safelapply(clusters, function(cluster) {
    # try({
    bulk.all <- so.2.bulk(so, root.name = paste0(cluster.ident, "_", cluster), assay = assay, slot= slot, 
                          ident = cluster.ident, sub.idents = cluster,
                            group.by = group.by, groups = groups, coldata.columns = coldata.columns)
    if (deseq2 == T) {
      bulk.all <- s2b.deseq(s2b.obj = bulk.all, design = design.formula , contrast = contrast,
                            filter.nz = filter.nz, filter.size = filter.size, sequential.filter = sequential.filter,
                            pca.ntop = pca.ntop, pca.groupings = pca.groupings,
                            force = T, sample.order = sample.order, try.hm = do.plot,
                            plot.dir = plot.dir, 
                            outlist = list(plot.dir = plot.dir, root.name.external = "cluster_" %>% paste0(cluster, "_")), 
                            p.adj.cutoff = p.adj.cutoff, p.cutoff = p.cutoff,
                            title = paste0("cluster_", cluster), single.sample = single.sample,
                            pivots = NULL)
    }

      # rm(list = ls() %>% .[!. == "bulk.all"])
      # gc()
      return(bulk.all)
    # })
  }, threads = threads)
  names(bulk.all.list) <- paste0(cluster.ident, "_", clusters)
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
  
  groupings <- groupings %>% .[names(.) %in% rownames(mat)]
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

