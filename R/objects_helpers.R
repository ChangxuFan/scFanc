count.stat <- function(so, metas, outfile = NULL) {
  df <- data.frame(cluster = so@meta.data$seurat_clusters %>% unique() %>% sort() %>% as.character())
  for (meta in metas) {
    meta.mn <- paste0(meta, "_mean")
    meta.md <- paste0(meta, "_median")
    df[, meta.mn] <- NA
    df[, meta.md] <- NA
    for (i in 1:nrow(df)) {
      df[i, meta.mn] <- so@meta.data %>% .[.$seurat_clusters == df[i, "cluster"], meta] %>% mean
      df[i, meta.md] <- so@meta.data %>% .[.$seurat_clusters == df[i, "cluster"], meta] %>% median
    }
  }
  if (!is.null(outfile)) {
    write.table(df, outfile, sep = "\t", col.names = T, row.names = F, quote = F)
  }
  return(df)
}

cluster.subset <- function(so, clusters) {
  if (!is.character(clusters))
    stop("cluster names must be characters")
  
  meta <- so@meta.data
  meta$seurat_clusters <- as.character(meta$seurat_clusters)
  cells <- rownames(meta)[meta$seurat_clusters %in% clusters]
  return(so[, cells])
}

count.by.cluster <- function(so,work.dir= NULL, root.name = NULL, meta="seurat_clusters", 
                             group.by = NULL, groups.list = NULL, mixedSort = T,
                             meta.df = NULL) {
  if (is.null(meta.df))
    meta.df <- so@meta.data %>% factor2character.fanc()
  if (is.null(group.by))
    group.by <- "sample"
  if (is.null(groups.list)) {
    groups.list <- meta.df[, group.by] %>% unique() %>% as.list()
  }
  groups.list <- c(list(unlist(groups.list)), groups.list)
  
  df <- lapply(seq_along(groups.list), function(i) {
    group.name <- paste0(groups.list[[i]], collapse = "|") %>% paste0(group.by, "_", .)
    if (i == 1) 
      group.name <- paste0(group.by, "_all")
    df <- meta.df %>% .[.[,group.by] %in% groups.list[[i]], meta] %>% table() %>% as.data.frame()
    colnames(df) <- c("cluster", "freq")
    df$frac <- df$freq/sum(df$freq)
    sum <- df
    sum$cluster <- NULL
    sum <- lapply(sum, sum) %>% as.data.frame()
    sum$cluster <- "sum"
    df <- rbind(df, sum)
    colnames(df) <- c("cluster", paste0(group.name, "_freq"), paste0(group.name, "_frac"))
    df <- df %>% factor2character.fanc()
    if (mixedSort == T) {
      df <- df[gtools::mixedorder(df$cluster), ]
    }
    return(df)
  }) %>% Reduce(left_join, .)
  
  if (!is.null(root.name))
    root.name <- paste0(root.name, "_")
  else
    root.name <- ""
  if (!is.null(work.dir))
    write.table(df, work.dir %>% paste0("/",root.name,"cluster_freq.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
  return(df)
}
rbind.fanc <- function(df1, df2) {
  cols <- colnames(df2) %>% .[!.%in% colnames(df1)]
  df1[,cols] <- NA
  
  cols <- colnames(df1) %>% .[!.%in% colnames(df2)]
  df2[,cols] <- NA
  
  return(rbind(df1, df2))
  
}
# int.corr <- function(q, s, meta) {
#   # sync sample
#   q <- q@meta.data
#   s <- s@meta.data
#   q.sample <- q$sample %>% as.character() %>% unique()
#   s.sample <- s$sample %>% as.character() %>% unique()
#   f.sample <- intersect(q.sample, s.sample)
#   
#   q <- q[q$sample %in% f.sample, ]
#   s <- s[s$sample %in% f.sample, ]
#   
#   q$cells <- rownames(q) %>% sub("_.+$", "",.)
#   s$cells <- rownames(s) %>% sub("_.+$", "",.)
#   
#   stats <- q %>% split(., f=factor(.[, meta])) %>% 
#     lapply(function(x) {
#       cells <- rownames(x) 
#       s.sub <- s[cells, ]
#       return(table(s.sub[, meta]))
#     })
#   return(stats)
# }

int.corr <- function(q.vec, s, meta, out.file = NULL) {
  s <- s@meta.data %>% factor2character.fanc()
  stats.full <- lapply(q.vec, function(q) {
    q <- q@meta.data %>% factor2character.fanc()
    q.sample <- q$sample %>% unique()
    s.sample <- s$sample %>% unique()
    f.sample <- intersect(q.sample, s.sample)
    
    q <- q[q$sample %in% f.sample, ]
    s <- s[s$sample %in% f.sample, ]
    
    q$cells <- rownames(q) %>% sub("_.+$", "",.)
    s$cells <- rownames(s) %>% sub("_.+$", "",.)
    # sync q and s so that pipe doesn't error out when some cells are missing:
    f.cells <- intersect(q$cells, s$cells)
    q <- q[q$cells %in% f.cells,]
    s <- s[s$cells %in% f.cells,]
    
    stats <- q %>% split(., f=factor(.[, meta])) %>%
      lapply(function(x) {
        #print(x$seurat_clusters[1])
        # if (x$seurat_clusters[1] == "20")
        #   browser()
        cells <- x$cells
        s.sub <- s[cells, ]
        stat <- table(s.sub[, meta])  %>% as.data.frame()
        rownames(stat) <- paste0("int__",stat$Var1)
        stat$Var1 <- NULL
        colnames(stat) <- x[, meta][1]
        stat <- t(stat) %>% as.data.frame()
        stat <- stat[, gtools::mixedsort(colnames(stat)), drop = F]
        return(stat)
      })  
    # print("miao")
    stats <- stats %>% Reduce(rbind.fanc, .)
    stats[is.na(stats)] <- 0
    rownames(stats) <- paste0(q.sample, "__", rownames(stats))
    return(stats)
  }) %>% Reduce(rbind.fanc, .)
  stats.full <- stats.full[gtools::mixedsort(rownames(stats.full)), gtools::mixedsort(colnames(stats.full))]
  stats.full$total <- apply(stats.full, 1, sum)
  stats.full["total", ] <- apply(stats.full, 2, sum)
  if (!is.null(out.file)) {
    write.table(stats.full, out.file, quote = F, row.names = T, col.names = T, sep = "\t")
  }
  return(stats.full)

}

get.cell.names.seurat <- function(so, filter.list = NULL, cells = NULL, style = "Seurat") {
  meta.df <- so@meta.data %>% factor2character.fanc()
  if (!is.null(cells))
    meta.df <- meta.df[cells, ]
  
  if (!is.null(filter.list)) {
    for (i in 1:length(filter.list)) {
      meta <- names(filter.list)[i]
      values <- filter.list[[i]]
      meta.df <- meta.df %>% .[.[, meta] %in% values,]
    }
  }

  if (style == "Seurat") {
    return(rownames(meta.df))
  }
  
  if (style == "ArchR") {
    return(paste0(meta.df$sample, "#", sub("_.+$", "", rownames(meta.df))))
  }
}

add.metrics.seurat <- function(so, bc.metrics.file.list, columns.to.add = NULL) {
  # bc.metrics.file.list should be a named list like this:
  ##
  meta.df <- so@meta.data %>% factor2character.fanc()
  meta.df$barcode <- rownames(meta.df) %>% sub("_.+$", "", .)
  metrics.df <- meta.df %>% split(., f = factor(.$sample, levels = unique(.$sample)) ) %>% 
    lapply(function(x) {
      sample <- x$sample[1]
      metrics <- read.csv(bc.metrics.file.list[[sample]])
      if (is.null(columns.to.add))
        columns.to.add <- colnames(metrics) %>% .[!grepl("barcode", .)]
      metrics <- metrics[, c("barcode", columns.to.add)]
      x[, columns.to.add] <- NULL
      x <- left_join(x, metrics)
      return(x[, c("barcode", columns.to.add)])
    }) %>% `names<-`(NULL) %>% Reduce(rbind, .)
  
  if (is.null(columns.to.add))
    columns.to.add <- colnames(metrics.df) %>% .[!grepl("barcode", .)]
  for (i in columns.to.add) {
    so[[i]] <- metrics.df[, i]
  }
  return(so)
}

vln.depth <- function(so, bc.metrics.file.list=NULL, plot.dir, metas = c("nCount_RNA", "atac_fragments"), force.re.add = F,
                      split.by = "sample", return.so=F) {
  p.list <- lapply(metas, function(meta) {
    if ((! meta %in% colnames(so@meta.data) || force.re.add == T) && meta != "nCount_RNA") {
      so <- add.metrics.seurat(so, bc.metrics.file.list = bc.metrics.file.list, columns.to.add = meta)
    }
    y.max <- NULL
    if (meta == "atac_fragments")
      y.max <- 50000
    p <- VlnPlot(object = so,features = meta, split.by = split.by, pt.size = 0.01, y.max = y.max)
    return(p)
  }) 
  
  n.splits <- so@meta.data[, split.by] %>% unique() %>% length()
  p <- wrap.plots.fanc(plot.list = p.list, n.col = 1, sub.width = 7.5*n.splits, sub.height = 3, plot.out = paste0(plot.dir, "/depth.png"))
  
  if (return.so == T) {
    return(so)
  }
  return()
}

vln.depth.2 <- function(so, ao = NULL, bc.metrics.file.list = NULL, plot.out, cellranger.metas = NULL, ao.embedding = NULL, embedding = "umap",
                        archr.metas=NULL, metas.plot = NULL,cells = NULL,order,ident = "seurat_clusters",
                        split.by = "sample", return.so=F, sub.width = NULL, numeric.only = T,
                        sub.height = NULL, n.col = 1, max.quantile = 0.999, violin = F, ...) {
  # task 1: add all metas
  if (!is.null(bc.metrics.file.list))
    so <- add.metrics.seurat(so = so, bc.metrics.file.list = BC.METRICS.FILE.LIST,
                           columns.to.add = cellranger.metas)
  if (!is.null(ao)) {
    so <- seurat.add.archr.meta(so = so, ao = ao,metas = archr.metas)
    if (!is.null(ao.embedding)) {
      so <- seurat.add.archr.embed(so = so, ao = ao, ao.embedding = ao.embedding)
      embedding <- paste0("ArchR", ao.embedding)
    }
      
  }
  # n.split <- length(unique(so@meta.data[, split.by]))
  if (!is.null(metas.plot))
    metas <- metas.plot
  else {
    if (violin == T || numeric.only == T) {
      metas <- colnames(so@meta.data)[sapply(as.list(so@meta.data), is.numeric)]
    } else {
      metas <- colnames(so@meta.data)
    }
  }
    
  
  n.splits <- so@meta.data[, split.by] %>% unique() %>% length()
  
  trash <- plot.panel.list(panel.list = metas, b2m = F, obj = so, split.by = split.by, assay = "SCT", order = order,
                            return.list = F, cells = cells, n.split = 1, plot.out = plot.out,
                           n.col = n.col, sub.width = sub.width*n.splits, sub.height = sub.height, reduction = embedding,
                           max.quantile = max.quantile, ident = ident, violin = violin ,...)
  return()
  
}

pos.neg.contrast <- function(so, cluster.ident=NULL, clusters=NULL, gene, archr.name = F) {
  # stop("untested function!!")
  if (!is.null(cluster.ident)) {
    cells <- so@meta.data %>% factor2character.fanc() %>% .[.[,cluster.ident] %in% clusters, ] %>% rownames()
  }
  exp.vec <- so@assays$RNA@counts[gene, cells]
  pos <- names(exp.vec)[exp.vec > 0]
  neg <- names(exp.vec)[exp.vec == 0]
  return.list <- list(pos, neg)
  names(return.list) <- paste0(cluster.ident, "__",paste0(clusters, collpase = "__"),"..", gene, "..", c("pos", "neg") ) %>% 
    sub("^_+", "",.)
  if (archr.name == T)
    return.list <- lapply(return.list, function(x) return(seurat.name.2.archr(so = so, cells = x)))
  return(return.list)
}

peak.watch.core <- function (df, id, stratify.var, stratify.quantiles = NULL, categorical.stratify = F,
                             use.smallest = F, use.smallest.quantile=NULL,
                             use.largest = F, use.largest.quantile = NULL,
                             out.file=NULL, out.bed = NULL,
                             min=NULL, max=NULL, n, seed = 42,
                             return.idx = F, return.idx.list = F, return.df.watch = F, zip = T, threads = 6, 
                             format.fun=NULL, format.fun.params = NULL) {
  if (!is.null(min))
    df <- df[df[,id] > min,]
  if (!is.null(max))
    df <- df[df[,id] < min,]
  
  ###### the code for stratification based peak watching:
  if (!is.null(stratify.quantiles)) {
    left <- right <- stratify.quantiles
    left <- left[-length(left)]
    right <- right[-1]
    stratify.df <- data.frame(left = left, right = right,
                              left.value = quantile(df[, stratify.var], left), 
                              right.value = quantile(df[, stratify.var], right),
                              id = 1:length(left))
    
    stratify.ids <- mclapply(df[, stratify.var], function(x) {
      id <- stratify.df %>% filter(left.value < x, right.value > x) %>% pull(id)
      # print(id)
      if (length(id) > 1) 
        stop("somehow 2 stratify ids returned for 1 line")
      if (length(id) == 0)
        id <- 0
      return(id)
    }, mc.cores = threads, mc.cleanup = T) %>% unlist()
    # browser()
    df <- df[stratify.ids != 0,]
    stratify.ids <- stratify.ids[stratify.ids != 0]
    stratify.vec <- paste0(stratify.df$left, "~", stratify.df$right)
    names(stratify.vec) <- stratify.df$id %>% as.character()
    df$stratify <- stratify.vec[as.character(stratify.ids)]
  } else {
    if (categorical.stratify == T) {
      stratify.ids <- df[, stratify.var]
      df$stratify <- df[, stratify.var]
    } else {
      stratify.ids <- rep(1, nrow(df))
      df$stratify <- "no_stratify"
    }
  }
  # print(stratify.ids)
  df$abs.id <- 1:nrow(df)
  # abs: absolute ID. 
  df.list <- df %>% split(f = stratify.ids)
  abs.id.list <- lapply(df.list, function(df.s) {
    if (n > nrow(df.s))
      n <- nrow(df.s)  
    ################# The sorting step!!!
    df.s <- df.s[order(df.s[, id]),]
    
    set.seed(seed = seed)
    if (use.smallest == T) {
      if (!is.null(use.smallest.quantile)) {
        idx <- which(df.s[, id] < quantile(df.s[, id], use.smallest.quantile))
      } else {
        idx <- 1:n
      }
    }  else {
      if (use.largest == T) {
        if (!is.null(use.largest.quantile)) {
          idx <- which(df.s[, id] > quantile(df.s[, id], use.largest.quantile))
        } else {
          idx <- (nrow(df.s) - n + 1) : nrow(df.s)
        }
      } else {
        idx <- sample(1:nrow(df.s), n, replace = F) %>% sort()
      }
    }
    abs.ids <- df.s[idx, "abs.id"]
    return(abs.ids)
  })
  
  abs.ids <- abs.id.list %>% unlist()  
  df.watch <- df[abs.ids, ] # %>% arrange(stratify, !!as.name(id))
  if (!is.null(out.file) ) {
    utilsFanc::write.zip.fanc(df = df.watch, out.file = out.file, zip = F, col.names = T)
  }
  
  if (!is.null(out.bed)) {
    if (!is.null(format.fun))
      df.bed <- do.call(format.fun, c(list(df.watch = df.watch), format.fun.params)) 
    else
      df.bed <- df.watch
    utilsFanc::write.zip.fanc(df = df.bed, out.file = out.bed, zip = T)
  }
  if (return.idx.list == T)
    return(abs.id.list)
  if (return.idx == T)
    return(abs.ids)
  if (return.df.watch == T)
    return(df.watch)
  return(NULL)
}

a2b.peakwatch.format <- function(df.watch) {
  df <- df.watch %>% mutate(chr = sub(":.+$", "", gene),
                            start = sub("(.+):(.+)-(.+)", "\\2", gene) %>% as.numeric(),
                            end = sub("(.+):(.+)-(.+)", "\\3", gene) %>% as.numeric(),
                            forth = paste0(stratify, "..", pvalue)) %>% 
    select(chr, start, end, forth)
  return(df)
}

pct.detect <- function(so, cells = NULL, genes = NULL, cluster.ident, clusters) {
  if (is.null(cells)) {
    get.cell.list <- list(clusters)
    names(get.cell.list) <- cluster.ident
    cells <- get.cell.names.seurat(so = so, filter.list = get.cell.list)
  }
  raw.mat <- so@assays$RNA@counts[, cells]
  if (!is.null(genes))
    raw.mat <- raw.mat[genes, ]
  pct <- (rowSums(raw.mat > 0)/ncol(raw.mat))
  return(pct)
}

parse.comp <- function(comp, pos) {
  pos1 <- sub(":.+$", "",comp)
  pos2 <- sub("^.+:", "",comp)
  if (pos == 1)
    return(pos1)
  if (pos == 2)
    return(pos2)
}

purity.check.seurat <- function(so, marker.genes, assay,
                                cluster.ident, clusters = NULL,
                                return.df = F, round = NULL) {
  if (is.null(clusters)) {
    clusters <- so@meta.data[, cluster.ident] %>% as.character() %>%  unique()
  }
  
  so <- so[rownames(so) %in% marker.genes, so@meta.data[, cluster.ident] %in% clusters]
  counts.mat <- GetAssayData(object = so, slot = "counts", assay = assay)
  pct.detect <- aggregate.fanc(mat = counts.mat, margin = 2, groupings = so$seurat_clusters,
                               binarize = T, take.mean = T, sort = T)
  # crosscheck:
  # Browse[2]> t <- so@assays$RNA@counts["Gfi1", so$seurat_clusters == "12"] 
  # Browse[2]> mean(as.numeric(t > 0))
  # [1] 0.032
  # matches
  
  norm.mat <- GetAssayData(object = so, slot = "data", assay = assay)
  norm.exp <- aggregate.fanc(mat = norm.mat, margin = 2, groupings = so$seurat_clusters,
                             binarize = F, take.mean = T, sort = T)
  # crosscheck:
  # Browse[2]> so@assays$RNA@data["Gfi1", so$seurat_clusters == "12"] %>% mean()
  # matches
  res <- c(list(pct.detect = pct.detect), sub.f.mat.order(df = pct.detect),
           list(norm.exp = norm.exp), sub.f.mat.order(df = norm.exp))
  if (!is.null(round)) {
    res <- lapply(res, function(x) return(round(x, digits = round)))
  }
  if (return.df == T) {
    res <- lapply(res, function(x) {
      x <- as.data.frame(x)
      x <- utilsFanc::add.column.fanc(df1 = x, df2 = data.frame(gene = rownames(x)), pos = 1)
      return(x)
    })
  }
  return(res)
  
}

sub.f.mat.order <- function(df) {
  df.order <- apply(df, MARGIN = 1, order) %>% t()
  dimnames(df.order) <- dimnames(df)
  df.sorted <- lapply(1:nrow(df.order), function(i) {
    return(df[i, df.order[i, ]])
  }) %>% Reduce(rbind, .)
  colnames(df.sorted) <- NULL
  rownames(df.sorted) <- rownames(df)
  df.sorted.clusters <- lapply(1:nrow(df.order), function(i) {
    return(colnames(df)[df.order[i, ]])
  }) %>% Reduce(rbind, .)
  
  colnames(df.sorted.clusters) <- NULL
  rownames(df.sorted.clusters) <- rownames(df)
  return(list(sorted = df.sorted, sorted.names = df.sorted.clusters, order =  df.order))
}

subset.by.embedding <- function(so, embedding, embedding.df = NULL, range.list) {
  # format of range.list: list(UMAP_1 = c(-10, 5), UMAP_2 = c(5, 7))
  if (is.null(embedding.df)) {
    embedding.df <- soi@reductions[[embedding]]@cell.embeddings %>% as.data.frame()
  }
  cells.pass <- df.range.filter(in.df = embedding.df, range.list = range.list) %>% rownames()

  so <- so[, colnames(so) %in% cells.pass]
  return(so)
}

df.range.filter <- function(in.df, range.list) {
  # format of range.list: list(UMAP_1 = c(-10, 5), UMAP_2 = c(5, 7))
  in.df <- as.data.frame(in.df)
  in.df$id <- 1:nrow(in.df)
  ids <- lapply(names(range.list), function(range.name) {
    range <- range.list[[range.name]]
    id <- in.df$id[in.df[, range.name] > range[1] & in.df[, range.name] < range[2]]
    return(id)
  }) %>% Reduce(intersect, .)

  out.df <- in.df %>% .[.$id %in% ids,] 
  out.df$id <- NULL
  return(out.df)
}
