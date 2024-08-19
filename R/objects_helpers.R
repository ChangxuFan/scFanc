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
                             meta.df = NULL, write.xls = F) {
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
    if (mixedSort == T) {
      df <- df[gtools::mixedorder(df$cluster), ]
    }
    colnames(df) <- c("cluster", "freq")
    df$frac <- df$freq/sum(df$freq)
    sum <- df
    sum$cluster <- NULL
    sum <- lapply(sum, sum) %>% as.data.frame()
    sum$cluster <- "sum"
    df <- rbind(df, sum)
    colnames(df) <- c("cluster", paste0(group.name, "_freq"), paste0(group.name, "_frac"))
    df <- df %>% factor2character.fanc()
    return(df)
  }) %>% Reduce(left_join, .)
  
  if (!is.null(root.name))
    root.name <- paste0(root.name, "_")
  else
    root.name <- ""
  if (!is.null(work.dir)) {
    if (write.xls) {
      xlsx::write.xlsx(x = df, file = work.dir %>% paste0("/",root.name,"cluster_freq.xlsx"),
                       sheetName = "cluster_frequency",
                       append = F, row.names = F, col.names = T)
    } else {
      write.table(df, work.dir %>% paste0("/",root.name,"cluster_freq.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
    }
  }
    
  return(df)
}
rbind.fanc <- function(df1, df2) {
  cols <- colnames(df2) %>% .[!.%in% colnames(df1)]
  df1[,cols] <- NA
  
  cols <- colnames(df1) %>% .[!.%in% colnames(df2)]
  df2[,cols] <- NA
  
  return(rbind(df1, df2))
  
}


int.corr <- function(q.vec, s, meta, out.file = NULL) {
  s <- s@meta.data %>% factor2character.fanc()
  stats.full <- lapply(q.vec, function(q) {
    q <- q@meta.data %>% factor2character.fanc()
    q.sample <- q$sample %>% unique()
    s.sample <- s$sample %>% unique()
    f.sample <- intersect(q.sample, s.sample)
    
    q <- q[q$sample %in% f.sample, ]
    s <- s[s$sample %in% f.sample, ]
    
    q$cells <- rownames(q)
    s$cells <- rownames(s)
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
  stats.full <- cbind(data.frame(samples = rownames(stats.full)), stats.full)
  if (!is.null(out.file)) {
    write.table(stats.full, out.file, quote = F, row.names = F, col.names = T, sep = "\t")
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
    if (!grepl("#", rownames(meta.df)[1]))
      return(paste0(meta.df$sample, "#", sub("_.+$", "", rownames(meta.df))))
    else
      return(rownames(meta.df))
  }
}

add.metrics.seurat <- function(so, bc.metrics.file.list, columns.to.add = NULL) {
  # bc.metrics.file.list should be a named list like this:
  ##
  meta.df <- so@meta.data # %>% factor2character.fanc()
  meta.df$barcode <- rownames(meta.df) %>% sub(".*#", "", .) %>% sub("_.+$", "", .) 
  
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

vln.depth.2 <- function(so = NULL, ao = NULL, bc.metrics.file.list = NULL, plot.out, cellranger.metas = NULL, ao.embedding = NULL, embedding = "umap",
                        archr.metas=NULL, metas.plot = NULL,cells = NULL,order,ident = "seurat_clusters", as.factor = T,
                        split.by = "sample", return.so=F, sub.width = NULL, numeric.only = T,
                        sub.height = NULL, n.col = 1, max.quantile = 0.999, violin = F, ...) {
  # task 1: add all metas
  if (is.null(so)) {
    so <- fake.so.gen(ao = ao, ao.embedding = "UMAP")
  }
  if (!is.null(bc.metrics.file.list))
    so <- add.metrics.seurat(so = so, bc.metrics.file.list = bc.metrics.file.list,
                           columns.to.add = cellranger.metas)
  if (!is.null(ao)) {
    if (!is.null(archr.metas))
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
  if (as.factor == T) {
    so@meta.data[, ident] <- so@meta.data[, ident] %>% factor(., levels = gtools::mixedsort(unique(.)))
  }
  
  trash <- plot.panel.list(panel.list = metas, b2m = F, obj = so, split.by = split.by, assay = "RNA", order = order,
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

pct.detect <- function(so, cells = NULL, genes = NULL,
                       assay = "RNA", cluster.ident, clusters) {
  if (is.null(cells)) {
    get.cell.list <- list(clusters)
    names(get.cell.list) <- cluster.ident
    cells <- get.cell.names.seurat(so = so, filter.list = get.cell.list)
  }
  raw.mat <- so@assays[[assay]]@counts[, cells]
  if (!is.null(genes)) {
    raw.mat <- raw.mat[intersect(rownames(raw.mat), genes), ]
  }
    
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

read.sol <- function(samples, master.dir) {
  sol <- mclapply(samples, function(sample) {
    so <- readRDS(paste0(master.dir, "/", sample ,"/soi.Rds"))
    return(so)
  }, mc.cores = 6)
  
  names(sol) <- samples
  return(sol)
}

metrics.distro.simple <- function(sol = NULL, cell.list, bc.metrics.file.vec, out.dir, 
                                  metrics = c("intergenic_ratio")) {
  # either sol or cell.list must be specified. they must be named by sample names
  if (!is.null(sol)) {
    cell.list <- lapply(sol, function(so) return(colnames(so)))
    names(cell.list) <- names(sol)
    rm(sol)
  }
  samples <- names(cell.list)
  if (is.null(samples)) {
    stop("lists must be named for input")
  }
  bc.stat.list <- utilsFanc::safelapply(samples, function(sample) {
    df <- read.csv(bc.metrics.file.vec[sample]) 
    if ("intergenic_ratio" %in% metrics) {
      df <- df %>% mutate(intergenic_ratio = gex_conf_intergenic_reads/gex_mapped_reads)
    }
    df <- df [, c(metrics, "is_cell" ,"barcode")]
    return(df)
  }, threads = length(samples))
  names(bc.stat.list) <- samples
  lapply(metrics, function(metric) {
    pl <- utilsFanc::safelapply(samples, function(sample) {
      cell.list.2 <- list()
      cell.list.2$all <- bc.stat.list[[sample]]$barcode
      cell.list.2$cell_10x <- bc.stat.list[[sample]] %>% filter(is_cell == 1) %>% pull(barcode)
      cell.list.2$cell_fanc <- cell.list[[sample]]
      pl <- lapply(names(cell.list.2), function(type) {
        p <- ggplot(bc.stat.list[[sample]] %>% filter(barcode %in% cell.list.2[[type]]), aes_string(x = metric)) +
          geom_density() +
          ggtitle(type)
        if (grepl("ratio", metric))
          p <- p + xlim(c(0,1))
        return(p)
      })
      return(pl)
    }, threads = length(samples)) %>% Reduce(c, .)
    wrap.plots.fanc(plot.list = pl, n.col = 3, plot.out = paste0(out.dir, "/", metric, ".png"))
    return()
  })
  return()
}

mt.interRatio.cor <- function(sol, bc.metrics.file.list, out.file) {
  pl <- utilsFanc::safelapply(names(sol), function(sample) {
    so <- sol[[sample]]
    so <- add.metrics.seurat(so = so, bc.metrics.file.list = bc.metrics.file.list,
                             columns.to.add = c("gex_conf_intergenic_reads","gex_mapped_reads"))
    so[["intergenic_ratio"]] <- (so$gex_conf_intergenic_reads/so$gex_mapped_reads) # %>% round(digits = 3)
    p <- xy.plot(df = so@meta.data, x = "percent.mt", y = "intergenic_ratio", 
                 x.limit = c(0, 25), y.limit = c(0,1), add.abline = F) + 
      ggtitle(sample)
    return(p)
  }, threads = length(sol))
  trash <- wrap.plots.fanc(plot.list = pl, plot.out = out.file)
  return()
}

int.corr.plot <- function(soi, sol, samples = NULL,  
                          cluster.ident, clusters = NULL, order = F,
                          plot.dir, root.name = NULL,
                          threads = 10) {
  if (is.character(soi))
    soi <- readRDS(soi)
  if (is.character(sol))
    sol <- readRDS(sol)
  if (!is.null(samples))
    names(sol) <- samples
  
  if (is.null(names(sol)))
    stop("samples should be provided or sol should be named")
  
  int.cluster.ident <- paste0("int.", cluster.ident)
  soi[[int.cluster.ident]] <- soi[[cluster.ident]]
  soi@meta.data <- soi@meta.data %>% factor2character.fanc(re.factor = T)
  soi$cells <- get.cell.names.seurat(so = soi, style = "ArchR")
  
  if (is.null(clusters))
    clusters <- soi@meta.data[, cluster.ident] %>% unique() %>% gtools::mixedsort()
  
  sol <- utilsFanc::safelapply(sol, function(so) {
    so$cells <- get.cell.names.seurat(so = so, style = "ArchR")
    df <- soi@meta.data[, c("cells", int.cluster.ident)]
    so@meta.data <- so@meta.data %>% left_join(df)
    if (nrow(so@meta.data) != ncol(so)) {
      stop("nrow(so@meta.data) != ncol(so)")
    }
    rownames(so@meta.data) <- colnames(so)
    for (cluster in clusters) {
      so[[paste0("IntClus_", cluster)]] <- 0
      so@meta.data[so@meta.data[, int.cluster.ident] == cluster, paste0("IntClus_", cluster)] <- 1
    }
    return(so)
  }, threads = threads)  
  
  for (cluster in clusters) {
    soi[[paste0("IntClus_", cluster)]] <- 0
    soi@meta.data[soi@meta.data[, int.cluster.ident] == cluster, paste0("IntClus_", cluster)] <- 1
  }
  
  plot.out <- "int_corr.png"
  
  if (!is.null(root.name))
    plot.out <- paste0(root.name, "..", plot.out)
  plot.out <- paste0(plot.dir, "/", plot.out)
  
  
  pl <- plot.panel.list.m(panel.list = paste0("IntClus_", clusters), obj = c(list(soi), sol), order = order, 
                          assay = "RNA", plot.out = plot.out, page.limit = 100, raster = T)
  
  # pl <- utilsFanc::safelapply(clusters, function(cluster) {
  #   # sol <- lapply(sol, function(so) {
  #   #   so[["hl"]] <- 0
  #   #   so$hl[so@meta.data[, int.cluster.ident] == cluster] <- 1
  #   #   return(so)
  #   # })
  #   # soi[["hl"]] <- 0
  #   # soi$hl[soi@meta.data[, int.cluster.ident] == cluster] <- 1
  #   
  #   browser()
  # 
  #   return(pl)
  #   
  # }, threads = threads) %>% Reduce(c, .)
  # 
  # 
  # trash <- wrap.plots.fanc(plot.list = pl, n.col = length(sol) + 1, plot.out = plot.out, page.limit = 100)
  return()
}

get.cell.list <- function(obj, is.ao = F, style = "ArchR", cells.include = NULL, n.cells.each = NULL,
                          split.by = NULL, splits = NULL, group.by, groups = NULL,
                          na.rm = T, return.named.vec = F) {
  print ("using get.cell.list. a similar function is archr.get.cells.grid()")
  if (is.ao == F) {
    # obj is so
    df <- obj@meta.data[, c(split.by, group.by), drop = F] %>% factor2character.fanc()
    df$cells <- get.cell.names.seurat(so = obj, style = style)
  } else {
    # obj is ao
    df <- getCellColData(ArchRProj = obj, select = c(split.by, group.by), drop = F) %>% as.data.frame() %>%
      mutate(., cells = rownames(.))
  }
  if (!is.null(cells.include))
    df <- df %>% filter(cells %in% cells.include)
  df <- df %>% na.omit()
  # note this works because df only contains relevant columns (split.by and group.by)
  if (!is.null(groups)) {
    df <- df[df[, group.by] %in% groups,]
  }
  
  if (na.rm == T) {
    df <- df[!is.na(df[, group.by]),]
  }
  df$key <- df[, group.by]
  if (!is.null(split.by)) {
    df$key <- paste0(df$key, "..", df[, split.by])
    if (!is.null(splits)) {
      df <- df[df[, split.by] %in% splits,]
    }
  }
    
  cell.list <- split(df$cells, f = df$key)
  if (!is.null(n.cells.each)) {
    cell.list <- lapply(cell.list, function(x) {
      if (length(x) > n.cells.each)
        x <- sample(x, size = n.cells.each, replace = F)
      return(x)
    })
  }
  if (return.named.vec == T) {
    res <- list.switch.name(cell.list)
  } else {
    res <- cell.list
  }
  return(res)
}

list.switch.name <- function(x) {
  # x must be a named list of vectors. such as a list of cell names. the list
  # is named using cluster names.
  # the goal of the function is to switch cell names and cluster names, generating 
  # a named vector of clusters, with cell names as element names
  if (is.null(names(x)))
    stop("is.null(names(x))")
  res <- lapply(names(x), function(name) {
    y <- rep(name, length(x[[name]]))
    names(y) <- x[[name]]
    return(y)
  }) %>% unlist()
  return(res)
}

binarize.columns <- function(df, cols, include.items = NULL) {
  for (col in cols) {
    items <- df[, col] %>% .[!is.na(.)] %>% unique()
    df[, col][is.na(df[, col])] <- "NA"
    if (!is.null(include.items))
      items <- items %>% .[.%in% include.items]
    for (i in items) {
      b.col <- paste0(col, "..", i)
      df[, b.col] <- 0
      df[ df[, col] == i, b.col] <- 1
      
    }
    df[,col][df[, col] == "NA"] <- NA
  }
  return(df)
}

select.cells.by.rect <- function(embed.df, rect.df) {
  # write rect.df just as if you are writing a polygon df.
  if (nrow(rect.df) != 4) {
    stop("only selection by rect is allowed")
  }
  embed.df <- as.data.frame(embed.df)
  cells <- embed.df %>% .[.[, 1] > min(rect.df$x) & .[, 1] < max(rect.df$x) &
                            .[, 2] > min(rect.df$y) & .[, 2] < max(rect.df$y),] %>% 
    rownames()
  return(cells)
}

select.cells.by.polygon <- function(embed.df, x = 1, y = 2, polygon.df) {
  if (nrow(polygon.df) != 4) {
    stop("only selection by rect is allowed")
  }
  embed.df <- as.data.frame(embed.df)
  a <- sp::point.in.polygon(point.x = embed.df[, x], point.y = embed.df[, y], 
                            pol.x = polygon.df$x, pol.y = polygon.df$y)
  cells <- embed.df[a == 1, ] %>% rownames()
  return(cells)
}


match.bg.cells.by.umap <- function(obj, embedding = NULL, fg.cells, bg.cells = NULL,
                                   n.bg.each, method = "euclidean") {
  if ("Seurat" %in% class(obj)) {
    if (is.null(embedding))
      embedding <- "umap"
    mat <- obj@reductions[[embedding]]@cell.embeddings %>% as.matrix()
  } else if ("ArchRProject" %in% class(obj)) {
    if (is.null(embedding))
      embedding <- "UMAP"
    mat <- obj@embeddings[[embedding]]$df  %>% as.matrix()
  } else {
    mat <- as.matrix(obj)
  }
  
  if (ncol(mat) != 2) {
    stop("ncol(mat) != 2")
  }

  colnames(mat) <- paste0("UMAP", 1:2)
  fg.not.found <- fg.cells %>% .[!. %in% rownames(mat)]
  if (length(fg.not.found) > 0) {
    stop(paste0("some fg.not.found \n",
                paste0(fg.not.found[1:5], collapse = "\n")))
  }
  
  if (!is.null(bg.cells)) {
    bg.not.found <- bg.cells %>% .[!. %in% rownames(mat)]
    if (length(bg.not.found) > 0) {
      stop(paste0("some bg.not.found \n",
                  paste0(bg.not.found[1:5], collapse = "\n")))
    }
    mat <- mat[rownames(mat) %in% c(fg.cells,bg.cells),]
  }
  
  if (length(fg.cells) * n.bg.each > nrow(mat)) {
    stop("length(fg.cells) * n.bg.each > nrow(mat)")
  }
  
  bg.vec <- bg.gen.2(mat = mat, fg.vec = fg.cells, n.bg.each = n.bg.each, method = method, 
                     scale = "none", no.replace = T)
  return(bg.vec)
}

get.metadata.df <- function(obj) {
  if ("Seurat" %in% class(obj)) {
    df <- obj@meta.data
  } else if ("ArchRProject" %in% class(obj)) {
    df <- getCellColData(obj, drop = F) %>% as.data.frame()
  } else {
    df <- obj
  }
  return(df)
}

scale.row.fanc <- function(mat, force.to.mat = F) {
  if (force.to.mat == T)
    mat <- as.matrix(mat)
  res <- mat %>% t () %>% scale() %>% t() 
  return(res)
}


select.by.quantile <- function(x, frac.cutoff, na.method = "keep", larger.than = T) {
  # NA behaviour: keep: NA elements in the input will give NA elements in output.
  # toT: NA changed to TRUE.
  # toF: NA changed to FALSE
  value.cutoff <- quantile(x, frac.cutoff, na.rm = T)
  res <- x > value.cutoff
  if (larger.than == F)
    res <- ! res
  if (na.method == "toT")
    res[is.na(res)] <- T
  if (na.method == "toF")
    res[is.na(res)] <- F
  return(res)
}

normalizePath.partial <- function(path) {
  # main difference from normalizePath:
  # it works as long as the directory exists.
  dir <- dirname(path)
  file <- basename(path)
  dir <- normalizePath(dir)
  return(paste0(dir, "/", file))
}

marker.overlap <- function(so, x, y, x.name = NULL, y.name = NULL,
                           cell.list = NULL, group.by, groups = NULL, split.by = NULL, splits = NULL,
                           assay, slot, exp.threshold = 0,
                           out.dir, root.name) {
  if (is.null(x.name))
    x.name <- paste0(x, collapse = "..")
  if (is.null(y.name))
    y.name <- paste0(y, collapse = "..") 
  # if (is.null(overlap.name))
  #   overlap.name <- "overlap"
  
  if (is.null(cell.list)) {
    cell.list <- get.cell.list(obj = so, is.ao = F, group.by = group.by, groups = groups,
                               split.by = split.by, splits = splits)
  }
  mat <- GetAssayData(object = so, assay = assay, slot = slot)[c(x, y), unlist(cell.list)]
  groupings <- c(rep("x", length(x)), rep("y", length(y)))
  names(groupings) <- c(x, y)
  
  df <- lapply(names(cell.list), function(name) {
    cells <- cell.list[[name]]
    mat.sub <- mat[, cells]
    mat.aggr <- aggregate.fanc(mat = mat.sub, margin = 1, groupings = groupings, sort = T)
    mat.aggr <- mat.aggr > exp.threshold
    
    df <- data.frame(cellGroup = name)
    df[, "total"] <- length(cells)
    df[, x.name] <- mat.aggr["x",] %>% sum()
    df[, y.name] <- mat.aggr["y",] %>% sum()
    df[, "overlap"] <- sum(mat.aggr["x",] & mat.aggr["y",])
    df[, "frac.in.x"] <- df[, "overlap"]/df[, x.name]
    df[, "frac.in.y"] <- df[, "overlap"]/df[, y.name]
    return(df)
  }) %>% Reduce(rbind, .)
  
  system(paste0("mkdir -p ", out.dir))
  write.table(df, file = paste0(out.dir, "/", root.name, ".tsv"))
  return(df)
}

so.label.pos.cells <- function(so, genes, meta.name, pos.name, neg.name = NA,
                               assay, slot,
                               threshold = 0) {
  mat <- GetAssayData(object = so, assay = assay, slot = slot)[genes, , drop = F]
  mat <- mat > threshold
  if (nrow(mat) == 1) {
    b.vec <- mat %>% as.vector()
  } else {
    groupings <- rep(meta.name, length(genes))
    names(groupings) <- genes
    b.vec <- aggregate.fanc(mat = mat, margin = 1, groupings = groupings)
  }
  named.vec <- rep(neg.name, length(b.vec))
  named.vec[b.vec] <- pos.name
  so[[meta.name]] <- named.vec
  return(so)
}

so.balance.cell.number <- function(so, ident = "seurat_clusters", n) {
  cells <- so@meta.data %>% mutate(., cellNames = rownames(.)) %>% 
    dplyr::group_by(!!as.name(ident)) %>% 
    slice_sample(n = n, replace = F) %>% ungroup() %>% pull(cellNames)
  so <- so[, colnames(so) %in% cells]
  return(so)
}

outer.paste <- function(string.list, sep = "..") {
  Reduce(function(x, y) outer(x, y, FUN = function(x, y) paste(x,y, sep = "..")),
         string.list) %>% as.character() %>% return()
}

gene.peak.df.by.dist <- function(genes = NULL, promoters.gr, peaks, distance.1side = 50000) {
  # promoters.gr: something like "~/genomes/mm10/TSS/gencode.vM24.TSS.bed"
  # peaks: chr:start-end
  if ("gene_name" %in% colnames(mcols(promoters.gr))) {
    promoters.gr$gene <- promoters.gr$gene_name
  }
  
  utilsFanc::check.intersect(
    c("gene"), "required fields",
    colnames(mcols(promoters.gr)), "colnames(mcols(promoters.gr))")
  
  if (is.null(genes)) {
    genes <- promoters.gr$gene %>% unique()
  }
  
  utilsFanc::check.intersect(
    genes, "genes", promoters.gr$gene, "promoter.gr$gene"
  )
  
  promoters.gr <- promoters.gr[promoters.gr$gene %in% genes]
  promoters.gr <- resize(promoters.gr, width = 2 * distance.1side, fix = "center")
  peaks <- utilsFanc::loci.2.df(loci.vec = peaks, remove.loci.col = F, return.gr = T)
  df <- df <- plyranges::join_overlap_left(promoters.gr, peaks) %>% 
    utilsFanc::gr2df() %>% dplyr::select(gene, loci) %>% unique() %>% na.omit()
  colnames(df) <- c("gene", "peak")
  return(df)
}

pct.detect.mat <- function(so, genes, group.by, groups = NULL, split.by = NULL, splits = NULL,
                           style = "ArchR", 
                           plot.out = NULL, width = NULL, height = NULL,
                           label.pct = F,
                           hm.colors = c("white", "green3"), hm.values = NULL, show_column_dend = F,
                           cluster_columns = T, clustering_distance_columns = "pearson",
                           cluster_rows = T, clustering_distance_rows = "pearson", ...) {
  utilsFanc::check.intersect(genes, "genes", rownames(so), "rownames(so)")
  cells <- get.cell.list(so, group.by = group.by, groups = groups, split.by = split.by, splits = splits, 
                         is.ao = F, return.named.vec = T, style = style)
  
  pct.mat <- aggregate.fanc(mat = soi@assays$RNA@counts[genes, ],
                            margin = 2, groupings = cells, binarize = T, take.mean = T)
  pct.mat <- as.matrix(pct.mat)
  
  if (!is.null(plot.out)) {
    data.out <- tools::file_path_sans_ext(plot.out) %>% paste0(".csv")
    system(paste0("mkdir -p ", dirname(plot.out)))
    write.csv(pct.mat, data.out, quote = F)
    
    if (is.null(width)) width <- 2 + 0.2 * ncol(pct.mat)
    if (is.null(height)) height <- 1 + 0.2 * nrow(pct.mat)
    
    suppressMessages(extrafont::loadfonts())
    suppressMessages(library(ComplexHeatmap))
    ht_opt$HEATMAP_LEGEND_PADDING <- unit(0.1, "in")
    ht_opt$DENDROGRAM_PADDING <- unit(0, "in")
    
    if (is.null(hm.values)) {
      hm.values <- c(min(pct.mat), max(pct.mat))
    }
    col_fun = circlize::colorRamp2(hm.values, hm.colors)
    
    if (label.pct) {
      cell.mat <- pct.mat %>% round(digits = 3)
      cell_fun <- function(j, i, x, y, width, height, fill) {
        grid.text(cell.mat[i, j], x, y,
                  gp = gpar(fontsize = 6))
      }
    } else {
      cell_fun <- NULL
    }
    
    hm.params <- list(matrix = pct.mat, col = col_fun,
                      show_column_names = T, show_row_names = T,
                      show_column_dend = show_column_dend, 
                      cluster_columns = cluster_columns, 
                      clustering_distance_columns = clustering_distance_columns, 
                      clustering_method_columns = 'complete',
                      cluster_rows = cluster_rows,
                      clustering_distance_rows = clustering_distance_rows,
                      row_names_gp = gpar(fontsize = 6, fontfamily = "Arial"),
                      column_names_gp = gpar(fontsize = 6, fontfamily = "Arial"),
                      column_dend_height = unit(0.1, "in"),
                      column_dend_gp = gpar(lwd = 0.5),
                      show_heatmap_legend = T,
                      heatmap_legend_param = list(
                        title = '', labels_gp = gpar(fontsize = 5), title_gp = gpar(fontsize = 6),
                        legend_height = unit(0.3, "in"), grid_width = unit(0.05, "in"), gap = unit(2, "in")
                      ),
                      cell_fun = cell_fun,
                      ...
                      
    )
    
    hm <- do.call(what = ComplexHeatmap::Heatmap, args = hm.params)
    
    if (!grepl(".pdf$", plot.out)) {
      stop("plot.out must be pdf")
    }

    cairo_pdf(filename = plot.out, width = width, height = height, family = "Arial")
    try({print(draw(hm))})
    dev.off()
    
    ht_opt(RESET = T)
    
    invisible(hm)
    
  }
  return(pct.mat)
}
