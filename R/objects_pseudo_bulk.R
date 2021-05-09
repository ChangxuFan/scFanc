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


so.2.bulk <- function(so, assay, slot, cells=NULL, genes = NULL, ident, sub.idents, 
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
  
  return.list <- list(bulk.mat = bulk.mat, coldata=coldata)
  if (debug == T) {
    return.list <- list(bulk.mat = bulk.mat, coldata=coldata, mat = exp.df)
  }
  return(return.list)
  
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

ao.2.bulk <- function(ao=NULL, cell.list = NULL, gr = NULL, gr.sample.cols =NULL, peakmat = NULL,
                      coldata.columns = NULL, coldata.df = NULL,
                      pseudo, seed = NULL, nprep = 2,
                      cluster.ident, cluster, group.ident, groups=NULL) {
  if (pseudo != T)
    stop("getting coldata from non-pseudo based setting is not yet developed")
  
  if (is.null(gr)) {
    if (is.null(cell.list)) {
      
      if (length(cluster) > 1)
        stop("ao.2.bulk currently only handle 1 cluster")
      cell.list <- archr.get.cells.grid(ao = ao, cluster.ident = cluster.ident, clusters = cluster,
                                        group.ident = group.ident, groups = groups, melt = T, 
                                        seed = seed, pseudo.gen = pseudo, nprep = nprep)
      names(cell.list) <- names(cell.list) %>% sub(paste0(cluster.ident, "_", cluster, ".."), "", .)
      # if (is.null(coldata.columns) && is.null(coldata.df) && pseudo == F) {
      #   getCellColData()
      # }
    }
    
    gr.sample.cols <- names(cell.list)
    if (is.null(peakmat))
      stop("peakmat must be specified when gr is not")
    gr <- vectorization.core(cell.list = cell.list, ao = ao, mat = peakmat, mat.name = "PeakMatrix", 
                             deseq2.norm = F)
  }
  
  if (pseudo == T && is.null(coldata.df)) {
    coldata.df <- data.frame(sample = sub("..prep.*$", "", gr.sample.cols))
    rownames(coldata.df) <- gr.sample.cols
  }
  
  # from gr to DESeq2 format:
  df <- gr %>% `names<-`(NULL) %>% as.data.frame() %>% .[, c("seqnames", "start", "end", gr.sample.cols)] %>% 
    mutate(peak = paste0(seqnames, ":", start, "-", end)) 
  rownames(df) <- df$peak
  mat <- df[, gr.sample.cols] %>% as.matrix()
  return(list(bulk.mat = mat, coldata = coldata.df, gr = gr))
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

ao.2.bulk.list <- function(ao, peakmat, pseudo,
                           cluster.ident, clusters = NULL, group.ident, groups = NULL, threads = 6,
                           design.formula, contrast, work.dir, plot.dir = NULL) {
  if (is.null(clusters)) {
    clusters <- getCellColData(ao, select = cluster.ident)[, cluster.ident] %>% unique() %>% 
      gtools::mixedsort()
  }
  a2b.list <- mclapply(clusters, function(cluster){
    a2b <- ao.2.bulk(ao = ao, peakmat = peakmat, pseudo = pseudo, cluster.ident = cluster.ident,
                     cluster = cluster, group.ident = group.ident, groups = groups)
    a2b <- s2b.deseq(s2b.obj = a2b, design = design.formula, contrast = contrast, try.hm = F,
                     force.hm = F, force = F)
    return(a2b)
  }, mc.cores = threads, mc.cleanup = T)
  # browser()
  names(a2b.list) <- paste0(cluster.ident, "_",clusters)
  saveRDS(a2b.list, paste0(work.dir, "/ao.bulk.list.Rds"))
  return(a2b.list)
}
  
  
s2b.deseq <- function(s2b.obj, design, contrast=NULL, try.hm = T, force.hm = F,
                      force = F, sample.order = NULL, outlist = NULL, single.sample = F,
                      p.adj.cutoff, p.cutoff, title = character(0), pivots = NULL) {
  if (!is.null(s2b.obj$dds) && force == F)
    stop("dds already present in s2b object. use force = T to overide")
  
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
  for (i in 2:3) {
    df[, contrast[i]] <- df[, grepl(paste0(contrast[i], "..prep"), colnames(df))] %>% rowSums()
  }
 
  df$fc <- (df[, contrast[2]]/df[, contrast[3]]) %>% format(digits = 3)
  
  colnames(df) <- colnames(df) %>% paste0("bulkNorm_", .)
  df$gene <- rownames(df)
  rownames(df) <- NULL
  s2b.obj[["bulkNorm"]] <- df
  
  if (single.sample == T)
    return(s2b.obj)
  
  if (!is.null(contrast)) {
    s2b.obj[["res"]] <- DESeq2::results(object = dds, contrast = contrast)
    # res.exp means "result.expanded"
    s2b.obj[["res.exp"]] <- deseq2.add.value(dds=s2b.obj$dds, res=s2b.obj$res, sample.order = sample.order)
    rownames(s2b.obj[["res.exp"]]) <- NULL
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


bulk.list <- function(so, cluster.ident = "seurat_clusters", clusters = NULL,
                      group.by = "sample", groups = NULL, threads = NULL,
                      assay = "RNA", slot = "counts", design.formula, contrast, 
                      work.dir, plot.dir=NULL, p.adj.cutoff = 0.05, p.cutoff =0.05,
                      single.sample = T) {
  system("mkdir -p " %>% paste0(work.dir, " "))
  meta.df <- so@meta.data %>% factor2character.fanc()
  if (is.null(clusters)) {
    clusters <- meta.df[, cluster.ident] %>% unique() %>% gtools::mixedsort()
  }
  
  if (is.null(threads)) {
    threads <- length(clusters)
  }
  threads <- min(12, threads)
  
  bulk.all.list <- mclapply(clusters, function(cluster) {
    bulk.all <- so.2.bulk(so, assay = assay, slot= slot, ident = cluster.ident, sub.idents = cluster,
                          group.by = group.by, groups = groups, coldata.columns = NULL)
    bulk.all <- s2b.deseq(s2b.obj = bulk.all, design = design.formula , 
                          contrast = contrast,
                          force = T, sample.order = NULL, 
                          outlist = list(plot.dir = plot.dir, root.name.external = "cluster_" %>% paste0(i, "_")), 
                          p.adj.cutoff = p.adj.cutoff, p.cutoff = p.cutoff,
                          title = paste0("cluster_", cluster), single.sample = single.sample,
                          pivots = NULL)
  }, mc.cores = threads, mc.cleanup = T)
  names(bulk.all.list) <- paste0(cluster.ident, "_", clusters)
  saveRDS(bulk.all.list, paste0(work.dir, "/bulk.list.Rds"))
  return(bulk.all.list)
}

aggregate.fanc <- function(mat, margin, groupings, binarize = F, take.mean = F, sort = F) {
  ### groupings: must be a named vector. the names are the colnames or rownames. mat will be 
  #rearranged to match the sequence of groupings
  if (is.factor(groupings)) {
    names <- names(groupings)
    groupings <- as.character(groupings)
    names(groupings) <- names
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
  
  mat <- mat[names(groupings), ]
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
  return(aggr)
}