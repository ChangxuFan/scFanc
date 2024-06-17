find.elbow.fanc <- function(df, spar=0.6, mask.rank.pct = c(80 ,100)) {
  # the 2 columns of df: x and y
  df$x.ori <- df$x
  df$y.ori <- df$y
  df <- df %>%  mutate(x= x/max(x) * 100, y = y/max(y) * 100)
  smooth <- smooth.spline(x = df$x, y = df$y, spar = spar)
  df$sd <- predict(smooth, x = df$x, deriv = 2)$y
  max.sd <- df[df$x %between% mask.rank.pct, "sd"]  %>% max()
  elbow.x <- df[df$sd == max.sd, "x.ori"][1]
  elbow.y <- df[df$sd == max.sd, "y.ori"][1]
  return(list(x=elbow.x, y = elbow.y))
  
}

sc.qc.construct.so <- function(sample.info, project.name, mt.pattern,
                               construct.FUN = NULL, construct.arg.list=NULL, 
                               min.gene=200, min.cell.w.gene=5, type = "Gene Expression",
                               threads = NULL) {
  if (is.character(sample.info))
    sample.info <- read.table(sample.info, header = T, as.is = T, sep = "\t")
  
  if (is.null(threads)) {
    threads <- nrow(sample.info) +1
  }
  sol <- split(sample.info, f=sample.info$sample %>% factor(.,levels=.)) %>%
    utilsFanc::safelapply(function(x) {
      add.assays <- list()
      if (!is.null(construct.FUN))
        so <- do.call(construct.FUN, c(list(x$dir), project.name = project.name,construct.arg.list))
      else {
        so <- Read10X(data.dir = x$dir)
        if (is.list(so)) {
          if ("Antibody Capture" %in% names(so)) {
            add.assays$ADT <- CreateAssayObject(counts = so[["Antibody Capture"]])
          }
          so <- so[[type]]
        }
        so <- CreateSeuratObject(so, project = project.name)
      }
      so <- so[so@assays$RNA@counts %>% as.matrix() %>% apply(1, function(x) return(sum(x>0) >= min.cell.w.gene)),]
      # have to add additional assays after filtering genes
      # when you filter genes, all the other assays will magically disappear...
      if (length(add.assays) > 0) {
        for (assay in names(add.assays)) {
          so[[assay]] <- add.assays[[assay]]
        }
      }
      so <- so[, so$nFeature_RNA >= min.gene]
      so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = mt.pattern)
      metadata <- colnames(x)[!colnames(x) %in% c("dir")]
      for (m in metadata) {
        so[[m]] <- x[1,m]
      }
      return(so)
    }, threads = threads)
  return(sol)
}

csv2seurat <- function(csv, project.name) {
  # this is adopted from vivier's dataset at ~/4dn/nk/fanc/single_cell_vivier_3_29_20.R, and hence the "nk.sc" nomenclature.
  nk.sc <- read.csv(file = csv,
                    header = T)
  row.names(nk.sc) <- nk.sc[,1]
  nk.sc <- nk.sc[,-1]
  nk.sc <- as.matrix(nk.sc)
  nk.scs <- CreateSeuratObject(counts = nk.sc, project= project.name)
  return(nk.scs)
}

tsv2seurat <- function(tsv, project.name, amit.clean = F) {
  df <- read.table(tsv, header = T, row.names = 1, sep = "\t")
  mat <- as.matrix(df)
  if (amit.clean) {
    # developed to deal with https://pubmed.ncbi.nlm.nih.gov/29915358/
    print("Removing gene names with ';' or 'Rik'")
    mat <- mat[!grepl(";|Rik", rownames(mat)),]
  }
  
  so <- CreateSeuratObject(counts = mat, project= project.name)
  return(so)
}

split.by.metadata <- function(so, sample.info,reduction, outdir=NULL, metas.include=NULL, ...) {
  if (is.character(sample.info))
    sample.info <- read.table(sample.info, as.is =T, header = T)
  metas <- colnames(sample.info)[colnames(sample.info) != "dir"]
  if (!is.null(metas.include))
    metas <- metas[metas %in% metas.include]
  n.plots <- length(metas)
  if (n.plots < 1)
    stop("nothing to plot")
  plot.list <- lapply(metas, function (meta) {
    p.s <- DimPlot(so, reduction = reduction, label = F, group.by = meta, split.by = meta, ...)
    
    n.sub.plots <- sample.info[,meta] %>% unique() %>% length()
    print(n.sub.plots)
    p.g <- lapply(sample.info[,meta] %>% unique(), function(x) {
      DimPlot(so, reduction = reduction, label = T, group.by = meta, order = x,repel = T, label.size = 5,...)
    })
    p.g <- cowplot::plot_grid(plotlist = p.g, nrow = 1)
    if (!is.null(outdir)) {
      ggsave(paste0(outdir, "/", reduction, "_split_", meta, ".png"), p.s, device = "png", width = 5*(n.sub.plots), height = 5,
             dpi = 100, units = "in", limitsize = F)
      ggsave(paste0(outdir, "/", reduction, "_group_", meta, ".png"), p.g, device = "png", width = 7*n.sub.plots, height = 5,
             dpi = 100, units = "in", limitsize = F)
    }
    return(list(ps=p.s, pg=p.g))
  })
  
  names(plot.list) <- metas
  return(plot.list)
}

split.by.metadata.2 <- function(so, reduction, outdir=NULL, metas, ...) {
  print("untested function!!")
  missing.metas <- ! metas %in% colnames(so@meta.data)
  if (sum(missing.metas) > 0)
    stop("these provided metas are not found in object: " %>% paste0(paste0(metas[missing.metas], collapse = ";")))
  plot.list <- lapply(metas, function (meta) {
    p.s <- DimPlot(so, reduction = reduction, label = T, group.by = meta, split.by = meta, ...)
    
    n.sub.plots <- sample.info[,meta] %>% unique() %>% length()
    print(n.sub.plots)
    p.g <- lapply(sample.info[,meta] %>% unique(), function(x) {
      DimPlot(so, reduction = reduction, label = T, group.by = meta, order = x,repel = T, label.size = 5,...)
    })
    p.g <- cowplot::plot_grid(plotlist = p.g, nrow = 1)
    if (!is.null(outdir)) {
      ggsave(paste0(outdir, "/", reduction, "_split_", meta, ".png"), p.s, device = "png", width = 5*(n.sub.plots), height = 5,
             dpi = 100, units = "in", limitsize = F)
      ggsave(paste0(outdir, "/", reduction, "_group_", meta, ".png"), p.g, device = "png", width = 7*n.sub.plots, height = 5,
             dpi = 100, units = "in", limitsize = F)
    }
    return(list(ps=p.s, pg=p.g))
  })
  
  names(plot.list) <- metas
  return(plot.list)
}

meta.cell.count <- function(so=NULL, so.list=NULL,x.meta, x.meta.order=NULL, group.meta,
                            group.meta.order=NULL,outdir=NULL, root.name=NULL, dodge=T) {
  if (!is.null(so.list)) {
    meta.df <- lapply(so.list, function(x) {
      meta.df.sub <- x$so@meta.data %>% factor2character.fanc()
      # meta.df.sub[, x.meta] <- paste0(x$prefix, meta.df.sub[, x.meta])
      meta.df.sub[, group.meta] <- paste0(x$prefix, meta.df.sub[, group.meta])
      return(meta.df.sub)
    }) %>% Reduce(rbind, .)
  } 
  else 
    meta.df <- so@meta.data %>% factor2character.fanc()

  if (is.null(root.name ))
    root.name <- meta.df$orig.ident %>% as.character() %>% unique() %>% paste0(collapse = "_")
  if (grepl("\\|",x.meta.order) %>% sum() > 0 ) {
    meta.df <- comp.meta(meta.df = meta.df, meta.name = x.meta, meta.sets = x.meta.order, out.comp.name = x.meta)
  }
  if (grepl("\\|",group.meta.order) %>% sum() > 0 ) {
    meta.df <- comp.meta(meta.df = meta.df, meta.name = group.meta, meta.sets = group.meta.order,
                    out.comp.name = group.meta)
  }
  
  
  df <- meta.df[, c(x.meta, group.meta)]
  colnames(df) <- c("x", "group")

  if (!is.null(x.meta.order)) {
    df <- df %>% filter(x %in% x.meta.order)
    df$x <- factor(df$x, levels = x.meta.order)
  }

  if (!is.null(group.meta.order)) {
    df <- df %>% filter(group %in% group.meta.order)
    df$group <- factor(df$group, levels = group.meta.order)
  }
  
  df <- df %>% group_by(x, group) %>% summarise(n.cells = n()) %>% ungroup() 
  if (dodge == T) {
    df <- df %>% 
      group_by(group) %>% mutate(pct.cells = n.cells/sum(n.cells))  %>% ungroup() %>% 
      group_by(x) %>% mutate(pct.cells.standardized = pct.cells/max(pct.cells)) %>% ungroup()
  } else {
    df <- df %>% 
      group_by(x) %>% mutate(pct.cells = n.cells/sum(n.cells)) %>% ungroup()
  }
  
  df <- df %>%
    reshape2::melt(id.vars = c("x", "group"))
  
  p <- df %>% ggplot(aes(x=x,
                         y=value,
                         fill = group ))
  if (dodge == T)
    p <- p + geom_bar(position=position_dodge(), stat="identity")
  else 
    p <- p + geom_bar(stat="identity")
  p <- p + 
    xlab(x.meta) +
    # scale_fill_manual(name=group.meta, limits = df$group %>% unique(), 
    #                   values = gg_color_hue( df$group %>% unique() %>% as.integer(), up.limit = )) +
    scale_fill_discrete(name=group.meta) +
    facet_grid(variable ~ ., scales='free')
  if (!is.null(outdir))
    ggsave(paste0(outdir, "/", root.name, "_n_cells_", x.meta, "_by_", group.meta, ".png"), p,
           units = "in", device = "png", width = df$x %>% unique() %>% length() %>% min(20) %>% max(5),
           height = 5, dpi = 200, limitsize = F)
  return(p)  
}





sc.rank.qc <- function(sol,feature,ymax ,plot.dir=NULL, project.name=NULL) {
  if (is.null(project.name))
    project.name <- sol[[1]]$orig.ident[1] %>% as.character()
  
  
  n.feature.plots <- lapply(seq_along(sol), function(i) {
    colnames(sol[[i]]@meta.data)[colnames(sol[[i]]@meta.data) == feature] <- "feature_to_rank"
    n.feature.df <- sol[[i]]@meta.data %>% mutate(rank.Feature = rank(feature_to_rank)) %>% select(rank.Feature, feature_to_rank) %>%
      arrange(rank.Feature) 
    p <- n.feature.df %>%
      ggplot(aes(x=rank.Feature, y=(feature_to_rank))) +
      geom_point(size = 0.1) +
      ylim(0, ymax) +
      ggtitle(sol[[i]]$sample[1] %>% as.character()) +
      xlab(paste0("rank_", feature)) +
      ylab(feature)
    elbow <- sol[[i]]@misc[[paste0("elbow_", feature)]]
    if (!is.null(elbow))
      p <- p+geom_hline(yintercept = elbow$y, color = "red")
    
    return(p)
  })
  plot.out <- NULL
  if(!is.null(plot.dir))
    plot.out <- paste0(plot.dir,"/",project.name,  "_", feature,"_rank_distro.png")
  
  p.all <- wrap.plots.fanc(plot.list = n.feature.plots, n.split = 1, 
                           plot.out = plot.out)
  return(p.all)
}




sc.qc.metadata.find.elbow.core <- function(sol, elbow.list=NULL, meta.term, override = T) {
  meta.sum.df <- list()
  meta.rank <- paste0("rank_", meta.term)
  meta.elbow <- paste0("elbow_", meta.term)
  for (i in seq_along(sol) ) {
    # meta.sum.df[[i]] <- sol[[i]]@meta.data %>% mutate(rank.Feature = rank(percent.mt))
    if (!is.null(sol[[i]]@misc[[meta.elbow]]) && override == F)
      next
    
    meta.sum.df[[i]] <- sol[[i]]@meta.data
    meta.sum.df[[i]][, meta.rank] <- rank(meta.sum.df[[i]][, meta.term])
    meta.sum.df[[i]] <- meta.sum.df[[i]][, c(meta.rank, meta.term)]
    meta.sum.df[[i]] <- meta.sum.df[[i]][order(meta.sum.df[[i]][, meta.rank]),]
    colnames(meta.sum.df[[i]]) <- c("x", "y")
    
    if (is.null(elbow.list[[sol[[i]]@meta.data$sample[1]]])) {
      sol[[i]]@misc[[meta.elbow]] <- meta.sum.df[[i]] %>%
        find.elbow.fanc()
    } else {
      sol[[i]]@misc[[meta.elbow]] <- list(y=elbow.list[[sol[[i]]@meta.data$sample[1]]])
    }
  }
  return(sol)
  
}

sc.qc.elbow.filter <- function(sol, metas, project.name, take.lower = T) {
  for (meta in metas) {
    sol <- sc.qc.metadata.find.elbow.core(sol, meta.term = meta, override = F)
  }

  
  sol <- lapply(sol, function(x) {
    for(meta in metas) {
      if (take.lower == T) {
        x <- x[, x@meta.data[, meta] < x@misc[[paste0("elbow_", meta)]]$y]
      } else {
        x <- x[, x@meta.data[, meta] >= x@misc[[paste0("elbow_", meta)]]$y]
      }
        
    }
      
    x[["orig.ident"]] <- project.name
    return(x)
  })
  return(sol)
}





sct.cycle.regress <- function(x, n.var.feature = 3000) {
  s.genes <- cc.genes$s.genes %>% human.mouse.convert(human.vec = .)
  g2m.genes <- cc.genes$g2m.genes %>% human.mouse.convert(human.vec = .)
  DefaultAssay(x) <- "SCT"
  x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  x <- SCTransform(x, assay = "RNA", do.correct.umi = T,
                   vars.to.regress = c("S.Score", "G2M.Score"),
                   verbose = T, variable.features.n = n.var.feature)
  return(x)
}



sct.int.prep <- function(sol, n.var.feature, regress.cycle = F, out.sol.Rds = NULL, sct=T) {
  if (is.character(sol))
    sol <- readRDS(sol)
  if (sct == T)
  sol <- mclapply(sol, function(x) SCTransform(x, assay = "RNA", do.correct.umi = T,
                                               verbose = T, variable.features.n = n.var.feature),
                  mc.cores = length(sol) + 1)
  if (regress.cycle == T) {
    s.genes <- cc.genes$s.genes %>% human.mouse.convert(human.vec = .)
    g2m.genes <- cc.genes$g2m.genes %>% human.mouse.convert(human.vec = .)
    
    sol <- mclapply(sol, function(x) {
      DefaultAssay(x) <- "SCT"
      x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
      x <- SCTransform(x, assay = "RNA", do.correct.umi = T,
                       vars.to.regress = c("S.Score", "G2M.Score"),
                       verbose = T, variable.features.n = n.var.feature)
      return(x)
    }, mc.cores = length(sol) + 1)
  }
  if (!is.null(out.sol.Rds))
    saveRDS(sol, out.sol.Rds)
  return(sol)
}

# t <- sct.int.prep(sol, 3000, T)

sct.int.pipe <- function(sol, n.int.features, n.threads, outfile= NULL, reference = NULL,
                         runPCA = F, reduction = "cca", dims = 1:30, k.anchor = 5) {
  if (is.character(sol))
    sol <- readRDS(sol)
  
  if (n.threads > 1) {
    plan("multiprocess", workers = n.threads)
  }

  options(future.globals.maxSize = 12000 * 1024^2)
  
  
  int.features <- SelectIntegrationFeatures(object.list = sol, nfeatures = n.int.features, 
                                            assay = rep("SCT", length(sol)))

  sol <- PrepSCTIntegration(object.list = sol, anchor.features = int.features, assay = rep("SCT", length(sol)),  
                            verbose = T)
  
  if (runPCA == T) {
    sol <- utilsFanc::safelapply(sol, FUN = RunPCA, fetures = int.features, assay = "SCT", 
                                 threads = n.threads)
  }
  
  anchors <- FindIntegrationAnchors(object.list = sol, normalization.method = "SCT",
                                    assay = rep("SCT", length(sol)),  reference = reference,
                                    anchor.features = int.features, verbose = T,
                                    reduction = reduction, dims = dims, k.anchor = k.anchor)
  soi<- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = dims,
                      verbose = T)
  if (!is.null(outfile)) {
    system(paste0("mkdir -p ", dirname(outfile)))
    saveRDS(soi, outfile)
  }
    
  plan("sequential")
  return(soi)
}

cluster.pipe <- function(soi, assay, assay.pre = "RNA", is.sct.int = T,
                         vars.to.regress = NULL,
                         pc.run = 60, pc.dim = 1:30, 
                         run.pca = T,
                         cluster.resolution = 0.7, cluster.seed = 0, umap.seed = 42,
                         work.dir, plot.dir,
                         project.name ,
                         metas.include = NULL,
                         add.all.normalization = F, 
                         sample.tsv, plot.only = F, split.by.metadata = T, 
                         save.rds = T, plot.common.markers = T, plot.metrics = T,
                         split.by = NULL, split.order = NULL,
                         to.human = F,
                         bc.metrics.file.list = NULL, 
                         do.sct = F, sct.n.vars=3000, hm.dims = 11:50,
                         re.umap.with.seed = NULL) {

  system(paste0("mkdir -p ", work.dir, " ", plot.dir))
  if (is.character(soi))
    soi <- readRDS(paste0(work.dir, "/soi.Rds"))
  
  if (plot.only == F) {
    if (do.sct == T) {
        soi <- SCTransform(soi, assay = assay.pre, do.correct.umi = T,
                           verbose = T, variable.features.n = sct.n.vars)
    } 
    
    if (assay == "integrated") {
      if (is.sct.int == F) {
        soi <- ScaleData(object = soi, assay = "integrated")
      } else {
        print("is.sct.int set to TRUE. Skipping ScaleData()")
      }
        
    } else if (!grepl("^SCT", assay)) {
      soi <- NormalizeData(object = soi, normalization.method = "LogNormalize", assay = assay)
      soi <- FindVariableFeatures(soi, assay = assay)
      soi <- ScaleData(object = soi, assay = assay, vars.to.regress = vars.to.regress)
    }
    
    DefaultAssay(soi) <- assay
    
    if (run.pca == T)
      soi <- RunPCA(soi, npcs = pc.run, assay = assay)
    
    soi <- RunUMAP(soi, reduction = "pca", dims = pc.dim, assay = assay, seed.use = umap.seed)
    soi <- FindNeighbors(soi, reduction = "pca", dims = pc.dim)
    soi <- FindClusters(soi, resolution = cluster.resolution, random.seed = cluster.seed)
    if (assay != assay.pre) {
      DefaultAssay(soi) <- assay.pre
      soi <- NormalizeData(object = soi, normalization.method = "LogNormalize", assay = assay.pre)
      DefaultAssay(soi) <- assay
    }
    if (add.all.normalization == T) {
      soi <- add.all.normalization(so = soi, methods = "LogNormalize", scale.factor = 10000)
      DefaultAssay(soi) <- assay
    }
    
    DefaultAssay(soi) <- assay
    
    # if (assay != "RNA")
    #   DefaultAssay(soi) <- "SCT"
    
    soi <- RenameCells(soi, new.names = get.cell.names.seurat(soi, style = "ArchR"))
    if (save.rds == T)
      saveRDS(soi, paste0(work.dir, "/soi.Rds"), compress = F)
    
    trash <- count.by.cluster(soi, work.dir, group.by = "sample")
  }
  
  
  
  # try({
  #   png(paste0(plot.dir, "/elbow.png"))
  #   print(ElbowPlot(soi, ndims = 60, reduction = "pca"))
  #   dev.off()
  # })
  # try({
  #   png(paste0(plot.dir, "/pchm.png"), units = "in", width = 40, height = 40, res=200)
  #   DimHeatmap(soi, dims = hm.dims, assays = assay, cells = 500, balanced = TRUE, combine = T, fast=T, ncol = 3)
  #   dev.off()
  # })

  if (!is.null(re.umap.with.seed)) {
    soi <- RunUMAP(soi, reduction = "pca", dims = pc.dim, assay = assay, seed.use = re.umap.with.seed)
    if (save.rds == T)
    saveRDS(soi, paste0(work.dir, "/soi.Rds"))
  }

  trash <- DimPlot(soi, reduction = "umap", label = T, pt.size = 0.5, label.size = 10, shuffle = T) +
    ggsave(paste0(plot.dir, "/umap_cluster.png"), device = "png", width = 10, height = 10, dpi= 200, units = "in")
  
  #soi[["age_skewing"]] <- paste0(soi[["age"]]$age, "_", soi[["skewing"]]$skewing)
  if (split.by.metadata == T && !is.null(metas.include)) {
    # trash <- split.by.metadata(so = soi, sample.info = sample.tsv, reduction = "umap", outdir = plot.dir,
    #                            metas.include = metas.include)
    plot.panel.list(panel.list = metas.include, obj = soi, binarize.panel = T, 
                    order = T, assay = assay, 
                    raster = T, publication = F,
                    root.name = paste0(plot.dir, "/umap_meta_hl"),
                    invisible = T)
  }
  
  
  if (plot.metrics ==T ) {
    try(trash <- vln.depth.2(so = soi, ao = NULL, bc.metrics.file.list = bc.metrics.file.list,metas.plot = NULL,
                         plot.out = paste0(plot.dir,"/metrics.pdf"), sub.width = 7.5, violin = T, page.limit = 4,
                         max.quantile = 0.99, n.col = 4, numeric.only = T, sub.height = 3))
    # 
    # try(trash <- vln.depth(so = soi, bc.metrics.file.list = bc.metrics.file.list, plot.dir = plot.dir, split.by = "sample",
    #                      return.so = F, force.re.add = T))
  }
  
  if (plot.common.markers == T) {
    trash <- hmtp.common.markers(obj = soi, plot.dir = plot.dir, to.human = to.human, threads = NULL, 
                                 split.by = split.by, split.order = split.order)
  }

  return(soi)
}

titrate.umap <- function(so, assay, reduction = "pca", dims.list, seeds, out.dir,
                         split.by = NULL,
                         threads = 1,
                         sub.width = 6, sub.height = 5, ...) {
  utilsFanc::safelapply(dims.list, function(dims) {
    pl <- utilsFanc::safelapply(seeds, function(seed) {
      so <- RunUMAP(so, reduction = reduction, dims = dims, assay = assay, 
                    seed.use = seed)
      p <- DimPlot(so, label = T, split.by = split.by, pt.size = 0.01)
      p <- p + ggtitle("seed: " %>% paste0(seed))
      return(p)
    }, threads = threads)
    plot.out <- paste0(out.dir, "/dims_",dims[length(dims)], ".png" )
    trash <- wrap.plots.fanc(plot.list = pl, plot.out = plot.out, 
                             sub.height = sub.height, sub.width = sub.width, ...)
    return()
  }, threads = 1)
  return()
}

# ooo means one on one
ooo.enrich <- function(so=NULL, so.list = NULL, comparison.df=NULL, debug =F) {
  if (is.null(comparison.df)) {
    comparison.df <- data.frame(sample.1 = c("PyMT_26w", "PyMT_25w", "PyMT_9w"), 
                                sample.2 = c("WT_24w", "PyMT_8w", "WT_9w"))
  }
  
  samples <- c(comparison.df$sample.1, comparison.df$sample.2) %>% unique()
  
  if (!is.null(so.list)) {
    meta.df <- lapply(so.list, function(x) {
      meta.df.sub <- x$so@meta.data %>% factor2character.fanc()
      meta.df.sub$seurat_clusters <- paste0(x$prefix, meta.df.sub$seurat_clusters)
      return(meta.df.sub)
    }) %>% Reduce(rbind, .)
  } else
    meta.df <- so@meta.data %>% factor2character.fanc()
  
  cell.num.mat <- so@meta.data %>% group_by(sample, seurat_clusters) %>% 
    summarise(n.cell = n()) %>% ungroup() %>% rename(cl = seurat_clusters) %>% 
    # mutate(cl = as.character(cl)) %>%
    reshape2::acast(formula = cl ~ sample)
  # return(cell.num.mat)
  df <- lapply(so$seurat_clusters %>% unique() %>% as.character(), function(cl1) {
    lapply(so$seurat_clusters %>% unique() %>% as.character(), function(cl2) {
      df <- data.frame(cl1 = cl1, cl2=cl2)
      # print(df)
      for (sample in samples) {
        df[1,sample] <- cell.num.mat[cl1,sample]/cell.num.mat[cl2, sample]
        
      }
      
      for (i in 1:nrow(comparison.df)) {
        df[1, paste0(comparison.df[i, "sample.1"], "/",comparison.df[i, "sample.2"] )] <- 
          df[1, comparison.df[i, "sample.1"]]/df[1, comparison.df[i, "sample.2"]]
      }
      
      if (debug == T) {
        for (sample in samples) {
          df[1,sample %>% paste0("_cl1_raw")] <- cell.num.mat[cl1,sample]
          df[1,sample %>% paste0("_cl2_raw")] <- cell.num.mat[cl2,sample]
        }
      }

      return(df)
    }) %>% Reduce(rbind,.)
  }) %>% Reduce(rbind, .) 
  return(list(mat = cell.num.mat, enrich = df))
}

# t <- ooo.enrich(prcr.q) 
# View(t$enrich)

comp.meta <- function(meta.df, meta.name,  meta.sets, out.comp.name) {
  # called by meta.cell.count and ooo.enrich. used to treat several metas as a set.
  if (meta.sets %>% strsplit("\\|") %>% unlist() %>% duplicated() %>% sum() > 0) {
    stop("there are overlapping metas in meta.sets")
  }
  
  for (set.name in meta.sets) {
    metas.in.set <- set.name %>% strsplit("\\|") %>% unlist()
    meta.df[meta.df[, meta.name] %in% metas.in.set, out.comp.name] <- set.name
  }
  
  return(meta.df)
}

factor2character.fanc <- function(df, re.factor = F) {
  df[] <- lapply(df, function(x)  {
    if (is.factor(x)) {
      x <- as.character(x)
      if (re.factor == T) {
        x <- factor(x, levels = unique(x) %>% gtools::mixedsort())
      }
    } 
    return(x)
  })
  return(df)
}

meta.cell.count.core <- function(so=NULL, so.list=NULL,x.meta, x.meta.order=NULL, group.meta,
                                 group.meta.order=NULL, prefix.on.x = F, prefix.on.group = T) { 
  if (!is.null(so.list)) {
    meta.df <- lapply(so.list, function(x) {
      meta.df.sub <- x$so@meta.data %>% factor2character.fanc()
      if (prefix.on.x == T)
        meta.df.sub[, x.meta] <- paste0(x$prefix, meta.df.sub[, x.meta])
      if (prefix.on.group == T)
        meta.df.sub[, group.meta] <- paste0(x$prefix, meta.df.sub[, group.meta])
      return(meta.df.sub)
    }) %>% Reduce(rbind, .)
  } 
  else 
    meta.df <- so@meta.data %>% factor2character.fanc()
  
  # if (is.null(root.name ))
  #   root.name <- meta.df$orig.ident %>% as.character() %>% unique() %>% paste0(collapse = "_")
  if (grepl("\\|",x.meta.order) %>% sum() > 0 ) {
    meta.df <- comp.meta(meta.df = meta.df, meta.name = x.meta, meta.sets = x.meta.order, out.comp.name = x.meta)
  }
  if (grepl("\\|",group.meta.order) %>% sum() > 0 ) {
    meta.df <- comp.meta(meta.df = meta.df, meta.name = group.meta, meta.sets = group.meta.order,
                         out.comp.name = group.meta)
  }
  
  
  df <- meta.df[, c(x.meta, group.meta)]
  colnames(df) <- c("x", "group")
  
  if (!is.null(x.meta.order)) {
    df <- df %>% filter(x %in% x.meta.order)
    df$x <- factor(df$x, levels = x.meta.order)
  }
  
  if (!is.null(group.meta.order)) {
    df <- df %>% filter(group %in% group.meta.order)
    df$group <- factor(df$group, levels = group.meta.order)
  }
  
  df <- df %>% group_by(x, group) %>% summarise(n.cells = n()) %>% ungroup() 
  return(df)
}

meta.cell.corr <- function(so=NULL, so.list = NULL, point.meta, point.meta.include=NULL, 
                           point.color.map = NULL,
                           axes.meta, axes.meta.include, root.name, out.dir,
                           sub.width = 5, sub.height = 5) {
  # the format of point.color.map: 2 column dataframe. first col: point.meta; 2nd col: color
  
  df <- meta.cell.count.core(so = so, so.list = so.list, prefix.on.x = T, prefix.on.group = F,
                             x.meta = point.meta, x.meta.order = point.meta.include,
                             group.meta = axes.meta, group.meta.order = axes.meta.include %>% unique()) %>% 
    reshape2::acast(formula = x~group) %>% as.data.frame() 
  df$point.meta <- rownames(df)
  plots <- lapply(axes.meta.include, function(x) {
    lapply(axes.meta.include, function(y) {
      if (length(axes.meta.include) ==2) {
        if (x == y || y == axes.meta.include[1])
          return(NULL)
      }
      
      df.sub <- df[, c(x, y, "point.meta")]
      colnames(df.sub) <- c("x", "y", "point.meta")
      df.sub$x <- df.sub$x/max(df.sub$x)
      df.sub$y <- df.sub$y/max(df.sub$y)
      if (!is.null(point.color.map)) {
        df.sub <- left_join(df.sub, point.color.map)
        df.sub[is.na(df.sub)] <- "other"
      }
      
      p <- df.sub %>% ggplot(., aes(x=x, y=y, label = point.meta))
      
      if (!is.null(point.color.map))
        p <- p + geom_point(aes(color = color))
      else 
        p <- p+ geom_point()
      
      p <- p+
        geom_text(nudge_x = 0.02) +
        geom_smooth(method = "lm", se = F)+
        xlab(x) +
        ylab(y) +
        xlim(0,1.2) +
        ylim(0,1.2) +
        theme_light()
      return(p)
    }) %>% return()
  }) %>% Reduce(c,.) 
  # return(plots)
  # plots
  plots <- plots[sapply(plots, function(x) return(!is.null(x)))]
  p <- wrap.plots.fanc(plot.list = plots, col.row.ratio = 1, dpi = 200, sub.width = sub.width, sub.height = sub.height,
                  plot.out = paste0(out.dir, "/", root.name,
                                    "_cell_count_corr_", point.meta, "_by_", axes.meta, ".png"))
  return(p)

}

# trash <- meta.cell.corr(so = prcr.q, point.meta = "seurat_clusters", axes.meta = "sample",
#                axes.meta.include = kc.sample.order.skew, root.name = "prcr_q", 
#                out.dir = "sct_6_17_20/prcr_q/plots/")
# 
# 
# t <- meta.cell.count.core(so.list = list(list(so=prcr.q, prefix = "q"), list(so=prcr.c, prefix = "c")), 
#                 x.meta = "sample", x.meta.order = kc.sample.order.skew,
#                 group.meta = "seurat_clusters", group.meta.order = c("c10","q3") %>% rev())

point.color.map.processing <- 
  function(map.tsv, color.order = c("myeloid", "lymphoid", 
                                    "ery/mega", "erythroid", "mega", "HSPC", "other")) {
  df <- read.table(map.tsv, header = T, sep = "\t", quote = "", as.is = T)
  df$point.meta <- as.character(df$point.meta)
  if (!is.null(color.order)) {
    df <- df %>% filter(color %in% color.order) %>% 
      mutate(color = factor(color, levels = color.order))
  }
  return(df)
  }

# hmtp.common.markers <- function(obj, plot.dir) {
#   plot.panel.list.2(panel.list = panel.list, obj = obj, assay = "RNA",
#                     plot.out = paste0(plot.dir, "/panel.png"))
#   plot.panel.list.2(cycle.panel, obj, assay = "RNA",
#                     plot.out = paste0(plot.dir, "/cycle_panel.png"))
#   plot.panel.list.2(lineage.panel, obj, assay = "RNA",
#                     plot.out = paste0(plot.dir, "/lineage_panel.png"))
  
#   plot.panel.list.2(lsk.panel, obj, assay = "RNA", 
#                     plot.out = paste0(plot.dir, "/lsk_panel.png"))
#   plot.panel.list.2(clp.panel, obj, assay = "RNA",
#                     plot.out = paste0(plot.dir, "/clp_panel.png"))
  
#   plot.panel.list.2(GMP.panel, obj, assay = "RNA", 
#                     plot.out = paste0(plot.dir, "/GMP_panel.png"))
  
#   trash <- plot.panel.list.2(mono.panel, obj, assay = "RNA", 
#                              plot.out = paste0(plot.dir, "/mono_panel.png"))
  
#   trash <- plot.panel.list.2(dc.panel, obj, assay = "RNA", 
#                              plot.out = paste0(plot.dir, "/dc_panel.png"))
  
#   trash <- plot.panel.list.2(ery.panel, obj, assay = "RNA", 
#                              plot.out = paste0(plot.dir, "/ery_panel.png"))
#   trash <- plot.panel.list.2(mega.panel, obj, assay = "RNA", 
#                              plot.out = paste0(plot.dir, "/mega_panel.png"))
#   trash <- plot.panel.list.2(me.panel, obj, assay = "RNA", 
#                              plot.out = paste0(plot.dir, "/me_panel.png"))
#   return()
# }


sct.re.noc <- function(soi, work.dir, n.threads, do.regress=T, do.sct=T) {
  system(paste0("mkdir -p ", work.dir))
  if (is.character(soi))
    soi <- readRDS(soi)
  samples <- soi@meta.data$sample %>% as.character() %>% unique()
  sol <- mclapply(samples, function(s) return(soi[,soi@meta.data$sample == s]), mc.cores = n.threads, mc.cleanup = T)
  sol <- sct.int.prep(sol = sol, n.var.feature = 3000, regress.cycle = do.regress, out.sol.Rds = paste0(work.dir, "/sol.Rds"), sct = do.sct )
  soi.c <- sct.int.pipe(sol = sol, n.int.features = 3000, n.threads = n.threads, outfile = paste0(work.dir, "/soi.Rds"))
  return(soi.c)
}

# TopFeatures(object = q419[['pca']], dim = "PC_1", nfeatures = 20)

add.all.normalization <- function(so, methods = c("LogNormalize", "CLR", "RC"), scale.factor = 1000000) {
  assay.bk <- so@active.assay
  for (method in methods) {
    so[[method]] <- CreateAssayObject(counts = so@assays$RNA@counts)
    so <- NormalizeData(object = so, assay = method, normalization.method = method, scale.factor = scale.factor)
  }
  return(so)
}

cluster.stability.test <- function(so, seeds, assay, resolution, 
                                   comp.n.hsc = F,
                                   threads = 1, out.dir) {
  ###
  # anchor <- colnames(so)[so$seurat_clusters == "21"][6]
  # DimPlot(so, cells.highlight = anchor)
  system(paste0("mkdir -p ", out.dir))
  l <- utilsFanc::safelapply(seeds, function(i) {
    
    DefaultAssay(so) <- assay
    so <- FindClusters(object = so, resolution = resolution, random.seed = i)
    p <- DimPlot(object = so, label = T, shuffle = T) # + 
      # ggsave("clean/UMAP/plots/per_sample/WT//FindCluster_titrate/seed_" %>% paste0(i, ".png"), 
      #        device = "png", width = 6, height = 6, units = "in", dpi = 150)
    if (comp.n.hsc == T) {
      hsc <- so@meta.data["AATAGCTGTGATTACG-1", "seurat_clusters"] %>% as.character()
      n.hsc <- so@meta.data %>% filter(seurat_clusters == hsc) %>% nrow()
      names(n.hsc) <- paste0("seed_", i, "_cluster_", hsc)
      p <- p + ggtitle(paste0(names(n.hsc), ": ", n.hsc))
    } else {
      p <- p + ggtitle(paste0("seed_", i))
      n.hsc <- NA
    }
    
    return(list(p=p, n.hsc = n.hsc))
  }, threads = threads)
  
  p <- lapply(l, function(x) return(x$p)) %>% 
    wrap.plots.fanc(plot.out = paste0(out.dir, "/", assay, "_", resolution, "_", "seeds.png"), sub.width = 6)
  stats <- sapply(l, function(x) return(x$n.hsc))
  return(stats)
}

titrate.clusters <- function(soi, sol = NULL, seeds, assay, resolutions, 
                                           threads = 1, out.dir) {
  system(paste0("mkdir -p ", out.dir))
  utilsFanc::safelapply(resolutions, function(res) {
    pl <- utilsFanc::safelapply(seeds, function(seed) {
      
      DefaultAssay(soi) <- assay
      soi <- FindClusters(object = soi, resolution = res, random.seed = seed)
      p <- DimPlot(object = soi, label = T, shuffle = T) 
      p <- p + ggtitle(paste0("seed_", seeds))
      if (!is.null(sol)) {
        int.corr.plot(soi = soi, sol = sol, cluster.ident = "seurat_clusters", order = F,
                 plot.dir = paste0(out.dir, "/res_", res, "/seed_", seed, "/"))
        int.corr(q.vec = sol, s = soi, meta = "seurat_clusters", 
                 out.file = paste0(out.dir, "/res_", res, "/seed_", seed, "/int_corr.tsv"))
      }
      return(p)
    }, threads = threads)
    trash <- pl %>% wrap.plots.fanc(plot.out = paste0(out.dir, "/", assay, "_", res, "_", "seeds.png"),
                                    sub.width = 6)
  }, threads = 1)
  return()
}

adjust.cluster.core <- function(cells.df,
                           clusters.rm.vec = NULL, clusters.combine.list = NULL, 
                           clusters.sep.list = NULL) {
  # cells.df: 2 columns: cells and cluster. both col names fixed
  # rules:
  # cluster.rm.vec: just simply remove these clusters
  # cluster.combine.list: each element of the list is a vector, containing the clusters to be 
  ##combined into 1.
  # cluster.sep.list: each element of the list is also a list, with each element being the names of the 
  #cells. cluster.sep.list = list(C9 = list(subclus1 = c("cell1", "cell2", "cell3"), subclus2 = c("cell4", "cell5", "cell6")))
  #subclus1/2 are optional. 
  
  # cells.list <- list()
  # if (!is.null(so))
  # cells.list$so <- get.cell.names.seurat(so, style = "ArchR")
  # cells.list$ao <- ao$cellNames
  # cells <- cells.list[[which.max(cells.so, cells.ao)]]
  # rm(cells.list)
  clusters.to.touch <- c(clusters.rm.vec, unlist(clusters.combine.list), names(clusters.sep.list))
  if (any(duplicated(clusters.to.touch)))
    stop("some clusters have more than 1 operation specified")
  if (!is.null(clusters.rm.vec)) {
    cells.df <- cells.df %>% filter(!cluster %in% clusters.rm.vec)
  }

  if (!is.null(clusters.combine.list)) {
    for (i in 1:length(clusters.combine.list)) {
      combined.name <- names(cluster.combine.list)[i]
      to.combine <- clusters.combine.list[[i]]
      if (is.null(combined.name)) {
        combined.name <- to.combine %>% gtools::mixedsort() %>% .[1]
      }
      cells.df$clusters[cells.df$clusters %in% to.combine ] <- combined.name
    }
  }
  
  if (!is.null(clusters.sep.list)) {
    stop("clusters.sep.list part has not been tested!")
    dup.cells <- unlist(clusters.sep.list) %>% .[duplicated(.)]
    if (length(dup.cells) > 0) {
      stop("duplicated cell names detected. they are: " %>% 
             paste0(dup.cells[1:max(length(dup.cells), 5)], sep = "\n"))
    }
    for (i in 1:length(clusters.sep.list)) {
      cluster.ori <- names(clusters.sep.list)[i]
      if (is.null(cluster.ori))
        stop("is.null(cluster.ori)")
      partition.list <- clusters.sep.list[[i]]
      if (is.null(names(partition.list)))
        names(partition.list) <- paste0(cluster.ori, ".", 1:length(partition.list))
      partition.cells <- unlist(partition.cells) %>% unique() 
      
      cells.in.cluster.ori <- cells.df[cells.df$cluster == cluster.ori,"cells"]
      partition.cells.not.in.ori <- partition.cells %>% .[! partition.cells %in% cells.in.cluster.ori]
      if (length(partition.cells.not.in.ori) > 0) {
        stop(paste0("some cells specified to separate cluster ", cluster.ori,
                    " are not found in this cluster. a few of them are: \n",
                    paste0(partition.cells.not.in.ori[1:5], collapse = "\n")))
      }
      cells.add <- cells.in.cluster.ori %>% .[!. %in% partition.cells] %>% list()
      names(cells.add) <- paste0(cluster.ori, ".0")
      partition.list <- c(partition.list, cells.add)
      for (j in 1:length(partition.list)) {
        cells.df[cells.df$cells %in% partition.list[[j]], "cluster"] <- names(partition.list)[j]
      }
      
    }
  }
  
  return(cells.df)
}


# int.titrate <- function(sol, soi.rpca = NULL, work.dir, k.anchor, int.assess.panels, plot.only = F) {
#   plot.dir <- paste0(work.dir, "/plots/")
#   dir.create(plot.dir, showWarnings = T, recursive = T)
#   if (plot.only == F) {
#     features <- SelectIntegrationFeatures(object.list = sol, nfeatures = 3000)
#     sol <- PrepSCTIntegration(object.list = sol, anchor.features = features)
#     sol <- lapply(X = sol, FUN = RunPCA, features = features)
#     
#     immune.anchors <- FindIntegrationAnchors(object.list = sol, normalization.method = "SCT",
#                                              anchor.features = features, dims = 1:30, 
#                                              reduction = "rpca", k.anchor = k.anchor)
#     soi.rpca <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)
#     
#     soi.rpca <- cluster.pipe(soi = soi.rpca, assay = "integrated",
#                              pc.dim = 1:30, cluster.resolution = 0.8,
#                              work.dir = work.dir, plot.dir = plot.dir,
#                              project.name = "intImprove_rep3_rpca", metas.include = "sample",
#                              sample.tsv = "rna_exonOnly/sample_info.tsv",
#                              save.rds = T, do.sct = F,
#                              add.all.normalization = T,
#                              plot.common.markers = F)
#   } else {
#     if (is.null(soi.rpca))
#       soi.rpca <- readRDS(paste0(work.dir, "/soi.Rds"))
#   }
#   
#   try({
#     soi.rpca <- add.bad.cells(soi.rpca)
#   })
#   try({
#     plot.panel.list(panel.list = int.assess.panels, obj = soi.rpca, order = F, assay = "RNA", 
#                     split.by = "sample", plot.out = paste0(plot.dir, "/int_assess.png"))
#   })
#   return(soi.rpca)
# }

add.bad.cells <- function(soi, bad.cells.rds = "sync/bad.cells.Rds") {
  bad.cells <- readRDS(bad.cells.rds)
  if (!grepl("#", colnames(soi)[1]))
    soi <- RenameCells(soi, new.names = get.cell.names.seurat(soi, style = "ArchR"))
  soi[["bad.cells"]] <- 0
  soi$bad.cells[colnames(soi) %in% bad.cells] <- 1
  return(soi)
}

add.highlight <- function(so, add.cells.list, is.Rds = F) {
  # add.cells.list: must be named. the names will be used as colnames in so@meta.data
  for (col.name in names(add.cells.list)) {
    add.cells <- add.cells.list[[col.name]]
    if (is.Rds == T)
      add.cells <- readRDS(add.cells)
    if (!grepl("#", colnames(so)[1]))
      so <- RenameCells(so, new.names = get.cell.names.seurat(so, style = "ArchR"))
    so[[col.name]] <- 0
    so@meta.data[colnames(so) %in% add.cells, col.name] <- 1
  }

  return(so)
}

merge.2.int <- function(som, samples, nfeatures = 2000, k = 5, work.dir, plot.dir = NULL,
                        assess.list = INT.ASSESS.LIST) {
  if (is.null(plot.dir))
    plot.dir <- paste0(work.dir, "/plots")
  system(paste0("mkdir -p ", work.dir, " ", plot.dir))
  features <- som@assays$RNA@var.features[1:nfeatures] %>% .[!is.na(.)]
  sol <- lapply(samples, function(sample) {
    so <- som[, som$sample == sample]
    return(so)
  })
  names(sol) <- samples
  rm(som)
  gc()

  immune.anchors <- FindIntegrationAnchors(object.list = sol, anchor.features = features, 
                                           reduction = "rpca", k.anchor = k, assay = rep("RNA", length(sol)))
  rm(sol)
  gc()
  soi <- IntegrateData(anchorset = immune.anchors)
  soi <- cluster.pipe(soi = soi, assay = "integrated", is.sct.int = F,
                      pc.run = 30, pc.dim = 1:30, cluster.resolution = 1.0,
                      work.dir = work.dir, 
                      plot.dir = plot.dir, 
                      project.name = "som2soi", metas.include = "sample",
                      sample.tsv = "rna_exonOnly/sample_info.tsv",
                      save.rds = T, do.sct = F,
                      plot.common.markers = F)

  try({
    soi <- add.highlight(so = soi, add.cells.list = assess.list, is.Rds = T)
  })
  try({
    trash <- plot.panel.list(panel.list = names(assess.list), obj = soi, order = F, assay = "RNA",
                             plot.out = paste0(plot.dir, "/int_assess.png"), split.by = "sample", 
                             split.order = samples)
  })
  try(hmtp.common.markers(obj = soi, plot.dir = plot.dir, threads = 4, to.human = F, assay = "RNA", 
                          panels.to.plot = "general"))
  return(soi)
}

so.titrate.non.int <- function(so, assay, vars.to.regress, var.features = NULL, n.var.features = 2000,
                               pca.dims, cluster.resolution, out.dir, root.name) {
  
  DefaultAssay(so) <- assay
  
  so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
  
  so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = n.var.features)
  
  so <- ScaleData(so, vars.to.regress = c("sample", "percent.mt"), features = var.features)
  
  so <- RunPCA(so, features = var.features)
  so <- FindNeighbors(so, dims = pca.dims)
  so <- FindClusters(so, resolution = cluster.resolution)
  
  so <- RunUMAP(so, dims = pca.dims)
  plot.panel.list(panel.list = "sample", obj = so, binarize.panel = T, assay = "desoup_v1", order = F, 
                  root.name = paste0(out.dir, "/",root.name,"_umap_sample_binary"), invisible = T)
  
  plot.panel.list(panel.list = "seurat_clusters", obj = so, binarize.panel = F, assay = "desoup_v1", order = F, 
                  plot.out = paste0(out.dir, "/", root.name,"_umap_cluster.png"), invisible = T)
  trash <- hmtp.common.markers(so, assay = "desoup_v1", plot.dir = out.dir, to.human = F, panels.to.plot = "general")
  return(so)
}

cluster.titrate <- function(so, assay, 
                            titrate.umap = T, umap.pc.dims.list,
                            threads = 1,
                            out.dir) {
  
}

hash.clustering <- function(sol, grep = "Hash", assay.name = "hash",
                            bRemove.hash = T, resolution = 0.4, 
                            plot.out.root) {
  dir.create(dirname(plot.out.root), showWarnings = F, recursive = T)
  names <- names(sol)
  sol <- lapply(names, function(name) {
    so <- sol[[name]]
    hash.mat <- so@assays$ADT@counts
    if (bRemove.hash)
      rownames(hash.mat) <- sub("\\-.+$", "", rownames(hash.mat))
    hash.mat <- hash.mat %>% .[grepl(grep, rownames(.)),]
    n.hash <- nrow(hash.mat)
    hashes <- rownames(hash.mat)
    # somehow after some operations hash.assay stopped being an assay...
    so[[assay.name]] <- CreateAssayObject(counts = hash.mat)
    DefaultAssay(so) <- assay.name
    so <- NormalizeData(so)
    so <- FindVariableFeatures(so)
    so <- ScaleData(so)
    npcs <- n.hash - 1
    so <- RunPCA(so, verbose = FALSE, npcs = npcs)
    so <- FindNeighbors(so, dims = 1:npcs)
    so <- FindClusters(so, resolution = resolution, verbose = FALSE)
    so <- RunUMAP(so, dims = 1:npcs)
    plot.out.root <- paste0(plot.out.root, "_", name)
    umap.file <- paste0(plot.out.root, "_umap.png")
    DimPlot(so, label = TRUE) +
      ggtitle(name) +
      ggsave(umap.file, width = 6, height = 5, units = "in", dpi = 100)
    plot.panel.list(panel.list = hashes, obj = so, assay = assay.name, order = F, 
                    root.name = plot.out.root)
    return(so)
  })
  names(sol) <- names
  return(sol)
}

extract.bam <- function(so, sample.map = NULL,n.cells.each = NULL, scAR.naming = F, ATAC.naming = F,
                        bam.subset.range = NULL, # eg. chr6:100000-200000
                        xf.filter = NULL, filter.thread = 1, # for ADT try 25. xf: the xf tag from cellranger
                        method = "sinto", 
                        cell.list = NULL, cluster.ident, clusters = NULL, 
                        group.ident = NULL, groups = NULL,
                        threads.each = 1, npar = 1,
                        out.dir,
                        samtools = "samtools",
                        sambamba = "sambamba",
                        run = F) {
  # sample.map: 2 columns: sample and dir. dir is the cellranger counts output dir
  if (is.null(sample.map)) {
    if (is.null(so@meta.data[, "dir"])) {
      stop("is.null(so@meta.data[, 'dir])")
    }
    sample.map <- so@meta.data[, c("sample", "dir")]
  }
  sample.map <- unique(sample.map)
  if (any(!c("sample", "dir") %in% colnames(sample.map))) {
    stop("colnames of sample.map must contain sample and dir")
  }
  if (any(duplicated(sample.map$sample))) {
    stop("any(duplicated(sample.map$sample))")
  }
  
  sample.map <- sample.map[, c("sample", "dir")]
  rownames(sample.map) <- sample.map$sample
  if (is.null(cell.list))
    cell.list <- get.cell.list(obj = so, is.ao = F, group.by = cluster.ident,
                             groups = clusters, split.by = group.ident, splits = groups,
                             na.rm = T, return.named.vec = F, n.cells.each = n.cells.each)
  df <- lapply(names(cell.list), function(group) {
    cells <- cell.list[[group]]
    df <- data.frame(group = group, bc = sub(".+#", "", cells),
                     sample = paste0(sub("#.+$", "", cells)))
    
    if (length(unique(df$sample)) !=1 ) {
      stop("length(unique(df$sample)) !=1 ")
    }
    return(df)
  }) %>% do.call(rbind, .)
  
  utilsFanc::check.intersect(x = unique(df$sample), "samples", 
                             y = sample.map$sample, "sample.map")
  file <- "possorted_genome_bam.bam"
  if (scAR.naming)
    file <- "gex_possorted_bam.bam"
  if (ATAC.naming)
    file <- "atac_possorted_bam.bam"
  # file <- ifelse(scAR.naming, "gex_possorted_bam.bam", "possorted_genome_bam.bam")
  sample.map$bam.full <- paste0(sample.map$dir, "/outs/", file) %>% 
    normalizePath(mustWork = T)
  cmd.files <- list()
  if (!is.null(bam.subset.range)) {
    dir.create(out.dir, showWarnings = F, recursive = T)
    sample.map$bam <- paste0(normalizePath(out.dir), "/", sample.map$sample, "_subset.bam")
    cmds <- paste0(samtools, " view -hbo ", sample.map$bam, " ",
                   sample.map$bam.full, " ", bam.subset.range, 
                   " && ", samtools, " index ", sample.map$bam)
    if (run) {
      utilsFanc::safelapply(cmds, system, threads = npar)
    }
    cmd.file <- paste0(out.dir, "/subset.sh")
    cmd.files$subset <- cmd.file
    write(cmds, cmd.file, sep = "\n")
    system(paste0("chmod 755 ", cmd.file))
  } else {
    sample.map$bam <- sample.map$bam.full
  }
  
  if (method == "sinto") {
    cmds <- df %>% split(., f = .$sample) %>% 
      lapply(function(df) {
        sample <- df$sample[1]
        dir.create(out.dir, showWarnings = F, recursive = T)
        bc.file <- paste0(normalizePath(out.dir, mustWork = T), "/bc_",sample, ".txt")
        write.table(df[, c("bc", "group")], bc.file,
                    sep = "\t", quote = F, col.names = F, row.names = F)
        bam <- sample.map[sample, "bam"]
        cmd <- paste0("sinto filterbarcodes -b ", bam, " -c ", bc.file, " -p ",
                      threads.each, " --outdir ", normalizePath(out.dir))
        return(cmd)
      }) %>% unlist()
    if (run) {
      utilsFanc::safelapply(cmds, system, threads = npar)
    }
    cmd.file <- paste0(out.dir, "/sinto.sh")
    cmd.files$sinto <- cmd.file
    write(cmds, cmd.file, sep = "\n")
    system(paste0("chmod 755 ", cmd.file))
  }
  if (!is.null(xf.filter)) {
    bams.bc <- paste0(normalizePath(out.dir), "/", names(cell.list), ".bam")
    bams.xf <- utilsFanc::insert.name.before.ext(
      name = bams.bc, insert = paste0("xf", xf.filter), delim = "_")
    
    filter.expression <- paste0("[xf] == ", xf.filter)
    cmds <- paste0(sambamba, " view -f bam -F ", "\"", 
                   filter.expression, "\"", " -t ", filter.thread, " ", bams.bc, 
                   " > ", bams.xf, 
                   " && ", samtools, " index ", bams.xf)
    if (run) {
      utilsFanc::safelapply(cmds, system, threads = npar)
    }
    cmd.file <- paste0(out.dir, "/xf", xf.filter, ".sh")
    cmd.files$xf <- cmd.file
    write(cmds, cmd.file, sep = "\n")
    system(paste0("chmod 755 ", cmd.file))
  }
  cmd.master <- paste0(out.dir, "/run.sh")
  write("#!/bin/bash", cmd.master)
  cmd.files <- unlist(cmd.files)
  cmds <- paste0("cat ", normalizePath(cmd.files,mustWork = T), " | parallel -j ", npar, " {} ")
  write(cmds, cmd.master, append = T, sep = "\n")
  system(paste0("chmod 755 ", cmd.master))
  return(cmd.master)
}

sc.bam.df <- function(bams, sample.names = NULL,
                      fields = c("qname", "flag", "rname", "strand", "pos", "mapq"),
                      tags= c("xf", "MM", "GN", "CB")) {
  not.found <- bams %>% .[!file.exists(.)]
  if (length(not.found) > 0) {
    stop(paste0("bams not found: \n", 
                paste0(not.found, collapse = "\n")))
  }
  if (!is.null(sample.names)) {
    if (length(sample.names) == 1) {
      sample.names <- rep(sample.names, length(bams))
    }
    if (length(bams) != length(sample.names)) {
      stop("length(bams) != length(sample.names)")
    }
  }

  bam.df <- lapply(seq_along(bams), function(i) {
    bam <- bams[i]
    bam.list <- Rsamtools::scanBam(
      file = bam, param = ScanBamParam(what = fields, tag = tags))[[1]]
    df1 <- bam.list %>% .[names(.) != "tag"] %>% as.data.frame()
    df2 <- bam.list$tag %>% as.data.frame()
    df <- cbind(df1, df2)
    if (!is.null(sample.names)) {
      df$cell <- paste0(sample.names[i], "#", df$CB)
      df$sample <- sample.names[i]
    }
    df$bam <- bam
    return(df)
  }) %>% do.call(rbind, .)
  return(bam.df)
}

sc.bam.qc <- function(bams, bam.sample.names, genes,
                      out.dir = NULL, root.name = NULL) {
  bam.df <- sc.bam.df(bams = bams, sample.names = bam.sample.names)
  if (nrow(bam.df) < 1) {
    stop("nrow(bam.df) < 1")
  }
  utilsFanc::check.intersect(c("sample", "cell","GN", "MM", "xf"),
                             "required fields or tags", 
                             y = colnames(bam.df), "colnames(bam.df)")

  genes.in <- bam.df$GN %>% unique() %>% .[. %in% genes]
  if (length(genes.in) < 1) {
    stop("no common genes between bam files and genes supplied")
  }
  bam.df <- bam.df %>% dplyr::filter(GN %in% genes)
  bam.df$MM[is.na(bam.df$MM)] <- 0
  qc <- list()
  qc <- lapply(c("sample", "cell"), function(field) {
    if (field == "cell") {
      df <- bam.df %>% dplyr::group_by(sample, cell, GN)
    } else if (field == "sample") {
      df <- bam.df %>% dplyr::group_by(sample, GN)
    } else {
      stop()
    }
    df <- df %>% dplyr::summarise(n.UMI = n(), n.MM = sum(MM)) %>% 
      dplyr::ungroup() %>% as.data.frame() %>% 
      dplyr::mutate(frac.MM = round(n.MM/n.UMI, digits = 3))
    if (!is.null(out.dir)) {
      dir.create(out.dir, showWarnings = F, recursive = T)
      if (is.null(root.name)) {
        root.name <- paste0(basename(out.dir), "_MM")
      }
      qc.file <- paste0(out.dir, "/", root.name, "_", field, ".tsv")
      write.table(df, qc.file, quote = F, sep = "\t", col.names = T,
                  row.names = F)
    }
    return(df)
  })
  names(qc) <- c("sample", "cell")
  invisible(qc)
}

sc.MM.correct <- function(so, MM.df) {
  # MM.df: required columns: cell, GN, n.MM. 
  # MM refers to the MM tag that indicates that cellranger changed STAR's MAPQ.
  # GN: the cellranger tag that contains gene name
  if (is.character(MM.df)) {
    MM.df <- read.table(MM.df, header = T, sep = "\t", comment.char = "")
  }
  utilsFanc::check.intersect(c("cell", "n.MM", "GN"), "required fields",
                             colnames(MM.df), "colnames(MM.df)")
  # if (!identical(sort(colnames(so)), unique(sort(MM.df$cell)))) {
  #   stop("!identical(sort(colnames(so)), sort(MM.df$cell))")
  #   # this part will never be identical. Because to generate MM.df, you almost
  #   # always look at a sinlge region of the genome and some cells might 
  #   # have no reads in that region.
  # }
  utilsFanc::check.intersect(unique(MM.df$cell), "MM.df$cell",
                             colnames(so), "colnames(so)", warning.only = F)

  utilsFanc::check.intersect(colnames(so), "colnames(so)",
                             unique(MM.df$cell), "MM.df$cell",
                             warning.only = T)
  utilsFanc::check.intersect(unique(MM.df$GN), "MM.df$GN",
                             rownames(so), "rownames(so)")
  exp.mat <- MM.df %>% 
    reshape2::acast( GN ~ cell, value.var = "n.MM")
  exp.mat[is.na(exp.mat)] <- 0
  for (gene in rownames(exp.mat)) {
    so@assays$RNA@counts[gene, colnames(exp.mat)] <- 
      so@assays$RNA@counts[gene, colnames(exp.mat)] - 
      exp.mat[gene, ]
  }
  
  if (any(so@assays$RNA@counts[, colnames(exp.mat)]) < 0) {
    browser()
  }
  so <- Seurat::NormalizeData(so, assay = "RNA")
  warning("you might need to re-do a lot of things after changing gene expression.")
  return(so)
}
# sc.mapping.qc <- function(so, sample.map,
#                           bam.subset.range, genes, 
#                           xf.filter,
#                           cluster.ident, clusters = NULL,
#                           group.ident, groups = NULL,
#                           threads = 1, npar = 1,
#                           out.dir) {
#   
# }

cluster.correspondence <- function(df, x, y, step.pct = 0.5, max.steps = 3,
                                   out.file = NULL) {
  # something like: df=soi@meta.data, x = "seurat_clusters", y = "Clusters"
  # initially written to find which ATAC clusters correspond to a given RNA cluster.
  df <- df[, c(x, y)]
  names(df) <- c("x", "y")
  mat <- table(df) %>% as.matrix()
  y.names <- colnames(mat)
  
  out.df <- lapply(1:nrow(mat), function(i) {
    x.name <- rownames(mat)[i]
    x <- mat[i, ] %>% as.vector()
    names(x) <- y.names
    found <- c()
    max.value.last <- 1
    
    while (length(x) > 0) {
      max.id <- which.max(x)
      max.value <- x[max.id]
      pct <- max.value/max.value.last
      if (pct < step.pct) break
      
      found <- c(found, names(x)[max.id])
      x <- x[-max.id]
      max.value.last <- max.value
    }
    raw <- paste0(found, collapse = ",")
    if (length(found) > max.steps) {
      corr <- NA
    } else {
      corr <- raw
    }
    df <- data.frame(x = x.name, y = corr, raw = raw)
    return(df)
  }) %>% do.call(rbind, .)
  colnames(out.df) <- c(x, y, "raw")
  if (!is.null(out.file)) {
    dir.create(dirname(out.file), showWarnings = F, recursive = T)
    write.table(out.df, out.file, quote = F, sep = "\t", col.names = T, row.names = F)
  }
  invisible(out.df)
} 

filter.by.cluster.correspondence <- function(
  so, cluster.ident = "seurat_clusters", filter.by = "Clusters",
  new.cluster.ident = "scc", 
  filter.df = NULL, step.pct = 0.5, max.steps = 3,
  out.dir) {
  
  if (is.null(filter.df)) {
    filter.df <- cluster.correspondence(df = so@meta.data, x = cluster.ident, y = filter.by,
                                        step.pct = step.pct, max.steps = max.steps, 
                                        out.file = paste0(out.dir, "/cluster_ident_correspondence.tsv"))
  }
  
  x <- filter.df[, cluster.ident]
  y <- filter.df[, filter.by] %>% strsplit(",")
  
  possible.clusters <- lapply(1:length(x), function(i) {
    paste0(x[i], "@", y[[i]])
  }) %>% unlist()
  
  so@meta.data$tmp <- paste0(so@meta.data[, cluster.ident], "@", so@meta.data[, filter.by])
  so@meta.data[, new.cluster.ident] <- so@meta.data[, cluster.ident]
  so@meta.data[, new.cluster.ident][! so@meta.data$tmp %in% possible.clusters] <- NA
  
  cluster.freq <- table(so@meta.data[, c(new.cluster.ident, "sample")]) %>% as.matrix()
  
  dir.create(out.dir, showWarnings = F, recursive = T)
  write.csv(cluster.freq, paste0(out.dir, "/cluster_freq.csv"))
  return(so)
}