fg.archr <- function(ao, cluster.ident, cluster, group.by, groups = NULL) {
  # note: cluster, not clusters. only support subsetting by one cluster
  ao <- ao[getCellColData(ao, cluster.ident, drop = T) == cluster,]
  cell.list <- get.cell.list(obj = ao, is.ao = T, group.by = group.by, groups = groups, na.rm = T)
  lapply(names(cell.list), function(group.name) {
    cells <- cell.list[[group.name]]
    ao <- ao[ao$cellNames %in% cells]

    t <- getMarkerFeatures(ao, groupBy = "Sample", useGroups = "tumor_rep1", bgdGroups = "ctrl_rep1", 
                           useMatrix = "PeakMatrix", testMethod = "wilcoxon",
                           bias = c("TSSEnrichment", "log10(nFrags)"))
  })
}


bg.gen.2 <- function(mat, fg.vec, n.bg.each, method = "euclidean", scale = "none", no.replace = F,
                     samples.use = NULL, return.details = F) {
  # mat has to be named
  if (!is.matrix(mat))
    mat <- as.matrix(mat)
  if (!is.null(samples.use)) {
    mat <- mat[, samples.use]
  }
  dup.fg.vec <- fg.vec %>% .[duplicated(.)]
  if (length(dup.fg.vec) > 0)
    stop(paste0("certain fg input are duplicated: \n",
                paste0(dup.fg.vec[1:5], collapse = "\n")))

  if (no.replace == T) {
    mat.fg <- mat[rownames(mat) %in% fg.vec, ]
    mat.bg <- mat[!rownames(mat) %in% fg.vec, ]
    bg.list <- list()
    for (fg in fg.vec) {
      mat.sub <- rbind(mat.fg[fg,], mat.bg %>% .[!rownames(.) %in% unlist(bg.list),])
      bg <- genefilter::genefinder(X = mat.sub, ilist = 1, numResults = n.bg.each, 
                                   method = method, scale = "none")
      bg <- rownames(mat.sub)[bg[[1]]$indices]
      bg.list <- c(bg.list, list(bg))
    }
    bg.vec <- unlist(bg.list)
    if (any(duplicated(bg.vec))) {
      stop("any(duplicated(bg.vec))")
    }
  } else {
    fg.id <- sapply(fg.vec, function(x) {
      id <- which(rownames(mat) == x)
      if (length(id) != 1)
        stop(paste0(x, "did not return one unique id"))
      return(id)
    })
    bg.list <- genefilter::genefinder(X = mat, ilist = fg.id, numResults = n.bg.each, 
                                      method = method, scale = "none")
    names(bg.list) <- fg.vec
    
    if (return.details == T)
      return(bg.list)
    bg.id <- lapply(bg.list, function(x) return(x$indices)) %>% Reduce(c, .) %>% 
      unique()
    bg.id <- bg.id %>% .[!.%in% fg.id]
    bg.vec <- rownames(mat)[bg.id]
  }

  return(bg.vec)
}

bg.gen.bin <- function(mat, n.bins, samples.use = NULL, return.df = F) {
  # mat has to be named
  if (!is.matrix(mat))
    mat <- as.matrix(mat)
  if (!is.null(samples.use)) {
    mat <- mat[, samples.use]
  }
  mean.rank <- rowMeans(mat) %>% rank(ties.method = "random")
  n.per.bin <- ceiling(nrow(mat)/n.bins)
  split.vec <- ceiling(mean.rank/n.per.bin)
  bg.list <- rownames(mat) %>% split(f = split.vec)
  if (return.df == T) {
    bg.list <- bg.list %>% lapply(function(x) return(utilsFanc::loci.2.df(loci.vec = x)))
  }
  return(bg.list)
}

bg.assess <- function(mat, fg.vec, bg.vec, scatter.xy.df, 
                      transformation = NULL, plot.out, return.pl = F, ...) {
  # scatter.xy.df: cols: x, y.
  df <- mat %>% as.data.frame()
  df$gene <- rownames(df)
  rownames(df) <- NULL
  gene.list <- list(fg = fg.vec, bg = bg.vec)
  pl <- scatter.xy.df %>% unique() %>% split(., f = 1:nrow(.)) %>% 
    lapply(function(comp) {
      # comp: comparison
      pl <- lapply(names(gene.list), function(type) {
        genes <- gene.list[[type]]
        p <- xy.plot(df = df, x = comp$x, y = comp$y, transformation = transformation,
                     highlight.var = "gene", highlight.values = genes, show.highlight.color.var = F,
                     plotly.var = "gene", color.density = T, ...) +
          ggtitle(paste0(comp$y, ":", comp$x, "..", type))
        return(p)
      })
      names(pl) <- names(gene.list)
      return(pl)
    }) %>% Reduce(c, .)
  names(pl) <- paste0(scatter.xy.df$y, "_", scatter.xy.df$x, "_", names(pl))
  if (return.pl == T)
    return(pl)
  trash <- wrap.plots.fanc(plot.list = pl, plot.out = plot.out, n.split = length(gene.list))
  return()
}

bg.gen <- function(fg, n.bg, bulkNorm, column, out.root.name=NULL) {
  # out.root.name should contain directory
  
  # some of the nomenclature doesn't exact make sense. This is because 
  # the code is mostly taken from https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#gene-ontology-enrichment-analysis
  # and I didn't bother to re-write everything
  overallBaseMean <- bulkNorm[, c("gene", column), drop = F]
  rownames(overallBaseMean) <- overallBaseMean$gene
  overallBaseMean$gene <- NULL
  overallBaseMean <- as.matrix(overallBaseMean)
  
  de.genes <- fg
  sig_idx <- match(de.genes, rownames(overallBaseMean))
  all.genes <- rownames(overallBaseMean)
  
  stats <- lapply(sig_idx, function(i) {
    ind <- genefilter::genefinder(overallBaseMean, i, 20, method = "manhattan")[[1]]$indices
    bg.genes <- all.genes[ind]
    bg.exp <- overallBaseMean[bg.genes, column]
    stat <- data.frame(fg.id = i,
                       fg.gene = all.genes[i],
                       bg.id = ind,
                       bg.genes = bg.genes,
                       fg.exp =  overallBaseMean[all.genes[i], column],
                       bg.exp = bg.exp,
                       bg.mean = mean(bg.exp),
                       bg.median = median(bg.exp),
                       bg.min = min(bg.exp),
                       bg.max = max(bg.exp))
    
    rownames(stat) <- NULL
    
    return(stat)
  }) %>% Reduce(rbind,.) 
  # length(which(base::duplicated(stats$bg.genes))) # 482 out of 920
  stats.sum <- stats %>% group_by(fg.id, fg.gene) %>% 
    summarise(bg.id = paste0(bg.id, collapse = ","),
              bg.genes = paste0(bg.genes, collapse = ","),
              fg.exp = fg.exp[1],
              bg.exp = paste0(bg.exp, collapse = ","),
              bg.mean = bg.mean[1],
              bg.median = bg.median[1],
              bg.min = bg.min[1],
              bg.max = bg.max[1]) %>% ungroup()
  
  backG <- unique(stats$bg.genes)

  backG <- setdiff(backG,  de.genes)
  
  if (!is.null(out.root.name)) {
    write.table(stats, paste0(out.root.name, "_stats.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
    write.table(stats.sum, paste0(out.root.name, "_sum.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
    write(backG, paste0(out.root.name, "_bg.txt"), sep = "\t")
  }
  
  png(filename = paste0(out.root.name, "_distro.png"))
  try(geneplotter::multidensity( list( 
    all= log2(overallBaseMean[,column] + 1),
    foreground =log2(overallBaseMean[de.genes, column]), 
    background =log2(overallBaseMean[backG, column])), 
    xlab="log2 mean normalized counts", main = "expression distro"))
  dev.off()
  
  
  return(list(bg = backG, stats = stats, sum = stats.sum))
}

gsea.fanc <- function(rnk.vec, gmt.vec, out.dir, thread.rnk=1, thread.gmt=1, n.plot = 20,
                      log.file = NULL, run = T, zip = T, force = F, parse = T, parse.only = F,
                      q.filter = NULL, n.filter = 50, 
                      gsea = "~/software/gsea/GSEA_Linux_4.1.0/gsea-cli.sh", gsea.params=GSEA.DEFAULT) {
  # rnk.vec and gmt.vec should both be named. Their names will be used to create output rpt_label.
  system(paste0("mkdir -p ", out.dir))
  # if (is.character(rnk.vec))
  #   rnk.vec <- readRDS(rnk.vec)
  # if (is.character(gmt.vec))
  #   gmt.vec <- readRDS(gmt.vec)
  
  res <- mclapply(seq_along(rnk.vec), function(i) {
    rnk <- rnk.vec[i]
    rnk.name <- names(rnk.vec)[i]
    if (is.null(rnk.name))
      rnk.name <- basename(rnk) %>% tools::file_path_sans_ext()
    mclapply(seq_along(gmt.vec), function(j) {
      gmt <- gmt.vec[j]
      gmt.name <- names(gmt.vec)[j]
      #print(gmt.name)
      if (is.null(gmt.name))
        gmt.name <- basename(gmt) %>% sub(".symbols.gmt$", "", .)
      label <- paste0(rnk.name, "..", gmt.name)
      out.dir.sub <- paste0(out.dir, "/", rnk.name)
      out.dir.base <- paste0(out.dir.sub, "/", label, ".GseaPreranked*")
      # print(out.dir.base)
      if (parse.only == F) {
        if (length(Sys.glob(out.dir.base)) > 0) {
          print(paste0("note: Gsea result for ", out.dir.base, " already exists"))
          if (force == T)
            system(paste0("rm -rf ",out.dir.base ))
          else
            stop("force = T is required to remove previous results and proceed")
        }
        # time.stamp <- as.character(floor(as.numeric(Sys.time())))
        cmd <- paste0(gsea, " GSEAPreranked", " -gmx ", gmt, " -rnk ", rnk, 
                      " -rpt_label ", label, " -plot_top_x ", n.plot,
                      " -out ", out.dir.sub, 
                      " ", gsea.params)
        utilsFanc::cmd.exec.fanc(cmd, stdout.file = log.file, run = run, intern = F)
      }
      
      if (parse == T) {
        out.dir.base <- Sys.glob(out.dir.base)
        if (length(out.dir.base) == 0)
          return(NULL)
        df <- lapply(c("pos", "neg"), function(x) {
          # browser()
          df.sub <- read.table(Sys.glob(paste0(out.dir.base, "/gsea_report_for_na_", x, "*.tsv")),
                               as.is = T, header = T, sep = "\t", quote = "")
          
          if (!is.null(n.filter))
            df.sub <- df.sub[1:n.filter,]
          if (!is.null(q.filter))
            df.sub <- df.sub %>% filter(FDR.q.val <= q.filter)
          n.row <- nrow(df.sub)
          add.df <- data.frame(label = rep(label, n.row), direction = rep(x, n.row),
                               rnk = rep(rnk.name, n.row), gmt = rep(gmt.name, n.row))
          df.sub <- utilsFanc::add.column.fanc(df1 = df.sub, df2 = add.df, pos = 1)
          return(df.sub)
        }) %>% Reduce(rbind, .)
        return(df)
      } else {
        return()
      }
      
    }, mc.cores = thread.gmt) %>% Reduce(rbind, .) %>% return()
  }, mc.cores = thread.rnk) %>% Reduce(rbind, .)
  if (!is.null(res)) {
    write.table(res, paste0(out.dir, "/parse.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
  }
  if (zip == T && parse.only == F) {
    cmd <- paste0("zip -r -q ", sub("/$", "", out.dir), ".zip", " " ,out.dir)
    print(cmd); try(system(cmd))
  }
  return(res)
}

de.2.rnk <- function(de.grid = NULL, pbl = NULL, pbl.topn = NULL, pbl.samples, pbl.clusters = NULL,
                     comp, pos.filter=NULL, neg.filter=NULL, out.dir, root.name = NULL) {
  # pbl: pseudobulk list
  # pbl.samples is really just used to figure out what are the columns corresponding to samples when
  #taking topn
  system(paste0("mkdir -p ", out.dir))
  if (!is.null(de.grid)) {
    comparison <- comp
    df <- de.grid$grid.sum %>% filter(comp == comparison) %>% 
      mutate(log_p = -1 * utilsFanc::log2pm(p, base = 10) * (logFC/abs(logFC)) , gene = toupper(gene)) %>% 
      filter(log_p != 0) %>% select(gene, log_p, master.ident)
  } else {
    if (is.character(pbl))
      pbl <- readRDS(pbl)
    pbl <- deseq2.summary(pbl = pbl, padj.cutoff = 0.05, force = T)
    # only used to remove invalid elements in the list
    if (!is.null(pbl.clusters)) {
      pbl <- pbl[names(pbl) %in% pbl.clusters]
    }
    df <- lapply(names(pbl), function(name) {
      s2b <- pbl[[name]]
      if (!is.null(pbl.topn)) {
        top.genes <- lapply(pbl.samples, function(sample) {
          genes <- s2b$res.exp$gene[rank(-1 * s2b$res.exp[, sample], ties.method = "random") <= pbl.topn]
          return(genes)
        }) %>% Reduce(union, .)
        s2b$res.exp <- s2b$res.exp %>% filter(gene %in% top.genes)
      }
      df <- s2b$res.exp %>% select(gene, pvalue, log2FoldChange) %>% 
        mutate(log_p = -log10(pvalue) * (log2FoldChange/abs(log2FoldChange)), gene = toupper(gene), master.ident = name) %>% 
        filter(log_p != 0) %>% select(gene, log_p, master.ident)
      return(df)
    }) %>% Reduce(rbind, .)
  }
  
  rnk.list <- df %>% split(f = df$master.ident) %>% 
    lapply(function(x) {
      cluster <- x$master.ident[1]
      x <- x %>% select(gene, log_p) %>% arrange(desc(log_p))
      if (!is.null(pos.filter))
        x <- x %>% filter(grepl(pos.filter, gene))
      if (!is.null(neg.filter))
        x <- x %>% filter(!grepl(neg.filter, gene))
      rnk.name <- paste0("cluster_", cluster)
      if (!is.null(root.name))
        rnk.name <- paste0(root.name, "_", rnk.name)
      rnk <- paste0(out.dir,"/", rnk.name, ".rnk")
      write.table(x, rnk, sep = "\t", row.names = F, col.names = F, quote = F)
      names(rnk) <- rnk.name
      return(rnk)
    })
  names(rnk.list) <- NULL
  rnk.vec <- unlist(rnk.list)
  saveRDS(rnk.vec, paste0(out.dir, "/rnk.vec.Rds"))
  return(rnk.vec)
}

gsea.translate <- function(genes.gsea, so=NULL, lookup = MOUSE.GENES, na.2.b2m) {
  if (!is.null(so))
    lookup <- rownames(so)
  names(lookup) <- toupper(lookup)
  out <- lookup[genes.gsea]
  if (na.2.b2m == T) {
    nas <- is.na(out)
    out[nas] <- "B2m"
    names(out)[nas] <- "B2M"
  }
  return(out)
}

# gsea.translate("miao", sol[[1]], na.2.b2m = T)

gsea.get.gs <- function(tsv = NULL, master.dir = NULL, rnk = NULL, gmt = NULL, gene.set = NULL) {
  if (is.null(tsv)) {
    q <- paste0(master.dir, "/", rnk, "/", rnk, "..", gmt, "*/", gene.set, ".tsv")
    tsv <- Sys.glob(q)
    if (length(tsv) != 1) {
      stop (paste0("the number of files found through the glob term ",q," is: ", length(tsv),
                   ". It should be exactly 1"))
    }
  }
  
  df <- read.table(tsv, header = T, as.is = T, sep = "\t", quote = "")
  return(df)
}


gmt.get.gs <- function(gmt = GMT.FILES["msigdb.v7.2"], gs.names) {
  gs.df <- lapply(gs.names, function(gs) {
    if (is.character(gmt)) {
      gmt <- clusterProfiler::read.gmt(gmtfile = gmt)
    }
    df <- gmt[grepl(gs, gmt$ont), ] 
    df <- df %>% group_by(ont) %>% summarise(gene = paste0(gene, collapse = ";"))
    return(df)
  }) %>% Reduce(rbind, .)
  gs.list <- gs.df$gene %>% strsplit(";")
  names(gs.list) <- gs.df$ont
  return(gs.list)
}

# gpo: gsea parse object
gpo.get.gs <- function(gpo, master.dir, label=NULL, direction = NULL, rnk = NULL, gmt = NULL,
                       genesets = NULL, NES.cutoff = NULL, p.cutoff = NULL, q.cutoff = NULL) {
  char.args <- c("label", "direction", "rnk", "gmt", "genesets")
  for (arg in char.args) {
    value <- get(arg)
    if (!is.null(value)) {
      if (arg %in% c("genesets", "gmt")) {
        if (arg == "genesets")
          arg <- "NAME"
        gpo <- lapply(value, function(x){
          return(gpo %>% .[grepl(x,.[,arg]), ])
        }) %>% Reduce(rbind, .)
      } else {
        gpo <- gpo[gpo[, arg] %in% value,]
      }
    }
  }
  
  if (!is.null(NES.cutoff))
    gpo <- gpo[gpo$NES > NES.cutoff, ]
  if (!is.null(p.cutoff))
    gpo <- gpo[gpo$NOM.p.val < p.cutoff, ]
  if (!is.null(q.cutoff))
    gpo <- gpo[gpo$FDR.q.val < q.cutoff, ]
  
  gs.list <- gpo %>% split(., f = paste0(.$rnk, "..", .$gmt, "..", .$NAME)) %>% 
    lapply(function(x) {
      gs <- gsea.get.gs(master.dir = master.dir, rnk = x$rnk, gmt = x$gmt, gene.set = x$NAME)
      return(gs$SYMBOL)
    })
  return(gs.list)
}

gsea.scatter <- function(gs.list=NULL, gmt= GMT.FILES["msigdb.v7.2"], gs.names = NULL, 
                         so = NULL, assay, group.by, comp, cluster.ident, clusters = NULL, 
                          bulk.list = NULL, quantile.limit = NULL, transformation = NULL,
                         plot.dir, threads = 4) {
  if (is.null(gs.list))  {
    gs.list <- gmt.get.gs(gmt = gmt, gs.names = gs.names)
  }
  
  if (is.null(names(gs.list)))
    stop(paste0("gs.list must be named!!"))
  # browser()
  
  if (is.null(bulk.list)) {
    if (is.null(clusters)) {
      clusters <- so@meta.data[, cluster.ident] %>% as.character() %>% unique() %>% gtools::mixedsort()
    }
    
    df.list <- mclapply(clusters, function(cluster) {
      mean.df <- so.2.bulk(so = so, assay = assay, slot = "data", sub.idents = cluster, 
                           ident = cluster.ident, group.by = group.by, 
                           take.mean = T, coldata.columns = NULL)$bulk.mat %>% as.data.frame()
      mean.df$gene <- rownames(mean.df)
      return(mean.df)
    }, mc.cleanup = T, mc.cores = threads)
    names(df.list) <- paste0(cluster.ident, "_", clusters)
  } else {
    df.list <- lapply(bulk.list, function(x) return(x$bulkNorm))
  }
  
  mclapply(seq_along(gs.list), function(i) {
    genes <- gs.list[[i]] %>% gsea.translate(na.2.b2m = T)
    gs.name <- names(gs.list)[i]
    sample.x <- sub("^.+:", "", comp)
    sample.y <- sub(":.+$", "", comp)
    p.list <- lapply(seq_along(df.list), function(j) {
      df <- df.list[[j]]
      cluster.name <- names(df.list)[j]
      p <- xy.plot(df = df, x = sample.x, y = sample.y, highlight.var = "gene", 
                   highlight.values = genes, transformation = transformation,
                   quantile.limit = quantile.limit, outfile = NULL)
      p <- p + ggtitle(cluster.name)
      return(p)
    })
    trash <- wrap.plots.fanc(plot.list = p.list, page.limit = 100, 
                             plot.out = paste0(plot.dir, "/", gs.name, ".png"))
    return()
  }, mc.cleanup = T, mc.cores = threads)
  return()
  
}


gsea.panel.list <- function(gs.list = NULL, gs.df.list = NULL, so, assay, order = F, split.by = NULL, 
                            plot.dir, threads = 6, max.quantile = NULL, plot.vln = T, plot.scatter = T,
                            max.plot = 200, ...) {
  if (is.null(gs.list)) {
    if (is.null(names(gs.df.list)))
      stop("gs.df.list must be named")
    gs.list <- lapply(gs.df.list, function(gs.df) {
      return(gs.df$SYMBOL)
    })
    names(gs.list) <- names(gs.df.list)
  }
  # browser()
  if (is.null(names(gs.list))) {
    stop("gs.list must be named")
  }
  mclapply(seq_along(gs.list), function(i) {
    genes <- gs.list[[i]] %>% gsea.translate(na.2.b2m = T)
    if (length(genes) > max.plot)
      genes <- sample(genes, size = max.plot, replace = F)
    gs.name <- names(gs.list)[i] 
    if (plot.vln == T) {
      plot.panel.list(panel.list = genes, obj = so, order = order, assay = assay,
                      split.by = split.by,
                      pt.size = 0.05, raster = T, page.limit = 4, n.col = 4, violin = T,  
                      plot.out = paste0(plot.dir, "/", gs.name, "..vln.pdf"), add.median = F, 
                      max.quantile = max.quantile, threads = 4, sub.width = 7.5,
                      sub.height = 3)
    }
    
    if (plot.scatter == T) {
      plot.panel.list(panel.list = genes, obj = so, order = order, assay = assay,
                      split.by = split.by,
                      pt.size = 1, raster = T, page.limit = 20, violin = F,  
                      plot.out = paste0(plot.dir, "/", gs.name, "..scatter.pdf"), 
                      max.quantile = max.quantile, threads = 4)
    }

  }, mc.cores = threads, mc.cleanup = T)
}

gsea.plot.batch <- function(gs.vec.list, so, gpo, 
                            assay, group.by, comp, cluster.ident, clusters,
                            master.dir, plot.dir) {
  lapply(seq_along(gs.vec.list), function(i) {
    gs.vec <- gs.vec.list[[i]]
    rnk <- names(gs.vec.list)[i]
    gs.list <- gpo.get.gs(gpo = gpo, master.dir = master.dir, 
                          genesets = gs.vec, rnk = rnk)
    gsea.scatter(gs.list = gs.list, so = so, assay = assay, group.by = group.by, comp = comp,
                 cluster.ident = cluster.ident, clusters = clusters, quantile.limit = 0.999, 
                 plot.dir =  plot.dir)
    gsea.panel.list(gs.list = gs.list, so = so, assay = assay, split.by = group.by, 
                    plot.dir =  plot.dir)
    return()
  })
}

chromVAR.pipe <- function(dev = NULL, ao = NULL, ao.annotation.name, peakmat, cells = NULL, peaks = NULL,
                          motifs.use = NULL,
                          compute.dev = T, compute.syn = F, compute.cor = F,
                          out.dir, root.name = "chromVAR",
                          do.ridge = T, do.density = T, group.by, groups  = NULL, threads = 1,
                          motif.base = "homer", bs.genome) {
  # bs.genome: something like BSgenome.Mmusculus.UCSC.mm10
  # if dev is specified, this basically becomes a ploting function
  if (is.null(dev)) {
    rownames(peakmat) <- rowRanges(peakmat) %>% utilsFanc::gr.get.loci()
    if (!is.null(peaks)) {
      peakmat <- peakmat[peaks,]
    }
    if (!is.null(cells)) {
      peakmat <- peakmat[, cells]
    }
    
    peakmat <- addGCBias(object = peakmat, genome = bs.genome)
    assayNames(peakmat) <- "counts" 
    if (is.null(ao)) {
      if (motif.base != "homer") {
        stop("only homer motif base is developed")
      }
      data("homer_pwms")
      motifs <- homer_pwms
      obj <- .summarizeChromVARMotifs(motifs)
      motifs <- obj$motifs
      motifSummary <- obj$motifSummary
      motif_ix <- motifmatchr::matchMotifs(motifs, peakmat, 
                                           genome = bs.genome)
    } else {
      motif_ix <- readRDS(getPeakAnnotation(ao, ao.annotation.name)$Matches)
      assayNames(motif_ix) <- "matches"
      rownames(motif_ix) <- utilsFanc::gr.get.loci(rowRanges(motif_ix))
      utilsFanc::check.intersect(x = rownames(peakmat), x.name = "rownames(peakmat)",
                                 y = rownames(motif_ix), y.name = "rownames(motif_ix")
      motif_ix <- motif_ix[rownames(motif_ix) %in% rownames(peakmat),]
    }
    
    if (compute.dev == T) {
      utilsFanc::t.stat("computing deviation")
      dev <- computeDeviations(object = peakmat, annotations = motif_ix)
      dir.create(out.dir, showWarnings = F, recursive = T)
      saveRDS(dev, paste0(out.dir, "/", root.name, "_dev.Rds"))
    }
    if (compute.cor == T) {
      utilsFanc::t.stat("computing correlation")
      if (!is.null(motifs.use))
        motif_ix <- motif_ix[, motifs.use]
      res <- getAnnotationCorrelation(object = peakmat, annotations = motif_ix)
      dir.create(out.dir, showWarnings = F, recursive = T)
      saveRDS(res, paste0(out.dir, "/", root.name, "_cor.Rds"))
    }
    if (compute.syn == T) {
      utilsFanc::t.stat("computing synergy")
      if (!is.null(motifs.use))
        motif_ix <- motif_ix[, motifs.use]
      res <- getAnnotationSynergy(object = peakmat, annotations = motif_ix)
      dir.create(out.dir, showWarnings = F, recursive = T)
      saveRDS(res, paste0(out.dir, "/", root.name, "_syn.Rds"))
    }
    
  } else {
    if (is.character(dev))
      dev <- readRDS(dev)
  }
  
  if (do.ridge == T) {
    try(chromVAR.ridge(obj = dev, motif.regex = ".+", plot.out = paste0(out.dir, "/", root.name, "_zScore_ridge.png"),
                   group.by = group.by, groups = groups, threads = threads))
    
    
    try(chromVAR.ridge(obj = dev, motif.regex = ".+", plot.out = paste0(out.dir, "/", root.name, "_raw_ridge.png"),
                   group.by = group.by, groups = groups, assay = "deviations", threads = threads))
  }  
  
  if (do.density == T) {
    try(chromVAR.ridge(obj = dev, motif.regex = ".+", plot.out = paste0(out.dir, "/", root.name, "_zScore_density.png"),
                       group.by = group.by, groups = groups, threads = threads, use.ridge = F))
    
    
    try(chromVAR.ridge(obj = dev, motif.regex = ".+", plot.out = paste0(out.dir, "/", root.name, "_raw_density.png"),
                       group.by = group.by, groups = groups, assay = "deviations", threads = threads, 
                       use.ridge = F))
  }
  
  
  return(dev)
}