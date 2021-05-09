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
    all= log2(overallBaseMean[,"bulkNorm_WT"] + 1),
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

de.2.rnk <- function(de.grid, comp, pos.filter=NULL, neg.filter=NULL, out.dir, root.name = NULL) {
  system(paste0("mkdir -p ", out.dir))
  df <- de.grid$grid.sum %>% filter(comp == comp) %>% 
    mutate(log_p = -log10(p) * (logFC/abs(logFC)) , gene = toupper(gene)) %>% 
    filter(log_p != 0) %>% select(gene, log_p, master.ident)
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
t.f.gsea.diag <- function() {
  print("miao")
  print("miao")
  print("mial")
}