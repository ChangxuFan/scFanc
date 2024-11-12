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
    mat <- mat[, samples.use, drop = F]
  }
  dup.fg.vec <- fg.vec %>% .[duplicated(.)]
  if (length(dup.fg.vec) > 0)
    stop(paste0("certain fg input are duplicated: \n",
                paste0(dup.fg.vec[1:5], collapse = "\n")))

  if (no.replace == T) {
    mat.fg <- mat[rownames(mat) %in% fg.vec, , drop = F]
    mat.bg <- mat[!rownames(mat) %in% fg.vec, , drop = F]

    n.bg.needed <- length(fg.vec) * n.bg.each

    if (n.bg.needed > nrow(mat.bg)) {
      stop(paste0("more bg requested than there is available:\n",
                  "needed: ", length(fg.vec), " x ", n.bg.each, " = ", n.bg.needed, "\n",
                  "available: ", nrow(mat.bg)))
    }

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

bg.assess <- function(mat, fg.vec, bg.vec = NULL,
                      scatter.xy.df, scatter.meta = NULL, col.data,
                      transformation = NULL, plot.out, return.pl = F, ...) {
  # scatter.xy.df: cols: x, y.
  df <- mat %>% as.data.frame()
  df$gene <- rownames(df)
  rownames(df) <- NULL
  gene.list <- list(fg = fg.vec, bg = bg.vec)
  gene.list <- gene.list[!sapply(gene.list, is.null)]
  pl <- scatter.xy.df %>% unique() %>% split(., f = 1:nrow(.)) %>%
    lapply(function(comp) {
      pl <- lapply(names(gene.list), function(type) {
        genes <- gene.list[[type]]

        comp.meta <- scatter.meta
        comp.l <- comp %>% as.list()
        if (!is.null(comp.meta)) {
          for (i in 1:length(comp.l)) {
            samples.in.comp <- col.data %>% .[.[, comp.meta] %in% comp.l[[i]],] %>% rownames()
            df[, comp.l[[i]]] <- df[, samples.in.comp] %>% utilsFanc::pmean()
          }
        }

        p <- xy.plot(df = df, x = comp.l$x, y = comp.l$y, transformation = transformation,
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

fgsea.fanc <- function(rnk, gmt, nperm = NULL) {
  rank.df <- read.table(rnk, header = F)
  ranks <- rank.df$V2
  names(ranks) <- rank.df$V1
  ranks <- sort(ranks)

  genesets <- clusterProfiler::read.gmt(gmtfile = gmt)
  genesets <- split(genesets$gene, f = genesets$ont)
  params <- list(pathways = genesets,
                 stats    = ranks,
                 minSize  = 15,
                 maxSize  = 500,
                 gseaParam = 1)

  if (!is.null(nperm))
    params$nperm <- nperm

  res <- do.call(fgsea::fgsea, params)

  return(list(ranks = ranks, genesets = genesets, res = res))
}

fgsea.top.genesets <- function(fgsea.res, n = 10, direction = "up", padj.cutoff = NULL) {
  if (!is.null(padj.cutoff)) {
    fgsea.res <- fgsea.res %>% filter(padj < padj.cutoff)
  }

  if (direction == "up") {
    fgsea.res <- fgsea.res %>% filter(ES > 0)
    pathways <- fgsea.res %>% arrange(desc(NES)) %>% pull(pathway) %>% .[1:n]
  } else if (direction == "down") {
    fgsea.res <- fgsea.res %>% filter(ES < 0)
    pathways <- fgsea.res %>% arrange(NES) %>% pull(pathway) %>% .[1:n]
  } else if (direction == "both"){
    pathways <- fgsea.res %>% arrange(desc(abs(NES))) %>% pull(pathway) %>% .[1:n]
  } else {
    stop(paste0("direction has to be up or down"))
  }
  return(pathways)
}

fgsea.plot.top.pathways <- function(fg, # returned by fgsea.fanc()
                                    n = 20, out.dir, root.name = NULL,
                                    width = 5, height = 3) {
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  top.pathways <- lapply(c("up", "down"), function(direction) {
    fgsea.top.genesets(fgsea.res = fg$res, n = 20, direction = direction)
  })

  names(top.pathways) <- c("up", "down")

  lapply(c("up", "down"), function(direction) {
    label.style <- list(size = 6)
    p <- fgsea::plotGseaTable(pathways = fg$genesets[top.pathways[[direction]]],
                              stats = fg$ranks,
                              fgseaRes = fg$res,
                              gseaParam=1,
                              colwidths = c(7, 3, 0.8, 1.2, 1.2),
                              pathwayLabelStyle = list(size = 5),
                              headerLabelStyle = label.style,
                              valueStyle = label.style, axisLabelStyle = label.style)
    dir.create(out.dir, showWarnings = F, recursive = T)
    ggsave(paste0(out.dir, "/", root.name, "_top20_", direction,".pdf"), p,
           device = cairo_pdf, height = height, width = width)
  return()
  })
  return()
}

de.2.gsea.txt <- function(pbl, pheno.col, to.upper = T, grep.exclude = NULL,
                          out.dir, root.name = NULL) {
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  dir.create(path = out.dir, recursive = T, showWarnings = F)
  lapply(pbl, function(s2b) {
    s2b$bulkNorm <- s2b$bulkNorm %>%
      dplyr::mutate(DESCRIPTION = "na") %>%
      dplyr::rename(NAME = "gene")
    if (to.upper) {
      s2b$bulkNorm$NAME <- toupper(s2b$bulkNorm$NAME)
    }

    exp.df <- cbind(s2b$bulkNorm[, c("NAME", "DESCRIPTION")],
                s2b$bulkNorm %>% .[, ! colnames(.) %in% c("NAME", "DESCRIPTION")])
    if (!is.null(grep.exclude)) {
      exp.df <- exp.df %>%
        dplyr::filter(!grepl(paste0(grep.exclude, collapse = "|"), NAME))
    }

    write.table(exp.df, paste0(out.dir, "/", root.name, "_", s2b$root.name, "_exp.txt"),
                sep = "\t", quote = F, row.names = F, col.names = T)
    pheno <- s2b$coldata[colnames(exp.df)[-c(1,2)], pheno.col]
    cls <- paste0(out.dir, "/", root.name, "_", s2b$root.name, "_pheno.cls")
    line1 <- c(length(pheno), length(unique(pheno)), 1)
    write(line1, sep = " ", cls)
    line2 <- c("#", names(table(pheno))) %>% paste0(collapse = " ")
    write(line2, cls, sep = " ", append = T)
    write(pheno %>% paste0(collapse = " "), cls, sep = " ", append = T)
    return()
  })
  return()
}

de.2.rnk <- function(de.grid = NULL, pbl = NULL, microarray.mode = F,
                     pbl.slot = "res",
                     rank.by = "log_p", logp.max = 8,
                     to.upper = T,
                     pbl.topn = NULL, pbl.samples, pbl.clusters = NULL,
                     comp, pos.filter=NULL, neg.filter=NULL, out.dir, root.name = NULL) {
  # rank.by: log_p or log2FoldChange
  # pbl: pseudobulk list
  # pbl.samples is really just used to figure out what are the columns corresponding to samples when
  #taking topn
  system(paste0("mkdir -p ", out.dir))
  if (!is.null(de.grid)) {
    comparison <- comp
    df <- de.grid$grid.sum %>% filter(comp == comparison) %>%
      mutate(log_p = -1 * pmin(utilsFanc::log2pm(p, base = 10), logp.max) * (logFC/abs(logFC))) %>%
      filter(log_p != 0) %>% dplyr::select(gene, log_p, master.ident)
    if (to.upper) {
      df$gene <- toupper(df$gene)
    }
  } else {
    if (is.character(pbl))
      pbl <- readRDS(pbl)
    if (!microarray.mode)
      pbl <- deseq2.summary(pbl = pbl, padj.cutoff = 0.05, force = T)
    # only used to remove invalid elements in the list
    if (!is.null(pbl.clusters)) {
      pbl <- pbl[names(pbl) %in% pbl.clusters]
    }
    df <- lapply(names(pbl), function(name) {
      s2b <- pbl[[name]]
      if (!pbl.slot %in% names(s2b)) {
        stop("!pbl.slot %in% names(s2b)")
      }
      df <- s2b[[pbl.slot]] %>% as.data.frame()
      if (is.null(df$gene))
        df <- df %>% dplyr::mutate(., gene =rownames(.))

      if (!is.null(pbl.topn)) {
        if (is.null(pbl.samples)) {
          pbl.samples <- colnames(s2b$bulkNorm) %>% .[. != "gene"]
        }
        top.genes <- lapply(pbl.samples, function(sample) {
          res.exp <- ifelse(microarray.mode, "bulkNorm", "res.exp")
          genes <- s2b[[res.exp]]$gene[rank(-1 * s2b[[res.exp]][, sample], ties.method = "random") <= pbl.topn]
          return(genes)
        }) %>% Reduce(union, .)
        df <- df %>% filter(gene %in% top.genes)
      }
      df <-  df %>% dplyr::select(gene, pvalue, log2FoldChange) %>%
        dplyr::mutate(log_p = pmin(-log10(pvalue), 8) * (log2FoldChange/abs(log2FoldChange)), gene = toupper(gene), master.ident = name) %>%
        dplyr::filter(log_p != 0) %>% dplyr::select(gene, log_p, log2FoldChange, master.ident)
      return(df)
    }) %>% Reduce(rbind, .)
  }
  if (nrow(df) < 1) {
    stop ("nrow(df) < 1")
  }
  rnk.list <- df %>% split(f = df$master.ident) %>%
    lapply(function(x) {
      cluster <- x$master.ident[1]
      if (rank.by == "log_p")
        x <- x %>% dplyr::select(gene, log_p) %>% arrange(desc(log_p))
      else if (rank.by == "log2FoldChange")
        x <- x %>% dplyr::select(gene, log2FoldChange) %>% arrange(desc(log2FoldChange))
      else
        stop("rank.by has to be log_p or log2FoldChange")

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

gmt.gen <- function(gene.list = NULL, # start from raw gene list
                     genes.df = NULL, gmt.vec.in = GMT.FILES["msigdb.v7.2"], geneset.names = NULL,
                    # subset a gmt file. regex enabled, genes.df: gmt already read in.
                    to.upper = T,
                    out.dir, gmt.name) {
  # gene.list format: list(type1 = c("miao", "wang"), type2 = c("aa", "huhu"))
  if (is.null(gene.list)) {
    if (is.null(genes.df)) {
      genes.df <- lapply(gmt.vec.in, function(gmt.in) {
        df <- clusterProfiler::read.gmt(gmtfile = gmt.in) %>%
          factor2character.fanc()
        return(df)
      }) %>% Reduce(rbind, .)
    }
    if (!is.null(geneset.names)) {
      geneset.names <- tolower(geneset.names)
      genes.df <- genes.df %>% mutate(key = tolower(ont)) %>%
        filter(grepl(paste0(geneset.names, collapse = "|"), key))
    }
    gene.list <- genes.df$gene %>% split(f = genes.df$ont)
  }
  if (is.null(names(gene.list)))
    stop("gene.list must be named")
  outs <- lapply(names(gene.list), function(name) {
    genes <- gene.list[[name]]
    if (to.upper) {
      genes <- genes %>% toupper()
    }
    out <- c(name, "http://miao.html", genes)
    out <- paste0(out, collapse = "\t")
    return(out)
  }) %>% unlist()
  out.gmt <- paste0(out.dir, "/", gmt.name, ".gmt")
  system(paste0("mkdir -p ", dirname(out.gmt)))
  write(outs, out.gmt, sep = "\n")
  names(out.gmt) <- gmt.name
  res <- list(gmt = out.gmt, stats = list(genesets = names(gene.list)))
  return(res)
}

de.2.gmt <- function(de, slot = "summary", directions = c("up", "down"),
                     to.upper = T,
                     out.dir, gmt.name) {
  gene.list <- lapply(de, function(s2b) {
    summ <- s2b[[slot]]
    if (is.null(summ)) {
      stop("slot not found")
    }

    gene.list <- lapply(directions, function(direction) {
      genes <- summ[[paste0(direction, ".genes")]]
      if (length(genes) < 1) {
        return()
      }
      res <- list()
      res[[paste0(s2b$root.name, "_", direction)]] <- genes
      return(res)
    }) %>% Reduce(c, .)
    return(gene.list)
  }) %>% Reduce(c, .)
  if (is.null(gene.list)) {
    stop("is.null(gene.list)")
  }
  gmt.gen(gene.list = gene.list, to.upper = to.upper,
          out.dir = out.dir, gmt.name = gmt.name)
}

de.reciprocal.gsea <- function(des, slot = "summary", top.n = 12000,
                               out.dir, directions = c("up", "down"),
                               fgsea.plot = T, fgsea.width = 1.5, fgsea.height = 1,
                               run = F, force = F, n.par.de = 1) {
  # des <- list(my = de1, geo = de2)
  # by default reciprocal: each de in des will be used to do gsea against the
  # DEGs of each of the other de's.
  if (length(des) < 2) {
    stop("length(des) < 2")
  }
  if (is.null(names(des))) {
    stop("des must be named")
  }
  utilsFanc::check.dups(names(des), "names(des)")

  if (fgsea.plot) {
    required.fields <- c("gene.set.name.display.up", "gene.set.name.display.down",
                         "rank.name.display",
                         "up.name.display", "down.name.display")

    lapply(des, function(de) {
      lapply(de, function(s2b) {
        if (is.null(s2b$fgsea.meta)) {
          stop(paste0("fgsea.meta field not found for ", s2b$root.name))
        }
        utilsFanc::check.intersect(required.fields, "required.fields",
                                   names(s2b$fgsea.meta), "names(s2b$fgsea.meta)")
      })
    })

  }

  utilsFanc::safelapply(names(des), function(query.name) {
    print(paste0("Query: ", query.name))
    de.query <- des[[query.name]]
    des.ref <- des[names(des) != query.name ]
    out.dir <- paste0(out.dir, "/query_", query.name, "/")
    gmts <- lapply(names(des.ref), function(ref.name) {
      de.2.gmt(de = des.ref[[ref.name]], slot = slot,
               directions = directions, out.dir = paste0(out.dir, "/gmts/"),
               gmt.name = ref.name)$gmt
    }) %>% unlist()
    names(gmts) <- names(des.ref)
    rnk.vec <- de.2.rnk(pbl = de.query, out.dir = paste0(out.dir, "/rnk/"),
                        pbl.topn = top.n, pbl.samples = NULL, microarray.mode = T)
    # use the plot engine from fgsea:
    if (fgsea.plot) {
      lapply(rnk.vec, function(rnk) {
        query.cluster <- basename(rnk) %>% sub("cluster_", "", .) %>%
          sub(".rnk$", "", .)
        s2b.query <- de.query[[query.cluster]]

        lapply(gmts, function(gmt) {
          ref.de.name <- basename(gmt) %>% sub(".gmt", "", .)
          de.ref <- des[[ref.de.name]]
          lapply(de.ref, function(s2b.ref) {
            lapply(directions, function(direction) {
              gmt.gene.set.name <- paste0(s2b.ref$root.name, "_", direction)
              if (length(gmt.get.gs(gmt = gmt, gs.names = gmt.gene.set.name)) < 1) {
                print(paste0("Skipping: ", gmt.gene.set.name, "; not found in gmt"))
                return()
              }
              plot.out <- paste0(
                out.dir, "/plot_curve/",
                "q",query.name, "_", query.cluster,
                "_r", ref.de.name, "_",
                s2b.ref$root.name, "_", direction, ".pdf")
              gsea.plot.curve(gmt.file = gmt,
                              gmt.gene.set.name = gmt.gene.set.name,
                              gene.set.name.display = s2b.ref$fgsea.meta[[paste0("gene.set.name.display.", direction)]],
                              rank.name.display = s2b.query$fgsea.meta$rank.name.display,
                              up.name.display = s2b.query$fgsea.meta$up.name.display,
                              down.name.display = s2b.query$fgsea.meta$down.name.display,
                              rnk = rnk,
                              plot.out = plot.out,
                              width = fgsea.width, height = fgsea.height, calculate.p = T)
            })
          })
        })
      })
    }

    gsea.fanc(rnk.vec = rnk.vec, gmt.vec = gmts,
              out.dir = paste0(out.dir, "/gsea/"),
              thread.rnk = min(3, length(rnk.vec)), thread.gmt = min(3, length(gmts)),
              n.plot = 80, run = run, force = force, parse = F, parse.only = F, zip = F)



    return()
  }, threads = n.par.de)
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

### 2022-12-24 GSEA code:
gsea.get.result <- function(tsv, top.n = NULL, translate = T) {
  tsv <- Sys.glob(tsv)
  if (length(tsv) != 1) {
    stop("length(tsv) != 1")
  }
  df <- read.table(tsv, header = T, sep = "\t")
  bDown <- ifelse(df$CORE.ENRICHMENT[1] == "Yes", F, T)
  genes <- df %>% dplyr::filter(CORE.ENRICHMENT == "Yes") %>%
    dplyr::pull(SYMBOL)
  if (!is.null(top.n)) {
    if (bDown) {
      genes <- rev(genes)
    }
    genes <- genes[1:min(length(genes), top.n)]
  }
  if (translate) {
    genes <- stringr::str_to_title(genes)
  }
  return(genes)
}
### END: 2022-12-24 GSEA code;

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

de.enrich.da <- function(de, da, gtf.gr, da.fc.filter = 1, plot.out = NULL) {
  stats <- lapply(de, function(s2b) {
    a2b <- da[[s2b$root.name]]
    s2b$summary$non_de.genes <- s2b$res.exp$gene %>% .[!.%in% s2b$summary$de.genes]
    s2b$summary$n.non_de <- length(s2b$summary$non_de.genes)
    stats <- lapply(c("up", "down", "non_de"), function(type) {

      genes <- s2b$summary[[paste0(type, ".genes")]]
      # note: gene_type == "protein_coding" is not good enough.
      # some coding genes could still have non-coding transcripts.
      # Should be transcript_type == "protein_coding"
      tss.gr <- gtf.gr %>% plyranges::filter(gene_name %in% genes, type == "transcript",
                                             gene_type == "protein_coding") %>%
        resize(width = 1, fix = "start", ignore.strand = F) %>% as.data.frame() %>%
        dplyr::rename(chr = seqnames, gene = gene_name) %>%
        dplyr::select(chr, start, end, gene) %>% unique() %>% makeGRangesFromDataFrame(keep.extra.columns = T)
      mcols(tss.gr) <- dplyr::left_join(as.data.frame(mcols(tss.gr)),
                                        s2b$res.exp[, c("gene", "log2FoldChange", "pvalue", "padj")])
      da.gr <- a2b$res.exp %>% dplyr::rename(locus = gene) %>%
        .[, c("locus", "log2FoldChange", "pvalue", "padj")] %>%
        utilsFanc::loci.2.df(loci.col.name = "locus", remove.loci.col = F, return.gr = T)
      tss.da <- plyranges::join_overlap_left(tss.gr, da.gr)
      tss.da <- tss.da[!is.na(tss.da$locus)]
      tss.da.df <- as.data.frame(tss.da)
      tss.da.df <- tss.da.df %>%
        dplyr::mutate(bSameTrend = (log2FoldChange.x * log2FoldChange.y > 0))
      if (!is.null(da.fc.filter)) {
        tss.da.df$bSameTrend[abs(tss.da.df$log2FoldChange.y) < da.fc.filter] <- F
      }
      tss.da.sum <- tss.da.df %>%
        dplyr::group_by(gene) %>% dplyr::summarise(bSameTrend = sum(bSameTrend) > 0) %>%
        dplyr::ungroup() %>% as.data.frame()
      stats <- data.frame(root.name = s2b$root.name, type = type, n = s2b$summary[[paste0("n.", type)]],
                          n.with.peak = nrow(tss.da.sum), n.same.trend = sum(tss.da.sum$bSameTrend))
      stats$pct.same.trend <- stats$n.same.trend/stats$n#.with.peak
      return(stats)
    }) %>% Reduce(rbind, .)
  }) %>% Reduce(rbind, .)
  p <- ggplot(stats, aes(x = root.name, y = pct.same.trend, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() + theme(text = element_text(size = 20), aspect.ratio = 1)
  trash <- wrap.plots.fanc(list(p), plot.out = plot.out, sub.width = 7, sub.height = 5)
  return(p)
}

de.enrich.da.check <- function(de, da, slot = "summary", tss.file, out.dir,
                               bedtools = "/bar/cfan/anaconda2/envs/jupyter/bin/bedtools") {
  stats <- lapply(names(de), function(name) {
    s2b <- de[[name]]
    a2b <- da[[name]]
    lapply(c("up", "down"), function(type) {
      dar <- a2b[[slot]][[paste0(type, ".genes")]] %>%
        utilsFanc::loci.2.df(loci.vec = ., remove.loci.col = T)
      dar.file <- paste0(out.dir, "/", name, "_dar_", type, ".bed")
      tss.in.dar <- paste0(out.dir, "/", name, "_dar_", type, "_tss.bed")
      utilsFanc::write.zip.fanc(df = dar, out.file = dar.file, bed.shift = T)
      cmd <- paste0(bedtools, " intersect -a ", tss.file, " -b ", dar.file, " -wa > ", tss.in.dar)
      system(cmd)
      genes.dar <- read.table(tss.in.dar)$V4 %>% unique()
      deg <- s2b[[slot]][[paste0(type, ".genes")]]
      intsct <- intersect(genes.dar, deg)
      res <- data.frame(cluster = name, direction = type,
                        n.intsct = length(intsct), n.deg = length(deg),
                        frac.deg = length(intsct)/length(deg),
                        intsct = paste0(intsct, collapse = ","),
                        non.intsct = paste0(deg[!deg %in% genes.dar], collapse = ","),
                        deg = paste0(deg, collapse = ","),
                        genes.dar = paste0(genes.dar, collapse = ","))
      return(res)
    }) %>% do.call(rbind, .) %>% return()
  }) %>% do.call(rbind, .)
  write.table(stats, paste0(out.dir, "/stats.tsv"), quote = F,
              sep = "\t", row.names = F, col.names = T)
  return(stats)
}

msigdb.m <- function(pbl, clusters = NULL, summary.slot = "summary",
                     universe = "expressed",
                     use.n.genes = NULL, seed = NULL,
                     cats = "C5..GO:BP",
                     species = "Mus musculus",
                     threads = 1) {
  # universe: expressed. anything else: consider as missing.
  for (cat in cats) {
    sub.cats <- strsplit(cat, split = "@") %>% unlist()
    go_df <- lapply(sub.cats, function(sub.cat) {
      category <- sub("\\.\\..*$", "", sub.cat)
      subcategory <- sub("^.+\\.\\.", "", sub.cat)
      if (subcategory == "") {
        subcategory <- NULL
      }
      go_df <- msigdbr::msigdbr(species = species, category = category,
                                subcategory = subcategory) %>%
        dplyr::select(gs_name, human_gene_symbol) %>%
        dplyr::rename(gene.set = gs_name, gene.name = human_gene_symbol)
      return(go_df)
    }) %>% Reduce(rbind, .)

    if (species == "Mus musculus") {
      go_df <- go_df %>% dplyr::mutate(gene.name = stringr::str_to_title(gene.name))
    } else if (species != "Homo sapiens") {
      warning("Didn't convert the case of gene names while species is not human or mouse")
    }
    
    # go_size <- go_df %>% dplyr::group_by(gene.set) %>% dplyr::summarise(n = n()) %>% 
    #   dplyr::ungroup() %>% as.data.frame()
    # names(go_size) <- c("ID", "GS.size")
    
    if (is.null(clusters)) {
      clusters <- names(pbl)
    }
    pbl <- utilsFanc::safelapply(pbl, function(s2b) {
      if (!s2b$root.name %in% clusters) {
        return(s2b)
      }
      if (is.null(s2b$msigdb)) {
        s2b$msigdb <- list()
      }

      enrich <- lapply(c("up", "down"), function(type) {
        genes <- s2b[[summary.slot]][[paste0(type, ".genes")]]
        if (!is.null(use.n.genes) && length(genes) > use.n.genes) {
          set.seed(seed = seed)
          genes <- sample(genes, size = use.n.genes, replace = F)
        }
        if (universe[1] == "expressed") {
          print("using all genes in s2b$bulkNorm$gene")
          res <- clusterProfiler::enricher(genes, TERM2GENE = go_df, universe = s2b$bulkNorm$gene)
        } else if (universe[1] == "all_mouse") {
          print("using all mouse protein coding genes")
          res <- clusterProfiler::enricher(genes, TERM2GENE = go_df,
                                           universe = readLines("~/genomes/mm10/gencode/gencode.vM24.protein.coding.txt"))
        } else if(universe[1] == "default") {
          print("using default universe")
          res <- clusterProfiler::enricher(genes, TERM2GENE = go_df)
        } else {
          res <- clusterProfiler::enricher(genes, TERM2GENE = go_df, universe = universe)
        }
        return(res)
      })
      names(enrich) <- c("up", "down")
      cat.name <- cat
      # if (!is.null(use.n.genes)) {
      #   cat.name <- paste0(cat, "_rand", use.n.genes, "_se", seed)
      # }
      # if (universe != "expressed") {
      #   cat.name <- paste0(cat.name, "_gs")
      # }
      s2b$msigdb[[cat.name]] <- enrich
      return(s2b)
    }, threads = threads)
  }
  #  cat <- paste(category, subcategory, sep = "..")
  return(pbl)
}

msigdb.sum <- function(pbl, cats = NULL, gene.sets = NULL, slot = "msigdb", out.file = NULL) {
  df <- lapply(pbl, function(s2b) {
    if (is.null(s2b[[slot]])) {
      return()
    }
    cats.use <- names(s2b[[slot]])
    if (!is.null(cats)) {
      cats.use <- cats.use %>% .[grepl(paste0(cats, collapse = "|"),.)]
    }
    if (length(cats.use) < 1) {
      return()
    }
    s2b[[slot]] <- s2b[[slot]][cats.use]
    lapply(names(s2b[[slot]]), function(enrich.name) {
      lapply(c("up", "down"), function(type) {
        # print(s2b$root.name)
        # print(enrich.name)
        # print(type)
        enrich.o <- s2b[[slot]][[enrich.name]][[type]]
        if (is.null(enrich.o)) {
          return()
        }

        if (!is.null(gene.sets)) {
          df <- enrich.o@result[gene.sets,] %>% na.omit()
        } else {
          df <- head(enrich.o, n = nrow(enrich.o))
        }
        
        sizes <- enrich.o@geneSets %>% sapply(length)
        size.df <- data.frame(ID = names(sizes), GS.size = sizes)
        df <- dplyr::left_join(df, size.df, by = "ID")
        rownames(df) <- NULL
        
        df <- utilsFanc::add.column.fanc(
          df1 = df,
          df2 = data.frame(cluster = s2b$root.name, cat = enrich.name, type = type),
          pos = 1)
        return(df)
      }) %>% Reduce(rbind, .) %>% return()

    }) %>% Reduce(rbind, .) %>% return()

  }) %>% Reduce(rbind, .) %>% return()
  if (!is.null(out.file)) {
    df$GeneRatio <- paste0("# ", df$GeneRatio)
    df$BgRatio <- paste0("# ", df$BgRatio)
    utilsFanc::write.zip.fanc(df = df, out.file = out.file, zip = F, col.names = T, row.names = F)
  }
  return(df)
}

# msigdb.test.geneset <- function(pbl, gene.set, slot = "msigdb", cats = NULL, out.file = NULL) {
#   df <- lapply(pbl, function(s2b) {
#     lapply(c("up", "down"), function(type) {
#       if (is.null(s2b[[slot]])) {
#         return()
#       }
#       cats.use <- names(s2b[[slot]])
#       if (!is.null(cats)) {
#         cats.use <- cats.use %>% .[grepl(paste0(cats, collapse = "|"),.)]
#       }
#       if (length(cats.use) < 1) {
#         return()
#       }
#
#       lapply(cat.use, function(enrich.name)) {
#         enrich.o <- s2b[[slot]][[enrich.name]][[type]]
#       }
#     })
#   })
# }

msigdb.dotplot <- function(pbl, cat.to.plot,
                           n = 15, font.size = 7, wrap.width = 40,
                           plot.out = NULL, sub.width = 6, sub.height = 4,
                           threads = 1) {
  pl <- utilsFanc::safelapply(pbl, function(s2b) {
    pl <- lapply(c("up", "down"), function(type) {
      enrich <- s2b$msigdb[[cat.to.plot]][[type]]
      if (is.null(enrich)) {
        stop("is.null(enrich)")
      }
      df <- head(enrich)
      if (!is.null(wrap.width)) {
        enrich@result$ID <- enrich@result$ID %>%
          gsub("_"," ", . ) %>% stringr::str_wrap(width = wrap.width)
        enrich@result$Description <- enrich@result$ID
        rownames(enrich@result) <- enrich@result$ID
      }

      p <- clusterProfiler::dotplot(enrich, showCategory = n,
                   font.size = font.size,
                   title = paste0(s2b$root.name, " ", type))
      return(p)
    })
    return(pl)
  }) %>% Reduce(c, .)
  p <- wrap.plots.fanc(pl, plot.out = plot.out,
                       sub.width = sub.width, sub.height = sub.height)
  invisible(p)
}

# msigdb.hm <- function(summary.tsv, clusters = NULL, directions = c("up", "down"), remove.up.down = F,
#                       geneset.size.min = 10, geneset.size.max = 500, 
#                       fdr.cutoff, gene.ratio.cutoff,
#                       top.n, rank.by, 
#                       color.by = NULL, color.map,
#                       plot.out, width = 2, height = 4) {
#   stop("Did not test for lack of priority")
#   summ <- read.table(summary.tsv, header = T, sep = "\t", comment.char = "")
#   required.cols <- c("cluster", "type", "ID", "GeneRatio", "p.adjust", "GS.size")
#   
#   utilsFanc::check.intersect(required.cols, "required columns", colnames(summ), "colnames(summ)")
#   
#   if (!is.null(clusters)) {
#     utilsFanc::check.intersect(clusters, "clusters", summ$cluster, "summ$cluster")
#   } else {
#     clusters <- summ$cluster %>% unique()
#   }
#   utilsFanc::check.intersect(directions, "directions", summ$type, "summ$type")
#   
#   summ <- summ %>% dplyr::filter(cluster %in% clusters, type %in% directions)
#   
#   summ <- summ %>% dplyr::filter(
#     GS.size > geneset.size.min, GS.size < geneset.size.max,
#     p.adjust < fdr.cutoff, GeneRatio > gene.ratio.cutoff)
#   
#   if (rank.by %in% c("GeneRatio")) {
#     summ <- summ %>% dplyr::group_by(cluster, type) %>% 
#       dplyr::top_n(top.n, !!as.name(rank.by)) %>% 
#       dplyr::ungroup() %>% as.data.frame()
#   } else if (rank.by %in% c("p.adjust")) {
#     summ <- summ %>% dplyr::group_by(cluster,type) %>% 
#       dplyr::top_n(top.n, -1 * !!as.name(rank.by)) %>% 
#       dplyr::ungroup() %>% as.data.frame()
#   }
#   
#   if (remove.up.down && length(directions) == 1) {
#     summ$col <- summ$cluster
#   } else {
#     summ$col <- paste0(summ$cluster, "_", summ$type)
#   }
#   if (is.null(color.by))
#     color.by <- rank.by
#   
#   if (nrow(summ) < 1) stop("nrow(summ) < 1")
#   
#   mat <- reshape2::dcast(summ, cat ~ col, value.var = color.by)
#   
#   rds <- tools::file_path_sans_ext(plot.out) %>% paste0("_mat.Rds")
#   saveRDS(mat, rds)
#   
#   browser()
#   
# }

msigdb.hm <- function(de, pathways, slot = "msigdb", color.by = c("GeneRatio"),
                      cat, direction, padj.cutoff = 0.05,
                      cluster.order = NULL, add.FDR = T,
                      width = 3, height = 5.5, hm.values = NULL,
                      wrap.width = 40, strip.first.field = T,
                      plot.out) {
  # you could also color by enrich.
  utilsFanc::check.dups(x = pathways, "pathways")
  if (!is.null(cluster.order)) {
    utilsFanc::check.intersect(cluster.order, "cluster.order", names(de), "names(de)")
  } else {
    cluster.order <- names(de)
  }
  dir.create(dirname(plot.out), showWarnings = F, recursive = T)
  mats <- lapply(c("GeneRatio", "BgRatio", "p.adjust"), function(mat.type) {
    mat <- lapply(de, function(s2b) {
      results <- s2b[[slot]][[cat]][[direction]]
      if (is.null(results)) {
        stop("is.null(results)")
      }
      results <- results@result
      # results only contain those that have at least 1 gene in the geneset
      # so if geneRatio is 0, you won't find it in the object.
      # utilsFanc::check.intersect(pathways, "pathways", names(results@geneSets), "rownames(results)")
      res <- data.frame(pathway = pathways, x = results[pathways, mat.type])
      colnames(res) <- c("pathway", s2b$root.name)
      return(res)
    }) %>% Reduce(dplyr::left_join, .)
    rownames(mat) <- mat$pathway
    mat$pathway <- NULL
    mat <- as.matrix(mat)
    
    if (mat.type %in% c("GeneRatio", "BgRatio")) {
      mat <- apply(mat, c(1,2), function(x) eval(parse(text = x)))
      mat[is.na(mat)] <- 0
    } else {
      mat[is.na(mat)] <- 1
    }
    mat <- mat[, cluster.order]
    return(mat)
  })
  names(mats) <- c("GeneRatio", "BgRatio" , "p.adjust")
  
  mats$enrich <- mats$GeneRatio/mats$BgRatio
  mats$enrich[is.nan(mats$enrich)] <- 0
  
  saveRDS(mats, paste0(tools::file_path_sans_ext(plot.out), ".Rds"))
  
  mats <- lapply(mats, function(mat) {
    if (strip.first.field) {
      rownames(mat) <- rownames(mat) %>% sub("^[^_]+_"," ", . )
    }
    rownames(mat) <- rownames(mat) %>% gsub("_"," ", . ) %>% 
      stringr::str_to_title() %>% 
      stringr::str_wrap(width = wrap.width) %>% 
      gsub("\n", "\n  ", .)

    return(mat)
  })
  
  symbol.mat <- NULL
  if (add.FDR) {
    symbol.mat <- mats$p.adjust
    symbol.mat[symbol.mat < padj.cutoff] <- "*"
    tmp <- dimnames(symbol.mat)
    symbol.mat <- stringr::str_extract(symbol.mat, "\\*+") %>%
      matrix(nrow = length(tmp[[1]]))
    dimnames(symbol.mat) <- tmp
    
    symbol.mat[is.na(symbol.mat)] <- ""
    
  }
  
  if (is.null(hm.values)) {
    hm.values <- c(0, quantile(mats[[color.by]], 0.95))
  }
  plot.mat.rank.row(mat = mats[[color.by]], symbol.mat = symbol.mat,
                    no.col.cluster = T, show_column_names = T,
                    hm.colors = c("white", "deeppink4"), hm.values = hm.values,
                    show_row_names = T, 
                    plot.out = plot.out, width = width, height = height,
                    row_name_fontSize = 5)
  invisible(mats)
  
}

homer.plot.heatmap <- function(homer.txt.df, 
                               motif.names = NULL, 
                               FDR.cutoff = 0.05, pct.fg.cutoff = 0, pct.bg.cutoff = 100, enrich.cutoff = 0,
                               column.color = "pct.fg",
                               cluster.columns = F, hm.column.order = NULL,
                               cluster.rows = F,
                               hm.colors = c("white", "deeppink4"), hm.values = NULL,
                               add.FDR = T,
                               plot.out, width = 2, height = 2,
                               ...) {
  # homer.txt.df <- data.frame(name = cluster.1, path = /path/to/knownResults.txt)
  # if motif.names is not offered, we will get motifs based on unbiased filters.
  
  columns <- c("motif", "consensus", "p.value", "log10p", "FDR", "n.fg", "pct.fg", "n.bg", "pct.bg", "enrich")
  utilsFanc::check.intersect(column.color, "column.color", columns, "available columns")
  
  if (is.character(homer.txt.df)) {
    homer.txt.df <- read.table(homer.txt.df, header = T)
  }
  
  res <- lapply(1:nrow(homer.txt.df), function(i) {
    df <- read.table(homer.txt.df[i, "path"], header = F, sep = "\t", quote = "", skip = 1)

    if (length(columns) - 1 != ncol(df))
      stop("length(columns) -1 != ncol(df)")
    
    colnames(df) <- columns[-length(columns)]
    
    df$pct.fg <- sub("\\%$", "", df$pct.fg) %>% as.numeric()
    df$pct.bg <- sub("\\%$", "", df$pct.bg) %>% as.numeric()
    
    df$enrich <- round(df$pct.fg/df$pct.bg, digits = 2)
    
    df <- cbind(data.frame(name = rep(homer.txt.df[i, "name"], n = nrow(df))),
                df)
    # res <- data.frame(name = homer.txt.df[i, "name"], motif = motif.names) %>% 
    #   dplyr::left_join(df, by = "motif") %>% dplyr::select(name, motif, !!as.name(column.color), FDR)
    return(df)
  }) %>% Reduce(rbind, .)
  
  if (is.null(motif.names)) {
    motif.names <- res %>% 
      dplyr::filter(FDR < FDR.cutoff, pct.fg >= pct.fg.cutoff,
                    pct.bg <= pct.bg.cutoff, enrich >= enrich.cutoff) %>% 
      dplyr::pull(motif) %>% unique()
    
    if (length(motif.names) < 1) {
      stop("length(motif.names) < 1")
    }
    
    black.list <- "p53(p53)/Saos-p53-ChIP-Seq(GSE15780)/Homer"
    # this is a duplicate of the motif "p53(p53)/Saos-p53-ChIP-Seq/Homer"
    motif.names <- motif.names %>% .[!.%in% black.list]
    
    suffix <- paste0("colorBy_", column.color, "_FDR", FDR.cutoff, '_pct.fg', 
                     pct.fg.cutoff, "_pct.bg", pct.bg.cutoff,
                     "_enrich", enrich.cutoff)
    plot.out <- utilsFanc::insert.name.before.ext(name = plot.out, insert = suffix, delim = "_")
  }
  
  utilsFanc::check.intersect(motif.names, "motif.names", res$motif, "res$motif")
  res <- res %>% dplyr::filter(motif %in% motif.names)
  # utilsFanc::check.dups(paste0(res$name, "_", res$motif), x.name = "paste0(res$name, res$motif)")
  # some of these are inherently duplicated, for example RORgt(NR)
  
  # colnames(res) <- c("name", "motif", "col", "FDR")
  res$motif.ori <- res$motif
  res$motif <- res$motif %>% sub("/.+$", "", .)
  
  summ.fun <- ifelse(column.color %in% c("fg.pct", "enrich"), max, min)
  
  tmp.df <- res %>% dplyr::group_by(name, motif) %>% 
    dplyr::summarise(!!as.name(column.color) := summ.fun(!!as.name(column.color))) %>% 
    dplyr::ungroup() %>% as.data.frame()  
  
  mat <- reshape2::acast(tmp.df, motif ~ name, value.var = column.color)
  
  if (!is.null(hm.column.order)) {
    utilsFanc::check.intersect(hm.column.order, "hm.column.order", colnames(mat), "colnames(mat)")
    mat <- mat[, hm.column.order]
  }
  
  
  if (add.FDR) {
    tmp.df <- res %>% dplyr::group_by(name, motif) %>% dplyr::summarise(FDR = min(FDR)) %>% 
      dplyr::ungroup() %>% as.data.frame()
    symbol.mat <- reshape2::acast(tmp.df, motif ~ name, value.var = "FDR")
    if (!is.null(hm.column.order)) {
      symbol.mat <- symbol.mat[, hm.column.order]
    }
    symbol.mat[symbol.mat < FDR.cutoff] <- "*"
    tmp <- dimnames(symbol.mat)
    symbol.mat <- stringr::str_extract(symbol.mat, "\\*+") %>%
      matrix(nrow = length(tmp[[1]]))
    dimnames(symbol.mat) <- tmp
    
    symbol.mat[is.na(symbol.mat)] <- ""
    
  } else {
    symbol.mat <- NULL
  }
  if (is.null(hm.values)) {
    hm.values <- c(0, max(mat))
  }
  plot.mat.rank.row(mat = mat, symbol.mat = symbol.mat, 
                    cluster_rows = cluster.rows, show_column_dend = cluster.rows,
                    no.col.cluster = !cluster.columns, show_column_names = T,
                    hm.colors = hm.colors, hm.values = hm.values,
                    show_row_names = T, 
                    plot.out = plot.out, width = width, height = height,
                    row_name_fontSize = 5, ...)
  
  write.table(res, paste0(tools::file_path_sans_ext(plot.out), ".tsv"), sep = "\t",
              row.names = F, col.names = T, quote = F)
  
  invisible(res)
  
}

homer.plot.bar <- function(da.homer.dir = NULL, homer.txt.df,
                           motif.map, use.regex = F, col.num = 3,
                           flip.up.down = F, plot.out, plot.logo = F,
                           width = 3, height = 3, font.size = 12,
                           logo.width = 200, logo.height = 300,
                           ...) {
  if (!is.null(da.homer.dir)) {
    homer.txt.df <- deseq2.homer.gather(da.homer.dir = da.homer.dir, ...)
  }
  if (is.character(homer.txt.df)) {
    homer.txt.df <- read.table(homer.txt.df, header = T)
  }
  if (file.exists(motif.map[1])) {
    motif.map <- read.table(motif.map, header = T, sep = "\t", quote = "")
  }
  j <- left_join(motif.map, homer.txt.df)
  j$root <- j$name %>% sub("_up|_down", "", .)
  j$type <- j$name %>% stringr::str_extract("up|down")
  pl <- split(j, f = j$root) %>%
    lapply(function(df) {
      df$pvalue <- homer.core.extract(path = df$path, motifs = df$motif,
                                      use.regex = use.regex, col.num = col.num)

      if (plot.logo) {
        logos <- homer.core.extract(paths = df$path, motifs = df$motif,
                                      use.regex = use.regex, get.logo.path = T)
        logos <- logos[order(df$pvalue)]
        if (any(!file.exists(logos))) {
          stop(paste0("these files do not exist: ",
                      paste0(logos[!file.exists(logos)]), collapse = "\n"))
        }
        logos <- lapply(logos, rsvg::rsvg)
        widths <- sapply(logos, dim)[2, ]
        widths <- round(widths/max(widths), 2)
        dir.create(dirname(plot.out), showWarnings = F, recursive = T)
        png(utilsFanc::insert.name.before.ext(plot.out, paste0("logo_", df$root[1]), delim = "_"),
            width = logo.width, height = logo.height, res = 100)
        try({
          par(mar=rep(0,4))
          layout(matrix(1:length(logos), ncol=1, byrow=TRUE))
          for(i in 1:length(logos)) {
            print(plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n"))
            print(rasterImage(logos[[i]],0,0, widths[i],1))
          }
        })
        dev.off()
      }

      df$motif.ez <- df$motif %>% sub("/.+$", "", .) %>% sub(",", "\n", .) %>%
        factor(., levels = .[rev(order(df$pvalue))])
      if (flip.up.down) {
        df$type[df$type == "up"] <- "miao"
        df$type[df$type == "down"] <- "up"
        df$type[df$type == "miao"] <- "down"
      }
      p <- ggplot(df, aes(x = motif.ez, y = -log10(pvalue), fill = type)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_classic() +
        ggtitle(df$root[1]) +
        theme(text = element_text(size = font.size))

      return(p)
    })
  p <- wrap.plots.fanc(plot.list = pl, plot.out = plot.out,
                  sub.height = height, sub.width = width)
  invisible(p)
}

homer.core.extract <- function(paths, motifs, use.regex = F, col.num = 3,
                               get.logo.path = F) {
  # path <- j$path[1]
  if (length(paths) == 1) {
    paths <- rep(paths, length(motifs))
  }
  if (length(paths) != length(motifs)) {
    stop("length(paths) != length(motifs)")
  }
  value <- sapply(seq_along(motifs), function(i) {
    motif <- motifs[i]
    path <- paths[i]
    df <- read.table(file = path, header = F, sep = "\t", quote = "", skip = 1)
    if (get.logo.path) {
      if (use.regex) {
        id <- which(grepl(motif, df$V1))
      } else {
        id <- which(df$V1 == motif)
      }
      value <- paste0(dirname(path), "/knownResults/known", id, ".logo.svg")
    } else {
      if (use.regex) {
        value <- df %>% .[grepl(motif, .$V1), col.num]
      } else {
        value <- df %>% .[.$V1 == motif, col.num]
      }
    }
    if (length(value) > 1 ) {
      stop(paste0("error in homer.core.extract at motif ", motif, ": length(value) > 1"))
    }
    if (length(value) == 0 ) {
      stop(paste0("error in homer.core.extract at motif ", motif, ": length(value) == 0"))
    }
    return(value)
  }) %>% `names<-`(NULL)
  return(value)
}

deseq2.homer.gather <- function(da.homer.dir,
                                regex.include = NULL, regex.exclude = NULL,
                                regex.from = NULL, regex.to = "",
                                out.file = NULL) {
  results <- Sys.glob(paste0(da.homer.dir, "/*/*/homer/knownResults.txt"))
  if (!is.null(regex.include)) {
    results <- results %>% .[grepl(paste0(regex.include, collapse = "|"), .)]
  }
  if (!is.null(regex.exclude)) {
    results <- results %>% .[!grepl(paste0(regex.exclude, collapse = "|"), .)]
  }
  if (length(results) < 1) {
    stop("length(results) < 1")
  }
  name <- results %>% sub(paste0(da.homer.dir, "/"), "", .)
  name <- sub("/.+$", "", name) %>% paste0("_", stringr::str_extract(name, "up|down"))
  if (!is.null(regex.from)) {
    name <- name %>% gsub(regex.from, regex.to, .)
  }
  df <- data.frame(name = name, path = results)
  if (!is.null(out.file)) {
    dir.create(dirname(out.file), showWarnings = F, recursive = T)
    write.table(df, out.file, quote = F, sep = "\t", row.names = F, col.names = T)
  }
  return(df)
}

homer.plot.pct <- function(homer.txt, motifs,  rm.legend = F, out.dir, root = NULL) {
  if(is.null(root)) root <- basename(out.dir)
  if (is.null(names(motifs)))
    stop("motifs must be named. homer motif names are too long!")
  df <- read.table(file = homer.txt, header = F, sep = "\t", quote = "", skip = 1)
  # foreground percentage: V7; bg percentage: V9

  utilsFanc::check.intersect(motifs, "motifs", df$V1, "df$V1")

  df <- df[df$V1 %in% motifs, c("V1", "V7", "V9")]
  colnames(df) <- c("motif.longname", "fg", "bg")
  map.df <- data.frame(motif.longname = motifs,
                       motif = names(motifs))
  df <- left_join(df, map.df, by = "motif.longname")
  df$motif.longname <- NULL

  df <- reshape2::melt(df, id.vars = "motif", variable.name = "type", value.name = "pct")
  df$pct <- df$pct %>% sub("\\.", "", .) %>% sub("%", "", .) %>% paste0("0.", .) %>% as.numeric()
  df$motif <- factor(df$motif, levels = names(motifs))
  p <- ggplot(df, aes(x = motif, y = pct)) +
    geom_bar(aes(fill = type, color = type), stat = "identity", position = "dodge", alpha = 0.5) +
    scale_y_continuous(labels = scales::percent)
  p <- utilsFanc::theme.fc.1(p, text.size = 6, rm.x.ticks = F, rotate.x.45 = F, italic.x = F)
  if (rm.legend)
    p <- p + theme(legend.position = "none")
  dir.create(out.dir, showWarnings = F, recursive = T)
  file <- paste0(out.dir, "/", root, "_pct_bar.pdf")
  ggsave(file, p, device = cairo_pdf, width = 1, height = 1, dpi = 300)
  invisible(p)
}

homer.plot.table <- function(da.homer.dir = NULL, homer.txt.df,
                           motif.map, use.regex = F, col.num = 3,
                           flip.up.down = F, out.dir, plot.logo = F,
                           width = 3, height = 3, font.size = 12,
                           logo.width = 200, logo.height = 300,
                           ...) {
  if (!is.null(da.homer.dir)) {
    homer.txt.df <- deseq2.homer.gather(da.homer.dir = da.homer.dir, ...)
  }
  if (is.character(homer.txt.df)) {
    homer.txt.df <- read.table(homer.txt.df, header = T)
  }
  if (file.exists(motif.map[1])) {
    motif.map <- read.table(motif.map, header = T, sep = "\t", quote = "")
  }
  j <- left_join(motif.map, homer.txt.df)
  j$root <- j$name %>% sub("_up|_down", "", .)
  j$type <- j$name %>% stringr::str_extract("up|down")
  j$type[is.na(j$type)] <- ""
  
  pl <- split(j, f = j$root) %>%
    lapply(function(df) {
      split(df, f = df$type) %>%
        lapply(function(df) {
          direction <- df$type[1] # up or  down
          columns <- c(3, 5, 7, 9)
          names(columns) <- c("pvalue", "FDR", "pct.fg", "pct.bg")

          df.extract <- lapply(columns, function(i) {
            homer.core.extract(
              path = df$path, motifs = df$motif,
              use.regex = use.regex, col.num = i)
          }) %>% as.data.frame()

          df <- cbind(df, df.extract)

          df$motif.ez <- df$motif %>% sub("/.+$", "", .) %>% sub(",", "\n", .)
          df$motif.ez <- df$motif.ez %>% sub("\\(.+$", "", .)
          df$motif.ez <- df$motif.ez %>% factor(., levels = .[rev(order(df$pvalue))])

          df <- df %>% dplyr::arrange(pvalue)

          df$motif.ez.copy <- df$motif.ez

          df$pct.fg <- sub("%", "", df$pct.fg) %>% as.numeric()
          df$pct.bg <- sub("%", "", df$pct.bg) %>% as.numeric()
          df$enrich <- df$pct.fg / df$pct.bg

          df$FDR[df$FDR < 0.0001] <- "< 0.0001"

          # df.melt <- reshape2::melt(
          #   df, id.vars = "motif.ez.copy", measure.vars = c("motif.ez", "FDR", "pct.fg", "pct.bg", "enrich"),
          #   variable.name = "column", value.name = "value")
          #
          # map.df <- data.frame(column = c("motif.ez", "FDR", "pct.fg", "pct.bg", "enrich"),
          #                      column.pub = c("Motif", "FDR", paste0("% ", direction, "DARs with Motif"),
          #                                     "% Background with Motif", "Enrichment"))
          # df.melt <- left_join(df.melt, map.df, by = "column")
          df.useful <- df[, c("motif.ez", "FDR", "pct.fg", "pct.bg", "enrich")]
          df.useful$enrich <- round(df.useful$enrich, digits = 2)
          colnames(df.useful) <- c("Motif", "FDR", paste0("% ", direction, "DARs\nwith Motif"),
                                   "% Background\nwith Motif", "Fold\nEnrichment")
          df.useful <- df.useful[, -ncol(df.useful)]
          fake.df <- data.frame(x = 1:3, y = 1:3)

          p <- ggplot(fake.df, aes(x = x, y = y)) + geom_point(color = "white", alpha = 0) +
            scale_x_continuous(expand = expansion(0, 0)) +
            scale_y_continuous(expand = expansion(0, 0)) +
            ggpp::annotate("table", x = 1, y = 1, label = df.useful,
                           table.theme = ggpp::ttheme_gtbw(base_size = 6, base_family = "Arial")) +
            theme_void() +
            theme(plot.margin = margin(unit = "in"))

          plot.out <- paste0(out.dir, "/", df$name[1], "_table.pdf")
          dir.create(out.dir, showWarnings = F, recursive = T)
          ggsave(plot.out, p, width = 2.3, height = 0.5 + 0.1 * nrow(df), device = cairo_pdf)

          if (plot.logo) {
            logos <- homer.core.extract(paths = df$path, motifs = df$motif,
                                        use.regex = use.regex, get.logo.path = T)
            logos <- logos[order(df$pvalue)]
            if (any(!file.exists(logos))) {
              stop(paste0("these files do not exist: ",
                          paste0(logos[!file.exists(logos)]), collapse = "\n"))
            }
            logos <- lapply(logos, rsvg::rsvg)
            widths <- sapply(logos, dim)[2, ]
            widths <- round(widths/max(widths), 2)
            dir.create(dirname(plot.out), showWarnings = F, recursive = T)
            png(paste0(tools::file_path_sans_ext(plot.out), "_logo.png"),
                width = logo.width, height = logo.height, res = 100)
            try({
              par(mar=rep(0,4))
              layout(matrix(1:length(logos), ncol=1, byrow=TRUE))
              for(i in 1:length(logos)) {
                print(plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n"))
                print(rasterImage(logos[[i]],0,0, widths[i],1))
              }
            })
            dev.off()
          }
          return()
        })
    })

}


t.f.titrate.gene.enrich <- function(de.list, lfc, summ, msigdb.cats, universes, root.name = NULL, out.dir) {
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  lapply(names(de.list), function(de.name) {
    de <- de.list[[de.name]]
    lapply(universes, function(universe) {
      lapply(lfc, function(i) {
        lapply(summ, function(summ) {
          if (summ == "summary") {
            res.slot <- "res"
          } else if (summ == "sum_shrink") {
            res.slot <- "shrink_ashr"
          } else {
            stop()
          }
          de <- deseq2.summary(pbl = de, res.slot = res.slot,
                               summary.slot = paste0(summ, "_lfc", i),
                               log2fc.cutoff = i,
                               rank.by = "log2FoldChange",
                               gene.out.dir = paste0(out.dir, "/", root.name, "_deGenes/", de.name, "_", summ, "_lfc/lfc", i, "/")
          )

          de <- msigdb.m(de, summary.slot = paste0(summ, "_lfc", i),
                         cats = msigdb.cats, threads = 4, universe = universe)
          out.root <- paste0(out.dir, "/", root.name, "_", de.name, "_", universe, "_", summ, "_lfc", i)

          trash <- msigdb.sum(pbl = de,
                              out.file = paste0(out.root, ".tsv"),
                              cats = msigdb.cats )

          msigdb.dotplot(de, cat.to.plot = msigdb.cats, n = 10,
                         font.size = 10,
                         plot.out = paste0(out.root, ".png"))
          return()
        })
      })

    })
  })
}



gsea.plot.curve <- function(genes = NULL, gmt.file, gmt.gene.set.name,
                            # 3 different entrys: you can directly supply a vector of genes as the gene set.
                            # or, you can give a gmt file and specify which gene set from the gmt
                            rnk, # rnk file, you know how it works...
                            gene.set.name.display = "gene_set", # the title of the plot
                            gseaParam = 1, # the exponent for score calculation
                            calculate.p = F, add.p.to.plot = T,
                            rank.name.display = "rank", # x axis of the plot
                            up.name.display = "up", # something like KO > WT
                            down.name.display = "down", # something like KO < WT
                            plot.out = NULL, width = 1.5, height = 1
) {
  if (is.null(genes)) {
    genes <- gmt.get.gs(gmt = gmt.file, gs.names = gmt.gene.set.name)[[gmt.gene.set.name]]
  }

  df <- read.table(rnk, header = F)
  stats <- df$V2
  names(stats) <- df$V1
  pathway <- genes

  # mostly copying code from fgsea::plotEnrichment ...
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                                 returnAllExtremes = TRUE)
  if (calculate.p) {
    gsea.stats <- fgsea::fgsea(pathways = list(pathway = genes),
                          stats = stats, gseaParam = gseaParam)
  }
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  # diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  ymin <- toPlot$y %>% min()
  ymax <- toPlot$y %>% max()
  diff <- ymax - ymin
  ticksSize <- 0.08
  tick.length <- diff * 0.2
  tick.min <- ymin - tick.length
  tick.max <- ymin
  tick.color <- "green3"

  #>>>>>>>>>>>> FANC:: add ribbon to plot
  names(df) <- c("gene", "metric")
  df <- df %>%
    dplyr::arrange(desc(metric)) %>%
    dplyr::mutate(
      rank = 1:nrow(df),
      trend = ifelse(metric > 0, "up", "down"))

  df <- df %>% split(., f= .$trend) %>% lapply(function(df) {
    r <- df$rank
    if (df$metric[1] > 0) {
      r <- -1 * r
    }
    df$cat <- cut(r, breaks = 3) %>% as.numeric()

    if (df$metric[1] < 0) {
      df$cat <- -1 * df$cat
    }

    df <- df %>% split(., f = .$cat) %>% lapply(function(df){
      data.frame(
        xmin = min(df$rank),
        xmax = max(df$rank),
        type = df$cat[1]
      )
    }) %>% do.call(rbind, .)
    return(df)
  }) %>% do.call(rbind, .)
  # df$type <- factor(df$type, levels = c(paste0("down", 1:3 ), paste0("up", 1:3))) %>% as.numeric()

  ribbon.size = 0.08 # as a fraction of diff, which is ymin - ymax
  ribbon.min <- tick.min - ribbon.size
  df$ymin <- ribbon.min
  df$ymax <- tick.min
  #<<<<<<<< FANC:: add ribbon to plot

  #>>>>>>>> FANC: manually add x.label to plot:
  #
  xmin = min(toPlot$x)
  xmax = max(toPlot$x)
  xrange = xmax - xmin
  xmid <- floor(xmin + 0.5 * xrange)
  x.lab.min <- ribbon.min - 0.3 * diff
  updown.min <- x.lab.min - 0.3 * diff
  # print(xmid)
  p <- ggplot(toPlot, aes(x = x, y = y)) +
    # geom_point(color = "green", size = 0.1) +
    ggrastr::rasterize(geom_line(color = "gold", size = 0.3), dpi = 300) +
    ggrastr::rasterize(geom_segment(data = data.frame(x = pathway),
                 mapping = aes(x = x,y = tick.min, xend = x, yend = tick.max),
                 size = ticksSize, alpha = 1, color = tick.color), dpi = 300) +
    geom_rect(data = df,
              mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type),
              size = 0, show.legend = F) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    # geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") +
    # geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
    # geom_hline(yintercept = 0, colour = "black") +
    # geom_text(data = text.df, aes(x = x, y = y, label = text), inherit.aes = F, size = 0.36 * 6) +
    annotate("text", x = xmin, y = ribbon.min, label = rank.name.display,
             hjust = 0, vjust = 1.4, size = 0.36 * 5.5,
             family = "Arial") +
    annotate("text", x = xmin, y = x.lab.min, label = up.name.display,
             hjust = 0, vjust = 1.3, color = "red", size = 0.36 * 5.5,
             family = "Arial") +
    annotate("text", x = Inf, y = x.lab.min, label = down.name.display,
             hjust = 1.0, vjust = 1.3, color = "blue", size = 0.36 * 5.5,
             family = "Arial") +
    annotate("segment", x = xmin + 0.02 * xrange, xend = xmin + 0.02 * xrange,
             y = updown.min - tick.length, yend = updown.min,
             size = 10* ticksSize, alpha = 1, color = tick.color) +
    annotate("text", x = xmin + 0.05 * xrange, y = updown.min, label = "Overlapping Gene Set",
             hjust = 0, vjust = 1.1, size = 0.36 * 5,
             family = "Arial")

  p <- p +  coord_cartesian(ylim = c(tick.min , ymax),clip="off") +
    labs(y = "Enrichment Score", title = gene.set.name.display) +
    theme_bw(base_size = 6, base_family = "Arial") +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 6, family = "Arial", hjust = 0),
          plot.margin = margin(t = 0.01, b = 0.3 * height, r = 0, l = 0.01, unit = "in"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 6, family = "Arial"),
          axis.title.y = element_text(size = 6, family = "Arial", vjust = 0))

  if (calculate.p && add.p.to.plot) {
    p.value <- gsea.stats$pval[1]

    if (p.value < 0.001) {
      p.value <- format(p.value, digits = 2, scientific = T)
      p.value <- as.character(p.value)
      if (!grepl("^\\d(\\.\\d+)*e\\-\\d+$", p.value)) {
        stop(paste0("Unexpected format for p value: ", p.value))
      }
      a <- sub("e.+$", "", p.value)
      b <- sub("^.+e", "", p.value)

      label <- substitute(italic(p) == A %*% 10^B, list(A = a, B = b))

    } else {
      p.value <- format(p.value, digits = 2, scientific = F)
      label <- substitute(italic(p) == A, list(A = p.value))
    }

    p <- p + annotate("text", x = Inf, y = Inf,
                      label = label,
                      hjust = 1.1, vjust = 1.1, family = "Arial", size = 0.36 * 5)
  }
  if (!is.null(plot.out)) {
    dir.create(dirname(plot.out), showWarnings = F, recursive = T)
    ggsave(plot.out, p, width = width, height = height, units = "in", device = cairo_pdf)
    if (calculate.p) {
      stat.file <- paste0(tools::file_path_sans_ext(plot.out), "_stats.txt")
      gsea.stats$leadingEdge <- gsea.stats$leadingEdge[[1]] %>% paste0(collapse = ",")
      write.table(as.data.frame(gsea.stats), stat.file, sep = "\t", quote = F,
                  col.names = T, row.names = F)
    }
  }

  invisible(p)
}

gsea.plot.curve.m <- function(job.df, out.dir) {
  stop("untested")
  func.fields <- c("gmt.file", "gmt.gene.set.name", "gene.set.name.display",
                   "rnk", "rank.name.display",
                   "up.name.display", "down.name.display")
  required.fields <- c(func.fields, "out.root")

  if (is.character(job.df)) {
    job.df <- read.table(job.df, header = T)
  }

  utilsFanc::check.intersect(required.fields, colnames(job.df))

  lapply(1:nrow(job.df), function(i) {
    print(paste0("Processing job ", i))
    df <- job.df[i, ]
    params <- as.list(df[, func.fields])
    params$plot.out <- paste0(out.dir, "/", df$out.root, ".pdf")
    do.call(gsea.plot.curve, params)
    return()
  })
}
