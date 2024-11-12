coexpr.grid <- function(so, assay, slot,
                        genes, gene.alias = NULL, use.gene.order = F,
                        method,
                        MI.use.bootstrap = F, MI.force.bootstrap = F,
                        MI.diag.na = F,
                        pos.ratio.log2 = F,
                        cluster.ident = NULL, clusters = NULL,
                        group.ident = "sample", groups = NULL,
                        hm.values = c(-0.5, 0, 0.5), hm.use.cell_fun = T,
                        hm.colors = c("green", "white", "red"),
                        format = "pdf",
                        gene.name.fontsize = 6,
                        out.dir, root = NULL,
                        pub = F,
                        ...) {
  if (is.null(root)) root <- paste0(basename(out.dir))
  dir.create(out.dir, showWarnings = F, recursive = T)

  sig.levels <- c(0.05, 0.01, 0.001)
  sig.stars <- c("*", "**", "***")

  suppressMessages(library(circlize))
  suppressMessages(library(ComplexHeatmap))

  if (length(sig.levels) != length(sig.stars)) {
    stop("length(sig.levels) != length(sig.stars)")
  }

  mat <- GetAssayData(object = so, assay = assay, slot = slot)
  utilsFanc::check.intersect(genes, "supplied genes", rownames(mat),
                             "genes in so")
  mat <- mat[genes, ]
  if (!is.null(gene.alias)) {
    if (gene.alias %in% c("translate", "translate_strip")) {
      temp <- gene.alias
      gene.alias <- utilsFanc::klra.ly49.translate(genes)
      gene.alias <- gene.alias %>% .[!is.na(.)]
      if (temp == "translate_strip") {
        gene.alias <- sub(".*Ly49", "", gene.alias)
      }
    }
    if (length(genes) != length(gene.alias)) {
      stop("length(genes) != length(gene.alias)")
    }
    rownames(mat) <- gene.alias
  }

  cell.list <- get.cell.list(obj = so, is.ao = F, group.by = group.ident,
                             groups = groups, split.by = cluster.ident, splits = clusters,
                             na.rm = T, return.named.vec = F, n.cells.each = NULL)

  master.stat <- paste0(out.dir, "/", root, "_stats.tsv")
  system(paste0("rm -rf ", master.stat))
  p.info <- lapply(names(cell.list), function(name) {
    cells <- cell.list[[name]]
    mat <- mat[, cells] %>% as.matrix()
    nFeature.df <- lapply(seq_along(genes), function(x) {
      gene <- genes[x]
      alias <- gene.alias[x]
      meta.df <- so@meta.data[cells, "nFeature_RNA", drop = F]
      meta.df$expr <- NA
      meta.df$expr[mat[alias,] > 0] <- "ip.ave.nFeature"
      meta.df$expr[mat[alias,] == 0] <- "in.ave.nFeature"
      if (any(is.na(meta.df))) {
        stop()
      }
      meta.df <- meta.df %>% dplyr::group_by(expr) %>%
        dplyr::summarise(ave.nFeature = round(mean(nFeature_RNA), digits = 0)) %>%
        dplyr::ungroup() %>% as.data.frame()
      meta.df$i <- alias
      return(meta.df)
    }) %>% do.call(rbind, .)
    nFeature.df <- reshape2::dcast(nFeature.df, i ~ expr, value.var = "ave.nFeature")

    prefix <- paste0(root, "_", name, "_", assay, "_", slot, "_", method)
    if (method %in% c("pearson", "MI", "enrich", "pos.ratio")) {
      # calculate pearson for both binarized and non-binarized.
      # remove purely zero genes.
      mat <- mat[rowSums(mat > 0) > 0,]
      mats <- list(num = mat, bin = (mat > 0) + 0)
      p.info <- lapply("bin", function(mat.type) {
        # decided to just use the binary matrix.
        prefix <- paste0(prefix, "_", mat.type)
        if (mat.type == "bin")
          enrich.res <- coexpr.enrich(bmat = mats[[mat.type]], log2.OR = pos.ratio.log2)

        if (method == "pearson") {
          mat.use <- as.matrix(t(mats[[mat.type]]))
          corr <- cor(mat.use, method = "pearson")
          corr <- round(corr, digits = 3)
          # get p value matrix:
          p.df <- lapply(1:(ncol(mat.use)- 1), function(i) {
            lapply((i + 1):(ncol(mat.use)), function(j) {
              p <- cor.test(mat.use[, i], mat.use[, j])$p.value
              df <- data.frame(i = colnames(mat.use)[i], j = colnames(mat.use)[j], p = p)
            }) %>% do.call(rbind, .)
          }) %>% do.call(rbind,. )
          write.table(p.df, paste0(out.dir, "/", prefix, "_p.tsv"), quote = F, sep = "\t",
                      col.names = T, row.names = F)

          p.df2 <- p.df
          colnames(p.df2) <- c("j", "i", "p")
          p.df <- rbind(p.df, p.df2)

          p.df$q <- p.df$p %>% p.adjust(method = "bonferroni")

          q <- p.df %>% reshape2::acast(i ~ j, value.var = "q")
        }
        else if (method == "MI") {
          mat.use <- as.matrix(t(mats[[mat.type]]))
          corr <- infotheo::mutinformation(as.data.frame(mat.use))
          corr <- round(corr, digits = 3)
          if (MI.use.bootstrap) {
            p.df.file <- paste0(out.dir, "/", prefix, "_p.tsv")
            if (!file.exists(p.df.file) || MI.force.bootstrap) {
              print(paste0("calculating MI bootstraping. will save to: ", p.df.file))
              p.df <- MI.bootstrap(mat = mat.use, n = 100, threads = 10, return.df = T)
              write.table(p.df, p.df.file, quote = F, sep = "\t",
                          col.names = T, row.names = F)
            } else {
              print(paste0("reading MI bootstraping p values from: ", p.df.file))
              p.df <- read.table(p.df, header = T)
            }
          } else {
            p.df <- lapply(1:(ncol(mat.use)- 1), function(i) {
              lapply((i + 1):(ncol(mat.use)), function(j) {
                G <- suppressWarnings(AMR::g.test(x = mat.use[, i], y = mat.use[, j]))
                p <- G$p.value
                df <- data.frame(i = colnames(mat.use)[i], j = colnames(mat.use)[j], p = p)
              }) %>% do.call(rbind, .)
            }) %>% do.call(rbind,. )
          }

          p.df2 <- p.df
          colnames(p.df2) <- c("j", "i", "p")
          p.df <- rbind(p.df, p.df2)

          p.df$q <- p.df$p %>% p.adjust(method = "bonferroni")

          q <- p.df %>% reshape2::acast(i ~ j, value.var = "q")
          diag(q) <- NA
        } else if (method %in% c("enrich", "pos.ratio")) {
          if (mat.type == "num") {
            return()
          }
          write(paste0("# ", prefix), master.stat, append = T)
          write.table(enrich.res$df, master.stat, quote = F, sep = "\t", append = T,
                      col.names = T, row.names = F)
          write("\n", master.stat, append = T)

          write.table(df, paste0(out.dir, "/", prefix, "_details.tsv"), quote = F, sep = "\t",
                      col.names = T, row.names = F)
          if (method == "enrich")
            corr <- enrich.res$mats$enrich
          else if (method == "pos.ratio")
            corr <- enrich.res$mats$pos.ratio
          else
            stop()
          corr <- round(corr, digits = 3)
          corr[!is.finite(corr)] <- NA
          q <- enrich.res$mats$q
        }


        df <- cbind(data.frame(gene = rownames(corr)), as.data.frame(corr))

        write.table(df, paste0(out.dir, "/", prefix, ".tsv"), quote = F, sep = "\t",
                    col.names = T, row.names = F)

        write(paste0("# ", prefix), master.stat, append = T)
        write.table(df, master.stat, quote = F, sep = "\t", append = T,
                    col.names = T, row.names = F)
        write("\n", master.stat, append = T)

        if (method %in% c("pearson", "MI", "enrich", "pos.ratio")) {
          symbol.mat <- corr

          for (i in 1:length(sig.levels)) {
            symbol.mat[which(q < sig.levels[i])] <- symbol.mat[which(q < sig.levels[i])] %>%
              sub("\n\\**", "", .)
            symbol.mat[which(q < sig.levels[i])] <-
              paste0(symbol.mat[which(q < sig.levels[i])], "\n", sig.stars[i])
          }
          cell.fun.mat <- symbol.mat
        } else {
          cell.fun.mat <- corr
        }

        col_fun = colorRamp2(hm.values, hm.colors)
        if (hm.use.cell_fun) {
          cell_fun <- function(j, i, x, y, width, height, fill) {
            grid.text(cell.fun.mat[i, j], x, y,
                      gp = gpar(fontsize = 10))
          }
        } else {
          cell_fun <- NULL
        }

        p <- Heatmap(matrix = corr, col = col_fun, cell_fun = cell_fun,
                     name = paste0("m", sample(1:100, 1)),
                     clustering_distance_rows = "pearson",
                     clustering_distance_columns = "pearson",
                     row_names_gp = gpar(fontface = "italic", fontsize = gene.name.fontsize),
                     column_names_gp = gpar(fontface = "italic", fontsize = gene.name.fontsize),
                     ...)
        save.base.plot(p = p, file = paste0(out.dir, "/", prefix, "_hm.", format),
                       width = 800, height = 700)
        if (mat.type == "num")
          return()

        # return information required for p value calculations
        enrich.res$df <- left_join(enrich.res$df, nFeature.df, by = "i")
        enrich.res$df$name <- name
        res <- list(enrich.df = enrich.res$df, corr = corr)
        if (method %in% c("pearson", "MI", "pos.ratio")) {
          res$q <- q
        }
        return(res)
      }) %>% Reduce(c, .)

      return(p.info)
    }
  })

  # mean corr
  common.genes <- lapply(p.info, function(x) return(rownames(x$corr))) %>%
    Reduce(intersect, .)

  # get linear regression based p values
  p.info <- lapply(p.info, function(x) {
    x$corr <- x$corr[common.genes, common.genes]
    x$enrich.df <- x$enrich.df %>% dplyr::filter(i %in% common.genes, j %in% common.genes)
    if (!is.null(x$q)) {
      x$q <- x$q[common.genes, common.genes]
    }
    return(x)
  })
  corr.list <- lapply(p.info, function(x) return(x$corr))
  corr <- round(Reduce(`+`, corr.list)/length(corr.list), digits = 3)
  # now we calculate p value
  formula <- jp ~ istat
  p.df <- lapply(gene.alias, function(gi) {
    lapply(gene.alias, function(gj) {
      if (gi == gj) {
        return()
      }
      df <- lapply(p.info, function(x) {
        x$enrich.df %>% dplyr::filter(i == gi, j == gj) %>% return()
      }) %>% do.call(rbind, .)
      df.jp <- df %>% reshape2::melt(measure.vars = c("jp.in.in", "jp.in.ip"),
                                  variable.name = "type", value.name = "jp")
      df.jp <- df.jp[, c("name", "type", "jp")]
      df.jp$istat <- df.jp$type %>% sub("jp.in.", "", .)
      df.nf <- df[, c("name", "in.ave.nFeature", "ip.ave.nFeature")]
      colnames(df.nf) <- c("name", "in", "ip")
      df.nf <- df.nf %>% reshape2::melt(measure.vars = c("in", "ip"),
                                        variable.name = "istat", value.name = "ave.nFeature")
      df.p <- left_join(df.jp, df.nf)
      # turned out that including ave.nFeature makes the p value hard to compute
      # when only 2 samples are available.

      fit <- lm(formula, data = df.p)
      coef <- fit$coefficients["istatip"]
      p <- summary(fit)$coefficients["istatip", 4]
      res <- data.frame(i = gi, j = gj, coef = coef, p = p)
      return(res)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
  p.df$q <- p.adjust(p = p.df$p, method = "BH")
  
  p.df <- utilsFanc::change.name.fanc(df = p.df, cols.from = c("p", "q"),
                                      cols.to = c("lm.p", "lm.q"))

  if (method %in% c("pearson", "MI", "pos.ratio")) {
    q.list <- lapply(p.info, function(x) return(x$q))
    max.q <- do.call(pmax, q.list)
    max.q <- reshape2::melt(max.q, varnames = c("i", "j"), value.name = paste0(method, ".q"))
    p.df <- p.df %>% dplyr::left_join(max.q)
  }

  p.tsv <- paste0(out.dir, "/", root, "_p.tsv")
  write(as.character(formula), p.tsv, sep = "\t")
  write.table(p.df, p.tsv, quote = F,
              sep = "\t", row.names = F, col.names = T)

  corr.coef <- reshape2::acast(p.df, i~j, value.var = "coef")

  if (method %in% c("pearson", "MI")) {
    p.types <- c("lm.p", "lm.q", "none")
  } else if (method == "pos.ratio") {
    p.types <- c()
  }
  if (method %in% c("pearson", "MI", "pos.ratio")) {
    p.types <- c(p.types, paste0(method, ".q"))
  }
  lapply(p.types, function(sta) {
    # sta: statistic
    if (sta != "none" )
      sta.mat <- reshape2::acast(p.df, i ~ j, value.var = sta)
    base.mats <- list(ave = corr, coef = corr.coef)
    if (method == "pos.ratio") {
      base.mats <- list(ave = corr)
    }
    lapply(names(base.mats), function(base.type) {
      base.mat <- round(base.mats[[base.type]], digits = 3)
      symbol.mat <- base.mat
      if (sta != "none") {
        for (i in 1:length(sig.levels)) {
          symbol.mat[which(sta.mat < sig.levels[i])] <- symbol.mat[which(sta.mat < sig.levels[i])] %>%
            sub("\n\\**", "", .)
          symbol.mat[which(sta.mat < sig.levels[i])] <-
            paste0(symbol.mat[which(sta.mat < sig.levels[i])], "\n", sig.stars[i])
        }

      }

      cell_fun <- function(j, i, x, y, width, height, fill) {
        grid.text(symbol.mat[i, j], x, y,
                  gp = gpar(fontsize = gene.name.fontsize - 1))
      }

      if (!hm.use.cell_fun) {
        cell_fun <- NULL
      }

      if (base.type == "ave") {
        col_fun = colorRamp2(hm.values, hm.colors)
      } else {
        col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
      }

      hm.params <- list(
        matrix = base.mat, col = col_fun, cell_fun = cell_fun,
        row_names_gp = gpar(fontface = "italic", fontsize = gene.name.fontsize),
        column_names_gp = gpar(fontface = "italic", fontsize = gene.name.fontsize),
        width = unit(2.5, "in"), height = unit(2.5, "in"),
        heatmap_legend_param = list(
          title = "", labels_gp = gpar(fontsize = gene.name.fontsize - 1),
          legend_height = unit(0.3, "in"), grid_width = unit(0.05, "in"), gap = unit(2, "in")
        )
      )

      if (use.gene.order) {
        gene.alias <- gene.alias %>% .[ .%in% colnames(base.mat)]
        order.params <- list(
          column_order = gene.alias,
          row_order = gene.alias
        )
      } else {
        order.params <- list(
          clustering_distance_rows = "pearson",
          clustering_distance_columns = "pearson",
        )
      }
      hm.params <- c(hm.params, order.params, list(...))
      p <- do.call(ComplexHeatmap::Heatmap, hm.params)
      save.base.plot(p = p, file = paste0(out.dir, "/", root, "_", base.type, "_w_", sta ,".", format),
                     width = 300, height = 300)
      # for making main figure:
      if (pub) {
        tmp <- dimnames(symbol.mat)
        symbol.mat <- stringr::str_extract(symbol.mat, "\\*+") %>%
          matrix(nrow = length(tmp[[1]]))
        ## previously written as "ncol". didn't matter in that case because the matrix was symmetrical.
        dimnames(symbol.mat) <- tmp
        rm(tmp)
        symbol.mat[is.na(symbol.mat)] <- ""
        cell_fun <- function(j, i, x, y, width, height, fill) {
          grid.text(symbol.mat[i, j], x, y,
                    gp = gpar(fontsize = 5))
        }

        suppressMessages(extrafont::loadfonts())
        ht_opt$HEATMAP_LEGEND_PADDING <- unit(0, "in")
        ht_opt$DENDROGRAM_PADDING <- unit(0, "in")

        hm.params <- list(
          matrix = base.mat, col = col_fun, cell_fun = cell_fun,
          row_names_gp = gpar(fontsize = 6, fontfamily = "Arial", fontface = "italic"),
          column_names_gp = gpar(fontsize = 6, fontfamily = "Arial", fontface = "italic"),
          width = unit(0.8, "in"), height = unit(0.8, "in"),
          show_heatmap_legend = T,
          heatmap_legend_param = list(
            title = "", labels_gp = gpar(fontsize = 5),
            legend_height = unit(0.3, "in"), grid_width = unit(0.05, "in"), gap = unit(2, "in")
          )
        )

        hm.params <- c(hm.params, order.params)
        p <- do.call(ComplexHeatmap::Heatmap, hm.params)
        save.base.plot(p = p, file = paste0(out.dir, "/", root, "_pub_", base.type, "_w_", sta ,".", format),
                       width = 150, height = 150)
        ht_opt(RESET = T)
      }
      return()
    })
  })

}

coexpr.enrich <- function(bmat, by.col = F, enhance.q = F, log2.OR = F) {
  # matrix must be binary
  if (!identical(as.numeric(sort(unique(as.vector(bmat)))), as.numeric(c(0, 1)))) {
    stop("bmat must be a matrix of 0 and 1s")
  }
  if (by.col) {
    bmat <- t(bmat)
  }

  df <- lapply(rownames(bmat), function(i) {
    lapply(rownames(bmat), function(j) {
      # n: number. ip and in: i positive and i negative (think about the Ly49s)
      # j is an urn containing white balls (positive) and black balls (negative)
      # ip represents the result of a draw from j.
      n <- ncol(bmat)
      n.ip = sum(bmat[i, ])
      n.in = n - n.ip

      n.jp = sum(bmat[j, ])
      n.jn = n - n.jp

      # dp: double positive
      n.dp.expected <- n * (n.ip/n) * (n.jp/n)
      n.dp.obs <- sum(bmat[i, ] > 0 & bmat[j, ] > 0)
      enrich <- round(n.dp.obs/n.dp.expected, digits = 3)
      # p <- phyper(q = n.dp.obs, m = n.jp, n = n.jn, k = n.ip, lower.tail = F)
      # p <- 2 * pmin(p, 1-p)
      # p <- round(p, digits = 3)
      #
      # pos ratio
      n.ipjn <- sum(bmat[i, ] > 0 & bmat[j, ] == 0)
      n.injp <- sum(bmat[i, ] == 0 & bmat[j, ] > 0)
      n.dn <- sum(bmat[i, ] == 0 & bmat[j, ] == 0)

      fisher.df <- data.frame(i.n = c(n.dn, n.injp), i.p = c(n.ipjn, n.dp.obs))
      rownames(fisher.df ) <- c("j.n", "j.p")
      fisher <- fisher.test(fisher.df)
      p <- fisher$p.value
      # p <- round(p, digits = 4)

      pos.ratio <- (n.dp.obs/n.ipjn)/(n.injp/n.dn)
      if (log2.OR) {
        pos.ratio <- log2(pos.ratio)
      }
      pos.ratio <- round(pos.ratio, digits = 3)
      jp.in.in <- round(n.injp/(n.injp + n.dn), digits = 3)
      jp.in.ip <- round(n.dp.obs/(n.dp.obs + n.ipjn), digits = 3)
      res <- data.frame(i = i, j = j,
                        dp.expected = n.dp.expected, dp.obs = n.dp.obs,
                        pos.ratio = pos.ratio,
                        enrich = enrich, p = p,
                        dn = n.dn, ipjn = n.ipjn, injp = n.injp, dp = n.dp.obs,
                        jp.in.in = jp.in.in, jp.in.ip = jp.in.ip)
      return(res)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)

  # calculate q:
  df$p[df$i == df$j] <- NA
  df$enrich[df$i == df$j] <- NA

  df$iF <- factor(df$i) %>% as.numeric()
  df$jF <- factor(df$j) %>% as.numeric()

  if (enhance.q) {
    df.q <- df %>% dplyr::filter(jF > iF) %>%
      dplyr::select(i, j, p) %>%
      dplyr::mutate(q = p.adjust(p, method = "fdr"))
      # dplyr::mutate(q = round(p.adjust(p, method = "fdr"), digits = 4))
    df.q2 <- df.q
    colnames(df.q2) <- c("j", "i", "p", "q")
    df.q <- rbind(df.q, df.q2)
  } else {
    
    df.q <- df %>% dplyr::filter(j != i) %>%
      dplyr::select(i, j, p) %>%
      dplyr::mutate(q = p.adjust(p, method = "fdr"))
      # dplyr::mutate(q = round(p.adjust(p, method = "fdr"), digits = 4))
  }

  df <- dplyr::left_join(df, df.q) %>%
    dplyr::mutate(iF = NULL, jF = NULL)
  mats <- lapply(c("enrich", "p", "q", "pos.ratio"), function(x) {
    mat <- reshape2::acast(df, i ~ j, value.var = x)
  })
  names(mats) <- c("enrich", "p", "q", "pos.ratio")
  res <- list(df = df, mats = mats)
  return(res)
}

MI.bootstrap <- function(mat, n = 100, threads = 1, return.df = T) {
  if (is.null(colnames(mat))) {
    colnames(mat) <- 1:ncol(mat)
  }
  if (is.null(rownames(mat))) {
    rownames(mat) <- 1:nrow(mat)
  }

  p.df <- lapply(1:(ncol(mat) - 1), function(i) {
    lapply((i + 1):(ncol(mat)), function(j) {
      print(paste0("i: ", i, "; j: ", j))
      MI <- infotheo::mutinformation(X = mat[, i], Y = mat[, j])
      # permutate:
      MI.per <- utilsFanc::safelapply(1:n, function(x) {
        i.per <- sample(mat[, i], size = length(mat[, i]), replace = T)
        j.per <- sample(mat[, j], size = length(mat[, j]), replace = T)
        infotheo::mutinformation(X = i.per, Y = j.per)
      }, threads = threads) %>% unlist()
      fun <- ecdf(MI.per)
      p <- 1 - fun(MI)
      df <- data.frame(i = colnames(mat)[i], j = colnames(mat)[j], p = p)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind,. )
  if (return.df) {
    return(p.df)
  } else {
    stop("only return.df has been developed.")
  }
}
