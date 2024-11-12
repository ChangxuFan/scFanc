de.da.intersect <- function(pbl1, pbl2, name1, name2,
                            directions = c("up", "down"),
                            slot = "summary", clusters = NULL,
                            rule = "promoter", promoter.ext.1side = 1000,
                            out.dir = NULL, root.name = NULL) {
  # pbl1 and pbl2 could be either de or da.
  # when both de: intersect by gene name
  # when de and da: convert gene name to loci, then perform loci intersection (GRanges based)
  # when both da: loci intersection (GRanges based)
  # rule is used to instruct how to tell if a DEG is intersecting a DAR.
  if (is.null(root.name) && !is.null(out.dir)) {
    root.name <- basename(out.dir)
  }
  all.clusters <- intersect(names(pbl1), names(pbl2))
  if (length(all.clusters) < 1) {
    stop("length(all.clusters) < 1")
  }

  if (is.null(clusters)) {
    clusters <- all.clusters
  }
  clusters <- clusters %>% .[.%in% all.clusters]
  if (length(clusters) < 1) {
    stop("length(clusters) < 1")
  }
  pbl1 <- pbl1[clusters]
  pbl2 <- pbl2[clusters]

  des <- list(pbl1, pbl2)
  names(des) <- c(name1, name2)
  type <- sapply(des, function(de) {
    gn <- rownames(de[[1]]$bulk.mat)[1]
    return(ifelse(grepl(":", gn), "da", "de"))
  })
  bG <- which(type == "de")
  bP <- which(type == "da")
  stats <- lapply(clusters, function(cluster) {
    if (sum(type == "de") > 0 && sum(type == "da") > 0) {
      if (rule == "promoter") {
        if (length(des) != 2) stop("length(des) != 2")

        stats <- lapply(directions, function(direc) {
          pros <- des[[bG]][[cluster]][[slot]]$promoters[[direc]]
          if (is.null(pros)) {
            stop(paste0("promoter slot empty for ", names(des)[bG],
                        ", cluster ", cluster,", slot ", slot))
          }
          peaks <- des[[bP]][[cluster]][[slot]][[paste0(direc, ".genes")]]
          if (is.null(peaks)) {
            stop(paste0("peak slot empty for ", names(des)[bP],
                        ", cluster ", cluster,", slot ", slot))
          }

          pros <- GenomicRanges::resize(pros, width = 2 * promoter.ext.1side,
                                        fix = "center")
          peaks <- peaks %>% utilsFanc::loci.2.df(loci.vec = ., remove.loci.col = F, return.gr = T)

          o <- plyranges::join_overlap_inner(pros, peaks)

          if (!is.null(out.dir)) {
            # merge rows that correspond to the same gene.
            utilsFanc::write.zip.fanc(
              o, out.file = paste0(out.dir, "/", root.name, "_", cluster, "_", direc, "_overlap.bed"),
              bed.shift = T)
          }

          stats <- data.frame(cluster = cluster, direction = direc,
                              n.de = des[[bG]][[cluster]][[slot]][[paste0("n.", direc)]],
                              n.pro = pros$gene %>% unique() %>% length(),
                              n.da = length(peaks),
                              n.de.w.da = o$gene %>% unique() %>% length()) %>%
            dplyr::mutate(frac.de.w.da = round(n.de.w.da/n.de, digits = 3))
          return(stats)
        }) %>% do.call(rbind, .)
        return(stats)

      } else if (rule == "nearestGene") {
        stats <- lapply(c("up", "down"), function(direc) {
          anno <- des[[bP]][[cluster]][[slot]]$anno[[direc]]
          if (is.null(anno)) {
            stop(paste0("anno slot empty for ", names(des)[bP],
                        ", cluster ", cluster,", slot ", slot))
          }
          deg <- des[[bG]][[cluster]][[slot]][[paste0(direc, ".genes")]]
          if (is.null(deg)) {
            stop(paste0("deg slot empty for ", names(des)[bG],
                        ", cluster ", cluster,", slot ", slot))
          }

          o <- anno[anno$gene %in% deg]
          if (!is.null(out.dir)) {
            utilsFanc::write.zip.fanc(
              o, out.file = paste0(out.dir, "/", root.name, "_", rule, "_", cluster, "_", direc, "_overlap.bed"),
              bed.shift = T)
          }

          stats <- data.frame(
            cluster = cluster, direction = direc,
            n.de = des[[bG]][[cluster]][[slot]][[paste0("n.", direc)]],
            n.da = length(anno),
            n.da.w.anno = sum(!is.na(anno$gene)),
            n.da.w.de = anno$gene %>% .[!is.na(.)] %>% .[. %in% deg] %>% unique() %>% length(),
            n.de.w.da = deg %>% .[. %in% anno$gene] %>% unique() %>% length()
          ) %>% dplyr::mutate(
            frac.da.w.de = round(n.da.w.de/n.da, digits = 3),
            frac.de.w.da = round(n.de.w.da/n.de, digits = 3)
          )

          return(stats)
        }) %>% do.call(rbind, .)
        return(stats)
      } else {
        stop("rule can be: promoter")
      }
    } else stop("only de vs da has been developed")
  }) %>% do.call(rbind, .)
  if (!is.null(out.dir)) {
    dir.create(out.dir, showWarnings = F, recursive = T)
    write.table(stats, paste0(out.dir, "/", root.name, "_", rule, "_stats.tsv"), sep = "\t", quote = F,
                row.names = F, col.names = T)
  }
  invisible(stats)
}

de.da.intersect.titrate <- function(de, da, slot = "summary", directions = c("up", "down"),
                                    mode = "promoter", promoter.ext.1side = 500, genome.size,
                                    use.background = F, n.folds.bg,
                                    padj.cutoffs = c(0.05, 0.1, 0.2),
                                    log2FC.cutoffs = c(0, 0.25, 0.5, 1),
                                    out.dir, root.name = NULL, print.details = T) {

  if(is.null(root.name)) root.name <- basename(out.dir)
  mode.name <- mode
  if (mode == "promoter") mode.name <- paste0("ProExt", promoter.ext.1side)
  root.name <- paste0(root.name, "_", mode.name)

  if (mode %in% c("promoter")) {
    des <- list(de = de, da = da)

    to.titrate <- ifelse(mode == "promoter", "da", "de")
    bT <- which(names(des) == to.titrate)
    bA <- which(names(des) != to.titrate) # A as in "anchor"
    de.slot.assert(des[[bA]], is.s2b = F, slot = slot)
    anchor.sets <- slot
    if (use.background) {
      # check if the bg slot is available
      # bg can be added by de.add.bg() -> de.bg.expand()
      lapply(des[[bA]], function(s2b) {
        utilsFanc::check.intersect(paste0("fold", 1:n.folds.bg), "folds", names(s2b), paste0("cluster ", s2b$root.name))
      })
      anchor.sets <- c(slot, paste0("fold", 1:n.folds.bg))
    }

  } else {
    stop("Only promoter has been developed")
  }

  res <- lapply(padj.cutoffs, function(padj) {
    lapply(log2FC.cutoffs, function(log2FC) {
      if (mode %in% c("promoter")) {
        DAR.dir <- NULL
        if (print.details) DAR.dir <- paste0(out.dir, "/DARs/DAR_padj", padj, "_log2FC", log2FC, "/")
        des[[bT]] <- deseq2.summary(des[[bT]], padj.cutoff = padj, log2fc.cutoff = log2FC, summary.slot = slot,
                                    gene.out.dir = DAR.dir)
        res <- lapply(anchor.sets, function(anchor.set) {
          # e.g. for mode == "promoter": da uses the same slot, whereas de uses a different slot for fg or bg
          des[[bT]] <- lapply(des[[bT]], function(s2b) {
            s2b[[anchor.set]] <- s2b[[slot]]
            return(s2b)
          })

          # if (grepl("fold", anchor.set)) {
          #   des[[bA]] <- lapply(des[[bA]], function(s2b) {
          #     bg <- s2b[[slot]]$bg
          #     s2b[[anchor.set]] <- list(up.genes = bg$up[[anchor.set]], down.genes = bg$down[[anchor.set]],
          #                               n.up = length(bg$up[[anchor.set]]), n.down = length(bg$down[[anchor.set]]))
          #     return(s2b)
          #   })
          # }
          details.dir <- NULL
          if (print.details && !grepl("^fold\\d+$", anchor.set)) {
            details.dir <- paste0(out.dir, "/details/")
          }
          root.name <- paste0(root.name, "_padj", padj, "_log2FC", log2FC, "_slot", anchor.set)

          res <- de.da.intersect(pbl1 = des$de, pbl2 = des$da, name1 = "de", name2 = "da", slot = anchor.set,
                                 rule = mode, promoter.ext.1side = promoter.ext.1side, out.dir = details.dir,
                                 root.name = root.name,
                                 directions = directions)
          frame <- data.frame(anchor.set = anchor.set, padj = padj, log2FC = log2FC)
          res <- cbind(frame, res)

          res <- lapply(1:nrow(res), function(i) {
            # Add stats:
            res <- res[i, ]
            stats <- de.da.intersect.p(n.deg = res$n.pro, n.ext.1side = promoter.ext.1side, n.dar = res$n.da,
                                       genome.size = genome.size, n.overlap = res$n.de.w.da)
            stats <- stats[, c("frac.of.genome", "expected.overlap", "binomial.p")]
            res <- cbind(res, stats)
            return(res)
          }) %>% do.call(rbind, .)
          return(res)
        }) %>% do.call(rbind, .)
        return(res)
      } else {
        stop("Only promoter has been developed")
      }

    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
  dir.create(out.dir, showWarnings = F, recursive = T)
  write.table(res, paste0(out.dir, "/", root.name, ".tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
  invisible(res)
}

de.da.intersect.titrate.add.p <- function(in.tsv, ext.1side, genome.size = 2.5 * 10^9,
                                          out.tsv = NULL) {
  # take in the results from de.da.intersect.titrate, and calculate a p value
  df <- read.table(in.tsv, header = T)
  out <- lapply(1:nrow(df), function(i) {
    out <- df[i, ]
    stats <- de.da.intersect.p(n.deg = out$n.pro, n.ext.1side = ext.1side, n.dar = out$n.da,
                             genome.size = genome.size, n.overlap = out$n.de.w.da)
    stats <- stats[, c("frac.of.genome", "expected.overlap", "binomial.p")]
    out <- cbind(out, stats)
  }) %>% do.call(rbind, .)
  if (is.null(out.tsv)) {
    out.tsv <- utilsFanc::insert.name.before.ext(name = in.tsv, insert = "addP", delim = "_")
  }
  dir.create(dirname(out.tsv), showWarnings = F, recursive = T)
  write.table(out, out.tsv, quote = F, sep = "\t", row.names = F, col.names = T)
  return(out)
}

de.da.intersect.p <- function(n.deg, n.ext.1side, n.dar, genome.size = 2.5*10^9,
                            n.overlap) {
  # first calculate the total fraction occupied by deg + extension:
  out <- list(n.deg = n.deg, n.ext.1side = n.ext.1side, n.dar = n.dar, n.overlap = n.overlap)
  out$frac.of.genome <- n.deg * n.ext.1side * 2/genome.size
  # we model this by a binomial distribution: every dar is flipping a coin
  # frac.of.genome is the change of it being heads
  out$expected.overlap <- out$frac.of.genome * n.dar
  out$binomial.p <- pbinom(n.overlap, n.dar, out$frac.of.genome, lower.tail = F) %>%
    format(scientific = 4)
  out <- as.data.frame(out)
  return(out)
}

de.da.intersect.titrate.plot.venn <- function(in.tsv, out.dir = NULL) {
  res <- read.table(in.tsv, header = T)
  root.name <- basename(in.tsv) %>% tools::file_path_sans_ext()

  if (is.null(out.dir)) {
    out.dir <- paste0(dirname(in.tsv))
  }

  lapply(1:nrow(res), function(i) {
    res <- res[i, ]
    anchor.set <- res$anchor.set
    root.name <- paste0(root.name, "_padj", res$padj, "_log2FC", res$log2FC, "_slot_", anchor.set)

    venn.list <- list(DEG = 1, DAR = 1)
    names(venn.list) <- paste0(res$direction, names(venn.list))
    venn.dir <- paste0(out.dir, "/venn/")
    if (grepl("^fold\\d+$", anchor.set))
      venn.dir <- paste0(venn.dir, "/bg/")

    dir.create(venn.dir, showWarnings = F, recursive = T)
    out.file <- paste0(venn.dir, "/", root.name, "_", res$cluster, "_", res$direction, ".pdf")
    area.vector <- c(res$n.de, res$n.da, res$n.de.w.da)

    grobs <- VennDiagram::venn.diagram(
      venn.list,
      fill = utilsFanc::color.hue.fc(n = 2, palette = "R4.fc1"),
      alpha = c(0.5, 0.5), lwd = 0, filename = NULL,
      imagetype = "svg", height = 1.2, width = 1.2, units = "in",
      cex = 0.5, fontfamily = "Arial",
      cat.cex = 0.5, cat.fontfamily = "Arial", cat.default.pos = "outer", margin = 0.06,
      rotation.degree = 180,
      cat.pos = c(27, -27),
      cat.dist = c(0.055, 0.055),
      direct.area = T, area.vector = area.vector)

    cairo_pdf(filename = out.file, width = 1.2, height = 1.2)
    try(grid.draw(grobs))
    dev.off()
    return()
  })
}
