de.slot.assert <- function(de, is.s2b = F, slot) {
  if (is.s2b) {
    de <- list(de)
    names(de) <- de[[1]]$root.name
  }

  bExists <- lapply(de, function(s2b) {
    slot %in% names(s2b)
  }) %>% unlist()

  if(sum(!bExists) > 0) {
    stop(paste0("The following clusters do not have the slot '", slot, "':\n",
                paste0(names(de)[!bExists], collapse = "\n")))
  }
  return(0)
}

de.sync.name <- function(de) {
  # change s2b$root.name to the corresponding name in names(de)
  names <- names(de)
  de <- lapply(names(de), function(name) {
    a2b <- de[[name]]
    a2b$root.name <- name
    return(a2b)
  })
  names(de) <- names
  return(de)
}

de.type <- function(de) {
  genes <- de[[1]]$bulkNorm$gene[1]
  if (grepl(":", genes)) return("da")
  return("de")
}

de.grid <- function(so, master.ident="seurat_clusters", assay, slot, master.ident.groups = NULL, group.by,
                    ident.df, outlist=NULL, thread = 12, test.use = "wilcox", latent.vars = NULL, equal.cells,
                    logfc.threshold = 0.25, min.pct = 0.1) {
  # ident.df: 2 columns, must be ident.1 and ident.2
  Idents(so) <- master.ident
  if (is.null(master.ident.groups))
    master.ident.groups <- so@active.ident %>% as.character() %>% unique()
  grid.de.and.sum <- utilsFanc::safelapply(master.ident.groups, function(i) {
    de.and.sum <- ident.df %>% split(., f= paste0(.$ident.1, ":", .$ident.2)) %>%
      lapply(function(x) {
        de.df <- FindMarkers.fanc(object = so, subset.ident = i,
                                     ident.1 = x$ident.1, ident.2 = x$ident.2,
                                     group.by = group.by, assay = assay, slot = slot, test.use=test.use,
                                     latent.vars = latent.vars, add.cdr = T, equal.cells = equal.cells,
                                  logfc.threshold = logfc.threshold, min.pct = min.pct)
        # de.df$p_val <- NULL
        de.df <- utilsFanc::add.column.fanc(df1 = de.df, df2 = data.frame(gene = rownames(de.df)), pos = 1)
        sum.df <- de.df
        rownames(de.df) <- NULL
        comp <- paste0(x$ident.1, ":", x$ident.2)
        colnames(de.df) <- c("gene", paste0("p_", comp),paste0("logFC_", comp), paste0("pct.", x$ident.1),
                             paste0("pct.", x$ident.2), paste0("p.adj_", comp))
        colnames(sum.df) <- c("gene", "p","logFC", "pct.1","pct.2", "p.adj")
        sum.df <- utilsFanc::add.column.fanc(sum.df, df2 = data.frame(comp = rep(comp, nrow(sum.df))), after = "gene")
        return(list(de.df=de.df, sum.df=sum.df))
      })
    # print("miao")
    de.df <- Reduce(full_join, lapply(de.and.sum, function(x) return(x$de.df)))
    de.df <- utilsFanc::add.column.fanc(de.df, df2 = data.frame(master.ident = rep(i,nrow(de.df))),
                                        after = "gene")
    sum.df <- Reduce(rbind, lapply(de.and.sum, function(x) return(x$sum.df)))
    sum.df <- utilsFanc::add.column.fanc(sum.df, df2 = data.frame(master.ident = rep(i,nrow(sum.df))),
                                        after = "gene")

    return(list(de.df = de.df,sum.df = sum.df))
  }, threads = thread)

  grid.de <- Reduce(rbind, lapply(grid.de.and.sum, function(x) return(x$de.df)))
  grid.sum <- Reduce(rbind, lapply(grid.de.and.sum, function(x) return(x$sum.df)))

  rownames(grid.de) <- NULL
  rownames(grid.sum) <- NULL
  if (!is.null(outlist)) {
    write.table(grid.de, file = utilsFanc::plot.name.construct(outlist = outlist, sub.name = "_grid_de.tsv"),
                col.names = T, row.names = F, sep = "\t", quote = F)
    write.table(grid.sum, file = utilsFanc::plot.name.construct(outlist = outlist, sub.name = "_grid_sum.tsv"),
                col.names = T, row.names = F, sep = "\t", quote = F)
  }
  return(list(grid.de = grid.de, grid.sum = grid.sum))
}

s2b.add.de.grid <- function(pbl, de.grid, de.grid.master.ident) {
  s2b.names <- names(pbl)
  for (s2b.name in s2b.names) {
    pbl[[s2b.name]]$res.exp <- de.grid$grid.sum %>%
      filter(master.ident == sub(paste0(de.grid.master.ident, "_"), "", s2b.name)) %>%
      rename(log2FoldChange = logFC, padj = p.adj, pvalue = p) %>%
      mutate(log2FoldChange = log2(2.71828^log2FoldChange))
  }
  return(pbl)
}

de.grid.sum <- function(de.grid, outlist.upset = NULL,outlist.df=NULL, p.adj.cutoff=1, p.cutoff = 1,
                        upset.x.max = NULL, upset.y.max = NULL, debug = F) {
  # de.grid: output from function de.grid
  # stop("untested function")
  order.df <- de.grid$grid.sum %>% select(master.ident) %>% unique()
  sum.df <- de.grid$grid.sum %>%
    filter(p.adj <= p.adj.cutoff, p <= p.cutoff) %>%
    group_by(master.ident, comp) %>%
    summarise(n.up = sum(logFC > 0), n.down = sum(logFC < 0), n.sig = n(),
              genes.up = paste0(gene[logFC > 0], collapse = ";"),
              genes.down = paste0(gene[logFC < 0], collapse = ";"),
              genes.sig = paste0(gene, collapse = ";"))
  sum.df <- left_join(order.df, sum.df)
  if (!is.null(outlist.df)) {
    write.table(sum.df %>% select(master.ident, comp, n.up, n.down, n.sig),
                utilsFanc::plot.name.construct(outlist.df, sub.name = "_grid_summed.tsv"),
                quote = F, col.names = T, row.names = F, sep = "\t")
  }
  if (!is.null(outlist.upset)) {
    sub.f.draw.venn <- function(df, genes.col,outlist) {
      sum.list <- df %>% split(., f= factor(.$master.ident, levels = unique(.$master.ident)))
      #print("br 1")
      set.list <- lapply(sum.list, function(x) return(as.data.frame(x)[,genes.col] %>%
                                                        strsplit(";")))
      #print("br 2")
      comp.list <- lapply(sum.list, function(x) return(as.data.frame(x)[,"comp"] ))
      #print("br 3")
      outlist$root.name.external <- ifelse(is.null(outlist$root.name.external), genes.col,
                                           paste0(outlist$root.name.external, "_", genes.col))
      #print("br 3.5")
      multi.venn(x.list.list = set.list, name.vec.list = comp.list,
                 outlist = outlist, use.upset = T, mainbar.y.max = upset.y.max,
                 set_size.scale_max = upset.x.max, debug=debug,
                 plot.titles.vec = paste0("cluster ",unique(df$master.ident)) )
      #print("br 4")
      return(NULL)
    }
    lapply(c("up","down", "sig"), function(x) {
      sub.f.draw.venn(df = sum.df, genes.col = paste0("genes.", x),
                      outlist = outlist.upset)
    })
  }


  return(sum.df)
}

de.grid.violin <- function(de.grid.sum, so, s2b.list = NULL, so.feature.split.by=NULL, so.feature.plot = NULL, out.dir= NULL, sub.height=2.5, assay, slot,
                           sub.width=2.5, master.ident, comp, p.adj.cutoff=1, p.cutoff=1,
                           add.feature.plot = F, add.genes = NULL, exclude.pattern = NULL,
                           feature.order, pt.size = 0.001, debug=F, cores = 1,  raster = F, format = "png",
                           ...) {
  # stop("plotting engine needs an upgrade to handle the new plot.panel.list functions")
  df <- de.grid.sum %>%
    filter(p <= p.cutoff, p.adj <= p.adj.cutoff)
  df %>% split(., f = factor(.$master.ident, unique(.$master.ident))) %>%
    mclapply(function(mi) {
      ident <- mi$master.ident[1]
      print(ident)
      so.mi <- so[,so@meta.data[,master.ident] == ident]
      mi %>% split(., f = factor(.$comp, unique(.$comp))) %>%
        lapply(function(cp) {
          # print("debug")
          so.mi.cp <- so.mi[, so.mi@meta.data[, comp] %in% unlist(strsplit(unique(cp$comp), ":"))]

          Idents(so.mi.cp) <- master.ident
          so.mi.cp$comp <- cp$comp[1]

          genes <- cp$gene %>% c(add.genes)

          if (!is.null(exclude.pattern)) {
            genes <- genes[!grepl(exclude.pattern, genes)]

          }

          if (debug == T)
            genes <- genes[1:2]
          plots <- Seurat::VlnPlot(object = so.mi.cp, features = genes, group.by = comp,
                                   idents = cp$master.ident %>% unique(), assay = assay,
                                   # split.by = comp, slot = slot, split.plot = T,
                                   combine = F)

          title <- paste0(genes, " ", formatC(cp$p, format = "e", digits = 2), " \n", cp$master.ident, " ", cp$comp)

          if (!is.null(s2b.list)) {
            bulkNorm.df <- s2b.list[[paste0(master.ident, "_", ident)]]$bulkNorm
            rownames(bulkNorm.df) <- bulkNorm.df$gene
            sample.1 <- cp$comp[1] %>% sub(":.+", "", .)
            sample.2 <- cp$comp[1] %>% sub(".+:", "", .)
            exp.1 <- bulkNorm.df[genes, paste0("bulkNorm_", sample.1)] %>% floor()
            exp.2 <- bulkNorm.df[genes, paste0("bulkNorm_", sample.2)] %>% floor()
            title <- paste0(title, " = ", exp.1, " : ", exp.2)
          }

          plots <- lapply(seq_along(plots), function(i) {
            plot <- plots[[i]] + xlab("") +
              theme(legend.position = "none",title =element_text(size=8, face='bold'),
                    axis.text.x = element_text(size = 8, angle = 0)) +
              ggtitle(title[i])
            return(plot)
          })
          n.split = 1
          if (add.feature.plot == T) {
            if (is.null(so.feature.plot))
              so.feature.plot <- so
            if (is.null(so.feature.split.by))
              so.feature.split.by <- comp
            genes.in.so.f <- rownames(so.feature.plot)
            genes.feature.plot <- genes
            genes.feature.plot[!genes.feature.plot %in% genes.in.so.f] <- "B2m"
            feature.plots <- plot.panel.list(panel.list = genes.feature.plot, obj = so.feature.plot, order = feature.order, assay = assay,
                                             split.by = so.feature.split.by,return.list = T,
                                             pt.size = pt.size, raster = raster, label.size = 3, ...) %>% Reduce(c, .)

            plots.final <- lapply(seq_along(plots), function(i) {
              sub.list <- list(plots[[i]], feature.plots[[2*i-1]], feature.plots[[2*i]])
              return(sub.list)
            }) %>% unlist(recursive = F)
            n.split = 3
          } else {
            plots.final <- plots
          }

          if (!is.null(out.dir)) {
            trash <- wrap.plots.fanc(plot.list = plots.final, sub.height = sub.height,
                                    sub.width = sub.width,n.split = n.split,
                                    plot.out = paste0(out.dir, "/", master.ident, "_", cp$master.ident[1],"/", cp$comp[1], ".",format))
          }

          return(NULL)
        }) %>% Reduce(c,.) %>% return()

    }, mc.cores = cores, mc.cleanup = T) %>% Reduce(c,.) %>% return()
}

de.pos.filter.marker <- function(so, ident, s2b.obj) {

}

FindMarkers.fanc <- function (object, set.ident = NULL, ident.1 = NULL, ident.2 = NULL, group.by = NULL, active.ident = "seurat_clusters",
                              subset.ident = NULL, assay = NULL, slot = "data", reduction = NULL,
                              features = NULL, logfc.threshold = 0.25, test.use = "wilcox",
                              min.pct = 0.1, min.diff.pct = -Inf, verbose = TRUE, only.pos = FALSE,
                              max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
                              min.cells.feature = 3, min.cells.group = 3, pseudocount.use = 1,
                              add.cdr = F, equal.cells = F,
                              ...) {
  if (!is.null(set.ident)) {
    Idents(object) <- set.ident
  }
  if (test.use == "MAST" || add.cdr == T) {
    raw.mat <- object@assays$RNA@counts
    object[["cdr"]] <-  (colSums(raw.mat > 0)/nrow(raw.mat)) %>% scale()
    remove(raw.mat)
    print("FANC: adding cdr")
    latent.vars <- c(latent.vars, "cdr")
  }

  if (equal.cells == T) {
    if (is.null(group.by))
      stop("you have to offer group.by to use equal.cells")
    filter.list <- list(ident.1, subset.ident)
    names(filter.list) <- c(group.by, active.ident)
    n.cells.ident.1 <- get.cell.names.seurat(so = object, filter.list = filter.list, style = "Seurat") %>% length()

    filter.list <- list(ident.2, subset.ident)
    names(filter.list) <- c(group.by, active.ident)
    n.cells.ident.2 <- get.cell.names.seurat(so = object, filter.list = filter.list, style = "Seurat") %>% length()

    max.cells.per.ident <- min(n.cells.ident.1, n.cells.ident.2)


    print(paste0("equal.cells triggered. max is ",  max(n.cells.ident.1, n.cells.ident.2), ", min is ",
          max.cells.per.ident, ", max.cells.per.ident is set to: ", max.cells.per.ident))
  }

  return(Seurat::FindMarkers(object = object, ident.1 = ident.1, ident.2 = ident.2, group.by = group.by,
                             subset.ident = subset.ident, assay = assay, slot = slot,
                             reduction = reduction, features = features,
                             logfc.threshold = logfc.threshold,
                             test.use = test.use, min.pct = min.pct,
                             min.diff.pct = min.diff.pct, verbose = verbose,
                             only.pos = only.pos, max.cells.per.ident = max.cells.per.ident,
                             random.seed = random.seed, latent.vars = latent.vars,
                             min.cells.feature = min.cells.feature,
                             min.cells.group = min.cells.group, pseudocount.use = pseudocount.use, ...))

}

# FindAllMarkers.fanc <- function(object,  assay = NULL,  features = NULL,  logfc.threshold = 0.25,
#                                 test.use = "wilcox",  slot = "data",  min.pct = 0.1,
#                                 min.diff.pct = -Inf,  node = NULL,  verbose = TRUE,
#                                 only.pos = FALSE,  max.cells.per.ident = Inf,
#                                 random.seed = 1,  latent.vars = NULL,
#                                 min.cells.feature = 3,  min.cells.group = 3,
#                                 pseudocount.use = 1,  return.thresh = 0.01,) {
#
# }

de.filter.pos <- function(grid, so, s2bo.list, master.ident, ident.1, comp, plot.out) {
  # grid: object generated by de.grid ()
  samples <- strsplit(comp, ":") %>% unlist()
  markers <- lapply(samples, function(x) {
    Idents(so) <- "sample"
    marker <- FindMarkers(object = so, ident.1 = ident.1, only.pos = T,
                          group.by = master.ident, subset.ident  = x) %>%
      filter(p_val_adj < 0.05) %>% rownames()
  }) %>% Reduce(intersect, .)

  ident.1.markers <- FindMarkers(object = so, ident.1 = ident.1, only.pos = T, max.cells.per.ident = 300) %>%
    filter(p_val_adj < 0.05) %>% rownames()
  trash <- plot.panel.list.m(panel.list = ident.1.markers, obj = so,
                    order = F, assay = "RNA", split.by = "sample",
                    plot.out = plot.out, raster = T, return.list = F)
  de.ident <- grid %>%
    filter(`p.adj_PyMT:WT` < 0.05, master.ident == ident.1)
  de.up <- de.ident %>% filter(`logFC_PyMT:WT` > 0) %>% pull(gene)
  de.down <- de.ident %>% filter(`logFC_PyMT:WT` < 0) %>% pull(gene)
  de.up %>% str()
  tt <- intersect(de.up, ident.1.markers)
  t <- intersect(de.down, ident.1.markers)
  trash <- plot.panel.list.m(panel.list = t, obj = so,
                             order = F, assay = "RNA", split.by = "sample",
                             plot.out = plot.out, raster = T, return.list = F, pt.size = 0.8)
  trash <- plot.panel.list.m(panel.list = tt, obj = so,
                             order = F, assay = "RNA", split.by = "sample",
                             plot.out = "~/test/mpage/miao_2.pdf", raster = T, return.list = F, pt.size = 0.8)

}

replicated.markers <- function(so, cluster.ident= "seurat_clusters", cluster.idents = NULL, assay, slot,
                               group.ident, group.idents = NULL, p.adj.cutoff, threads = NULL,
                               plot.dir = NULL, debug = F) {
  meta.df <- so@meta.data %>% factor2character.fanc()
  if (is.null(cluster.idents))
    cluster.idents <- meta.df[, cluster.ident] %>% unique()
  if (is.null(group.idents))
    group.idents <- meta.df[, group.ident] %>% unique()
  if (is.null(threads))
    threads <- min(length(cluster.idents), 18)


  marker.list <- mclapply(cluster.idents, function(cl) {
    marker.intersect <- mclapply(group.idents, function(x) {
      Idents(so) <- group.ident
      marker <- FindMarkers(object = so, ident.1 = cl, only.pos = T, group.by = cluster.ident,
                            assay = assay, slot = slot, subset.ident = x, test.use = "wilcox")
      if (debug == T)
        return(marker)
      marker <- marker %>% filter(p_val_adj < p.adj.cutoff) %>% rownames()
      if (!is.null(plot.dir))
      trash <- plot.panel.list.m(panel.list = marker, obj = so, order = F, assay = assay, slot = slot,
                                 split.by = "sample", raster = T, pt.size = 0.8,
                                 plot.out = paste0(plot.dir,"/markers_intersect_", cluster.ident, "_", cl, ".pdf") )
      return(marker)
    }, mc.cores = 2)
    if (debug == F)
      marker.intersect <- marker.intersect %>% Reduce(intersect, .)
    return(marker.intersect)
  }, mc.cores = threads)

  print("miao")
  names(marker.list) <- paste0(cluster.ident, "_",cluster.idents)
  return(marker.list)

}


grid.diagnostics <- function(grids, so, norms, tests,  master.ident = "seurat_clusters", comp,
                            marker.list = NULL , n.col = NULL,
                            sub.idents=NULL, quantiles = c(0.95, 0.999, 1),
                            s2b.list = NULL, p.adj.cutoff, master.dir, threads = NULL) {
  # each cluster is a directory
  # each normalization is 2 files (one on itself the other one using s2b.list).
  Rpl.genes <- rownames(so) %>% .[grepl("Rp[ls]", .)]

  sample.1 <- comp %>% sub(":.+$", "",.)
  sample.2 <- comp %>% sub("^.+:", "",.)

  if (is.character(grids))
    grids <- readRDS(grids)
  if (is.null(sub.idents)) {
    sub.idents <- so@meta.data %>% factor2character.fanc() %>% pull(!!master.ident) %>% unique() %>%
      gtools::mixedsort()
  }

  if (is.null(threads))
    threads <- length(sub.idents)
  utilsFanc::safelapply(sub.idents, function(sub.ident) {

    # marker.genes <- marker.list[[paste0(master.ident, "_", sub.ident)]]

    p.list <- lapply(norms, function(norm) {
      # first we prepare a dataframe for scatter plots:
      mean.df <- so.2.bulk(so = so, assay = norm, slot = "data", sub.idents = sub.ident,
                           ident = master.ident, group.by = "sample",
                           take.mean = T, coldata.columns = NULL)$bulk.mat %>%
        as.data.frame()
      mean.df$gene <- rownames(mean.df)
      df.list <- list(mean.df = mean.df)
      if (!is.null(s2b.list)) {
        bulk.df <- s2b.list[[paste0(master.ident, "_", sub.ident)]]$bulkNorm
        colnames(bulk.df) <- colnames(bulk.df) %>% sub("bulkNorm_", "", .)
        df.list <- list(mean.df = mean.df, bulk.df = bulk.df)
      }


      lapply(seq_along(df.list), function(i) {
        df <- df.list[[i]]
        df.name <- names(df.list)[i]
        lapply(seq_along(quantiles), function(j) {
          quantile <- quantiles[j]
          p.list.tests <- lapply(tests,function(test) {
            de.genes <- grids[[norm]][[test]]$grid.de %>% filter(master.ident == sub.ident) %>%
              .[.[, paste0("p.adj_", comp)] < p.adj.cutoff, "gene"]
            # de.markers <- intersect(de.genes, marker.genes)
            de.Rpl <- intersect(de.genes, Rpl.genes)

            gene.list <- list(# all.genes = NULL,
                              # marker.genes = marker.genes,
                              # Rpl.genes = Rpl.genes,
                              de.genes = de.genes
                              # de.markers = de.markers,
                              # de.Rpl = de.Rpl
                              )

            p.list <- lapply(seq_along(gene.list), function(i) {
              p <- xy.plot(df = df, x = sample.2, y = sample.1, highlight.var = "gene", highlight.values = gene.list[[i]],
                           quantile.limit = quantile, outfile = NULL)
              p <- p + ggtitle(paste0(sub.ident, " ",norm, " q:",quantile, " ",  test, " " ,names(gene.list)[i]))
              return(p)
            })

            return(p.list)

          }) %>% Reduce(c, .)

          return(p.list.tests)
        }) %>% Reduce(c, .) %>% return()

      }) %>% Reduce(c, .) %>% return()
    }) %>% Reduce(c,.)
    trash <- wrap.plots.fanc(plot.list = p.list, n.col = n.col,
                             plot.out = paste0(master.dir, "/", master.ident, "_", sub.ident, ".png"), page.limit = 2000000)
    return()
  }, threads = threads)

  return()
}

grid.diagnostics.2 <- function(grids, so, norms, tests,  master.ident = "seurat_clusters", comp,
                               gene.list = NULL , n.col = NULL, n.split = 1,
                               slide.column = NULL, slide.window = 10, add.direction = T,
                               genes.label = NULL,
                               sub.idents=NULL, max.quantile = NULL,
                               s2b.list = NULL, p.adj.cutoff, master.dir, root.name = NULL,
                               threads = 1,
                               use.plotly = F) {
  # different from grid.diagnostics: it tries to put all clusters into 1 file
  # the format of gene.list: list(up = c("gene1", "Gene2"), down = c("miao", "wang"))
  label.var <- NULL
  if (!is.null(genes.label))
    label.var <- "gene"

  sample.1 <- comp %>% sub(":.+$", "",.)
  sample.2 <- comp %>% sub("^.+:", "",.)

  if (is.character(grids))
    grids <- readRDS(grids)
  if (is.null(sub.idents)) {
    sub.idents <- so@meta.data %>% factor2character.fanc() %>% pull(!!master.ident) %>% unique() %>%
      gtools::mixedsort()
  }

  lapply(norms, function(norm) {
    mat <- so.2.bulk.m(so = so, root.name = "miao", work.dir = NULL,
                       assay = norm, slot = "data", cluster.ident = master.ident,
                       clusters = sub.idents,
                       group.ident = "sample", groups = c(sample.1, sample.2),
                       take.mean = T)$bulk.mat

    # df.list <- utilsFanc::safelapply(sub.idents, function(sub.ident) {
    #   mean.df <- so.2.bulk(so = so, assay = norm, slot = "data", sub.idents = sub.ident,
    #                        ident = master.ident, group.by = "sample",
    #                        take.mean = T, coldata.columns = NULL)$bulk.mat %>%
    #     as.data.frame()
    #   mean.df$gene <- rownames(mean.df)
    #   return(mean.df)
    # }, threads = threads)

    # names(df.list) <- sub.idents
    lapply(tests, function(test) {
      p.list <- utilsFanc::safelapply(sub.idents, function(sub.ident) {
        # df <- df.list[[sub.ident]]
        df <- mat[, grepl(paste0(sub.ident, "\\.\\."), colnames(mat))] %>% as.data.frame()
        colnames(df) <- sub(paste0(sub.ident, "\\.\\."), "", colnames(df))
        df$gene <- rownames(df)
        if (is.null(gene.list)) {
          if (!is.null(slide.column)) {
            grid <- grids[[norm]][[test]]$grid.sum %>%
              .[.$comp == comp & .$master.ident == sub.ident,]
            if (! slide.column %in% colnames(grid)) {
              stop(paste0("slide.column: ", slide.column, " not found in grid"))
            }
            if (add.direction == T) {
              # grid[, slide.column] <- grid[, slide.column] * (grid$logFC/abs(grid$logFC))
              grid$direction <- grid$logFC/abs(grid$logFC)
            } else {
              grid$direction <- 1
            }
            grid <- grid %>% group_by(direction)
            grid <- grid %>% mutate(rank = rank(!!as.name(slide.column))) %>%
              mutate(pct = 100*rank/max(rank)) %>%
              mutate(group = ceiling(pct/slide.window) * direction) %>%
              ungroup() %>% as.data.frame()
            gene.list <- grid$gene %>% split(., f = grid$group)
            names(gene.list)[as.numeric(names(gene.list)) < 0] <- as.numeric(names(gene.list)[as.numeric(names(gene.list)) < 0] ) - 100
            gene.list <- gene.list[order(abs(as.numeric(names(gene.list))))]
            names(gene.list) <- paste0("g", names(gene.list))
          } else {
            de.genes <- grids[[norm]][[test]]$grid.de %>% filter(master.ident == sub.ident) %>%
              .[.[, paste0("p.adj_", comp)] < p.adj.cutoff, "gene"] %>% .[!is.na(.)]
            gene.list <- list(# all.genes = NULL,

              de.genes = de.genes

            )
          }
        }

        p.list <- lapply(seq_along(gene.list), function(i) {
          p <- xy.plot(df = df, x = sample.2, y = sample.1,
                       highlight.var = "gene",
                       plotly.var = "gene",
                       highlight.values = gene.list[[i]],
                       highlight.only = use.plotly,
                       label.var = label.var, label.values = unlist(genes.label),
                       text.size = 3, use.repel = T,
                       quantile.limit = max.quantile,
                       outfile = NULL, title = paste0(master.ident, " ", sub.ident, " ", length(gene.list[[i]])))
          # p <- p + ggtitle(paste0(sub.ident, " ",norm, " q:", max.quantile, " ",  test, " " ,names(gene.list)[i]))
          return(p)
        })
        return(p.list)
      }, threads = threads) %>% Reduce(c, .)

      if (is.null(root.name))
        root.name <- master.ident
      root <- paste0(root.name, "_", comp, "..", norm, "..", test, "..q", max.quantile)
      if (use.plotly == T) {
        plot.out <- paste0(master.dir, "/",root, ".html")
      } else {
        plot.out <- paste0(master.dir, "/",root, ".png")
      }
      trash <- wrap.plots.fanc(p.list, tooltip = "gene", plot.out = plot.out, n.col = n.col,
                               n.split = n.split)
      return()
    })
  })

}

grids.plot.panel <- function(grids, norms, tests, comp,
                            so,
                            padj.cutoff, max.plots = 50,
                            plot.type = c("de.genes", "up.genes", "down.genes"),
                            order = F, assay, reduction = "umap",
                            cluster.ident, clusters = NULL,
                            sample.order = NULL,
                            plot.dir, threads = 6,
                            violin = F, violin.clusters = NULL) {
  # honest you should try to plot just 1 norm and 1 test.
  # because plot.panel.list can be a slow and dangerous function
  # it might choke R. It also hurts your eye to go through too many plots
  lapply(norms, function(norm) {
    lapply(tests, function(test) {
      grid <- grids[[norm]][[test]]
      comparison <- comp
      fake.pbl <- grid$grid.sum %>% filter(comp == comparison) %>%
        dplyr::rename(log2FoldChange = logFC, padj = p.adj, pvalue = p) %>%
        mutate(log2FoldChange = log2(2.71828^log2FoldChange)) %>%
        split(., f= .$master.ident)
      fake.pbl <- fake.pbl[gtools::mixedorder(names(fake.pbl))]
      if (!is.null(clusters))
        fake.pbl <- fake.pbl[names(fake.pbl) %in% clusters]
      if (length(fake.pbl) < 1) {
        stop("length(fake.pbl) < 1")
      }
      names(fake.pbl) <- paste0(cluster.ident, "_", names(fake.pbl))
      fake.pbl <- lapply(fake.pbl, function(x) {
        y = list(res.exp = x)
        return(y)
      })
      trash <- deseq2.plot.panel(pbl = fake.pbl, so = so, padj.cutoff = padj.cutoff,
                                 max.plots = max.plots, plot.type = plot.type, order = order,
                                 assay = assay, reduction = reduction, cluster.ident = cluster.ident,
                                 sample.order = sample.order, threads = threads, plot.dir = plot.dir,
                                 root.name = paste0(norm, "_", test),
                                 violin = violin, violin.clusters = violin.clusters)
      return()

    })
  })
}

cluster.marker.fanc <- function(marker.df,so=NULL, clusters, pct.1.cutoff = 0.70,
                                pct.2.cutoff = 0.2, n.plot = 20, plot.dir=NULL) {
  marker.df$gene <- rownames(marker.df) %>% sub("\\.\\d+$", "", .)
  marker.df <- marker.df %>% mutate(cluster = as.character(cluster)) %>%
    filter(cluster %in% clusters, pct.1 > pct.1.cutoff, pct.2 < pct.2.cutoff) %>%
    dplyr::group_by(cluster) %>% arrange(p_val) %>% slice_head(n = 20) %>% ungroup()
  if (!is.null(plot.dir)) {
    marker.df %>% split(f = marker.df$cluster) %>%
      lapply(function(x) {
        trash <- plot.panel.list(panel.list = x$gene, obj = so, order = F, assay = "SCT",
                                 plot.out = paste0(plot.dir, "/markers_cluster_", x$cluster, ".png"))
        return()
      })
  }

  return(marker.df)

}

assay.mat.gen <- function(so, out.Rds, cluster.ident, clusters = NULL, group.by, groups = NULL,
                            assay, slot, thread = 6 ) {
  if (is.character(so))
    so <- readRDS(so)
  if (is.null(clusters))
    clusters <- so@meta.data[, cluster.ident] %>% unique() %>% as.character() %>% gtools::mixedsort()
  s2b.list <- mclapply(clusters, function(cluster) {
    utilsFanc::t.stat(paste0(cluster.ident, "_", cluster))
    s2b <- so.2.bulk(so = so, assay = "RNA", slot = "data", ident = cluster.ident,
                     sub.idents = cluster, group.by = group.by, groups = groups,
                     take.mean = T, coldata.columns = NULL)
    return(s2b)
  }, mc.cleanup = T, mc.cores = thread)
  names(s2b.list) <- paste0(cluster.ident, "_", clusters)
  saveRDS(s2b.list, out.Rds)
  return(s2b.list)
}

deseq2.xyplot <- function(pbl, publication = F,
                          pbl.slot = "res.exp", summ.slot = "summary", # behaviour change: now prefer to use the summary slot.
                          # padj.cutoff = 0.05, use.p = F,
                          pbl.soup = NULL, de.grid = NULL, de.grid.master.ident, # mostly not used anymore
                          comp.vec, comp.meta = NULL,
                          transformation = NULL, quantile.limit = NULL,
                          no.highlight = F, highlight.by.overlap.gr = NULL,
                          highlight.genes = NULL, highlight.color = "red",
                          density.filter = NULL,
                          add.label = F, label.list = NULL,
                          plotly.label.all = F, plotly.highlight.only = T,
                          device = c("png"), plot.dir, root.name = NULL, # device could also be html
                          n.col = NULL, threads = 8, save.rds = F, plot.each.comp = F,
                          sub.width = NULL, sub.height = NULL,
                          theme.text.size = NULL, highlight.ptsize = NULL,
                          ...) {
  # pbl: pseudobulk list object generated by bulk.list()
  # de.grid: extend support to objects generated by de.grid()
  # p.flag <- "padj"
  # if (use.p == T)
  #   p.flag <- "p"

  if (publication) {
    if (is.null(sub.width)) sub.width <- 1.6
    if (is.null(sub.height)) sub.height <- 1.6
    pt.size = 0.008
    pt.shape = 18
    if (is.null(highlight.ptsize))
      highlight.ptsize = 0.05
    theme.text.size <- 6
  } else {
    if (is.null(sub.width)) sub.width <- 4
    if (is.null(sub.height)) sub.height <- 4
    pt.size = 0.05
    pt.shape = 19
    # highlight.ptsize = NULL
    # theme.text.size <- NULL
  }

  if (is.null(root.name)) {
    root.name <- summ.slot
  }
  if (add.label == T) {
    label.var <- "gene"
    if (is.null(label.list)) {
      if (!is.null(highlight.genes)) {
        label.list <- rep(list(highlight.genes), length(pbl))
        names(label.list) <- names(pbl)
      }
    } else if (!is.list(label.list)) {
      label.list <- rep(list(label.list), length(pbl))
      names(label.list) <- names(pbl)
    }
  } else {
    label.var <- NULL
  }
  if (is.character(pbl))
    pbl <- readRDS(pbl)
  pbl <- pbl[sapply(pbl, is.list)]
  if (length(pbl) < 1) {
    stop("length(pbl) < 1")
  }
  s2b.names <- names(pbl)
  if (is.null(s2b.names)) {
    stop("names(pbl) is null")
  }

  n.split <- length(pbl.soup) + 1

  pl <- lapply(comp.vec, function(comp) {
    lapply(device, function(dev) {
      if (dev == "html") {
        highlight.only <- plotly.highlight.only
      } else {
        highlight.only <- F
      }
      pl <- utilsFanc::safelapply(s2b.names, function(s2b.name) {
        # s2b <- pbl[[s2b.name]]
        cat(paste0("Processing cluster: ", s2b.name, "\n"))
        if (no.highlight == F) {
          if (!is.null(highlight.by.overlap.gr)) {
            loci <- pbl[[s2b.name]][[pbl.slot]]$gene %>% utilsFanc::loci.2.gr()
            loci <- loci %>% subsetByOverlaps(highlight.by.overlap.gr)
            de.df <- data.frame(gene = utilsFanc::gr.get.loci(loci),
                                direction = "gr")
            highlight.var <- "gene"
          } else if (!is.null(highlight.genes)) {
            de.df <- data.frame(gene = highlight.genes, direction = "geneset")
            highlight.var <- "gene"
          } else if (!is.null(summ.slot)) {
            warning("defaulted to use summ.slot (behavior change since 2022-07-07.")
            summ <- pbl[[s2b.name]][[summ.slot]]
            if (is.null(summ)) {
              stop(paste0("slot '", summ.slot, "' not found"))
            }
            
            de.df <- pbl[[s2b.name]][[pbl.slot]] %>% filter(gene %in% summ$de.genes)
            if (nrow(de.df) == 0) {
              highlight.var <- NULL
              no.highlight <- T
            } else {
              de.df$direction <- NA
              de.df$direction[de.df$gene %in% summ$up.genes] <- "up"
              de.df$direction[de.df$gene %in% summ$down.genes] <- "down"
              highlight.var <- "gene"
              de.df <- de.df[, c("gene", "direction")]
            }
          } else {
            stop("old code that hasn't been examined forever")
            # if (!is.null(de.grid)) {
            #   pbl[[s2b.name]][[pbl.slot]] <- de.grid$grid.sum %>%
            #     filter(master.ident == sub(paste0(de.grid.master.ident, "_"), "", s2b.name)) %>%
            #     rename(log2FoldChange = logFC, padj = p.adj, pvalue = p) %>%
            #     mutate(log2FoldChange = log2(2.71828^log2FoldChange))
            # }
            # if (!is.null(highlight.by.overlap.gr)) {
            #   stop("this part of the function is not working. because sample.y and sample.x is updated to incoporate multiple samples")
            #   gr <- highlight.by.overlap.gr
            #   de.df <- pbl[[s2b.name]][[pbl.slot]] %>% utilsFanc::loci.2.df(loci.col.name = "gene", return.gr = T)
            #   de.df <- subsetByOverlaps(de.df, gr, ignore.strand = T) %>% as.data.frame()
            #   de.df$log2FoldChange <- log2pm(de.df[, sample.y]/de.df[, sample.x])
            # } else {
            #   if (use.p == F) {
            #     de.df <- pbl[[s2b.name]][[pbl.slot]] %>% filter(padj < padj.cutoff) %>% dplyr::select(gene, log2FoldChange)
            #   } else {
            #     de.df <- pbl[[s2b.name]][[pbl.slot]] %>% filter(pvalue < padj.cutoff) %>% dplyr::select(gene, log2FoldChange)
            #   }
            # }
            #
            # de.df$log2FoldChange[de.df$log2FoldChange > 0] <- "up"
            # de.df$log2FoldChange[de.df$log2FoldChange < 0] <- "down"
            # de.df <- de.df[, c("gene", "log2FoldChange")]
            # highlight.var <- "gene"
          }
        } else {
          highlight.var <- NULL
        }

        pbl.2 <- c(pbl[s2b.name], pbl.soup)
        pl <- lapply(names(pbl.2), function(s2b.name.2) {
          s2b <- pbl.2[[s2b.name.2]]
          df <- s2b$bulkNorm
          colnames(df) <- sub("^bulkNorm_", "", colnames(df))
          if (no.highlight == F)
            df <- left_join(df, de.df, by = "gene")

          comp.l <- list (y = comp %>% sub(":.+$", "",.),
                          x =  comp %>% sub("^.+:", "",.))

          if (!is.null(comp.meta)) {
            for (i in 1:length(comp.l)) {
              samples.in.comp <- s2b$coldata %>% .[.[, comp.meta] %in% comp.l[[i]], , drop = F] %>% rownames()
              df[, comp.l[[i]]] <- df[, samples.in.comp, drop = F] %>% utilsFanc::pmean()
            }
          }
          df <- cbind(data.frame(id = 1:nrow(df)), df)
          if (dev == "html") {
            ref.file <- paste0(c(root.name, "_", s2b.name.2, "_", paste0(comp, "_", summ.slot), ".tsv"), collapse = ".") %>%
              paste0(plot.dir, "/", .)
            dir.create(plot.dir, showWarnings = F, recursive = T)
            write.table(df[!is.na(df$direction),], ref.file, row.names = F, col.names = T, sep = "\t", quote = F)
          }
          if (add.label) {
            if (!is.null(label.list[[s2b.name]])) {
              label.values <- label.list[[s2b.name]]
            } else {
              label.values <- s2b[[summ.slot]]$de.genes
              if (length(label.values) < 1) {
                label.var <- NULL
              }
            }
          } else {
            label.values <- NULL
          }

          color.map <- utilsFanc::gg_color_hue(2)
          names(color.map) <- c("down", "up")
          if (!is.null(highlight.color) && !is.null(highlight.genes)) {
            df$direction <- "hl"
            color.map <- highlight.color
            names(color.map) <- "hl"
          }

          p <- xy.plot(df = df, x = comp.l$x, y = comp.l$y, transformation = transformation,
                       quantile.limit = quantile.limit, density.filter = density.filter,
                       highlight.var = highlight.var, highlight.values = de.df$gene,
                       label.var = label.var, label.values = label.values,
                       use.geom_text = T,
                       plotly.var = "gene", highlight.color.var = "direction",
                       color.map = color.map,
                       highlight.only = highlight.only,
                       plotly.label.all = plotly.label.all,
                       show.highlight.color.var = F, raster = T,
                       pt.size = pt.size, pt.shape = pt.shape,
                       highlight.ptsize = highlight.ptsize,
                       theme.text.size = theme.text.size,
                       ...) +
            ggtitle(paste0(root.name, "  ", s2b.name.2))
          if (publication) {
            p <- p %>% utilsFanc::theme.fc.1(italic.x = F, font = NULL)
            # p <- p + theme(aspect.ratio = 1, plot.title = element_text(
            #   size = 6, hjust = 0.5, margin = margin(b = 0, unit = "in")))
          }


          return(p)
        })
        return(pl)
      }, threads = threads) %>% Reduce(c, .)
      system(paste0("mkdir -p ", plot.dir))

      plot.out <- paste0(root.name, "_", comp, "_", summ.slot, ".", dev) %>%
        paste0(plot.dir, "/", .)
      if (plot.each.comp)
        wrap.plots.fanc(plot.list = pl, plot.out = plot.out, n.split = n.split, n.col = n.col, tooltip = "gene",
                        sub.width = sub.width, sub.height = sub.height,
                        pdf.use.cairo = F, pdf.by.row = T)
      if (save.rds == T)
        saveRDS(pl, paste0(tools::file_path_sans_ext(plot.out), ".Rds"))
      return(pl)
    }) %>% Reduce(c, .) %>% return()

  }) %>% Reduce(c, .) %>% return()
  plot.out <- paste0(plot.dir, "/", root.name, ".", device)
  wrap.plots.fanc(plot.list = pl, plot.out = plot.out, tooltip = "gene", n.col = n.col,
                  sub.width = sub.width, sub.height = sub.height, pdf.use.cairo = F, pdf.by.row = T)
  if (length(pl) == 1) {
    pl <- pl[[1]]
  }
  invisible(pl)
}

# pbl.annotate.peaks <- function(pbl, genome) {
#   pbl <- lapply(pbl, function(a2b) {
#     a2b$bulkNorm %>% head()
#     gr <- a2b$bulkNorm %>% utilsFanc::loci.2.df(loci.col.name = "gene", return.gr = T)
#     anno <- utilsFanc::annotate.fanc(gr = gr, genome = genome)
#     anno
#   })
#   return(pbl)
# }

deseq2.de.from.density <- function(pbl, comp,
                                   slot.name = NULL,
                                   transformation = NULL, two.sided = F,
                                   density.filter = 0.01,
                                   genome,
                                   plot.dir, plot.root.name = NULL,
                                   do.html = F,
                                   threads = 8) {
  sample.y <- comp %>% sub(":.+$", "",.)
  sample.x <- comp %>% sub("^.+:", "",.)
  if (is.null(slot.name))
    slot.name <- comp
  res.list <- utilsFanc::safelapply(pbl, function(s2b) {
    if (is.null(s2b$de.dens))
      s2b$de.dens <- list()

    # dens.df <- s2b$bulkNorm

    bulkNorm.transform <- s2b$bulkNorm
    bulkNorm.transform[, sample.x] <- bulkNorm.transform[, sample.x] %>% transformation()
    bulkNorm.transform[, sample.y] <- bulkNorm.transform[, sample.y] %>% transformation()

    dens.df <- xy.plot(df = bulkNorm.transform, x = sample.x, y = sample.y,
                       density.filter = density.filter, density.filter.2sided = two.sided,
                       return.df = T)

    p.png <- xy.plot(df = bulkNorm.transform, x = sample.x, y = sample.y,
                     # transformation = transformation,
                     highlight.var = "gene", highlight.values = dens.df$gene,
                     return.df = F, title = s2b$root.name,
                     density.filter.2sided = two.sided)
    if (do.html == T) {
      p.html <- xy.plot(df = bulkNorm.transform, x = sample.x, y = sample.y,
                        # transformation = transformation,
                        highlight.var = "gene", highlight.values = dens.df$gene,
                        highlight.only = T, plotly.var = "gene", density.filter.2sided = two.sided,
                        return.df = F, title = s2b$root.name)

    } else {
      p.html <- NULL
    }

    dens.list <- list()
    dens.list$dens.df <- dens.df
    dens.list$de <- s2b$bulkNorm %>%
      filter(gene %in% dens.df$gene) %>%
      mutate(dist = !!as.name(sample.y) - !!as.name(sample.x)) %>%
      mutate(sign = if_else(dist > 0, "+", "-")) %>%
      arrange(sign, !!as.name(sample.x))

    dens.list$de.transformed <- bulkNorm.transform %>%
      filter(gene %in% dens.df$gene) %>%
      mutate(dist = !!as.name(sample.y) - !!as.name(sample.x)) %>%
      mutate(sign = if_else(dist > 0, "+", "-")) %>%
      arrange(sign, !!as.name(sample.x))

    dens.list$filter <- density.filter
    dens.list$transformed <- bulkNorm.transform

    s2b$de.dens[[slot.name]] <- dens.list
    dir.create(paste0(plot.dir, "/xlsx/"), showWarnings = F, recursive = T)

    anno <- utilsFanc::annotate.fanc(gr = dens.list$de.transformed %>%
                                       utilsFanc::loci.2.df(loci.col.name = "gene",
                                                            return.gr = T,
                                                            remove.loci.col = T),
                                     genome = genome)
    anno <- anno %>% as.data.frame()
    xlsx::write.xlsx(x = anno,
                     file = paste0(plot.dir, "/xlsx/",
                                   plot.root.name, "_", s2b$root.name, ".xlsx"),
                     row.names = F, col.names = T)

    mat <- dens.list$de.transformed %>% `rownames<-`(., .$gene) %>%
      .[, c(sample.x, sample.y)] %>% as.matrix()

    hm <- ComplexHeatmap::Heatmap(mat,
                                  cluster_rows = T, show_row_dend = F, show_row_names = F,
                                  cluster_columns = F)
    trash <- save.base.plot(p = hm, file = paste0(plot.dir, "/hm/",
                                         plot.root.name, "_", s2b$root.name, ".png"),
                            # height = 15*nrow(dens.list$de.transformed)
                            height = 800)
    res <- list(s2b = s2b, p.png = p.png, p.html = p.html)
    return(res)
  })
  pbl <- res.list %>% lapply(function(x) return(x$s2b))
  if (!is.null(plot.root.name))
    root.name <- paste0(plot.dir, "/xy/", plot.root.name, "..")
  else
    root.name <- paste0(plot.dir, "/xy/")

  plot.out.png <- paste0(root.name, "_dens_", density.filter, ".png")
  p.png <- res.list %>% lapply(function(x) return(x$p.png))
  wrap.plots.fanc(p.png, plot.out = plot.out.png)

  if (do.html == T) {
    plot.out.html <- paste0(root.name, "_dens_", density.filter, ".html")
    p.html <- res.list %>% lapply(function(x) return(x$p.html))
    wrap.plots.fanc(p.png, plot.out = plot.out.html)
  }

  return(pbl)
}

deseq2.summary <- function(pbl, res.slot = "res", summary.slot = "summary",
                           rank.by = "padj",
                           padj.cutoff = 0.05, log2fc.cutoff = 1,
                           min.log2.exp = NULL, min.log2.exp.samples,
                           annot = F, annot.genome,
                           exclude = NULL, exclude.by.overlap = F,
                           force = T, save = F,
                           stats.file = NULL, gene.out.dir = NULL,
                           use.topn.p = NULL) {
  # rank: padj or log2FoldChange
  if (save) {
    if (is.null(gene.out.dir)) {
      work.dir <- attr(pbl, "work.dir")
      if (is.null(work.dir)) {
        stop("attr(pbl, 'work.dir') returned null")
      } else {
        gene.out.dir <- paste0(work.dir, "/", summary.slot, "/")
      }
    }
    if (is.null(stats.file)) {
      stats.file <- paste0(gene.out.dir, "/stats.tsv")
    }
  }

  pbl <- pbl[sapply(pbl, is.list)]
  s2b.names <- names(pbl)

  if (!is.null(gene.out.dir)) {
    system(paste0("rm -rf ", gene.out.dir))
    system(paste0("mkdir -p ", gene.out.dir))
  }

  pbl <- lapply(s2b.names, function(s2b.name) {
    s2b <- pbl[[s2b.name]]
    if (is.null(s2b[[res.slot]])) {
      return(s2b)
    }
    if (is.null(s2b[[summary.slot]]) || force == T) {
      de.df <- s2b[[res.slot]] %>% as.data.frame()
      if (is.null(de.df$gene))
        de.df <- de.df %>%  dplyr::mutate(., gene = rownames(.))

      if (!is.null(min.log2.exp)) {
        utilsFanc::check.intersect(min.log2.exp.samples, "min.log2.exp.samples",
                                   colnames(de.df), "colnames(de.df)")
        bPass <- rowSums(de.df[, min.log2.exp.samples] + 1 > 2^min.log2.exp) > 0
        de.df <- de.df[bPass, ]
      }

      if (!is.null(exclude)) {
        if (exclude.by.overlap)
          stop("exclude by overlap has not been developped yet")
        de.df <- de.df %>% filter(!gene %in% exclude)
      }

      de.df <- de.df %>% filter(abs(log2FoldChange) > log2fc.cutoff)
      if (is.null(use.topn.p)) {
        de.df <- de.df %>% filter(padj < padj.cutoff)
      } else {
        de.df <- de.df %>% split(f = de.df$log2FoldChange > 0) %>%
          lapply(function(de.df) {
            de.df <- de.df %>% dplyr::arrange(pvalue)
            de.df <- de.df[1:min(use.topn.p, nrow(de.df)),]
          }) %>% do.call(rbind, .)
        rownames(de.df) <- NULL
        if (rank.by == "padj") rank.by <- "pvalue"
        padj.cutoff <- NA
      }

      if (rank.by %in% c("padj", "pvalue"))
        de.df <- de.df %>% arrange(!!as.name(rank.by))
      else if (rank.by == "log2FoldChange")
        de.df <- de.df %>% arrange(desc(abs(!!as.name(rank.by))))

      de.df <- de.df %>% dplyr::select(gene, log2FoldChange)
      s2b[[summary.slot]] <- list(padj.cutoff = padj.cutoff,
                          log2fc.cutoff = log2fc.cutoff,
                          n.de = nrow(de.df),
                          n.up = nrow(de.df %>% filter(log2FoldChange > 0)),
                          n.down = nrow(de.df %>% filter(log2FoldChange < 0)),
                          de.genes = de.df$gene,
                          up.genes = de.df %>% filter(log2FoldChange > 0) %>% pull(gene),
                          down.genes = de.df %>% filter(log2FoldChange < 0) %>% pull(gene))

      if (!is.null(gene.out.dir)) {
        lapply(c("de.genes", "up.genes", "down.genes"), function(type) {
          write(s2b[[summary.slot]][[type]], paste0(gene.out.dir, "/", s2b.name, "_", type, ".txt"), sep = "\n")
          if (grepl(":", s2b[[summary.slot]][[type]][1])) {
            # if ATAC
            gr <- s2b[[res.slot]] %>% as.data.frame()
            if (is.null(gr$gene)) {
              gr <- gr %>% dplyr::mutate(., gene = rownames(.))
            }
            gr <- gr %>% dplyr::filter(gene %in% s2b[[summary.slot]][[type]]) %>%
              utilsFanc::loci.2.df(loci.col.name = "gene", remove.loci.col = T, return.gr = T)
            if(annot) {
              gr <- utilsFanc::gr.fast.annot(obj = gr, genome = annot.genome, use.klraps = F)
            }
            out.bed <- paste0(gene.out.dir, "/", s2b.name, "_", type, ".bed")
            utilsFanc::write.zip.fanc(
              gr, bed.shift = T,
              out.file = out.bed)
            write.table(utilsFanc::gr2df(gr), out.bed, quote = F, sep = "\t",
                        col.names = T, row.names = F)
          }
          return()
        })
      }
    }
    return(s2b)
  })
  names(pbl) <- s2b.names
  if (!is.null(gene.out.dir)) {
    if (is.null(stats.file))
      stats.file <- paste0(gene.out.dir, "/stats.tsv")
  }
  if (!is.null(stats.file)) {
    stats <- lapply(names(pbl), function(name) {
      if (is.null(pbl[[name]][[summary.slot]])) {
        return()
      }
      stats <- pbl[[name]][[summary.slot]][c("n.de", "n.up", "n.down")] %>% as.data.frame()
      stats <- utilsFanc::add.column.fanc(df1 = stats, df2 = data.frame(name = name), pos = 1)
      return(stats)
    }) %>% Reduce(rbind, .)
    system(paste0("mkdir -p ", dirname(stats.file)))
    write.table(stats, stats.file, row.names = F, col.names = T, sep = "\t", quote = F)
    deseq2.stats.barplot(stats.df = stats.file, padj.cutoff = padj.cutoff,
                         log2fc.cutoff = log2fc.cutoff, topn.cutoff = use.topn.p)
  }
  return(pbl)
}



deseq2.stats.barplot <- function(stats.df, cluster.ident.rm = c("seurat_clusters", "Clusters"),
                                 padj.cutoff = NULL, log2fc.cutoff = NULL, topn.cutoff = NULL,
                                 plot.out = NULL) {
  if (is.character(stats.df)) {
    if (is.null(plot.out)) {
      plot.out <- paste0(tools::file_path_sans_ext(stats.df), ".pdf")
    }
    stats.df <- read.table(stats.df, header = T)
  }
  required.cols <- c("name", "n.up", "n.down")
  utilsFanc::check.intersect(required.cols, "required columns",
                             colnames(stats.df), "colnames(stats.df)")
  for (rm in cluster.ident.rm) {
    stats.df$name <- sub(paste0(rm, "_*"), "", stats.df$name)
  }
  stats.df$name <- factor(stats.df$name, levels = gtools::mixedsort(unique(stats.df$name)))
  stats.df <- stats.df %>%
    reshape2::melt(measure.vars = c("n.up", "n.down"),
                   variable.name = "type", value.name = "n.genes")
  p <- ggplot(stats.df, aes(x = name, y = n.genes, fill = type)) +
    geom_bar(stat = 'identity', position = "dodge") +
    geom_text(aes(label = n.genes), position = position_dodge(width = 0.8))

  title <- paste("padj", padj.cutoff, "log2fc", log2fc.cutoff, "topn", topn.cutoff)
  p <- p + ggtitle(title)

  if (!is.null(plot.out)) {
    dir.create(dirname(plot.out), showWarnings = F, recursive = T)
    ggsave(plot.out, p, width = 0.5 * length(unique(stats.df$name)) + 2, height = 5, dpi = 100)
  }

  invisible(p)
}
deseq2.filter.summ.by.overlap <- function(pbl, slot = "summary", slot.new,
                                          filter.by, center.pbl = F, center.filter = F,
                                          out.dir = NULL) {
  # filter.by: a list of loci that matches the root.names of each cluster.
  # eg: filter.by = list(GMP = c("chr1:100-200", "chr2:200-300"),
  #                      MEP = c("chr1:300-400", "chr3:400-500"))
  clusters <- sapply(pbl, function(x) return(x$root.name)) %>%
    `names<-`(NULL)
  if (!is.list(filter.by)) {
    filter.by <- rep(list(filter.by), length(clusters))
    names(filter.by) <- clusters
  }

  if (!identical(sort(clusters), sort(names(filter.by)))) {
    stop("!identical(sort(clusters), sort(names(filter.by)))")
  }
  pbl <- lapply(pbl, function(a2b) {
    new.summ <- a2b[[slot]]
    for (type in c("de", "up", "down")) {
      g <- paste0(type, ".genes")
      gf <- paste0(type, ".genes.failed")
      bPass <- new.summ[[g]] %>%
        in.loci(x = ., y = filter.by[[a2b$root.name]],
                  center.x = center.pbl, center.y = center.filter)
      new.summ[[gf]] <- new.summ[[g]][!bPass]
      new.summ[[g]] <- new.summ[[g]][bPass]

      n <- paste0("n.", type)
      nf <- paste0("n.", type, ".failed")
      new.summ[[n]] <- length(new.summ[[g]])
      new.summ[[nf]] <- length(new.summ[[gf]])

      if (!is.null(out.dir)) {
        dir.create(out.dir, showWarnings = F, recursive = T)
        lapply(c("filtered", "failed"), function(pass.type) {
          suffix <- if_else(pass.type == "failed", ".failed", "")
          loci <- new.summ[[paste0(type, ".genes", suffix)]]
          out.file <- paste0(out.dir, "/",a2b$root.name, "_", type, "_", pass.type, ".genes.txt")
          write(loci, out.file, sep = "\n")
          out.bed <- sub(".txt$", ".bed", out.file)
          utilsFanc::write.zip.fanc(df = utilsFanc::loci.2.gr(loci), out.file = out.bed, bed.shift = T)
          return()
        })
      }
    }
    new.summ$filtered.by.overlap <- T
    a2b[[slot.new]] <- new.summ
    return(a2b)
  })
  if (!is.null(out.dir)) {
    stats.file <- paste0(out.dir, "/stats_filtered.tsv")
    stats <- lapply(names(pbl), function(name) {
      if (is.null(pbl[[name]][[slot.new]])) {
        return()
      }
      stats <- pbl[[name]][[slot.new]][c("n.de", "n.up", "n.down",
                                         "n.de.failed", "n.up.failed", "n.down.failed")] %>%
        as.data.frame()
      stats <- utilsFanc::add.column.fanc(df1 = stats, df2 = data.frame(name = name), pos = 1)
      return(stats)
    }) %>% Reduce(rbind, .)
    system(paste0("mkdir -p ", dirname(stats.file)))
    write.table(stats, stats.file, row.names = F, col.names = T, sep = "\t", quote = F)
  }
  return(pbl)
}

deseq2.plot.panel <- function(pbl, so, max.plots = 50, summary.slot = "summary", padj.cutoff,
                              plot.type = c("up.genes", "down.genes"),
                              order = F, assay, reduction = "umap", cluster.ident =NULL,
                              sample.order = NULL,
                              n.col = NULL, page.limit = 20, panel.label = T,
                              pt.size = 1.2,
                              violin.swap = F,
                              plot.dir, root.name = NULL, threads = 6,
                              violin = F, violin.clusters = NULL,
                              violin.single.cluster = F, ident.root.name, ...) {
  if (!is.null(sample.order)) {
    so <- so[, so$sample %in% sample.order]
    so$sample <- factor(so$sample, levels = sample.order)
  }
  n.split <- so$sample %>% unique() %>% length()
  trash <- utilsFanc::safelapply(names(pbl), function(name) {
    s2b <- pbl[[name]]
    if (is.null(summary.slot)) {
      df <- s2b$res.exp %>% filter(padj < padj.cutoff) %>%
        select(gene, padj, log2FoldChange) %>% arrange(padj)
      gene.list <- list()
      gene.list$de.genes <- df$gene
      gene.list$up.genes <- df %>% filter(log2FoldChange > 0) %>% pull(gene)
      gene.list$down.genes <- df %>% filter(log2FoldChange < 0) %>% pull(gene)
    } else {
      if (!is.null(s2b[[summary.slot]])) {
        gene.list <- s2b[[summary.slot]][c("de.genes", "up.genes", "down.genes")]
        gene.list <- lapply(gene.list, function(genes) {
          genes <- s2b$res %>% as.data.frame() %>% .[genes, ] %>%
            dplyr::arrange(padj) %>% rownames()
          return(genes)
        })
      } else {
        stop("!is.null(s2b[[summary.slot]])")
      }
    }
    lapply(names(gene.list), function(type) {
      if (! type %in% plot.type)
        return()
      genes <- gene.list[[type]]
      if (length(genes) > max.plots)
        genes <- genes[1:max.plots]
      if (violin == F) {
        if (!is.null(root.name)) {
          name <- paste0(root.name, "..", name)
        }
        plot.out <- paste0(plot.dir, "/", name, "..", assay, "..", type, "..", order, ".pdf")

        trash <- plot.panel.list(panel.list = genes, obj = so, order = order, assay = assay, split.by = "sample",
                                 plot.out = plot.out, raster = T, label = panel.label,
                                 pt.size = pt.size, auto.adjust.raster.pt.size = F, page.limit = page.limit,
                                 max.quantile = 0.99,
                                 n.split = n.split, reduction = reduction, ident = cluster.ident, n.col = n.col,
                                 pdf.by.row = F, pdf.use.cairo = F,
                                 ...)

      } else {
        print(genes)
        if (!is.null(root.name)) {
          name <- paste0(root.name, "..", name)
        }
        plot.out <- paste0(plot.dir, "/", name, "..", assay, "..", type, "..", order, ".vln.pdf")

        if (violin.single.cluster == T) {
          violin.clusters <- name %>% sub(ident.root.name, "", .) %>% sub("_", "", .)
        }
        if (violin.swap == T) {
          so <- so[, so@meta.data[, cluster.ident] %in% violin.clusters]
          split.by <- cluster.ident
          split.order <- violin.clusters
          cluster.ident <- "sample"
          violin.clusters <- sample.order
        } else {
          split.by <- "sample"
          split.order <- sample.order
        }
        trash <- plot.panel.list(panel.list = genes, obj = so, order = order, assay = assay,
                                 split.by = split.by, split.order = split.order,
                                 ident = cluster.ident, cluster = violin.clusters,
                                 n.col = 6,
                                 plot.out = plot.out,
                                 raster = T,
                                 pt.size = 0.5, auto.adjust.raster.pt.size = F,
                                 page.limit = 6, max.quantile = 0.99, violin = T,add.median = F,
                                 sub.width = 15, sub.height = 3.5, ...)
      }
      return()
    })
    return()
  }, threads = threads)
  return()
}

de.soup.stability <- function(main.pbl, clusters = NULL, rhos, titrate.dir, rank.cutoff.vec,
                              do.summary.main = F, do.summary.rhos = F, padj.cutoff = 0.05,
                              out.dir, threads = 4) {
  if (is.character(main.pbl)) {
    main.pbl <- readRDS(main.pbl)
  }
  if (!is.null(clusters))
    main.pbl <- main.pbl[names(main.pbl) %in% clusters]
  if (do.summary.main == T)
    main.pbl <- deseq2.summary(pbl = main.pbl, padj.cutoff = padj.cutoff, force = T, stats.file = NULL, gene.out.dir = NULL)

  bl.rho.list <- utilsFanc::safelapply(rhos, function(rho) {
    bl.rho <- readRDS(paste0(titrate.dir, "/rho_", rho,"/bulk.list.Rds"))
    if (do.summary.rhos == T) {
      bl.rho <- deseq2.summary(pbl = bl.rho, padj.cutoff = padj.cutoff, force = T, stats.file = NULL, gene.out.dir = NULL)
    }
    return(bl.rho)
  }, threads = threads)
  names(bl.rho.list) <- paste0("rho_", rhos)

  s2b.names <- names(main.pbl)
  lapply(s2b.names, function(s2b.name) {
    s2b <- main.pbl[[s2b.name]]
    if (is.null(s2b$summary)) {
      stop("s2b$summary is NULL")
    }

    lapply(c("up.genes", "down.genes"), function(type) {
      if (type == "up.genes")
        trend.fun <- `>`
      else
        trend.fun <- `<`
      genes <- s2b$summary[[type]]
      if (is.null(genes)) {
        warning(paste0("s2b.name: ", s2b.name, "; type: ", type, " returned NULL"))
        return()
      }
      n.genes <- length(genes)
      rank.cutoff.vec <- c(rank.cutoff.vec, n.genes) %>% sort() %>% unique()
      for (rank.cutoff in rank.cutoff.vec) {
        if (rank.cutoff > n.genes)
          break
        genes.2 <- genes[1:rank.cutoff]
        n.genes.2 <- length(genes.2)
        df <- utilsFanc::safelapply(rhos, function(rho) {
          s2b.rho <- bl.rho.list[[paste0("rho_",rho)]][[s2b.name]]
          if (is.null(s2b.rho$summary[[type]])) {
            warning(paste0("s2b.name: ", s2b.name, "; rho: ",rho,  "; type: ", type, " returned NULL"))
            return()
          }
          still.sig <- genes.2 %>% .[. %in% s2b.rho$summary[[type]]]
          same.trend <- s2b.rho$res.exp %>% filter(gene %in% genes.2) %>%
            filter(trend.fun(log2FoldChange, 0)) %>% pull(gene)
          same.trend.nosig <- same.trend %>% .[! . %in% still.sig]
          opposite.trend <- s2b.rho$res.exp %>% filter(gene %in% genes.2) %>%
            filter(!trend.fun(log2FoldChange, 0)) %>% pull(gene)

          df <- data.frame(rho = rho, n.still.sig = length(still.sig), n.same.trend = length(same.trend),
                           n.same.trend.nosig = length(same.trend.nosig),
                           n.opposite.trend = length(opposite.trend),
                           still.sig = still.sig %>% paste0(collapse = ","),
                           same.trend = same.trend %>% paste0(collapse = ","),
                           same.trend.nosig = same.trend.nosig %>% paste0(collapse = ","),
                           opposite.trend = opposite.trend %>% paste0(collapse = ","))
          return(df)
        },threads = threads) %>% Reduce(rbind, .)
        df.head <- data.frame(rho = "ori", n.still.sig = n.genes.2,
                              n.same.trend = n.genes.2, n.same.trend.nosig = 0,
                              n.opposite.trend = 0,
                              still.sig = genes.2 %>% paste0(collapse = ","), same.trend.nosig = "",
                              same.trend = genes.2 %>% paste0(collapse = ","), opposite.trend = "")
        df <- rbind(df.head, df)
        file.name <- paste0(out.dir, "/", s2b.name, "/top", rank.cutoff, "..", type, ".tsv")
        dir.create(dirname(file.name), showWarnings = F, recursive = T)
        write.table(df, file.name , sep = "\t",
                    row.names = F, col.names = T, quote = F)
      }
      return()
    })
    return()

  })
  return()

}

deseq2.write.results <- function(pb.m, write.raw = F,
                                 out.dir, root.name = NULL) {
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  system(paste0("mkdir -p ", out.dir))
  trash <- lapply(names(pb.m), function(name) {
    s2b <- pb.m[[name]]
    write.table(s2b$res.exp, paste0(out.dir, "/", root.name, "_", name, "_deseq2.tsv"),
                quote = F, row.names = F, col.names = T, sep = "\t")
    if (write.raw) {
      df <- cbind(data.frame(gene = rownames(s2b$bulk.mat)), s2b$bulk.mat)
      write.table(df, paste0(out.dir, "/", root.name, "_", name, "_raw.tsv"),
                  quote = F, row.names = F, col.names = T, sep = "\t")
    }
  })
}

deseq2.pipe.ez <- function(soi, plot.only = T, pb.m = NULL, samples, pbl.soup = NULL,
                           assay, order = F, reduction, comp.vec, comp.meta = NULL,
                           summary.slot = "summary",
                           cluster.ident, clusters,
                           panel.page.limit = 20, n.col = NULL, panel.label = T,
                           violin.cluster.ident = NULL,violin.single.cluster = F, violin.clusters = NULL, violin.swap = F,
                           work.dir,
                           plot.xy = T, plot.panel = T, max.plots = 50, plot.violin = T,
                           padj.cutoff = 0.05, threads = 1, ...) {
  if (plot.only == F) {
    stop("the deseq2 part has not been generalized")
    pb.m <- bulk.list(so = soi, assay = assay, slot = "counts",
                      cluster.ident = cluster.ident, clusters = clusters,
                      coldata.columns = c("rep", "type"),
                      design.formula = ~ rep + type,
                      contrast = c("type", "tumor", "ctrl"),
                      work.dir = work.dir,
                      plot.dir = paste0(work.dir, "/plots"), threads = threads)
  } else {
    if (is.null(pb.m))
      pb.m <- readRDS(paste0(work.dir, "/bulk.list.Rds"))
  }
  if (plot.xy == T) {
    try({
      deseq2.xyplot(pbl = pb.m, pbl.soup = pbl.soup, summ.slot = summary.slot,
                    plot.dir = paste0(work.dir, "/plots_xy_w_soup/"),
                    comp.vec = comp.vec, comp.meta = comp.meta,
                    root.name = paste0("linear_0.99"), quantile.limit = 0.99,
                    padj.cutoff = padj.cutoff, threads = threads)
    })

    try({
      deseq2.xyplot(pbl = pb.m, pbl.soup = pbl.soup, summ.slot = summary.slot,
                    plot.dir = paste0(work.dir, "/plots_xy_w_soup/") ,
                    comp.vec = comp.vec, comp.meta = comp.meta,
                    root.name = paste0("log2"), quantile.limit = NULL,
                    transformation = log2,
                    padj.cutoff = padj.cutoff, threads = threads)
    })
  }

  if (plot.panel == T) {
    try({
      deseq2.plot.panel(pbl = pb.m, so = soi, padj.cutoff = padj.cutoff, plot.type = c("up.genes", "down.genes"),
                        sample.order = samples, cluster.ident = cluster.ident, summary.slot = summary.slot,
                        max.plots = max.plots, order = order, assay = assay, reduction = reduction,
                        threads = threads, plot.dir = paste0(work.dir, "/panel/"), root.name = summary.slot,
                        n.col = n.col, page.limit = panel.page.limit, panel.label = panel.label, ...)
    })
  }
  if (plot.violin == T) {
    try({
      # we need to specify violin.clusters because if we plot all clusters,
      # it will be too wide (too many violins), especially if
      # you have multiple samples.
      if (is.null(violin.cluster.ident))
        violin.cluster.ident <- cluster.ident
      deseq2.plot.panel(pbl = pb.m, so = soi, padj.cutoff = padj.cutoff, plot.type = c("up.genes", "down.genes"),
                        sample.order = samples,
                        max.plots = max.plots, order = F, assay = assay, reduction = reduction,
                        plot.dir = paste0(work.dir, "/violin/"), violin = T, violin.swap = violin.swap,
                        cluster.ident = violin.cluster.ident, violin.single.cluster = violin.single.cluster,
                        violin.clusters = violin.clusters, ident.root.name = cluster.ident,
                        threads = threads)
    })
  }
  return()
}

deseq2.assess.core <- function(dds, plot.out, pca.ntop = 10000, 
                               pca.groupings = c("type", "rep"),
                               return.plot.list = F, blind = T,
                               add.hm = F, add.3d = F,
                               publication = F, pub.pt.size = 1.2,
                               pub.width = 1.5, pub.height = 1.2,
                               add.jitter = F, jitter.seed = 42) {
  
  dir.create(dirname(plot.out), showWarnings = F, recursive = T)
  if (is.null(pca.ntop))
    pca.ntop <- 10000
  pl <- list()
  ntd <- normTransform(dds)
  rld <- rlog(dds,blind = blind, fitType = "local")
  
  if (!publication) {
    pl$ntd <- mean.sd.plot(mat = assay(ntd))
    pl$rld <- mean.sd.plot(mat = assay(rld))
  }

  data <- plotPCA.fanc(rld, intgroup = pca.groupings,
                     ntop = pca.ntop, returnData = T, dims.return = 1:3)
  pcaData <- data$d
  genes <- data$genes
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  saveRDS(list(d = data, pct.var = percentVar),
          tools::file_path_sans_ext(plot.out) %>% paste0(".Rds"))
  flag <- F
  if (sum(grepl("^PC", colnames(pcaData))) == 2) {
    flag <- T
  }
  rm(data)
  if (add.hm == T) {
    mat <- DESeq2::counts(object = dds, normalized = T)

    mat <- mat[genes, ]
    mat <- scale.row.fanc(mat)
    p <- Heatmap(matrix = mat, show_row_names = F, show_row_dend = F, show_column_names = T,
                 column_title = paste0("ntop: ", pca.ntop))
    save.base.plot(p = p, file = utilsFanc::insert.name.before.ext(name = plot.out, insert = "hm"),
                   width = 700, height = 1000)
  }
  if (!is.null(pcaData$rep)) {
    pcaData[, "rep"] <- pcaData$rep %>% sub("rep", "", .)
  }
  if (add.3d == T) {
    p.plotly <- plot_ly(x = pcaData[, "PC1"], y = pcaData[, "PC2"], z = pcaData[, "PC3"], color = pcaData[, "de_v3"], symbol = pcaData[, "type"],
                        symbols = c('circle','x'), text = pcaData[, "rep"]) %>% add_markers() %>% add_text() %>%
      layout(scene = list(xaxis = list(title = paste0("PC1: ",percentVar[1],"% variance")),
                          yaxis = list(title = paste0("PC2: ",percentVar[2],"% variance")),
                          zaxis = list(title = paste0("PC3: ",percentVar[3],"% variance"))))
    htmlwidgets::saveWidget(as_widget(p.plotly),
                            normalizePath.partial(path = tools::file_path_sans_ext(plot.out) %>% paste0(".html")))
  }
  pc.list <- list(c(1,2), c(1,3), c(2,3))
  if (flag == T || publication == T)
    pc.list <- list(c(1,2))

  for (i in pc.list) {
    comp <- paste0("PC", i)
    p.name <- paste0(comp[1], ".", comp[2])
    if (add.jitter) {
      stop("jitter not validated")
      set.seed(seed = jitter.seed)
      for (i in 1:2) {
        vec <- pcaData[, comp[i]]
        span <- max(vec) - min(vec)
        noise <- rnorm(length(vec), mean = 0, sd = 0.05 * span)
        pcaData[, comp[i]] <- pcaData[, comp[i]] + noise
      }
    }
    
    if (length(pca.groupings) >= 1) {
      if (length(pca.groupings) == 1) {
        pl[[p.name]] <- ggplot(pcaData, aes_string(comp[1], comp[2], color=pca.groupings[1]))
      } else {
        pl[[p.name]] <- ggplot(pcaData, aes_string(comp[1], comp[2], color=pca.groupings[1],
                                     shape=pca.groupings[2]))
      }
      if (publication) {
        pl[[p.name]] <- pl[[p.name]] +
          geom_point(size = pub.pt.size) +
          xlab(paste0(comp[1], ": ",percentVar[i[1]],"% variance")) +
          ylab(paste0(comp[2], ": ",percentVar[i[2]],"% variance"))

        pl[[p.name]] <- utilsFanc::theme.fc.1(pl[[p.name]], italic.x = F, remove.axis.titles = F) +
          theme(legend.position = "right", legend.direction = "vertical", aspect.ratio = 1,
                plot.margin = unit(c(0, 0.15, 0, 0.05), "in"))
        
        if (length(pca.groupings) == 2) {
          groups <- pcaData[, pca.groupings[2]] %>% unique() %>% sort()
          shape.map <- 1:length(groups)
          if (length(groups) == 2) 
            shape.map <- c(1, 4)
          
          names(shape.map) <- groups
          pl[[p.name]] <- pl[[p.name]] + scale_shape_manual(values = shape.map)
        }
        
        # Hyung Joo's palette.
        groups <- pcaData[, pca.groupings[1]] %>% unique() %>% sort()

        if (length(groups) == 4) {
          color.map <- c("#440077", "#00AAFF", "#AAAAAA", "#FF00AA")
          names(color.map) <- groups
          pl[[p.name]] <- pl[[p.name]] + scale_color_manual(values = color.map)
        }
        
      }  else {
        pl[[p.name]] <- pl[[p.name]] +
          geom_point(size=3) +
          xlab(paste0(comp[1], ": ",percentVar[i[1]],"% variance")) +
          ylab(paste0(comp[2], ": ",percentVar[i[2]],"% variance")) +
          theme_bw() +
          theme(aspect.ratio = 1) +
          ggtitle("ntop: " %>% paste0(pca.ntop))
      }
      
    }
    if (!is.na(pca.groupings[3])) {
      y.span <- max(pcaData[, comp[2]]) - min(pcaData[, comp[2]])
      pl[[p.name]] <-  pl[[p.name]] +
        geom_text(aes_string(label = pca.groupings[3]), color = "black", nudge_y = 0.05 * y.span)
    }
  }


  # for (grouping in pca.groupings[-1]) {
  #   pl[[paste0("pca.", grouping)]] <-
  #     ggplot(pcaData, aes_string("PC1", "PC2", color=pca.groupings[1], shape=grouping)) +
  #     geom_point(size=3) +
  #     xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #     ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #     coord_fixed() +
  #     theme(aspect.ratio = 1)
  # }
  #
  if (publication) {
    pl[[1]] <- pl[[1]] + theme(legend.text = element_text(margin = margin(b = 0.02, unit = "in")))
    ggsave(plot.out, pl[[1]], width = pub.width, height = pub.height, device = cairo_pdf)
  } else {
    trash <- wrap.plots.fanc(pl, plot.out = plot.out)
  }
  

  if (return.plot.list == T)
    return(pl)
  else
    return()
}

deseq2.homer <- function(dds, fg.vec, n.bg.each = 25, no.replace = F,
                         samples.for.bg, scatter.xy.df, scatter.meta = NULL,
                      work.dir, write.log = T,
                      genome, size = "given", thread, denovo = F, other.params = "",
                      find.motifs.path = "/bar/cfan/software/homer/bin/findMotifsGenome.pl",
                      homer.path = "/bar/cfan/software/homer/bin/", run = T) {
  norm.mat <- DESeq2::counts(dds, normalized = T)
  bg.dir <- paste0(work.dir, "/bg/")
  prep.dir <- paste0(work.dir, "/prep/")
  homer.dir <- paste0(work.dir, "/homer/")
  system(paste0("mkdir -p ", bg.dir, " ", prep.dir, " ", homer.dir))
  fg.bed <- paste0(prep.dir, "/fg.bed")
  trash <- utilsFanc::write.zip.fanc(df = utilsFanc::loci.2.df(loci.vec = fg.vec)[, 1:3], out.file = fg.bed)

  if (n.bg.each > 0) {
    bg.vec <- bg.gen.2(mat = norm.mat, fg.vec = fg.vec, n.bg.each = n.bg.each, samples.use = samples.for.bg,
                       no.replace = no.replace)
    write(bg.vec, paste0(bg.dir, "/bg.txt"), sep = "\n")
    bg.bed <- paste0(prep.dir, "/bg.bed")
    trash <- utilsFanc::write.zip.fanc(df = utilsFanc::loci.2.df(loci.vec = bg.vec)[, 1:3], out.file = bg.bed)
  } else {
    bg.vec <- NULL
    bg.bed <- NULL
  }
  bg.assess(mat = norm.mat, fg.vec = fg.vec, bg.vec = bg.vec,
            scatter.xy.df = scatter.xy.df, scatter.meta = scatter.meta, col.data = dds@colData,
            plot.out = paste0(bg.dir, "/bg_assess.png" ))
  bg.assess(mat = norm.mat, fg.vec = fg.vec, bg.vec = bg.vec,
            scatter.xy.df = scatter.xy.df, scatter.meta = scatter.meta, col.data = dds@colData,
            transformation = function(x) return(log2(x + 1)),
            plot.out = paste0(bg.dir, "/bg_assess_log2.png" ))

  if (write.log == T)
    log.file <- paste0(homer.dir, "/stdout.log")
  else
    log.file <- NULL
  liteRnaSeqFanc::homer.motif.enrich(fg.bed = fg.bed, genome = genome, outdir = homer.dir, bg.bed = bg.bed,
                                     size = size, thread = thread, denovo = denovo, other.params = other.params,
                                     homer.path = homer.path, find.motifs.path = find.motifs.path, run = run,
                                     stdout.file = log.file)
  return()
}

a2bl.homer.pipe <- function(a2bl, work.dir, slot, summ.slot = "summary", # behaviour change
                            rank.key, desc = F, top.n.vec = c(100, 200, 300, 500),  # foreground choice
                            updown.key = "log2FoldChange", trends = c("up", "down"),
                            quantile.interval = NULL, quantile.samples,
                            n.bg.each = 25, no.replace = F, samples.for.bg,# background choice
                            scatter.xy.df, scatter.meta = NULL,
                            threads.a2bl = 1, threads.titrate = 4,  threads.homer = 1,
                            genome, size = "given", denovo = F, other.params = "", # homer params
                            find.motifs.path = "/bar/cfan/software/homer/bin/findMotifsGenome.pl", # homer params
                            homer.path = "/bar/cfan/software/homer/bin/", homer.run = T)  {# homer params
  # slot could be either res.exp or res.shrink.exp
  trash <- utilsFanc::safelapply(names(a2bl), function(a2b.name) {
    a2b <- a2bl[[a2b.name]]
    utilsFanc::safelapply(top.n.vec, function(top.n) {
      lapply(trends, function(trend) {
        if (trend == "up")
          trend.f <- `>`
        else
          trend.f <- `<`

        fg.df <- a2b[[slot]]
        if (!is.null(summ.slot)) {
          warning("using de.genes in summ.slot. behaviour change since 2022-07-07")
          fg.df <- fg.df %>% dplyr::filter(gene %in% a2b[[summ.slot]]$de.genes)
        }
        if (!is.null(quantile.interval)) {
          filter.vec <- fg.df[, quantile.samples] %>% rowMeans()
          limits <- quantile(filter.vec, quantile.interval)
          fg.df <- fg.df[filter.vec > limits[1] & filter.vec < limits[2] , ]
        }
        fg.df <- fg.df %>% filter(trend.f(!!as.name(updown.key), 0)) %>% .[order(.[, rank.key]),]
        if (nrow(fg.df) < 1) {
          return()
        }
        if (desc == T)
          fg.df <- fg.df[nrow(fg.df):1,]
        if (top.n == 0) {
          # 0 means take all
          top.n <- nrow(fg.df)
        }
        fg.df <- fg.df[1:min(nrow(fg.df), top.n), ]

        work.dir <- paste0(work.dir, "/", a2b.name, "/", slot, "_",rank.key, "_top_", top.n,"_", trend,
                           paste0(quantile.interval, collapse = ""),"/")
        prep.dir <- paste0(work.dir, "/prep/")
        system(paste0("mkdir -p ", work.dir, " ", prep.dir))

        write.table(fg.df, paste0(prep.dir, "/fg.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
        fg.vec <- fg.df$gene
        deseq2.homer(dds = a2b$dds, fg.vec = fg.vec, n.bg.each = n.bg.each, no.replace = no.replace,
                     samples.for.bg = samples.for.bg,
                     scatter.xy.df = scatter.xy.df, scatter.meta = scatter.meta,
                     work.dir = work.dir, write.log = T, genome = genome,
                     size = size, denovo = denovo, thread = threads.homer, other.params = other.params,
                     find.motifs.path = find.motifs.path, homer.path = homer.path, run = homer.run)
        return()
      })
    }, threads = threads.titrate)
  }, threads = threads.a2bl)

  return()
}


deseq2.test.shrinkage <- function(dds, out.dir, shrink.cmd, top.n = 100) {
  res <- lfcShrink(dds, type="ashr", contrast = c("type", "tumor", "ctrl"))
  res$log2FoldChange

}

a2bl.add.shrinkage <- function(a2bl, method = "ashr", coef = NULL, contrast = NULL,
                               sample.order, threads = 1, shrink.name = NULL) {
  if (is.null(shrink.name))
    shrink.name <- method
  slot.name <- paste0("res.shrink.", shrink.name)
  slot.name.exp <- paste0(slot.name, ".exp")
  a2bl <- utilsFanc::safelapply(a2bl, function(s2b.obj) {
    params <- list(dds = s2b.obj$dds, type = method, coef = coef, contrast = contrast) %>%
      .[! sapply(., is.null)]
    s2b.obj[[slot.name]] <- do.call(DESeq2::lfcShrink, params)

    s2b.obj[[slot.name.exp]] <- deseq2.add.value(dds=s2b.obj$dds, res=s2b.obj[[slot.name]],
                                                    sample.order = sample.order)
    rownames(s2b.obj[[slot.name.exp]]) <- NULL
    return(s2b.obj)
  }, threads = threads)
  return(a2bl)
}

plotPCA.fanc <- function(object, intgroup = "condition", ntop = 500, returnData = FALSE,
                         dims.return = 1:2) {
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  # if (is.null(intgroup) || !all(intgroup %in% names(colData(object)))) {
  #   stop("the argument 'intgroup' should specify columns of colData(dds)")
  # }
  intgroup.df <- as.data.frame(colData(object)[, intgroup,
                                               drop = FALSE])

  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  #########FANC edits
  d1 <- pca$x %>% as.data.frame()
  dims.return <- dims.return[paste0("PC", dims.return) %in% colnames(d1)]
  d1 <- d1[, paste0("PC", dims.return)]
  d <- data.frame(group = group,
                  intgroup.df, name = colnames(object)) %>%
    cbind(d1, .)

  if (returnData) {
    attr(d, "percentVar") <- percentVar[dims.return]
    return(list(d = d, genes = rownames(object)[select]))
  }
  ##################
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
    geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] *
                                                        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] *
                                                                                                            100), "% variance")) + coord_fixed()
}

up.down.table <- function(pbl, slot, query.gr.list, use.bulkNorm = F,
                          trends = c("up", "down"), up.down.key,
                          write.query.overlap = T, out.dir,
                          threads.pbl = 1, threads.query = 1,
                          genome) {
  # sample.y <- comp %>% sub(":.+$", "",.)
  # sample.x <- comp %>% sub("^.+:", "",.)
  stats <- utilsFanc::safelapply(names(query.gr.list), function(gr.name) {
    query.gr <- query.gr.list[[gr.name]]
    if (is.character(query.gr)) {
      query.gr <- utilsFanc::loci.2.gr(query.gr)
    }
    stats <- utilsFanc::safelapply(pbl, function(a2b) {
      lapply(trends, function(trend) {
        if (trend == "up")
          trend.f <- `>`
        else
          trend.f <- `<`

        if (use.bulkNorm == T)
          all.df <- a2b$bulkNorm
        else
          all.df <- a2b[[slot]]

        de.df.bg <-  a2b[[slot]] %>% filter(trend.f(!!as.name(up.down.key), 0))
        n.bg <- nrow(all.df)
        n.de.bg <- de.df.bg %>% nrow()
        pct.de.bg <- n.de.bg/n.bg

        bg.gr <- utilsFanc::loci.2.df(df = all.df, loci.col.name = "gene", return.gr = T)
        fg.gr <- bg.gr %>% subsetByOverlaps(query.gr, ignore.strand = T)
        fg.df <- fg.gr %>% as.data.frame()

        bg.gr.de <- utilsFanc::loci.2.df(df = de.df.bg, loci.col.name = "gene", return.gr = T)
        fg.gr.de <- bg.gr.de %>% subsetByOverlaps(query.gr, ignore.strand = T)
        fg.df.de <- fg.gr.de %>% as.data.frame()

        n.fg.de <- nrow(fg.df.de)
        n.fg <- nrow(fg.df)

        stat <- data.frame(query.name = gr.name, cluster = a2b$root.name, trend = trend,
                           n.fg.de = n.fg.de, n.fg = n.fg,
                           pct.de.fg = n.fg.de/n.fg, pct.de.bg = pct.de.bg,
                           p = pbinom(n.fg.de, n.fg, prob = pct.de.bg, lower.tail = F))

        if (write.query.overlap == T && n.fg.de > 0) {
          system(paste0("mkdir -p ", out.dir))
          anno <- fg.df.de %>% utilsFanc::annotate.fanc(gr = fg.df.de %>%
                                                          makeGRangesFromDataFrame(keep.extra.columns = T),
                                                        genome = genome)
          anno$gene <- NULL
          anno <- anno %>% as.data.frame()

          xlsx::write.xlsx(anno, file = paste0(out.dir, "/", gr.name, "..",a2b$root.name, "_fg_de.xlsx"),
                           col.names = T, row.names = F)
          write.table(anno, paste0(out.dir, "/", gr.name, "..", a2b$root.name, "_fg_de.tsv"),
                      col.names = T, row.names = F)
        }

        return(stat)
      }) %>% Reduce(rbind, .)
    }, threads = threads.pbl) %>% Reduce(rbind, .)
    return(stats)
  }, threads = threads.query) %>% Reduce(rbind, .)

   utilsFanc::write.zip.fanc(stats, out.file = paste0(out.dir, "/stats.tsv"), row.names = F,
                             col.names = T, zip = F)
   return(stats)
}

comp.vec.gen <- function(group.var, comp.var, return.df = F) {
  if (length(comp.var) != 2) {
    stop("length(comp.var) != 2")
  }
  res <- lapply(group.var, function(group) {
    return(paste0(paste0(group, "..", comp.var), collapse = ":"))
  }) %>% unlist()
  if (return.df == T) {
    res <- data.frame(x = res %>% sub(".+:", "", .),
                      y = res %>% sub(":.+", "", .))
  }
  return(res)
}

deseq2.hm <- function(pbl, use.bubble = F,
                      genes = NULL, summ.slot = "summary",
                      plot.dir, root.name = NULL, threads = 1,
                      dense.hm = T, dense.cutoff = 30,
                      width = NULL, height = NULL, res = 100,
                      ...) {
  if (is.null(root.name)) {
    root.name <- basename(plot.dir)
  }
  utilsFanc::safelapply(pbl, function(s2b) {
    plot.out <- paste0(plot.dir, "/", root.name, "_", s2b$root.name, ".pdf")
    # update 2023-03-17: ("_q", s2b$summary$padj.cutoff) used to be in the file names as well

    hm.core.2(s2b = s2b, use.bubble = use.bubble,
              genes = genes, summ.slot = summ.slot,
              add.basemean = !use.bubble,
              dense.hm = dense.hm, dense.cutoff = dense.cutoff,
              plot.out = plot.out,
              width = width, height = height, res = res, ...)
    return(NULL)
  }, threads = threads)
  return()
}


deseq2.volcano.plot <- function(pbl, comp, label.top.n = 15,
                                res.slot = "res", summary.slot = "summary",
                                out.file, xintercept = c(-1, 1)) {
  pl <- lapply(pbl, function(s2b) {
    plot.df <- s2b[[res.slot]] %>% as.data.frame() %>%  mutate(., gene = rownames(.)) %>%
      mutate(mlog10pvalue = -log10(pvalue))
    genes.label <- c(s2b[[summary.slot]]$up.genes[1:label.top.n],
                     s2b[[summary.slot]]$down.genes[1:label.top.n]) %>% .[!is.na(.)]
    p <- xy.plot(df = plot.df, x = "log2FoldChange", y = "mlog10pvalue",
            add.abline = F, pt.color = "grey50", highlight.var = "gene",
            highlight.values = s2b[[summary.slot]]$de.genes, highlight.ptsize = 0.2,
            label.var = "gene", label.values = genes.label, use.repel = T) +
      geom_vline(xintercept = xintercept, linetype = "dashed", color = "blue", alpha = 0.5) +
      ggtitle(s2b$root.name) +
      theme(aspect.ratio = 1) +
      theme(text = element_text(size = 15))

    type.y <- comp %>% sub(":.+", "", .)
    type.x <- comp %>% sub(".+:", "", .)

    grob <- grobTree(textGrob(paste0(type.x, " > ",type.y,", n=", s2b[[summary.slot]]$n.down),
                              x=0.05,  y=0.05, hjust=0,
                              gp=gpar(col="blue", fontsize=13)))
    p <- p + annotation_custom(grob)

    grob <- grobTree(textGrob(paste0(type.y, " > ",type.x,", n=", s2b[[summary.slot]]$n.up),
                              x=0.7,  y=0.05, hjust=0,
                              gp=gpar(col="blue", fontsize=13)))
    p <- p + annotation_custom(grob)

    return(p)

  })

  trash <- wrap.plots.fanc(plot.list = pl, plot.out = out.file, sub.height = 8,
                  sub.width = 8)
  return()
}

a2bl.distro <- function(a2bl, summary.slot = "summary",
                        genome,
                        mode = "pro_vs_enh", promoter.boundaries = c(-1000, 500),
                        plot.out = NULL) {
  if (any(sapply(a2bl, function(x) is.null(x[[summary.slot]])))) {
    stop("summary.slot not found in at least 1 of the a2bl elements")
  }
  if (genome == "mm10") {
    genome.name <- "BSgenome.Mmusculus.UCSC.mm10"
    annoDb <- "org.Mm.eg.db"
    TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  } else {
    stop("only mm10 has been developed!")
  }
  a2bl <- lapply(a2bl, function(a2b) {
    if (mode == "pro_vs_enh") {
      gr <- rownames(a2b$dds) %>% utilsFanc::loci.2.gr()
      gr <- ChIPseeker::annotatePeak(gr, tssRegion = promoter.boundaries, TxDb = TxDb, annoDb = annoDb)@anno
      gr$proEnh <- NA
      gr$proEnh[gr$annotation == "Promoter"] <- "promoter"
      gr$proEnh[gr$annotation != "Promoter"] <- "enhancer"

      # it's definitely not a gene. But using the "gene" nomenclature to make left join easier
      df <- data.frame(gene = utilsFanc::gr.get.loci(gr), proEnh = gr$proEnh)

      summ <- a2b[[summary.slot]]
      de.df <- lapply(c("up", "down"), function(type) {
        sub.df <- data.frame(gene = summ[[paste0(type, ".genes")]], deStat = type)
        return(sub.df)
      }) %>% Reduce(rbind, .)

      df <- left_join(df, de.df)
      df$deStat[is.na(df$deStat)] <- "non_de"
      stats <- df %>% group_by(proEnh, deStat) %>%
        summarise(n = n()) %>% ungroup() %>% as.data.frame()
      summ[[mode]] <- stats
      a2b[[summary.slot]] <- summ
      return(a2b)
    } else if (mode == "location") {
      gr <- rownames(a2b$dds) %>% utilsFanc::loci.2.gr()
      gr <- ChIPseeker::annotatePeak(gr, tssRegion = promoter.boundaries, TxDb = TxDb, annoDb = annoDb)@anno
      gr$annotation <- gr$annotation %>% sub("Intron.+", "Intron", .)
      gr$annotation <- gr$annotation %>% sub("Exon.+", "Exon", .)
      df <- data.frame(gene = utilsFanc::gr.get.loci(gr), annotation = gr$annotation)

      summ <- a2b[[summary.slot]]
      de.df <- lapply(c("up", "down"), function(type) {
        sub.df <- data.frame(gene = summ[[paste0(type, ".genes")]], deStat = type)
        return(sub.df)
      }) %>% Reduce(rbind, .)

      df <- left_join(df, de.df)
      df$deStat[is.na(df$deStat)] <- "non_de"
      stats <- df %>% group_by(annotation, deStat) %>%
        summarise(n = n()) %>% ungroup() %>% as.data.frame()
      summ[[mode]] <- stats
      a2b[[summary.slot]] <- summ
      return(a2b)
    } else {
      stop("only pro_vs_enh and location has been developed")
    }
  })
  if (!is.null(plot.out)) {
    pl <- lapply(a2bl, function(a2b) {
      stats <- a2b[[summary.slot]][[mode]]
      if (mode == "pro_vs_enh") {
        pl <- lapply(c("proEnh", "deStat"), function(x) {
          group <- c("proEnh", "deStat") %>% .[!. %in% x]
          stats <- stats %>% group_by(!!as.name(x)) %>%
            dplyr::mutate(pct = round(n/sum(n), digits = 3)) %>%
            as.data.frame()
          pl <- lapply(c("n", "pct"), function(type) {
            p <- ggplot(stats, aes_string(x = x, y = type, fill = group)) +
              geom_bar(stat = "identity") +
              theme_bw() +
              theme(aspect.ratio = 1) +
              theme(text = element_text(size = 15)) +
              ggtitle(a2b$root.name)
            return(p)
          })
          return(pl)
        }) %>% Reduce(c, .)
        return(pl)
      } else if (mode == "pro_vs_enh") {
        stop("only pro_vs_enh has been developed for the plotting part")
      } else {
        stop("only pro_vs_enh has been developed")
      }

    }) %>% Reduce(rbind, .)
    trash <- wrap.plots.fanc(plot.list = pl, n.split = 4, plot.out = plot.out)
  }
  return(a2bl)
}

a2bl.distro.titrate.log2FC <- function(a2bl, df = NULL, mode = "pro_vs_enh",
                                       log2FC.vec = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2),
                                       de.type = "down",
                                       text.size = 20, cluster.levels = NULL,
                                       width = 8, height = 4,
                                       plot.out) {
  root <- tools::file_path_sans_ext(plot.out)
  plot.dir <- dirname(root)
  dir.create(plot.dir, recursive = T, showWarnings = F)

  if (is.null(df)) {

    if (mode == "pro_vs_enh") {
      group.col <- "proEnh"
    } else if (mode == "location") {
      group.col <- "annotation"
    }

    df <- lapply(log2FC.vec, function(log2fc) {
      print(paste0("Processing log2fc: ", log2fc))
      a2bl <- deseq2.summary(pbl = a2bl, log2fc.cutoff = log2fc)
      a2bl <- a2bl.distro(a2bl = a2bl, genome = "mm10", mode = mode)
      saveRDS(a2bl, paste0(root, "_a2bl_log2fc", log2fc, ".Rds"))

      stats <- lapply(a2bl, function(a2b) {
        stats <- a2b$summary[[mode]] %>%
          dplyr::group_by(!!as.name(group.col)) %>%
          dplyr::mutate(log2fc = log2fc, cluster = a2b$root.name, pct = n/sum(n)) %>%
          dplyr::ungroup() %>% as.data.frame()
        return(stats)
      }) %>% Reduce(rbind, .)
      return(stats)
    }) %>% Reduce(rbind, .)

    rds <- paste0(root, "_stats.Rds")
    saveRDS(df, rds)
  } else {
    if (is.character(df)) {
      df <- readRDS(df)
    }
  }
  df <- df %>% dplyr::filter(deStat %in% de.type, log2fc %in% log2FC.vec)
  pl <- df %>% split(f = df$deStat) %>%
    lapply(function(df) {
      if (mode == "pro_vs_enh") {
        df$lineGroup <- paste0(df$cluster, "_",  df$proEnh)
        if (!is.null(cluster.levels)) {
          if (grepl("seurat", cluster.levels[1])) {
            df$cluster <- sub("seurat_clusters_", "sc", df$cluster)
            cluster.levels <- sub("seurat_clusters_", "sc", cluster.levels)
          }
          df$cluster <- factor(df$cluster, levels = cluster.levels)
        }
        p <- ggplot(df, aes(x = log2fc, y = pct, group = lineGroup, color = cluster)) +
          geom_point(shape = 17) +
          geom_line(aes(linetype = proEnh)) +
          theme_classic() +
          ggtitle(df$deStat[1]) +
          scale_y_continuous(labels = scales::percent) +
          theme(text = element_text(size = text.size))
        # theme(aspect.ratio = 1)
      } else if (mode == "location") {
        if (!is.null(cluster.levels)) {
          if (grepl("seurat", cluster.levels[1])) {
            df$cluster <- sub("seurat_clusters_", "sc", df$cluster)
            cluster.levels <- sub("seurat_clusters_", "sc", cluster.levels)
          }
          df$cluster <- factor(df$cluster, levels = rev(cluster.levels))
        }
        df$log2fc <- factor(df$log2fc, levels = unique(sort(df$log2fc)))

        # df$annotation[grepl("Downstream", df$annotation)] <- "Downstream (<3kb)"
        # anno.levels <- c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream (<3kb)", "Distal Intergenic")
        anno.levels <- c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR",
                         "Downstream (1-2kb)", "Downstream (<1kb)", "Downstream (2-3kb)",
                         "Distal Intergenic")
        anno.levels <- anno.levels %>% .[. %in% df$annotation]

        utilsFanc::check.intersect(df$annotation, "df$annotation", anno.levels, "anno.levels")

        df$Annotation <- factor(df$annotation, levels = anno.levels)

        p <- ggplot(df, aes(x = Annotation, y = pct)) +
          facet_grid(cluster ~ log2fc , switch = c("both")) +
          geom_bar(stat = "identity", aes(fill = Annotation)) +
          scale_fill_brewer(palette = "Set1") +
          scale_y_continuous(labels = scales::percent)
        # utilsFanc::gg_color_hue()
        p <- utilsFanc::theme.fc.1(p)
        p <- p + theme(legend.position = "right", plot.margin = unit(c(0, 0.12, 0, 0), "in"),
                       strip.placement = "outside", axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(), strip.background = element_blank())
        ggsave(paste0(root, "_", df$deStat[1], ".pdf"), p, device = cairo_pdf, height = 2, width = 5)
      }

      return(p)
    })
  wrap.plots.fanc(plot.list = pl, sub.width = width, sub.height = height, plot.out = plot.out)

  return(df)
}

a2bl.distro.check <- function(da, slot = "summary", tss.file, ext.1side = 0,
                              out.dir, bedtools = "/bar/cfan/anaconda2/envs/jupyter/bin/bedtools") {
  print("written as a sanity check for a2bl.distro()")
  sub.f <- function(a, b, ext.b, bV, out) {
    if (bV) {
      v.flag <- "-v"
    } else {
      v.flag <- ""
    }
    cmd <- paste0("cat ", b,
                  " | awk -F '\\t' 'BEGIN {OFS = FS} {$2 = $2 - ", ext.b, "; $3 = $3 + ", ext.b, "; print $0}'", " | ",
                  bedtools, " intersect ",v.flag, " -a ", a, " -b stdin -wa -wb ",
                  " > ", out)
    print(cmd); system(cmd)
  }
  stats <- lapply(da, function(a2b) {
    name <- a2b$root.name
    peak.file <- paste0(out.dir, "/", name, "_peaks.bed")
    utilsFanc::write.zip.fanc(utilsFanc::loci.2.gr(a2b$bulkNorm$gene),
                              out.file = peak.file, bed.shift = T)
    types <- c("pro", "enh")
    files <- lapply(types, function(type) {
      file <- paste0(out.dir, "/", name, "_", type, "_ext", ext.1side, ".bed")
      sub.f(a = peak.file, b = tss.file, ext.b = ext.1side,
            bV = (type != "pro"), out = file)
      return(file)
    })
    names(files) <- types
    peaks <- lapply(files, function(file) {
      peaks <- read.table(file)
      peaks <- paste0(peaks$V1, ":", peaks$V2+1, "-", peaks$V3)
      return(peaks)
    })
    stats <- lapply(types, function(type) {
      stats <- data.frame(n = length(peaks[[type]]),
                          n.down = length(peaks[[type]][peaks[[type]] %in% a2b[[slot]]$down.genes]))
      stats$pct.down <- stats$n.down/stats$n
      colnames(stats) <- paste0(colnames(stats), ".", type)
      return(stats)
    }) %>% Reduce(cbind, .)
    stats <- cbind(data.frame(cluster = name), stats)
    return(stats)
  }) %>% do.call(rbind,.)
  write.table(stats, file = paste0(out.dir,"/ext_", ext.1side, "_stats.tsv"),
              quote = F, row.names = F, col.names = T, sep = "\t")
  return(stats)
}

de.non.de <- function(pbl, out.dir, fc.th = 1.1, slot = "res.exp",
                      n.rand = 1000, n.top = 1000,
                      root.name = NULL, rank.by = "baseMean") {
  root.name <- root.name %||% basename(out.dir)
  dir.create(out.dir, showWarnings = F, recursive = T)
  df <- lapply(pbl, function(s2b) {
    log2fc.th <- log2(c(1/fc.th, fc.th))
    df <- s2b[[slot]] %>%
      dplyr::filter(log2FoldChange > log2fc.th[1], log2FoldChange < log2fc.th[2]) %>%
      dplyr::arrange(desc(!!as.name(rank.by))) %>%
      dplyr::select(gene, !!as.name(rank.by))
    out.file <- paste0(out.dir, "/", root.name, "_", s2b$root.name,
                       "_", fc.th, "_by_", rank.by, ".tsv")
    if (!grepl(":", df$gene[1])) {
      write.table(df, out.file, sep = "\t", col.names = T, row.names = F, quote = F)
    } else {
      out.file <- sub(".tsv$", ".bed", out.file)
      df <- utilsFanc::loci.2.df(df, loci.col.name = "gene", remove.loci.col = T)
      df$start <- df$start - 1
      write.table(df, out.file, sep = "\t", col.names = F, row.names = F, quote = F)
    }
    df$cluster <- s2b$root.name
    return(df)
  }) %>% do.call(rbind, .)

  if ("gene" %in% colnames(df)) {
    suffix <- ".tsv"
    bed.flag <- F
  } else {
    suffix <- ".bed"
    bed.flag <- T
  }
  df.rand <- df %>% split(., f = .$cluster) %>%
    lapply(function(x) return(x[sample(1:nrow(x), size = min(n.rand, nrow(x)), replace = F),])) %>%
    do.call(rbind, .)
  out.file <- paste0(out.dir, "/", root.name, "_", "all",
                     "_", fc.th, "_by_", rank.by, "_rand",n.rand, suffix)
  write.table(df.rand, out.file, sep = "\t", col.names = !bed.flag, row.names = F, quote = F)

  df.top <- df %>% split(., f = .$cluster) %>%
    lapply(function(x) return(x[1:min(n.top, nrow(x)),])) %>%
    do.call(rbind, .)
  out.file <- paste0(out.dir, "/", root.name, "_", "all",
                     "_", fc.th, "_by_", rank.by, "_top",n.top, suffix)
  write.table(df.top, out.file, sep = "\t", col.names = !bed.flag, row.names = F, quote = F)
  return()
}

da.pileup.DAR <- function(summ.dir, clusters,
                          bws = NULL, bw.dir, use.all.bws = F,
                          n.peaks = 1000, seed = 42,
                          threads.clusters = 1, threads.each = 6,
                          out.dir) {
  bed.dir <- paste0(out.dir, "/bed/")
  peaks <- lapply(clusters, function(cluster) {
    lapply(c("up", "down"), function(type) {
      peaks <- readLines(paste0(summ.dir, "/", cluster, "_", type, ".genes.txt"))
      set.seed(seed = seed)
      rand <- peaks %>% .[sort(sample(1:length(.), min(n.peaks, length(.)), replace = F))]
      gr <- utilsFanc::loci.2.gr(rand)
      out.file <- paste0(bed.dir, "/", cluster, "_", type, "_rand", n.peaks, "_seed", seed, ".bed")
      utilsFanc::write.zip.fanc(df = gr, bed.shift = T, zip = T, out.file = out.file)
      return(out.file)
    }) %>% unlist() %>% return()
  })

  names(peaks) <- clusters

  if (!is.null(bws)) {
    bws <- lapply(clusters, function(x) return(bws))
  } else {
    bws <- lapply(clusters, function(cluster) {
      if (use.all.bws)
        glob <- paste0(bw.dir, "/*.bw")
      else
        glob <- paste0(bw.dir, "/", cluster, "*")

      bws <- Sys.glob(glob)
      if (length(bws) < 1) {
        stop(paste0("glob failed for: ", glob))
      }
      return(bws)
    })
  }

  names(bws) <- clusters

  utilsFanc::safelapply(clusters, function(cluster) {
    deeptools.refpoint(bw.vec = bws[[cluster]], regions.vec = peaks[[cluster]],
                       threads = threads.each, out.dir = out.dir, root.name = cluster)
  }, threads = threads.clusters)
  return()
}

### wrote for Jun
da.cmp <- function(query, ref, query.name, ref.name,
                   clusters = NULL, slot = "res.exp",
                   col = "log2FoldChange", abs.diff = F,
                   other.cols = NULL, bSimple = F,
                   query.summ.slot,
                   add.trend = F, color.density = F,
                   threads = 1, out.dir, root = NULL,
                   bScatter = F, use.plotly = F) {
  # query and ref should both be pbl.
  if (is.null(root)) {
    root <- basename(out.dir)
  }
  clusters.avail <- intersect(names(query), names(ref))
  if (is.null(clusters))
    clusters <- clusters.avail
  clusters <- clusters %>% .[. %in% clusters.avail]
  if (length(clusters) < 1) {
    stop("length(clusters) < 1")
  }

  utilsFanc::safelapply(clusters, function(cluster) {
    id.df <- data.frame(gene.x = rownames(query[[cluster]]$bulk.mat))
    id.df$id <- 1:nrow(id.df)
    gr.q <- query[[cluster]][[slot]][, c("gene", col, other.cols)] %>%
      utilsFanc::loci.2.df(loci.col.name = "gene", remove.loci.col = F, return.gr = T)
    summ <- query[[cluster]][[query.summ.slot]]
    if (! query.summ.slot %in% names(query[[cluster]])) {
      stop("! query.summ.slot %in% names(query[[cluster]])")
    }
    gr.q$type <- NA
    for (type in c("up", "down")) {
      gr.q$type[gr.q$gene %in% summ[[paste0(type, ".genes")]]] <- type
    }
    gr.q <- gr.q[!is.na(gr.q$type)]
    if (length(gr.q) < 1) {
      stop("length(gr.q) < 1")
    }
    gr.r <- ref[[cluster]][[slot]][, c("gene", col, other.cols)] %>%
      utilsFanc::loci.2.df(loci.col.name = "gene", remove.loci.col = F, return.gr = T)

    j <- plyranges::join_overlap_inner(gr.q, gr.r)
    df <- mcols(j) %>% as.data.frame()
    col.q <- paste0(col, ".x")
    col.r <- paste0(col, ".y")
    df$fc <- df[[col.q]]/df[[col.r]]
    df$diff <- df[[col.q]] - df[[col.r]]
    if (abs.diff) {
      df$diff <- abs(df$diff)
    }
    dir.create(out.dir, showWarnings = F, recursive = T)
    lapply(c("up", "down"), function(type) {
      df <- df[df$type == type, ]
      df <- left_join(df, id.df)
      df <- utilsFanc::df.rearrange.cols(df = df ,cols = "id", pos = 1)
      out.file <- paste0(out.dir, "/", root, "_", cluster, "_", query.name, "_vs_", ref.name,"_", slot,
                         "_", col, "_", query.summ.slot, "_", type, ".tsv")
      write.table(df, out.file, quote = F, sep = "\t", col.names = T, row.names = F)
      if (bScatter) {
        out.plot <- sub(".tsv", if_else(use.plotly, ".html" ,".png"), out.file)
        pl <- lapply(c(col, other.cols), function(v) {
          colx <- paste0(v, ".x")
          coly <- paste0(v, ".y")
          if (bSimple) {
            to.plot <- data.frame(x = colx, y = coly)
          } else {
            to.plot <- data.frame(x = c(colx, colx, coly),
                                  y = c(coly, "diff", "diff"))
          }

          pl <- lapply(1:nrow(to.plot), function(i) {
            x <- to.plot[i, "x"]
            y <- to.plot[i, "y"]
            type.x <- sub(".[xy]$", "", x)
            type.y <- sub(".[xy]$", "", y)
            y.lim <- x.lim <- NULL
            transformation.x <- transformation.y <- NULL
            if (type.x == "baseMean") {
              transformation.x <- log2
            } else if (type.x == "pvalue") {
              transformation.x <- function(x) return(-log10(x))
              x.lim <- c(0, 25)
            } else if (grepl("\\.\\.mean", type.x)){
              transformation.x <- log2
            } else if (!type.x %in% c("log2FoldChange", "diff")) {
              stop(paste0("column ", type.x, " not recognized"))
            }

            if (type.y == "baseMean") {
              transformation.y <- log2
            } else if (type.y == "pvalue") {
              transformation.y <- function(y) return(-log10(y))
              y.lim <- c(0, 25)
            } else if (grepl("\\.\\.mean", type.y)){
              transformation.y <- log2
            } else if (!type.y %in% c("log2FoldChange", "diff")) {
              stop(paste0("column ", type.y, " not recognized"))
            }
            df <- df[!is.na(df[, x]) & !is.na(df[, y]),]
            p1 <- xy.plot(df, x = x, y = y, plotly.var = "id", plotly.label.all = T,
                         transformation.x = transformation.x, transformation.y = transformation.y,
                         color.density = color.density, add.abline = F)
            if (!is.null(x.lim)) p1 <- p1 + xlim(x.lim)
            if (!is.null(y.lim)) p1 <- p1 + ylim(y.lim)
            if (add.trend)
              p1 <- p1 + geom_smooth(show.legend = F)

            pl <- list(p1)
            if (!use.plotly) {
              p2 <- grid.sampler(df, x = x, y = y,
                                 transf.x = transformation.x, transf.y = transformation.y,
                                 range.x = x.lim, range.y = y.lim, color.density = F, seed = 42,
                                 out.dir = paste0(out.dir, "/grid_sampler/"), return.plot = T,
                                 print.plot = F, add.abline = F,
                                 root.name = paste0(root, "_", cluster, "_", x, "_vs_", y))
              pl <- c(pl, list(p2))
            }

            return(pl)
          }) %>% do.call(c, .)
          return(pl)
        }) %>% do.call(c, .)
        wrap.plots.fanc(plot.list = pl, plot.out = out.plot, n.col = min(6, length(pl)),
                        sub.width = 6, sub.height = 6, tooltip = "id")
      }
    })
    return()
  }, threads = threads)
}

da.cmp.plot <- function(tsvs, fill.var, use.trend = F, brewer.palette = "Greens",
                        plot.out, text.size = 20, width = 6, height = 4) {
  if (is.null(names(tsvs))) {
    stop("is.null(names(tsvs))")
  }
  df <- lapply(names(tsvs), function(name) {
    df <- read.table(tsvs[[name]], header = T)
    df <- df %>% dplyr::mutate(direction = if_else(!!as.name(fill.var) > 0, "up", "down"),
                               name = name) %>%
      dplyr::select(name, direction) %>%
      dplyr::group_by(name, direction) %>%
      dplyr::summarise(n = n()) %>% as.data.frame()
    if (use.trend) {
      if (grepl("up.tsv$", tsvs[[name]])) {
        trend <- "up"
      } else if (grepl("down.tsv$", tsvs[[name]])) {
        trend <- "down"
      } else {
        stop("couldn't find grepping patterns (up.tsv$ or down.tsv$)")
      }
      df$old.d <- df$direction
      df$direction[df$old.d == trend] <- "same"
      df$direction[df$old.d != trend] <- "oppo"
    }
    return(df)
  }) %>% do.call(rbind, .)
  df <- df %>% dplyr::group_by(name) %>%
    dplyr::mutate(frac = n/sum(n)) %>% as.data.frame()
  df$name <- factor(df$name, levels = rev(unique(df$name)))
  p <- ggplot(df, aes(x = name, y = frac, fill = direction)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_brewer(palette = brewer.palette, direction = -1, type = "div") +
    theme_classic() +
    theme(text = element_text(size = text.size)) +
    coord_flip()
  wrap.plots.fanc(plot.list = list(p), sub.width = width, sub.height = height, plot.out = plot.out)
  return(p)
}

## wrote for comparing the differential results from RNA or ATAC based clustering
deseq2.cmp <- function(pbl.list, cmp.df, plot.out, summ.slot = "summary") {
  # pbl.list <- list(x = de.seurat_clusters, y = de.clusters)
  # cmp.df format: data.frame(x = c(paste0("seurat_clusters_", 6:7),
  #                           y = c(paste0("Clusters_C", 1:2)))
  # note the names of pbl.list needs to match the names of cmp.df
  utilsFanc::check.intersect(colnames(cmp.df), "colnames(cmp.df)",
                             names(pbl.list), "names(pbl.list)")
  lapply(colnames(cmp.df), function(pbl.name) {
    pbl <- pbl.list[[pbl.name]]
    utilsFanc::check.intersect(cmp.df[,pbl.name], paste0("cmp.df[,", pbl.name, "]"),
                               names(pbl), paste0("names(pbl.list[",pbl.name,"])"))
  })
  x.list.list <- cmp.df %>% split(., f = 1:nrow(.)) %>%
    lapply(function(cmp) {
      cmp.name <- unlist(cmp) %>% paste0(collapse = ":")
      x.list.list <- lapply(c("up", "down"), function(direction) {
        x.list <- lapply(colnames(cmp), function(pbl.name) {
          pbl <- pbl.list[[pbl.name]]
          s2b <- pbl[[cmp[, pbl.name]]]
          if (is.null(s2b[[summ.slot]])) {
            stop("is.null(s2b[[summ.slot]])")
          }
          genes <- s2b[[summ.slot]][[paste0(direction,".genes")]]
          return(genes)
        })
        names(x.list) <- colnames(cmp)
        return(x.list)
      })
      names(x.list.list) <- paste0(cmp.name, "_", c("up", "down"))
      return(x.list.list)
    }) %>% Reduce(c, .)
  multi.venn(x.list.list = x.list.list,
             plot.titles.vec = names(x.list.list), use.upset = T,
             height.sub = 800, width.sub = 800,
             out.file = plot.out)
  return(x.list.list)
}

de.summary.cmp <- function(des, slot = "summary", res.slot = "res.exp", clusters = NULL,
                           trend = "de",
                           out.dir, root.name = NULL,
                           bXy = F, comp.vec, comp.meta = NULL,
                           bHm = F) {
  # des <- list(de1 = de1, de2 = de2)
  if (length(des) < 2) {
    stop("length(des) < 2")
  } else {
    des <- des[1:2]
  }
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  if (is.null(names(des))) {
    stop("is.null(names(des))")
  }
  all.clusters <- intersect(names(des[[1]]), names(des[[2]]))
  if (is.null(clusters)) {
    clusters <- all.clusters
  }
  print(all.clusters)
  clusters <- clusters %>% .[. %in% all.clusters]
  if (length(clusters) < 1) {
    stop("length(clusters) < 1")
  }
  x <- names(des)[1]
  y <- names(des)[2]
  name <- paste0(root.name, "_",x, "_vs_", y, "_", trend, ".genes_")
  n.df <- lapply(clusters, function(cluster) {
    name <- paste0(name, "_", cluster)
    des <- lapply(des, function(de) return(de[cluster]))
    degs <- lapply(des, function(de) {
      s2b <- de[[cluster]]
      de.genes <- s2b[[slot]][[paste0(trend, ".genes")]]
      if (is.null(de.genes)) {
        stop(paste0("slot ", slot, " not found in at least one of the de's for cluster ", cluster))
      }
      return(de.genes)
    })
    summ <- utilsFanc::intersect.summary(x = degs[[1]], x.name = names(degs)[1],
                                         y = degs[[2]], y.name = names(degs)[2],
                                         ret.elements = T)
    degs.all <- unlist(degs) %>% unique()

    ress <- lapply(names(des), function(de.name) {
      de <- des[[de.name]]
      res <- de[[cluster]][[res.slot]]
      if (is.null(res)) {
        stop("de[[cluster]][[res.slot]]")
      }
      if (!is.data.frame(res)) {
        stop("!is.data.frame(res)")
      }
      res <- res %>% dplyr::filter(gene %in% degs.all)
      res$pass <- F
      res$pass[res$gene %in% degs[[de.name]]] <- T
      colnames(res) <- paste0(colnames(res), "_", de.name)
      colnames(res)[colnames(res) == paste0("gene_", de.name)] <- "gene"
      return(res)
    })

    j <- dplyr::full_join(ress[[1]], ress[[2]])
    for (i in c("lfcSE", "log2FoldChange", "pvalue", "padj")) {
      if (any(grepl(i, colnames(j)))) {
        j[[paste0(i, "_diff")]] <- abs(j[[paste0(i, "_", y)]]) - abs(j[[paste0(i, "_", x)]])
      }
    }
    rm(i)
    j <- j[, sort(colnames(j))]
    tsv <- paste0(out.dir, "/", name, ".tsv")
    dir.create(out.dir, showWarnings = F, recursive = T)
    write.table(j, tsv, sep = "\t", col.names = T, row.names = F, quote = F)

    reasons <- c("Filter", "Log2fc", "Fdr")
    summ.reason <- lapply(reasons, function(type) {
      res <- lapply(c(x, y), function(this) {
        that <- ifelse(this == x, y, x)
        only <- summ[[paste0(this, ".only")]]
        that.s2b <- des[[that]][[cluster]]
        res.that <- that.s2b[[res.slot]]
        if (!is.data.frame(res.that)) {
          stop("!is.data.frame(res.that)")
        }
        res.that <- res.that %>% dplyr::filter(gene %in% only)
        if (type == "Filter") {
          res <- list(only %>% .[!. %in% res.that$gene])
        } else if (type == "Log2fc") {
          res <- list(res.that$gene[res.that$log2FoldChange <= that.s2b[[slot]]$log2fc.cutoff])
        } else if (type == "Fdr") {
          res <- list(res.that$gene[res.that$padj >= that.s2b[[slot]]$padj.cutoff])
        }

        names(res) <- paste0(this, ".only.", type)
        return(res)
      }) %>% Reduce(c, .)
      return(res)
    }) %>% Reduce(c, .)

    summ <- c(summ, summ.reason)
    n.df <- lapply(summ, length) %>% as.data.frame()
    n.df <- cbind(data.frame(cluster = cluster), n.df)

    lapply(names(summ[1:3]), function(gene.type) {
      name <- paste0(name, "_", gene.type)
      genes <- summ[[gene.type]]
      lapply(names(des), function(id) {
        de <- des[[id]]
        if (bXy) {
          lapply(c(T, F), function(bLabel) {
            label.list <- list(genes)
            names(label.list) <- cluster
            deseq2.xyplot(
              pbl = de, comp.vec = comp.vec, comp.meta = comp.meta,
              transformation = function(x) {
                return(log2(x + 1))
              }, highlight.genes = genes,
              add.label = bLabel, label.list = label.list,
              hjust = 0.5, vjust = 0.5, use.repel = T, italic.label = T,
              text.size = 1.5, highlight.ptsize = 0.1,
              plot.dir = out.dir,
              root.name = paste0(name, "_xy_on_", id, "_", bLabel),
              device = "png", theme.text.size = 15,
              plot.each.comp = F
            )
          })
        }
        if (bHm) {
          sample.order <- de[[cluster]]$coldata %>% rownames()
          deseq2.hm(de,
                    genes = genes, use.bubble = T,
                    plot.dir = out.dir,
                    root.name = paste0(name, "_hm"),
                    dense.hm = F, sample.order = sample.order
          )
        }
      })
      return()
    })
    return(n.df)
  }) %>% do.call(rbind, .)
  dir.create(out.dir, showWarnings = F, recursive = T)
  write.table(n.df, paste0(name, "_summary.tsv"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  return(n.df)
}
##

de.add.mean <- function(pbl, slot = "res.exp", meta) {
  pbl <- lapply(pbl, function(s2b) {
    if (! meta %in% colnames(s2b$coldata))
      stop("! meta %in% colnames(a2b$coldata)")
    df.add <- s2b$coldata %>% split(., f = .[, meta]) %>%
      lapply(function(df) {
        samples <- rownames(df)
        res <- s2b[[slot]][, samples] %>% apply(1, mean)
        return(res)
      }) %>% as.data.frame()
    colnames(df.add) <- paste0(colnames(df.add), "..mean")
    s2b[[slot]] <- utilsFanc::add.column.fanc(s2b[[slot]], df.add)
    return(s2b)
  })
  return(pbl)
}


de.get.joint.mat <- function(de.list) {
  if (is.null(names(de.list))) {
    stop("is.null(names(de.list))")
  }
  de.list <- lapply(names(de.list), function(name) {
    de <- de.list[[name]]
    de <- lapply(de, function(s2b) {
      s2b$root.name <- paste0(s2b$root.name, "_", name)
      return(s2b)
    })
    return(de)
  })

  s2bs <- do.call(c, de.list)
  genes <- lapply(s2bs, function(s2b) return(rownames(s2b$bulk.mat))) %>%
    unlist() %>% unique() %>% gtools::mixedsort()

  mat.j <- lapply(s2bs, function(s2b) {
    mat <- s2b$bulk.mat
    colnames(mat) <- paste0(s2b$root.name, "_", colnames(mat))
    genes.not.in <- genes %>% .[!. %in% rownames(mat)]
    mat.add <- matrix(NA, nrow = length(genes.not.in), ncol = ncol(mat))
    rownames(mat.add) <- genes.not.in; colnames(mat.add) <- colnames(mat)
    mat <- rbind(mat, mat.add)
    mat <- mat[genes, ]
  }) %>% do.call(cbind, .)
  return(mat.j)
}

de.transfer.slot <- function(de.from, de.to, slot.from, slot.to, strict = T) {
  if (!identical(sort(names(de.from)), sort(names(de.to)))) {
    stop("!identical(sort(names(de.from)), sort(names(de.to)))")
  }
  de.to <- lapply(de.to, function(a2b) {
    cluster <- a2b$root.name
    a2b[[slot.to]] <- de.from[[cluster]][[slot.from]]
    if (strict && "de.genes" %in% names(a2b[[slot.to]])) {
      for (i in c("de", "up", "down")) {
        slot.genes <- paste0(i, ".genes")
        slot.n <- paste0("n.", i)
        a2b[[slot.to]][[slot.genes]] <- a2b[[slot.to]][[slot.genes]] %>% .[. %in% a2b$res.exp$gene]
        a2b[[slot.to]][[slot.n]] <- length(a2b[[slot.to]][[slot.genes]])
      }
    }
    return(a2b)
  })
  return(de.to)
}

de.transfer.slot.by.overlap <- function(de.from, de.to, expand.de.to = F, slot.from, slot.to,
                                        clusters.from, clusters.to, require.overlap = T) {
  if (expand.de.to) {
    if (length(de.to) != 1) {
      stop("when expand.de.to is specified, de.to has to be of length of 1")
    }
    de.from <- de.from[unique(clusters.from)]
    de.to <- rep(de.to, length(de.from))
    names(de.to) <- names(de.from)
    for ( i in names(de.from)) {
      de.to[[i]]$root.name <- i
    }
    rm(i)
    clusters.to <- clusters.from
  }
  if (length(clusters.from) == 1) {
    clusters.from <- rep(clusters.from, length(clusters.to))
  }
  if (length(clusters.from) != length(clusters.to)) {
    stop("length(clusters.from) != length(clusters.to)")
  }

  for (i in 1:length(clusters.from)) {
    c.from <- clusters.from[i]
    c.to <- clusters.to[i]
    prep <- de.from[[c.from]][[slot.from]]
    if (is.null(prep)) {
      stop(paste0("de.from[[", c.from, "]][[", slot.from, "]] not found"))
    }
    prep$transfer.stats <- "overlap not required"

    if (require.overlap) {
      prep$transfer.stats <- list()

      filter.slot <- "res.exp"
      if (is.null(de.to[[c.to]][[filter.slot]])) {
        filter.slot <- "bulkNorm"
      }

      for (j in c("de", "up", "down")) {
        slot.genes <- paste0(j, ".genes")
        slot.genes.no <- paste0(j, ".genes.noOverlap")

        slot.n <- paste0("n.", j)
        slot.n.no <- paste0("n.", j, ".noOverlap")
        slot.pct <- paste0("pct.", j, ".Overlap")

        genes.ori <- prep[[slot.genes]]
        if (!grepl(":", de.from[[c.from]]$bulkNorm$gene[1])) {
          # RNA-seq data
          bO <- genes.ori %in% de.to[[c.to]][[filter.slot]]$gene
          prep[[slot.genes]] <- genes.ori[bO]
        } else {
          # ATAC data
          print("ATAC-seq data; determining overlap by genomic overlap")
          bO <- utilsFanc::loci.in(genes.ori, de.to[[c.to]][[filter.slot]]$gene)
          ids <- utilsFanc::loci.in(de.to[[c.to]][[filter.slot]]$gene, genes.ori, return.id.y = T)
          prep[[slot.genes]] <- de.to[[c.to]][[filter.slot]]$gene[order(ids)][1:sum(!is.na(ids))]
        }

        prep[[slot.genes.no]] <- genes.ori[!bO]
        prep[[slot.n]] <- sum(bO)
        prep$transfer.stats[[slot.n]] <- sum(bO)
        prep$transfer.stats[[slot.n.no]] <- sum(!bO)
        prep$transfer.stats[[slot.pct]] <- round(sum(bO)/length(bO), digits = 3)
      }
      prep$transfer.stats <- as.data.frame(prep$transfer.stats)
    }
    de.to[[c.to]][[slot.to]] <- prep
  }
  return(de.to)
}

da.footprint <- function(da, slot = "summary", direction = c("up", "down"),
                         ao, motif.pos = NULL, motif.mat.name,
                         cluster.ident, type.ident, group.by = NULL,
                         norms, smooth.win = 20,
                         out.dir, root.name = "footprint") {
  if (is.null(motif.pos)) {
    motif.pos <- getPositions(ArchRProj = ao, name = motif.mat.name)
  }

  motif.pos.l <- footprint.motif.dar(motifPos = motif.pos, da = da, slot = slot)
  lapply(da, function(a2b) {
    cluster <- sub(paste0("^", cluster.ident, "_"), "", a2b$root.name)
    types <- getCellColData(aoi, c(cluster.ident, type.ident)) %>% as.data.frame() %>%
      dplyr::filter(!!as.name(cluster.ident) %in% cluster) %>%
      .[, type.ident] %>% unique() %>% gtools::mixedsort()
    if (is.null(group.by)) {
      group.by <- paste0(cluster.ident, "_", type.ident)
    }
    lapply(norms, function(norm) {
      lapply(direction, function(direc) {
        se.foot <- getFootprints(ArchRProj = ao, positions = motif.pos.l[[a2b$root.name]][[direc]],
                                 groupBy = group.by, useGroups = paste0(cluster, "..", types))
        plot.name <-  paste0(root.name, "_", a2b$root.name, "_", slot, "_", direction, "_", group.by,
                             "_norm", stringr::str_to_title(norm), "_win", smooth.win)
        plotFootprints(
          seFoot = se.foot,
          ArchRProj = ao,
          normMethod = norm,
          plotName = plot.name,
          addDOC = FALSE,
          smoothWindow = smooth.win
        )
        dir.create(out.dir, showWarnings = F, recursive = T)
        system(paste0("mv ", getOutputDirectory(aoi), "/Plots/", plot.name, ".pdf ",
                      out.dir, "/"))
        return()
      })

    })
  })
}

de.write.xls <- function(de, cluster.alias = NULL,
                         summ = "summary", slot = "res.exp",
                         cols = "all",
                         out.file,
                         append = F, prefix = NULL) {
  # cluster.alias: data.frame(name = "seurat_cluster_0", alias = "Neu")
  if (!append)
    system(paste0("rm -rf ", out.file))
  dir.create(dirname(out.file), showWarnings = F, recursive = T)
  lapply(de, function(s2b) {
    if (!is.null(cluster.alias)) {
      cluster <- cluster.alias %>% dplyr::filter(name == s2b$root.name) %>%
        dplyr::pull(alias)
      if (length(cluster) != 1) {
        stop("length(cluster) != 1")
      }
    } else {
      cluster <- s2b$root.name
    }
    utilsFanc::t.stat(paste0("Processing Cluster: ", cluster))

    df <- s2b[[slot]]

    if (cols != "all")
      df <- df[, unique(c("gene", cols))]

    lapply(c("up", "down"), function(direction) {
      genes <- s2b[[summ]][[paste0(direction, ".genes")]]
      df <- df %>% dplyr::filter(gene %in% genes) %>%
        dplyr::arrange(padj)
      sheet <- paste0(cluster, "_", direction)
      if (grepl(":", df$gene[1])) {
        colnames(df)[1] <- "peak"
      }

      if (!is.null(prefix)) {
        sheet <- paste0(prefix, "_", sheet)
      }

      if (nrow(df) == 0)
        df[1, ] <- rep("", ncol(df))
      xlsx::write.xlsx(x = df, file = out.file, sheetName = sheet,
                       append = T, row.names = F, col.names = T)
      return()
    })
  })
  return()
}

de.write.xls.all.genes <- function(de, cluster.alias = NULL,
                         slot = "res.exp", padj.cutoff = NULL,
                         out.file,
                         append = F,
                         prefix = NULL) {
  # cluster.alias: data.frame(name = "seurat_cluster_0", alias = "Neu")
  if (!append)
    system(paste0("rm -rf ", out.file))

  dir.create(dirname(out.file), showWarnings = F, recursive = T)

  lapply(de, function(s2b) {
    if (!is.null(cluster.alias)) {
      cluster <- cluster.alias %>% dplyr::filter(name == s2b$root.name) %>%
        dplyr::pull(alias)
      if (length(cluster) != 1) {
        stop("length(cluster) != 1")
      }
    } else {
      cluster <- s2b$root.name
    }
    df <- s2b[[slot]]
    sheetName <- cluster
    if (!is.null(prefix)) {
      sheetName <- paste0(prefix, "_", sheetName)
    }
    xlsx::write.xlsx(x = df, file = out.file, sheetName = sheetName,
                     append = T, row.names = F, col.names = T)
  })
}


de.write.summary <- function(de, slot, out.dir) {
  bDa <- grepl(":", rownames(de[[1]]$bulk.mat)[1])
  stats <- lapply(de, function(s2b) {
    if (is.null(s2b[[slot]]))
      return()
    lapply(c("de", "up", "down"), function(type) {
      genes <- s2b[[slot]][[paste0(type, ".genes")]]
      if (length(genes) == 0) {
        return()
      }
      dir.create(out.dir, showWarnings = F, recursive = T)
      out.file <- paste0(out.dir, "/", s2b$root.name, "_", slot, "_", type,
                         ifelse(bDa, ".bed", ".txt"))
      if (bDa) {
        utilsFanc::write.zip.fanc( utilsFanc::loci.2.gr(genes), out.file, bed.shift = T)
      } else {
        write(genes, out.file)
      }
      return()
    })
    stats <- s2b[[slot]][paste0("n.", c("de", "up", "down"))]
    stats <- cbind(data.frame(root.name = s2b$root.name), stats)
    if (!is.null(s2b[[slot]]$transfer.stats) && is.data.frame(s2b[[slot]]$transfer.stats)) {
      stats <- cbind(stats, s2b[[slot]]$transfer.stats)
    }
    return(stats)
  }) %>% do.call(rbind, .)
  dir.create(out.dir, showWarnings = F, recursive = T)
  write.table(stats, paste0(out.dir, "/", "all", "_", slot, "_stats.tsv"), quote = F,
              row.names = F, col.names = T, sep = "\t")
  return()
}

de.plot.bar <- function(de, genes, Ly49.strip = T,
                        clusters = NULL,
                        group.by, shape.by = NULL,
                        pt.size = 0.2, text.size = 6,
                        width = NULL, height = 1,
                        plot.out = NULL, return.df = F) {
  if (is.null(clusters)) {
    clusters <- names(de)
  }
  de <- de[names(de) %in% clusters]
  if (length(de) < 1) {
    stop("none of the clusters in de is in the 'clusters' parameter specified")
  }
  pl <- lapply(de, function(s2b) {
    cluster <- s2b$root.name
    df <- s2b$bulkNorm %>% dplyr::filter(gene %in% genes)
    if (nrow(df) < 1) {
      stop(paste0("none of the genes were found in the matrix of cluster: ", cluster))
    }

    df <- df %>% reshape2::melt(id.vars = "gene", variable.name = "sample", value.name = "expr")
    df <- dplyr::left_join(df, s2b$coldata)
    suppressMessages(library(ggpubr))
    if (Ly49.strip) {
      df$gene <- sub(".*Ly49", "", df$gene)
    }
    if (return.df) {
      return(df)
    }
    p <- ggbarplot(data = df, x = "gene", y = "expr",
                   fill = group.by, position = position_dodge(width = 0.75),
              color = "black", alpha = 0.5, add = "mean_sd",
              add.params = list(size = 0.2), size = 0.2) +
      geom_point(aes_string(fill = group.by, shape = shape.by),
                 position = position_jitterdodge(
                   dodge.width = 0.8, jitter.width = 0.25, jitter.height = 0),
                 size = pt.size, color = "grey25") +
      theme_classic() +
      theme(text = element_text(size = text.size, family = "Arial", color = "black"),
            # legend.position = "none",
            legend.position = "bottom",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"),
            legend.box.margin = margin(t = -0.12, r = -0.1, b = 0, l = -0.1, unit = "in"),
            legend.background = element_blank(),
            legend.spacing = unit(0.02, "in"),
            legend.key.size = unit(0.05, "in"),
            legend.box = "vertical",
            legend.title = element_text(size = text.size, family = "Arial"),
            axis.title = element_blank(),
            axis.text.x = element_text(
              face = "italic", size = text.size, family = "Arial", color = "black"),
            axis.text.y = element_text(
              size = text.size, family = "Arial", color = "black"),
            axis.line.x.bottom = element_line(size = 0.2),
            axis.line.y.left = element_line(size = 0.2),
            axis.ticks = element_line(size = 0.2),
            axis.ticks.length = unit(0.02, "in"),
            plot.margin = unit(c(0, 0.04, 0.01, 0), units = "in")) +
      guides(fill = guide_legend(order=1),
             shape = guide_legend(order=2))
    # ggsave("~/test/test.png", p, width = 25, height = 4, dpi = 100)
    return(p)
  })
  if (return.df) {
    if (length(pl) == 1) {
      pl <- pl[[1]]
    }
    return(pl)
  }
  if (is.null(width)) {
    n.genes <- sapply(de, function(s2b) {
      return(length(intersect(genes, s2b$bulkNorm$gene)))
    }) %>% max()
    width <- floor(n.genes / 4)
  }
  if (length(pl) == 1) {
    p <- pl[[1]]
    dir.create(dirname(plot.out), showWarnings = F, recursive = T)
    device <- NULL
    if (grep("pdf$", plot.out)) {
      device <- cairo_pdf
    }
    ggsave(plot.out, p, height = height, width = width, dpi= 300, device = device)
  } else {
    p <- wrap.plots.fanc(plot.list = pl, sub.height = height, sub.width = width,
                         plot.out = plot.out)
  }
  invisible(p)
}

de.plot.bar.3 <- function(
  de, genes, Ly49.strip = F, cluster = 1, color.by,
  spread.bin.size = 100,
  bar.width = 0.6, dodge.width = 0.6,
  add.pval = T, pval.group.1, pval.group.2,
  plot.out, width = 4.5, height = 1.5) {
  stop("use de.plot.bar.2. This function is quite useless..")
  genes <- genes %>% .[. %in% de[[cluster]]$res.exp$gene]
  if (length(genes) < 1) {
    stop("length(genes) < 1")
  }
  bar.df <- de.plot.bar(
    de = de, genes = genes, clusters = cluster,
    return.df = T, Ly49.strip = Ly49.strip
  )
  if (Ly49.strip) {
    levels <- sub("Ly49", "", genes)
  } else {
    levels <- genes
  }
  bar.df$gene <- factor(bar.df$gene, levels = levels)
  if (add.pval) {
    q.df <- de[[cluster]]$res.exp %>%
      dplyr::filter(gene %in% genes) %>%
      dplyr::select(padj, gene) %>%
      dplyr::mutate(group.1 = pval.group.1, group.2 = pval.group.2) %>%
      dplyr::rename(p = padj)
    if (Ly49.strip) {
      q.df$gene <- q.df$gene %>% sub("Ly49", "", .)
    }
  }
  dir.create(dirname(plot.out), showWarnings = F, recursive = T)
  p <- utilsFanc::barplot.pub.3(df = bar.df, x = "gene", y = "expr",
                              color.by = color.by,
                              spread.bin.size = spread.bin.size,
                              add.pval = add.pval, pval.group.1 = pval.group.1, pval.group.2 = pval.group.2,
                              external.p.df = q.df,
                              bar.width = bar.width, dodge.width = dodge.width) %>%
  utilsFanc::theme.fc.1() +
  ggsave(plot.out, device = cairo_pdf,
         width = 4.5, height = 1.5)
  invisible(p)
}

da.motif.foldChange <- function(da, peak.motif.df, motifs.use = NULL,
                                bViolin = F, bBoxplot = F, bDot = F,
                                publication = F, publication.small = F,
                                pt.size = 0.25, alpha = 0.4,
                                slot = "res.exp", comps,
                                # motif.gr.list,
                                motif.min.score = NULL,
                                n.levels = 5, level.fun = "max", bLog.fc = T,
                                violin.ceiling = 5, violin.floor = -5,
                                motifs.label = NULL, dot.n.label = 5,
                                height = 4, width = 4.5, n.col = 1,
                                write.df = F, write.each.df = F,
                                # threads.motif = 1,
                                plot.dir, root.name) {
  # motif.gr.list: the same format as the motifPositions you will get from ArchR::getPositions()
  # warning("this function calculates fold changes de novo based on normalized counts,
  #         unless 'beta' is comp")
  if (!is.null(motif.min.score)) {
    root.name <- paste0(root.name, "_score", motif.min.score)
  }
  if (bLog.fc) {
    value <- "log2FoldChange"
  } else {
    stop("only bLog.fc = T has been developed")
  }

  if (is.character(peak.motif.df)) {
    peak.motif.df <- readRDS(peak.motif.df)
  }
  df <- lapply(da, function(a2b) {
    # establish peak-motif association
    peak.by.motif <- peak.motif.df %>% dplyr::filter(gene %in% a2b[[slot]]$gene)
    if (!is.null(motifs.use))
      peak.by.motif <- peak.by.motif %>% dplyr::filter(motif %in% motifs.use)
    if (!is.null(motif.min.score)) {
      peak.by.motif <- peak.by.motif %>% dplyr::filter(score >= motif.min.score)
    }
    # peak.by.motif <- utilsFanc::safelapply(names(motif.gr.list), function(motif.name) {
    #   motif.gr <- motif.gr.list[[motif.name]]
    #   # mcols(motif.gr) <- mcols(motif.gr)[, "score", drop = F]
    #   if (!is.null(motif.min.score)) {
    #     motif.gr <- motif.gr[motif.gr$score >= motif.min.score]
    #   }
    #   gr <- a2b[[slot]]$gene %>% utilsFanc::loci.2.gr()
    #   gr <- gr %>% subsetByOverlaps(motif.gr)
    #   df <- data.frame(gene = utilsFanc::gr.get.loci(gr), motif = motif.name)
    #   return(df)
    # }, threads = threads.motif) %>% do.call(rbind, .)

    df <- lapply(comps, function(comp) {
      if (comp == "beta") {
        if (level.fun == "max") {
          max.df <- data.frame(gene = a2b$bulkNorm$gene,
                               max = rowMax(as.matrix(a2b$bulkNorm[, colnames(a2b$bulkNorm) != "gene"])))
        } else {
          stop("level.fun: only max has been developed")
        }

        df <- left_join(peak.by.motif, a2b[[slot]][,c("gene", "log2FoldChange")])
        df <- left_join(df, max.df)
      } else {
        x <- comp %>% sub("^.+:", "", .)
        y <- comp %>% sub(":.+$", "", .)
        log2fc <- log2((a2b[[slot]][, y] + 1)/(a2b[[slot]][, x] + 1))
        df <- left_join(peak.by.motif, data.frame(gene = a2b[[slot]]$gene, log2FoldChange = log2fc))
        if (level.fun == "max") {
          max.df <- data.frame(gene = a2b$bulkNorm$gene,
                               max = rowMax(as.matrix(a2b$bulkNorm[, c(x, y)])))
        } else {
          stop("level.fun: only max has been developed")
        }
        df <- left_join(df, max.df)
      }
      # divide into levels by quantile:
      df$level <- gtools::quantcut(df[, level.fun], n.levels)
      df$level_id <- df$level %>% as.integer()
      if (write.each.df) {
        dir.create(paste0(plot.dir, "/source/"), showWarnings = F, recursive = T)
        out.Rds <- paste0(plot.dir, "/source/", root.name, "_", a2b$root.name, "_", comp, ".Rds")
        saveRDS(df, out.Rds)
      }

      df$comp <- comp
      return(df)
    }) %>% do.call(rbind, .)
    df$cluster <- a2b$root.name
    return(df)
  }) %>% do.call(rbind, .)

  if (write.df) {
    dir.create(plot.dir, showWarnings = F, recursive = T)
    saveRDS(df, paste0(plot.dir, "/", root.name, ".Rds"))
  }

  # tsv files are too huge
  # write.table(df, paste0(plot.dir, "/", root.name, ".tsv"),
  #             sep = "\t", col.names = T, row.names = F, quote = F)
  # violin plot:
  if (bViolin || bBoxplot) {
    pl.violin <- df %>% split(., f = .$cluster) %>% lapply(function(df.cluster) {
      df.cluster %>% split(., f=.$comp) %>%
        lapply(function(df) {
          if (length(unique(df$motif)) > 4) {
            stop("violin shouldn't be used for more than 4 motifs")
          }
          df$level <- droplevels(df$level)
          nlevels <- max(df$level_id)
          df$quant <- paste0("(", round((df$level_id - 1) /nlevels, digits = 2), ",",
                             round(df$level_id / nlevels, digits = 2), "]")
          df$quant <- sub("\\(0\\.0,", "[0.0,", df$quant)
          df[, value][df[, value] > violin.ceiling] <- violin.ceiling
          df[, value][df[, value] < violin.floor] <- violin.floor

          n.stats <- df[, c("quant", "motif")] %>% table()
          n.out <- paste0(plot.dir, "/", root.name, "_", df$cluster[1], "_", df$comp[1], "_n_numbers.csv")
          dir.create(dirname(n.out), showWarnings = F, recursive = T)
          write.csv(n.stats, n.out)

          p <- ggplot(df, aes_string(x = "quant", y = value))
          if (publication) {
            if (bViolin) {
              p <- p + geom_violin()
            } else{
              stat.test <- df %>%
                group_by(quant) %>%
                rstatix::t_test(as.formula(paste0(value, " ~ motif"))) %>%
                rstatix::adjust_pvalue(method = "fdr") %>%
                rstatix::add_significance("p.adj")
              stat.test <- stat.test %>%
                rstatix::add_xy_position(x = "quant", dodge = 0.8)

              p <- p + geom_boxplot(
                aes(color = motif), size = 0.2,
                outlier.size = 0, outlier.stroke = 0, outlier.shape = 18) +
                geom_jitter(aes(color = motif, fill = motif), size = pt.size, stroke = 0.03, shape = 18,
                            position = position_jitterdodge(jitter.width = 0.25), alpha = alpha) +
                ggpubr::stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0,
                                           size = 10, label.size = 2) +
                scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
              p <- p +  ggtitle(paste0(df$cluster[1], "  ", df$comp[1]))
              p <- utilsFanc::theme.fc.1(p, rm.x.ticks = F, italic.x = F, rotate.x.45 = F)
              dir.create(plot.dir, showWarnings = F, recursive = T)
              file.out <- paste0(plot.dir, "/", root.name, "_", df$cluster[1], "_", df$comp[1], ".pdf")
              ggsave(file.out, p, width = 2, height = 1.2, device = cairo_pdf, dpi = 150)
             }
          } else if (publication.small) {
            stat.test <- df %>%
              group_by(quant) %>%
              rstatix::t_test(as.formula(paste0(value, " ~ motif"))) %>%
              rstatix::adjust_pvalue(method = "fdr") %>%
              rstatix::add_significance("p.adj")
            stat.test <- stat.test %>%
              rstatix::add_xy_position(x = "quant", dodge = 0.8)

            p <- p + geom_boxplot(
              aes(color = motif), size = 0.2,
              outlier.size = 0, outlier.stroke = 0, outlier.shape = 18) +
              geom_jitter(aes(color = motif, fill = motif), size = 0.5, stroke = 0.03, shape = 18,
                          position = position_jitterdodge(jitter.width = 0.25), alpha = 0.4, show.legend = F) +
              ggpubr::stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0,
                                         size = 10, label.size = 2) +
              scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
            p <- p +  ggtitle(paste0(df$cluster[1], "  ", df$comp[1]))
            p <- utilsFanc::theme.fc.1(p, rm.x.ticks = T, italic.x = F, rotate.x.45 = F)
            p <- p + theme(aspect.ratio = 1, axis.text.x = element_blank())
            dir.create(plot.dir, showWarnings = F, recursive = T)
            file.out <- paste0(plot.dir, "/", root.name, "_", df$cluster[1], "_", df$comp[1], ".pdf")
            ggsave(file.out, p, width = 1, height = 1, device = cairo_pdf, dpi = 150)
          }
          else {
            # p <- p + theme(plot.title = element_text(size=8))
            stop("Only publication branch is still maintained")
          }
          return(p)
        })
    }) %>% Reduce(c, .)
    wrap.plots.fanc(plot.list = pl.violin,
                    plot.out = paste0(plot.dir, "/", root.name, "_", ifelse(bViolin, "violin", "boxPlot"), ".png"),
                    sub.width = width, sub.height = height, n.col = n.col)
  }

  if(bDot) {
    pl.dot <- df %>% split(., f = .$cluster) %>% lapply(function(df.cluster) {
      df.cluster %>% split(., f=.$comp) %>%
        lapply(function(df) {
          df.sum <- df %>% dplyr::group_by(motif, level_id) %>%
            dplyr::summarise(median = median(!!as.name(value))) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(rank = rank(median, ties.method = "random"),
                          motif_level = paste0(motif, "_L", level_id)) %>%
            as.data.frame()

          if (write.each.df) {
            dir.create(paste0(plot.dir, "/source/"), showWarnings = F,recursive = T)
            write.table(df.sum, paste0(plot.dir, "/source/", root.name, "_",
                                       df$cluster[1], "_", df$comp[1], "_median.tsv"),
                        sep = "\t", col.names = T, row.names = F, quote = F)
          }

          if (is.null(motifs.label)) {
            if (nrow(df.sum) <= 2 * dot.n.label) {
              motifs.label <- df.sum$motif_level
            } else {
              motifs.label <- df.sum$motif_level[df.sum$rank %in% c(1:dot.n.label, (max(df.sum$rank) - 0:(dot.n.label - 1)))]
            }
          } else{
            motifs.label <- paste0(motifs.label, "_L", 1:n.levels)
          }
          p <- xy.plot(df = df.sum, x = "rank", y = "median", pt.size = 0.2,
                       label.var = "motif_level", label.values = motifs.label, use.repel = T, text.size = 2,
                       add.abline = F, highlight.var = "motif_level", highlight.values = motifs.label,
                       plotly.var = "motif_level", plotly.label.all = T,
                       title = paste0(df$cluster[1], "  ", df$comp[1])) +
            theme(plot.title = element_text(size=8))
          return(p)
        })
    }) %>% Reduce(c, .)
    lapply(c("png", "html"), function(dev) {
      wrap.plots.fanc(plot.list = pl.dot, plot.out = paste0(plot.dir, "/", root.name, "_dot.", dev),
                      sub.width = width, sub.height = height, n.col = NULL, n.split = length(comps),
                      tooltip = "text")
    })

  }
  invisible(df)
}


da.motif.pileup <- function(da, subset.n = 1000, seed = 42,
                            peak.motif.df, motifs.use,
                            samples = NULL, bw.dir,
                            motif.min.score = NULL,
                            n.par = 1, threads.each = 1,
                            plot.dir, root.name = NULL) {
  if (is.null(root.name)) root.name <- basename(plot.dir)
  if (!is.null(motif.min.score)) root.name <- paste0(root.name, "_score", motif.min.score)

  if (is.character(peak.motif.df)) {
    peak.motif.df <- readRDS(peak.motif.df)
    if (!is.null(motif.min.score))
      peak.motif.df <- peak.motif.df %>% dplyr::filter(score >= motif.min.score)
  }

  utilsFanc::check.intersect(motifs.use, "motifs.use",
                             peak.motif.df$motif, "peak.motif.df$gene")

  utilsFanc::safelapply(da, function(a2b) {
    loci <- a2b$bulkNorm$gene
    utilsFanc::check.intersect(loci, "a2b$bulk.norm$gene",
                               peak.motif.df$gene, "peak.motif.df$gene")

    all.samples <- a2b$bulk.mat %>% colnames()
    if (is.null(samples)) samples <- all.samples
    samples <- samples[samples %in% all.samples]
    if (length(samples) < 1) {
      stop("length(samples) < 1")
    }

    glob <- paste0(bw.dir, "/", ifelse(grepl("^\\d", samples[1]), "X", ""), samples, "*")
    bws <- utilsFanc::glob.single.files(x = glob)
    names(bws) <- samples

    beds <- sapply(motifs.use, function(motif.use) {
      df <- peak.motif.df %>% dplyr::filter(gene %in% loci, motif == motif.use)
      beds <- sapply(c("w", "wo"), function(type) {
        bed.dir <- paste0(plot.dir, "/", root.name, "_beds/")
        bed <- paste0(bed.dir, "/", root.name, "_", a2b$root.name, "_", motif.use, "_", type, ".bed")
        if (type == "w")
          out <- df$gene
        else
          out <- loci[!loci %in% df$gene]
        if (!is.null(subset.n) && length(out) > subset.n) {
          set.seed(seed = seed)
          out <- out[sort(sample(1:length(out), size = subset.n, replace = F))]
        }
        utilsFanc::write.zip.fanc(utilsFanc::loci.2.gr(out), out.file = bed, bed.shift = T)
        return(bed)
      })
      return(beds)
    }) %>% as.character()
    names(beds) <- basename(beds) %>%
      sub(paste0(a2b$root.name, "_"), "", ., fixed = T) %>%
      sub(".bed", "", .)

    liteRnaSeqFanc::deeptools.refpoint(
      bw.vec = bws, regions.vec = beds, threads = threads.each, out.dir = plot.dir,
      root.name = paste0(root.name, "_", a2b$root.name))
    return()
  }, threads = n.par)
  return()
}


de.plot.bar.2 <- function(de, genes, cluster = 1, color.by, color.levels = NULL,
                          shape.by = NULL,
                          spread.width = 0.3, spread.bin.size = 100,
                          pt.size = 0.8,
                          wrap.n = NULL,
                          plot.out = NULL,
                          add.pval = T, pval.group.1 = NULL, pval.group.2,
                          external.p.df = NULL,
                          no.legend = F,
                          ...) {

  genes <- genes %>% .[. %in% de[[cluster]]$res.exp$gene]
  if (length(genes) < 1) {
    stop("length(genes) < 1")
  }

  if (is.null(wrap.n)) {
    f <- rep(1, length(genes))
  } else {
    f <- ceiling((1:length(genes))/wrap.n)
  }
  pl <- genes %>% split(., f = f) %>%
    lapply(function(genes) {
      bar.df <- de.plot.bar(de = de, genes = genes, clusters = cluster, return.df = T)
      if (!is.null(color.levels)) {
        bar.df[, color.by] <- factor(bar.df[, color.by], levels = color.levels)
      }
      n.genes <- bar.df$gene %>% unique() %>% length()
      
      p <- utilsFanc::barplot.pub.3(df = bar.df, x = "gene", y = "expr",
                         color.by = color.by,
                         shape.by = shape.by,
                         pt.size = pt.size,
                         spread.bin.size = spread.bin.size,
                         bar.width = 0.6, dodge.width = 0.8, spread.width = spread.width,
                         add.pval = add.pval, pval.group.1 = pval.group.1, pval.group.2 = pval.group.2, 
                         external.p.df = external.p.df,
                         ...) %>%
        utilsFanc::theme.fc.1()
      
      if (no.legend) {
        p <- p + theme(legend.position = "none")
      }
      return(p)
    })
  if (length(pl) > 1) {
    p <- wrap.plots.fanc(pl, n.col = 1)
  } else {
    p <- pl[[1]]
  }
  if (!is.null(plot.out)) {
    dir.create(dirname(plot.out), showWarnings = F, recursive = T)
    ggsave(plot.out, p, device = cairo_pdf,
           width = 0.5 * length(genes), height = 3 * length(unique(f)))

  }
  invisible(p)
}

de.ggVenn <- function(de.1, de.2, name.1, name.2,
                      slot = "summary", plot.out = NULL) {
  clusters <- intersect(names(de.1), names(de.2))
  if (length(clusters) < 1) {
    stop(length(clusters) < 1)
  }
  de.1 <- de.1[clusters]
  de.2 <- de.2[clusters]
  pl <- lapply(clusters, function(cluster) {
    pl <- lapply(c("up", "down"), function(type) {
      del <- list(de.1, de.2)
      names(del) <- c(name.1, name.2)
      gl <- lapply(del, function(de) {
        genes <- de[[cluster]][[slot]][[paste0(type, ".genes")]]
        if (is.null(genes)) {
          stop(paste0("is.null(de$",cluster,"$",slot, "$", type, ".genes)"))
        }
        return(genes)
      })
      names(gl) <- paste0(names(gl), "_", type)
      p <- ggVennDiagram::ggVennDiagram(gl, set_size = 3) +
        scale_fill_gradient(low="#2297E6", high = "#F5C710") +
        ggtitle(paste0(cluster, " ", type))
      return(p)
    })
    return(pl)
  }) %>% Reduce(c, .)
  p <- wrap.plots.fanc(plot.list = pl, n.split = 2, plot.out = plot.out)
  invisible(p)
}

de.distro <- function(de, use.samples = NULL,
                      transforms = c("linear", "log2"),
                      rm.zero = F, # if F, just use the bulkNorm matrix, which carries the filtering from s2b.
                      rm.zero.samples = NULL, # only use genes where all these samples are not zero
                      rm.zero.n.samples = 0.5, # only use genes where at least this # samples are not zero
                      # note: if (0,1): used as fraction of n.samples; if >= 1, used as # samples
                      # to indicate all samples, use something like 0.99999
                      out.dir, root.name = NULL) {
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  pl <- lapply(de, function(s2b) {
    mat <- s2b$bulkNorm %>% dplyr::mutate(gene = NULL) %>% as.matrix()
    if (!is.null(use.samples)) {
      if (length(intersect(colnames(mat), use.samples)) == 0) {
        stop("none of the samples offered in use.samples were found")
      }
      mat <- mat[, colnames(mat) %in% use.samples, drop = F]
    }
    if (rm.zero) {
      if (!is.null(rm.zero.samples)) {
        utilsFanc::check.intersect(rm.zero.samples, "rm.zero.samples",
                                   colnames(mat), "colnames(s2b$bulkNorm)")
        bRm <- rowSums(mat[, rm.zero.samples, drop = F] > 0) == length(rm.zero.samples)
      } else {
        if (rm.zero.n.samples < 1) {
          n.nz <- round(ncol(mat) * rm.zero.n.samples, digits = 0)
        } else {
          n.nz <- rm.zero.n.samples
          if (n.nz > ncol(mat)) stop(paste0("rm.zero.n.samples is higher than the number of samples"))
        }
        bRm <- rowSums(mat > 0) >= n.nz

      }
      mat <- mat[bRm, ]
    }

    pl <- lapply(transforms, function(transf) {
      trf <- ifelse(transf == "linear", function(x) x,
                       ifelse(transf == "log2", function(x) log2(x + 1), NA))
      if (is.na(trf)) stop(paste0("transform ", transf, " not recognized"))

      mat <- trf(mat)
      df <- as.data.frame(mat)
      df.melt <- reshape2::melt(df, variable.name = "sample", value.name = "expr")
      title <- paste0(s2b$root.name, " ", transf)
      pl <- list()
      pl$rank <- rank.plot(df = df, vars = colnames(df)) +
        ggtitle(label = title)
      pl$distro <- ggplot(df.melt, aes(x = expr, color = sample, fill = sample)) +
        geom_density(alpha = 0.3) +
        theme(aspect.ratio = 1) +
        ggtitle(label = title)
      return(pl)
    }) %>% Reduce(c, .)
    return(pl)
  }) %>% Reduce(c, .)
  suffix <- paste0(rm.zero, "_", paste0(rm.zero.samples, collapse = "."), "_",
                   ifelse(!is.null(rm.zero.samples), "", rm.zero.n.samples))
  plot.out <- paste0(out.dir, "/", root.name, "_", suffix, ".png")
  p <- wrap.plots.fanc(plot.list = pl, n.split = 4, plot.out = plot.out)
  invisible(p)
}

de.add.promoter <- function(de, gtf.gr, slot = "summary",
                            out.dir = NULL, root.name = NULL) {
  if (!any(sapply(de, function(x) !is.null(x[[slot]])))) {
    stop("slot not found for any clusters within de")
  }
  if (!is.null(out.dir)) {
    if (is.null(root.name)) root.name <- basename(out.dir)
  }
  if (is.character(gtf.gr)) {
    gtf <- rtracklayer::import(gtf.gr)
  }
  de <- lapply(de, function(s2b) {
    if (is.null(s2b[[slot]])) {
      return(s2b)
    }
    pros <- list()
    useful.genes <- unname(unlist(s2b[[slot]][c("up.genes", "down.genes", "de.genes")]))
    pros$de <- gtf.gr %>%
      plyranges::filter(gene_name %in% useful.genes,
                        type == "transcript",
                        transcript_type == "protein_coding") %>%
      resize(width = 1, fix = "start", ignore.strand = F) %>% as.data.frame() %>%
      dplyr::rename(chr = seqnames, gene = gene_name) %>%
      dplyr::select(chr, start, end, gene) %>% unique() %>%
      makeGRangesFromDataFrame(keep.extra.columns = T)
    pros$up <- pros$de[pros$de$gene %in% s2b[[slot]]$up.genes]
    pros$down <- pros$de[pros$de$gene %in% s2b[[slot]]$down.genes]
    s2b[[slot]]$promoters <- pros

    if (!is.null(out.dir)) {
      utilsFanc::write.zip.fanc(
        df = pros$de, out.file = paste0(out.dir, "/", root.name, "_", s2b$root.name, "_pros.bed"),
        bed.shift = T)
    }

    return(s2b)
  })
  return(de)
}

de.add.promoter.all.genes <- function(de, gtf.gr, promoter.slot.name = "promoter",
                                      out.dir = NULL, root.name = NULL) {
  stop("Not finished/tested")
  if (!is.null(out.dir)) {
    if (is.null(root.name)) root.name <- basename(out.dir)
  }
  if (is.character(gtf.gr)) {
    gtf <- rtracklayer::import(gtf.gr)
  }
  de <- lapply(de, function(s2b) {
    s2b[[promoter.slot.name]] <- gtf.gr %>%
      plyranges::filter(gene_name %in% s2b$bulkNorm$gene, type == "transcript",
                        transcript_type == "protein_coding") %>%
      resize(width = 1, fix = "start", ignore.strand = F) %>% as.data.frame() %>%
      dplyr::rename(chr = seqnames, gene = gene_name) %>%
      dplyr::select(chr, start, end, gene) %>% unique() %>%
      makeGRangesFromDataFrame(keep.extra.columns = T)
    return(s2b)
  })
}

gene.promoter.map <- function(genes, gtf.gr) {
  if (is.character(gtf.gr))
    gtf.gr <- rtracklayer::import(gtf.gr)
  pros <- gtf.gr %>%
    plyranges::filter(gene_name %in% genes, type == "transcript",
                      transcript_type == "protein_coding") %>%
    resize(width = 1, fix = "start", ignore.strand = F) %>% as.data.frame() %>%
    dplyr::rename(chr = seqnames, gene = gene_name) %>%
    dplyr::select(chr, start, end, gene) %>% unique() %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  # if (collpase.each.gene) {
  #   pros <- pros %>% split(f = pros$gene) %>% lapply(function(gr) {
  #     gene <- gr$gene[1]
  #     gr <- GenomicRanges::reduce()
  ##### it is impractical to collapse for each gene. What if the 2 intervals for each gene do not
  # overlap to begin with?
  #   })
  # }
}

da.add.anno <- function(da, genome, slot = "summary",
                        out.dir = NULL, root.name = NULL) {
  if (!is.null(out.dir)) {
    if (is.null(root.name)) root.name <- basename(out.dir)
  }
  da <- lapply(da, function(a2b) {
    if (is.null(a2b[[slot]])) {
      return(a2b)
    }
    gr <- a2b[[slot]]$de.genes %>%
      utilsFanc::loci.2.df(loci.vec = ., return.gr = T, remove.loci.col = F)
    gr <- utilsFanc::gr.fast.annot(gr, genome = genome)
    gr$gene <- gr$SYMBOL
    mcols(gr)$SYMBOL <- NULL

    anno <- list(de = gr, up = gr[gr$loci %in% a2b[[slot]]$up.genes],
                 down = gr[gr$loci %in% a2b[[slot]]$down.genes])
    a2b[[slot]]$anno <- anno

    if (!is.null(out.dir)) {
      utilsFanc::write.zip.fanc(
        anno$de %>% as.data.frame() %>% dplyr::select(seqnames, start, end, gene),
        out.file = paste0(out.dir, "/", root.name, "_", a2b$root.name, "_anno.bed"),
        bed.shift = T)
    }

    return(a2b)
  })
}

de.add.bg <- function(de, slot = "summary",
                      n.bg.each = 25, split.into.folds = F, split.into.folds.random.seed = 42,
                      no.replace = F,
                      samples.for.bg, scatter.xy.df, scatter.meta = NULL,
                      write.bed = F,
                      work.dir, root.name = NULL) {
  if (is.null(root.name)) root.name <- basename(work.dir)
  de <- lapply(de, function(s2b) {
    if (is.null(s2b[[slot]])) {
      return(s2b)
    }
    norm.mat <- s2b$bulkNorm %>% .[, colnames(.) != "gene"] %>% as.matrix()
    rownames(norm.mat) <- s2b$bulkNorm$gene
    bg <- lapply(c("up", "down"), function(type) {
      fg.vec <- s2b[[slot]][[paste0(type, ".genes")]]
      if (length(fg.vec) < 1) {
        return(c())
      }
      prefix <- paste0(work.dir, "/", root.name, "_", s2b$root.name, "_", type)
      utilsFanc::check.intersect(x = samples.for.bg, "samples.for.bg",
                                 colnames(norm.mat), "colnames(norm.mat)")

      bg.vec <- bg.gen.2(mat = norm.mat, fg.vec = fg.vec, n.bg.each = n.bg.each, samples.use = samples.for.bg,
                         no.replace = no.replace)

      if (write.bed && grepl(":", fg.vec[1])) {
        fg.bed <- paste0(prefix, "_fg.bed")
        trash <- utilsFanc::write.zip.fanc(df = utilsFanc::loci.2.df(loci.vec = fg.vec)[, 1:3], out.file = fg.bed)
        bg.bed <- paste0(prefix, "_bg.bed")
        trash <- utilsFanc::write.zip.fanc(df = utilsFanc::loci.2.df(loci.vec = bg.vec)[, 1:3], out.file = bg.bed)
      }
      stats <- data.frame(
        cluster = s2b$root.name,
        direc = type,
        n.fg = length(fg.vec),
        n.bg.expected = length(fg.vec) * n.bg.each,
        n.bg = length(bg.vec)
      ) %>% dplyr::mutate(
        frac = round(n.bg/n.bg.expected, digits = 3)
      )
      stats.file <- paste0(prefix, "_stats.tsv")
      dir.create(dirname(stats.file), showWarnings = F, recursive = T)
      write.table(stats, stats.file, append = F, quote = F, col.names = T,
                  row.names = F, sep = "\t")
      bg.assess(mat = norm.mat, fg.vec = fg.vec, bg.vec = bg.vec,
                scatter.xy.df = scatter.xy.df, scatter.meta = scatter.meta, col.data = s2b$coldata,
                transformation = function(x) return(log2(x + 1)),
                plot.out = paste0(prefix, "_asssess.png" ))
      return(bg.vec)
    })
    names(bg) <- c("up", "down")
    if (split.into.folds) {
      set.seed(seed = split.into.folds.random.seed)
      bg <- lapply(bg, function(x) {
        if(length(x) > 1) {
          x <- sample(x, size = length(x), replace = F) # shuffle x
          f <- cut(1:length(x), breaks = n.bg.each) %>% as.numeric()
          x <- split(x, f = f)
          names(x) <- paste0("fold", names(x))
          return(x)
        } else {
          x <- rep(list(character()), n.bg.each)
          names(x) <- paste0("fold", 1:length(x))
          return(x)
        }

      })
    }
    s2b[[slot]]$bg <- bg
    return(s2b)
  })
  return(de)
}

de.add.bg.validate <- function(de, slot = "summary", n.folds,
                               cols.check = c("no.overlap", "split.into.folds", "number.match.fg", "n.folds.match")) {
  # this is to validate that a de object generated by de.add.bg() has valid slots
  de.slot.assert(de, is.s2b = F, slot = slot)
  assess.df <- lapply(de, function(s2b) {
    lapply(c("up", "down"), function(direction) {
      out <- data.frame(cluster = s2b$root.name, direction = direction,
                        bg.exists = NA, not.empty = NA,
                        no.overlap = NA, split.into.folds = NA, number.match.fg = NA,
                        n.folds.match = NA)
      bg <- s2b[[slot]]$bg[[direction]]

      if (is.null(bg)) {
        out$bg.exists <- F
      } else {
        out$bg.exists <- T

        if (length(unlist(bg)) < 1) {
          out$not.empty <- F
        } else {
          out$not.empty <- T

          if (any(duplicated(unlist(bg)))) {
            out$no.overlap <- F
          } else {
            out$no.overlap <- T
          }

          if (any(!grepl("^fold\\d+$", names(bg)))) {
            out$split.into.folds <- F
          } else {
            out$split.into.folds <- T
            fg <- s2b[[slot]][[paste0(direction, ".genes")]]
            if (any(sapply(bg, length) != length(fg))) {
              out$number.match.fg <- F
            } else {
              out$number.match.fg <- T
            }
            folds <- as.numeric(stringr::str_extract(names(bg), "\\d+"))
            if (identical(folds, as.numeric(1:n.folds))) {
              out$n.folds.match <- T
            } else {
              out$n.folds.match <- F
            }

          }
        }
      }
      return(out)
    }) %>% do.call(rbind, .) %>% return()
  }) %>% do.call(rbind, .)
  rownames(assess.df) <- NULL
  utilsFanc::check.intersect(cols.check, "cols.check", colnames(assess.df), "colnames(assess.df)")

  lapply(cols.check, function(col.check) {
    names <- paste0(assess.df$cluster, "_", assess.df$direction)
    vec <- assess.df[, col.check]
    vec[is.na(vec)] <- T
    failed <- names[vec == F]
    if (length(failed) > 0) {
      stop(paste0("Failed check: ", col.check, "\n",
                  "clusters that failed: \n",
                  paste0(failed, collapse = "\n")))
    }
  })

  return(assess.df)
}

de.bg.expand <- function(de, slot = "summary", n.folds,
                         cols.check = c("no.overlap", "split.into.folds", "number.match.fg", "n.folds.match")) {
  # change from s2b$summary$bg$up$fold1 to s2b$fold1$up.genes. Basically changing it
  # to the same format as summary.
  de.add.bg.validate(de = de, slot = slot, n.folds = n.folds, cols.check = cols.check)
  de <- lapply(de, function(s2b) {
    bg <- s2b[[slot]]$bg
    for (i in 1:n.folds) {
      fold.id <- paste0("fold", i)
      bg.summ <- list(up.genes = bg$up[[fold.id]],
                      n.up = length(bg$up[[fold.id]]),
                      down.genes = bg$down[[fold.id]],
                      n.down = length(bg$down[[fold.id]]))
      bg.summ$de.genes <- c(bg.summ$up.genes, bg.summ$down.genes) %>% unique()
      bg.summ$n.de <- length(bg.summ$de.genes)
      s2b[[paste0("fold", i)]] <- bg.summ
    }
    return(s2b)
  })
  # double check to make sure we are looking at the same things:
  lapply(de, function(s2b) {
    lapply(c("up", "down"), function(direction) {
      new <- lapply(1:n.folds, function(i) {
        s2b[[paste0("fold", i)]][[paste0(direction, ".genes")]]
      })
      names(new) <- paste0("fold", 1:n.folds)
      old <- s2b[[slot]]$bg[[direction]]
      if (!identical(new, old)) {
        stop(paste0("Failed checking at: ", s2b$root.name))
      }

    })
  })
  return(de)

}

de.motif.n <- function(de, peak.motif.df,
                       gene.peak.df = NULL, gene.motif.df = NULL,
                       peaks, promoters.gr,
                       # peaks and promoters.gr might be needed if gene.motif.df or peak.motif.df are not supplied.
                       binarize = F, normalize = F,
                       slot = "summary", directions = c("up", "down"),
                     motifs.use, motifs.write.plot = NULL,
                     n.bg.sets = 50, with.replacement = T, seed = 42,
                     method = "distance", distance.1side = 50000,
                     write.bed = F, write.details = F,
                     plot.only = F,
                     plot.each = F, plot.z = T, motifs.label = NULL, n.label.1side = 8,
                     work.dir, root.name = NULL,
                     return.gene.motif.df = F) {
  # binarize: for each gene, instead of counting how many AP1 motifs it has, we only
  # count whether it has an AP1 motif at all. We then normalize through dividing by the total number of genes in the group
  # normalize: sum(# of peaks in the group with the motif)/sum(# of peaks in the group)

  if (is.null(root.name)) root.name <- basename(work.dir)
  if (binarize && normalize) {
    stop("binarize and normalize can't be used together!")
  }
  if (binarize) root.name <- paste0(root.name, "_bin")
  if (normalize) root.name <- paste0(root.name, "_norm")

  # check.list <- list(motifs.use = motifs.use, method = method,
  #                    distance.1side = distance.1side)
  # lapply(names(check.list), function(name) {
  #   param <- check.list[[name]]
  #   if (name == "motifs.use") {
  #     # this check is very necessary. If you don't see a motif, there is
  #     # a chance that you just didn't include it into the calculation of gene.motif.df
  #     utilsFanc::check.intersect(param, name, attr(gene.motif.df, name),
  #                                paste0("attr(gene.motif.df, ", name, ")"))
  #   } else {
  #     if (param != attr(gene.motif.df, name)) {
  #       stop(paste0(name, " is not consistent with attr(gene.motif.df, ", name, ")"))
  #     }
  #   }
  # })

  if (normalize && is.null(gene.peak.df)) {
    utilsFanc::t.stat("Generating gene.peak.df")
    gene.peak.df <- gene.peak.df.by.dist(genes = NULL, promoters.gr = promoters.gr,
                                         peaks = peaks, distance.1side = distance.1side)
  }

  if (is.null(gene.motif.df)) {
    # > head(gene.motif.df %>% filter(gene == "Socs3"))
    #    gene            motif n
    # 1 Socs3     AMYB.HTH_167 3
    # 2 Socs3      AP.1.bZIP_1 2
    # 3 Socs3  AP.2alpha.AP2_3 8
    # 4 Socs3  AP.2gamma.AP2_2 8
    # 5 Socs3       Ap4.bHLH_4 5
    # 6 Socs3 AR.halfsite.NR_7 9
    #                                                                                                                                                                                                                                        peak
    # 1                                                                                                                                                             chr11:117891507-117892007,chr11:117946575-117947075,chr11:117958131-117958631
    # 2                                                                                                                                                                                       chr11:117989715-117990215,chr11:118000470-118000970
    # 3                           chr11:117871891-117872391,chr11:117878757-117879257,chr11:117891507-117892007,chr11:117935630-117936130,chr11:117954045-117954545,chr11:117974674-117975174,chr11:117979355-117979855,chr11:117991666-117992166
    # 4                           chr11:117871891-117872391,chr11:117878226-117878726,chr11:117878757-117879257,chr11:117891507-117892007,chr11:117935630-117936130,chr11:117954045-117954545,chr11:117974674-117975174,chr11:117979355-117979855
    # 5                                                                                                         chr11:117952126-117952626,chr11:117968120-117968620,chr11:117991164-117991664,chr11:118000470-118000970,chr11:118016157-118016657
    # 6 chr11:117873234-117873734,chr11:117881790-117882290,chr11:117894727-117895227,chr11:117946575-117947075,chr11:117952126-117952626,chr11:117954045-117954545,chr11:117974674-117975174,chr11:118014533-118015033,chr11:118021728-118022228
    utilsFanc::t.stat("Generating gene.motif.df")
    gene.motif.df <- archr.gene.motif.df.gen(
      promoters.gr = promoters.gr,
      peak.motif.df = peak.motif.df[peak.motif.df$gene %in% peaks, ],
      method = method, distance.1side = distance.1side)
    if (return.gene.motif.df) return(gene.motif.df)
  }

  if (is.null(motifs.use)) motifs.use <- gene.motif.df$motif %>% unique()

  j <- lapply(de, function(s2b) {
    if (is.null(s2b[[slot]])) return()
    j <- lapply(directions, function(direc) {
      prefix <- paste0(work.dir, "/", root.name, "_", s2b$root.name, "_", direc)

      if (is.null(s2b[[slot]]$bg[[direc]])) {
        stop(paste0("is.null(s2b[[slot]]$bg[[direc]]) for ", s2b$root.name))
      }
      fg <- s2b[[slot]][[paste0(direc, ".genes")]]
      bg <- s2b[[slot]]$bg[[direc]]
      if (length(bg) < 1) return()
      if (with.replacement) {
        bg.df <- lapply(1:n.bg.sets, function(bg.id) {
          set.seed(seed = seed + bg.id)
          bg.selected <- sample(bg, length(fg), replace = F)
          df <- data.frame(group = paste0("bg", bg.id),
                           gene = bg.selected)
        }) %>% do.call(rbind, .)
      } else {
        stop("only with.replacement = T has been developed")
      }

      df <- rbind(data.frame(group = "fg", gene = fg), bg.df)
      #> head(df)
      #   group    gene
      # 1    fg   Tgfbi
      # 2    fg  Ifitm1
      # 3    fg    Pim1
      # 4    fg Emilin2
      # 5    fg Slc30a4
      # 6    fg    Klf9
      if (normalize) {
        gene.peak.n <- gene.peak.df %>% dplyr::select(gene, peak) %>%
          unique() %>% dplyr::group_by(gene) %>%
          dplyr::summarise(n.enh = n()) %>% dplyr::ungroup() %>%
          as.data.frame()
        group.peak.n <- dplyr::left_join(df, gene.peak.n) %>%
          dplyr::mutate(n.enh = ifelse(is.na(n.enh), 0, n.enh)) %>%
          dplyr::group_by(group) %>%
          dplyr::summarise(n.enh = sum(n.enh))
        dir.create(dirname(prefix), showWarnings = F, recursive = T)
        tempt.list <- list(gene.peak.df = gene.peak.df, gene.peak.n = gene.peak.n, group.peak.n = group.peak.n)
        lapply(names(tempt.list), function(x) {
          to.write <- tempt.list[[x]]
          write.table(to.write, paste0(prefix, "_normFactors_", x, ".tsv"),
                      quote = F, sep = "\t", row.names = F, col.names = T)
        })

      }

      # Inner join seems wrong. However, the lost genes will be later "reclaimed"
      # using tidyr::comlete.
      j <- dplyr::inner_join(df, gene.motif.df, by = "gene", relationship = "many-to-many")

      # df1 <- data.frame(group = c("a", "b", "c"), gene = c("gene1", "gene2", "gene2"))
      # df2 <- data.frame(gene = c("gene1","gene2", "gene2"), motif = c("motif1", "motif1", "motif2"))
      # Browse[1]>       inner_join(df1, df2, by = "gene")
      # group  gene  motif
      # 1     a gene1 motif1
      # 2     b gene2 motif1
      # 3     b gene2 motif2
      # 4     c gene2 motif1
      # 5     c gene2 motif2
      if (write.details) {
        dir.create(dirname(prefix), showWarnings = F, recursive = T)
        write.table(j %>% dplyr::arrange(motif), file = paste0(prefix, "_details.tsv"),
                    sep = "\t", row.names = F, col.names = T, quote = F)
      }
# Browse[1]> head(j[j$n > 1,])
#    group    gene    motif                                                                    peak n
# 7     fg    Pim1    Stat3 chr17:29491880-29492380,chr17:29497995-29498495,chr17:29500708-29501208 3
# 9     fg Emilin2 CEBP.AP1                         chr17:71264679-71265179,chr17:71299833-71300333 2
# 10    fg Emilin2    Stat3                         chr17:71264679-71265179,chr17:71299833-71300333 2
# 13    fg    Klf9 CEBP.AP1                         chr19:23134768-23135268,chr19:23139947-23140447 2
# 15    fg   Ddit4 CEBP.AP1                         chr10:59943420-59943920,chr10:59952456-59952956 2
# 21    fg    Lifr CEBP.AP1                             chr15:7092935-7093435,chr15:7128667-7129167 2

      if (write.bed) {
        if (is.null(motifs.write.plot)) {
          motifs.write.plot <- j$motif %>% unique()
        } else {
          utilsFanc::check.intersect(motifs.write.plot, "motifs.write.plot",
                                     j$motif, "j$motif")
        }

        j %>% dplyr::filter(motif %in% motifs.write.plot) %>% split(., f = .$motif) %>% lapply(function(j) {
          motif <- j$motif[1]
          j <- tidyr::separate_rows(j, peak, sep = ",") %>% as.data.frame()
          j <- utilsFanc::loci.2.df(df = j, loci.col.name = "peak", remove.loci.col = T, return.gr = F)
          j$forth <- paste0(j$gene, "_", j$group)
          j <- j[, c("chr", "start", "end", "forth")]
          utilsFanc::write.zip.fanc(
            df = j, out.file = paste0(
              prefix, "_beds/", basename(prefix), "_", motif, ".bed"), bed.shift = T)
          return()
        })
      }

      j$group <- factor(j$group, levels = unique(df$group))
      j$motif <- factor(j$motif, levels = motifs.use)
      j <- tidyr::complete(j, group, motif)
# > df <- tibble(
# +   group = c(1:2, 1),
# +   item_id = c(1:2, 2),
# +   item_name = c("a", "b", "b"),
# +   value1 = 1:3,
# +   value2 = 4:6
# + )
# > df
# # A tibble: 3 x 5
#   group item_id item_name value1 value2
#   <dbl>   <dbl> <chr>      <int>  <int>
# 1     1       1 a              1      4
# 2     2       2 b              2      5
# 3     1       2 b              3      6
# > df$item_name <- factor(df$item_name, levels = c("a", "b", "c"))
# > df %>% tidyr::complete(group, item_id, item_name)
# # A tibble: 12 x 5
#    group item_id item_name value1 value2
#    <dbl>   <dbl> <fct>      <int>  <int>
#  1     1       1 a              1      4
#  2     1       1 b             NA     NA
#  3     1       1 c             NA     NA
#  4     1       2 a             NA     NA
#  5     1       2 b              3      6
#  6     1       2 c             NA     NA
#  7     2       1 a             NA     NA
#  8     2       1 b             NA     NA
#  9     2       1 c             NA     NA
# 10     2       2 a             NA     NA
# 11     2       2 b              2      5
# 12     2       2 c             NA     NA
      j$n[is.na(j$n)] <- 0
      # > j %>% filter(motif == "CDX4.Homeobox_31") %>% as.data.frame() %>% head()
      #   group            motif    gene n                     peak
      # 1    fg CDX4.Homeobox_31    <NA> 0                     <NA>
      # 2   bg1 CDX4.Homeobox_31    <NA> 0                     <NA>
      # 3   bg2 CDX4.Homeobox_31  Zfp758 1  chr17:22361457-22361957
      # 4   bg3 CDX4.Homeobox_31    <NA> 0                     <NA>
      # 5   bg4 CDX4.Homeobox_31 Cysltr1 1 chrX:106603313-106603813
      # 6   bg5 CDX4.Homeobox_31    <NA> 0                     <NA>

      j <- j %>% dplyr::group_by(group, motif)
      #>>>>>>>>>>>>>>>>>>>>>> Summarizing starts here. j no longer holds peak-level info.

      if (binarize) {
        j <- j %>% dplyr::summarise(n = sum(n!=0))
      } else {
        j <- j %>% dplyr::summarise(n = sum(n))
      }

      j <- j %>% dplyr::ungroup() %>% as.data.frame()
      j <- utilsFanc::factor2character.fanc(j)

      n.combs <- length(motifs.use) * length(unique(df$group))
      if (nrow(j) != n.combs) {
        stop(paste0("There are ", length(unique(df$group)), " groups and ",
                    length(motifs.use), " motifs. Therefore, we expect ", n.combs,
                    " lines in j. However, we got ", nrow(j)))
      }

      if (normalize) {
        j <- dplyr::left_join(j, group.peak.n, by = "group") %>%
          dplyr::mutate(mean = round(n/n.enh, digits = 3))
      } else {
        j$mean <- round(j$n/length(fg), digits = 3)
      }

      j <- cbind(data.frame(cluster = s2b$root.name, direc = direc), j)

      stats <- j %>% split(., f = .$motif) %>% lapply(function(j) {
        bg <- j %>% dplyr::filter(grepl("bg", group)) %>% dplyr::pull(mean)
        fg <- j %>% dplyr::filter(group == "fg") %>% dplyr::pull(mean)
        m <- mean(bg)
        sd <- sd(bg)
        ECDF <- ecdf(bg)
        stats <- data.frame(motif = j$motif[1], fg = fg, bg.mean = m, bg.sd = sd,
                            z = round((fg - m)/sd, digits = 2), p = 1 - ECDF(fg))
        return(stats)
      }) %>% do.call(rbind, .)
      stats$rank <- rank(stats$z, ties.method = "random")
      stats$rank.rev <- max(stats$rank) - stats$rank + 1
      stats$rank.pct <- round(stats$rank/max(stats$rank), digits = 3)
      write.table(stats, paste0(prefix, "_z.tsv"), sep = "\t", quote = F,
                  row.names = F, col.names = T)

      if (plot.z) {
        if (is.null(motifs.label)) {
          if (nrow(stats) <= 2 * n.label.1side) {
            motifs.label <- stats$motif
          } else {
            motifs.label <- stats$motif[
              stats$rank %in% c(1:n.label.1side, (max(stats$rank) - 0:(n.label.1side - 1)))]
          }
        }
        utilsFanc::check.intersect(motifs.label, "motifs.label",
                                   stats$motif, "stats$motif")
        motifs.w.na.zscore <- stats$motif[is.na(stats$z)] %>% unique()
        motifs.w.na.zscore <- paste0(motifs.w.na.zscore, collapse = " ")
        stats.plot <- stats[!is.na(stats$z),]

        p <- xy.plot(df = stats.plot, x = "rank", y = "z", pt.size = 0.2,
                     label.var = "motif", label.values = motifs.label, use.repel = T, text.size = 2,
                     add.abline = F, highlight.var = "motif", highlight.values = motifs.label,
                     plotly.var = "motif", plotly.label.all = T, title = motifs.w.na.zscore) +
          theme(plot.title = element_text(size=3))
        wrap.plots.fanc(plot.list = list(p), plot.out = paste0(prefix, "_z.png"))
      }

      if (plot.each) {
        if (is.null(motifs.write.plot)) {
          motifs.write.plot <- j$motif %>% unique()
        } else {
          utilsFanc::check.intersect(motifs.write.plot, "motifs.write.plot",
                                     j$motif, "j$motif")
        }
        pl <- lapply(motifs.write.plot, function(motif.use) {
          j <- j %>% dplyr::filter(motif == motif.use)
          stat <- stats %>% dplyr::filter(motif == motif.use)
          p <- ggplot(j[grepl("bg", j$group),], aes(x = mean)) +
            geom_density() +
            geom_vline(xintercept = j[j$group == "fg","mean"]) +
            ggtitle(paste0(motif.use),
                    subtitle = paste0("z: ", stat$z, "\n",
                      "Top ", stat$rank.rev, "; ", stat$rank.pct))
          return(p)
        })
        wrap.plots.fanc(plot.list = pl, plot.out = paste0(prefix, ".png"))
      }
      return(j)
    }) %>% do.call(rbind, .)
    return(j)
  }) %>% do.call(rbind, .)

  dir.create(work.dir, showWarnings = F, recursive = T)
  write.table(arrange(j, motif), paste0(work.dir, "/", root.name, "_j.tsv"),
              sep = "\t", quote = F, row.names = F, col.names = T)
  invisible(j)
}

de.motif.n.rank.plot <- function(z.tsv, motifs.label, title = NULL,
                                 plot.out, height = 1, width = 1,
                                 use.repel = T,
                                 ...) {
  zdf <- read.table(z.tsv, header = T, sep = "\t")
  p <- xy.plot(df = zdf, x = "rank", y = "z", y.lab = "z-score",
               pt.size = 0.1, highlight.ptsize = 0.1, 
               label.var = "motif", label.values = motifs.label, use.repel = use.repel, 
               text.size = 2, 
               add.abline = F, highlight.var = "motif", highlight.values = motifs.label,
               plotly.var = "motif", plotly.label.all = T, title = title, ...) %>% 
    utilsFanc::theme.fc.1(rm.x.ticks = F, italic.x = F, font = NULL)

  dir.create(dirname(plot.out), showWarnings = F, recursive = T)
  ggsave(plot.out, p, device = pdf, height = height, width = width)
  invisible(p)
}

de.filter.close.genes <- function(de, gtf, slot = "summary", cutoff = 100000) {
  if (is.character(gtf)) {
    gtf <- rtracklayer::import(gtf)
  }
  de <- lapply(de, function(s2b) {
    if (is.null(s2b[[slot]])) return(s2b)
    if (!is.null(s2b[[slot]]$close.filter)) {
      stop("close.filter already exists")
    }
    not.in.gtf <- s2b[[slot]]$de.genes %>% .[!. %in% gtf$gene_name]
    items <- c("up", "down")
    res <- lapply(items, function(x) {
      genes <- s2b[[slot]][[paste0(x, ".genes")]]
      genes <- genes %>% .[!. %in% not.in.gtf]
      too.close <- utilsFanc::gene.dist.mat(
        genes = genes, gtf = gtf, return.too.close = T, cutoff = cutoff)
      too.close <- c(too.close$x, too.close$y) %>% unique()
      genes <- genes %>% .[!. %in% too.close]
      return(genes)
    })
    names(res) <- paste0(items, ".genes")
    res$de.genes <- c(res$up.genes, res$down.genes)
    n.l <- lapply(res, length)
    names(n.l) <- paste0("n.", sub(".genes", "", names(res)))
    res <- c(res, n.l)
    copy <- c("padj.cutoff", "log2fc.cutoff")
    res <- c(res, s2b[[slot]][copy])
    f.l <- lapply(c("up", "down", "de"), function(x) {
      filtered <- s2b[[slot]][[paste0(x, ".genes")]] %>%
        .[! . %in% res[[paste0(x, ".genes")]]]
      return(filtered)
    })
    names(f.l) <- c("up", "down", "de")
    res$close.filter <- f.l
    res$close.filter$cutoff <- cutoff
    res$close.filter$not.in.gtf <- not.in.gtf
    s2b[[slot]] <- res
    return(s2b)
  })
  return(de)
}

de.fake.sc.from.bulk <- function(de.bulk, de.sc) {
  # for de objects made for bulk data (de.bulk), usually it's a list of 1.
  # the function will duplicate this to match the clusters from de.sc
  if (length(de.bulk) != 1) {
    stop("length(de.bulk) != 1")
  }
  de.bulk <- rep(de.bulk, length(de.sc))
  names(de.bulk) <- names(de.sc)
  for ( i in names(de.sc)) {
    de.bulk[[i]]$root.name <- i
  }
  rm(i)
  return(de.bulk)
}

# de.trend.cmp <- function(de1, de2, clusters.use = NULL, slot = "summary", out.file = NULL) {
#   # we take up and down genes from de1, and see their trends in de2
#   clusters <- intersect(names(de1), names(de2))
#   if (length(clusters) == 0)
#     stop("No shared clusters found between de1 and de2")
#   if (!is.null(clusters.use)) {
#     clusters <- clusters %>% .[.%in% clusters.use]
#     if (length(clusters) == 0) {
#       stop("No shared clusters found that are specified by clusters.use")
#     }
#   }
#
#   stat <- lapply(clusters, function(cluster) {
#     s1 <- de1[[cluster]]
#     if (!slot %in% names(s1)) {
#       stop("slot not found")
#     }
#     s2 <- de2[[cluster]]
#
#     stat <- lapply(c("up", "down"), function(type) {
#       genes <- s1[[slot]][[paste0(type, ".genes")]]
#       n <- length(genes)
#       genes.intsct <- genes %>% .[. %in% s2$bulkNorm$gene]
#       n.intsct <- length(genes.intsct)
#
#       trend.df1 <- s1$res.exp %>% dplyr::filter(gene %in% genes.intsct) %>%
#         dplyr::select(gene, log2FoldChange)
#       names(trend.df1) <- c("gene", "log2fc1")
#       trend.df2 <- s2$res.exp %>% dplyr::filter(gene %in% genes.intsct) %>%
#         dplyr::select(gene, log2FoldChange)
#       names(trend.df2) <- c("gene", "log2fc2")
#
#       trend.df <- trend.df1 %>% dplyr::left_join(trend.df2, by = "gene")
#
#       trend.df$is.same.trend <- (trend.df$log2fc1 * trend.df$log2fc2) > 0
#
#       n.trend <- sum(trend.df$is.same.trend)
#       pct.trend <- round(n.trend/n.intsct, digits = 3)
#
#       stat <- data.frame(cluster = cluster, type = type, n = n, n.intsct = n.intsct, n.trend = n.trend,
#                          pct.trend = pct.trend,
#                          genes.trend = paste0(trend.df$gene[trend.df$is.same.trend], collapse = ","),
#                          genes.oppo = paste0(trend.df$gene[!trend.df$is.same.trend], collapse = ","))
#       return(stat)
#     }) %>% do.call(rbind, .)
#     return(stat)
#   }) %>% do.call(rbind, .)
#
#   if (!is.null(out.file)) {
#     dir.create(dirname(out.file), showWarnings = F, recursive = T)
#     write.table(stat, out.file, quote = F, sep = "\t",
#                 row.names = F, col.names = T)
#   }
#   return(stat)
# }

de.trend.cmp <- function(de1, de2, clusters.1 = NULL, clusters.2 = NULL, slot = "summary", out.file = NULL) {
  # we take up and down genes from de1, and see their trends in de2
  # compared to the previous version: de1 and de2 don't have to have the same cluster names.
  if (is.null(clusters.1))
    clusters.1 <- names(de1)
  if (is.null(clusters.2))
    clusters.2 <- names(de2)

  utilsFanc::check.intersect(clusters.1, "clusters.1", names(de1), "names(de1)")
  utilsFanc::check.intersect(clusters.2, "clusters.2", names(de2), "names(de2)")


  stat <- lapply(clusters.1, function(cluster.1) {
    lapply(clusters.2, function(cluster.2) {
      s1 <- de1[[cluster.1]]
      if (!slot %in% names(s1)) {
        stop("slot not found")
      }
      s2 <- de2[[cluster.2]]

      stat <- lapply(c("up", "down"), function(type) {
        genes <- s1[[slot]][[paste0(type, ".genes")]]
        n <- length(genes)
        genes.intsct <- genes %>% .[. %in% s2$bulkNorm$gene]
        n.intsct <- length(genes.intsct)

        trend.df1 <- s1$res.exp %>% dplyr::filter(gene %in% genes.intsct) %>%
          dplyr::select(gene, log2FoldChange)
        names(trend.df1) <- c("gene", "log2fc1")
        trend.df2 <- s2$res.exp %>% dplyr::filter(gene %in% genes.intsct) %>%
          dplyr::select(gene, log2FoldChange)
        names(trend.df2) <- c("gene", "log2fc2")

        trend.df <- trend.df1 %>% dplyr::left_join(trend.df2, by = "gene")

        trend.df$is.same.trend <- (trend.df$log2fc1 * trend.df$log2fc2) > 0

        n.trend <- sum(trend.df$is.same.trend)
        pct.trend <- round(n.trend/n.intsct, digits = 3)

        stat <- data.frame(cluster.1 = cluster.1, cluster.2 = cluster.2,
                           type = type, n = n, n.intsct = n.intsct, n.trend = n.trend,
                           pct.trend = pct.trend,
                           genes.trend = paste0(trend.df$gene[trend.df$is.same.trend], collapse = ","),
                           genes.oppo = paste0(trend.df$gene[!trend.df$is.same.trend], collapse = ","))
        return(stat)
      }) %>% do.call(rbind, .)
      return(stat)
    }) %>% do.call(rbind, .) %>% return()
  }) %>% do.call(rbind, .)

  if (!is.null(out.file)) {
    dir.create(dirname(out.file), showWarnings = F, recursive = T)
    write.table(stat, out.file, quote = F, sep = "\t",
                row.names = F, col.names = T)
  }
  return(stat)
}

de.violin <- function(de, geneset.list, plot.out, quantile.filter = NULL, ...) {
  # each dot is a gene
  # x axis is sample, y axis is expression level
  # this function assesses the expression level of a gene set
  # geneset.list: a list of vectors. THe list must be named
  if (is.null(names(geneset.list))) {
    stop("geneset.list must be named")
  }
  pl <- lapply(de, function(s2b) {
    cluster <- s2b$root.name
    pl <- lapply(names(geneset.list), function(geneset.name) {
      genes <- geneset.list[[geneset.name]]
      exp <- s2b$bulkNorm %>% filter(gene %in% genes)
      if (nrow(exp) < length(genes)) {
        frac <- round(nrow(exp)/length(genes), digits = 3)
        warning(paste0(frac, " of all genes in the geneset '",
                       geneset.name, "' are found in cluster '", cluster, "'."))
        if (frac == 0) {
          stop()
        }
      }
      exp <- exp[, c("gene", colnames(s2b$bulk.mat))]
      if (!is.null(quantile.filter)) {
        expr <- exp
        expr$gene <- NULL
        expr <- as.matrix(expr)
        q <- quantile(rowMaxs(expr), quantile.filter)
        exp <- exp[rowMaxs(expr) <= q, ]
      }

      exp <- reshape2::melt(exp, id.vars = "gene", variable.name = "sample", value.name = "expr")

      p <- ggplot(exp, aes(x = sample, y = expr)) +
        geom_jitter(size = 0.1) +
        geom_violin(alpha = 0.5) +
        ggtitle(paste0(cluster, " ", geneset.name)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      return(p)
    })
    return(pl)
  }) %>% Reduce(c, .)
  trash <- wrap.plots.fanc(pl, n.split = length(geneset.list), plot.out = plot.out, ...)
  return()
}


de.chip.candidate <- function(de, da, clusters = NULL, genome,
                              summary.slot = "summary", res.slot = "res.exp",
                              nonDA.threshold = 0.1,
                              fg.peaks.include = NULL,
                              fg.peaks.exclude = NULL,
                              bg.peaks.exclude = NULL,
                              motif.beds = NULL, motif.beds.join.mode = "any",
                              chip.beds = NULL, chip.beds.join.mode = "any",
                              out.dir, root = NULL) {
  # written for the Baf155 KO project. So the foreground are down DARs/down DEGs.
  # if you have multiple motif bed files, "any" means that a motif instance is
  # considered as long as it appears in at least one of the motif bed files
  # "all" means that it needs to appear at all bed files.
  if (is.null(root)) root <- basename(out.dir)
  if (is.null(clusters)) {
    clusters <- intersect(names(de), names(da))
  }
  clusters.avail <- intersect(names(de), names(da))
  if (length(clusters.avail) == 0) {
    stop("length(clusters.avail) == 0")
  }

  utilsFanc::check.intersect(clusters, "clusters", clusters.avail, "clusters.avail")
  xlsx <- paste0(out.dir, "/", root, "_", paste0(clusters, collapse = "_") ,".xlsx")
  system(paste0("rm -rf ", xlsx))

  stats <- lapply(clusters,function(cluster) {
    s2b.e <- de[[cluster]]
    s2b.a <- da[[cluster]]

    if (is.null(s2b.a[[summary.slot]])) {
      stop("is.null(s2b.a[[summary.slot]])")
    }
    if (is.null(s2b.e[[summary.slot]])) {
      stop("is.null(s2b.e[[summary.slot]])")
    }
    stats <- list()
    stats$cluster <- cluster
    # annotate with closest gene
    gr.all <- s2b.a[[res.slot]] %>% as.data.frame()
    stats$n.total <- nrow(gr.all)
    if (is.null(gr.all$gene)) {
      gr.all <- gr.all %>% dplyr::mutate(., gene = rownames(.))
    }
    gr.all <- gr.all %>% utilsFanc::loci.2.df(loci.col.name = "gene", remove.loci.col = T, return.gr = T)
    gr.all$peak <- utilsFanc::gr.get.loci(gr.all)
    gr.all <- utilsFanc::gr.fast.annot(obj = gr.all, genome = genome, use.klraps = F)


    # filter for those overlapping chip
    if (!is.null(chip.beds)) {
      if (chip.beds.join.mode == "any") {
        chip.gr <- lapply(chip.beds, utilsFanc::import.bed.fanc, return.gr = T) %>%
          Reduce(c, .)
        gr.all <- gr.all %>% subsetByOverlaps(chip.gr)
      } else if (chip.beds.join.mode == "all") {
        chip.gr <- lapply(chip.beds, utilsFanc::import.bed.fanc, return.gr = T) %>%
          Reduce(plyranges::join_overlap_intersect, .)
        gr.all <- gr.all %>% subsetByOverlaps(chip.gr)
      } else {
        stop("only chip.beds.join.mode == any or all has been developed")
      }
    }
    stats$n.chip <- length(gr.all)
    # filter for those overlapping motif
    if (!is.null(motif.beds)) {
      if (motif.beds.join.mode == "any") {
        motif.gr <- lapply(motif.beds, utilsFanc::import.bed.fanc, return.gr = T) %>%
          Reduce(c, .)
        gr.all <- gr.all %>% subsetByOverlaps(motif.gr)
      } else {
        stop("only motif.beds.join.mode == any has been developed")
      }
    }
    stats$n.motif <- length(gr.all)

    # fg: filter by downDAR
    gr.fg <- gr.all[gr.all$peak %in% s2b.a[[summary.slot]][["down.genes"]]]
    stats$n.fg.downDAR <- length(gr.fg)
    # filter for down DARs that are assigned to downDEGs
    gr.fg <- gr.fg[gr.fg$SYMBOL %in% s2b.e[[summary.slot]][["down.genes"]]]
    stats$n.fg.downDEG <- length(gr.fg)

    if (!is.null(fg.peaks.include)) {
      gr.fg <- gr.fg %>% subsetByOverlaps(utilsFanc::loci.2.gr(fg.peaks.include), invert = F)
      stats$n.fg.peaks.include <- length(gr.fg)
    }
    if (!is.null(fg.peaks.exclude)) {
      gr.fg <- gr.fg %>% subsetByOverlaps(utilsFanc::loci.2.gr(fg.peaks.exclude), invert = T)
      stats$n.fg.peaks.exclude <- length(gr.fg)
    }

    stats$n.fg <- length(gr.fg)

    gr.fg <- gr.fg[utilsFanc::sort.by(gr.fg$SYMBOL, s2b.e[[summary.slot]][["down.genes"]], return.order = T)]

    utilsFanc::write.zip.fanc(gr.fg, out.file = paste0(out.dir, "/", root, "_", cluster, "_fg.bed"),
                              bed.shift = T)
    xlsx::write.xlsx2(utilsFanc::gr2df(gr.fg), xlsx,
                      sheetName = paste0(cluster, "_DAR"), append = T, col.names = T, row.names = T)

    ### background (those not changing)
    gr.bg <- gr.all[abs(gr.all$log2FoldChange) < nonDA.threshold]
    stats$n.bg.nonDA <- length(gr.bg)
    # assigned to upDEGs instead:
    gr.bg <- gr.bg[gr.bg$SYMBOL %in% s2b.e[[summary.slot]][["up.genes"]]]

    stats$n.bg.upDEG <- length(gr.bg)

    if (!is.null(bg.peaks.exclude)) {
      gr.bg <- gr.bg %>% subsetByOverlaps(utilsFanc::loci.2.gr(bg.peaks.exclude), invert = T)
      stats$n.bg.peaks.exclude <- length(gr.bg)
    }

    stats$n.bg <- length(gr.bg)

    utilsFanc::write.zip.fanc(gr.bg, out.file = paste0(out.dir, "/", root, "_", cluster, "_bg.bed"),
                              bed.shift = T)
    xlsx::write.xlsx2(utilsFanc::gr2df(gr.bg), xlsx,
                      sheetName = paste0(cluster, "_ctrl"), append = T, col.names = T, row.names = T)

    stats <- as.data.frame(stats)
  }) %>% do.call(rbind, .)
  write.table(stats, paste0(out.dir, "/", root, "_stats.tsv"), row.names = F, col.names = T,
              sep = "\t", quote = F)
  return(stats)
}

da.write.peakset <- function(da, out.dir, root = NULL, threads = 1) {
  if (is.null(root)) root <- basename(out.dir)
  utilsFanc::safelapply(da, function(a2b) {
    utilsFanc::write.zip.fanc(rownames(a2b$bulk.mat) %>% utilsFanc::loci.2.gr(),
                              out.file = paste0(out.dir, "/", root, "_", a2b$root.name, ".bed"))
  }, threads = threads)
}

da.motif.intersect <- function(da, motifs, motif.dir = "~/motifs/homer/sth/each_motif/",
                               out.dir, root = NULL, threads = 1) {
  if (is.null(root)) root <- basename(out.dir)
  if (is.null(names(motifs))) {
    stop("motifs must be named. set names(motifs)")
  }
  motifs.grl <- lapply(names(motifs), function(motif.name) {
    motif <- motifs[motif.name]
    motif.file <- paste0(motif.dir, "/", motif, ".bed")
    if (!file.exists(motif.file)) stop(paste0("motif file not found: ", motif.file))
    rtracklayer::import(motif.file)
  })
  names(motifs.grl) <- names(motifs)

  utilsFanc::safelapply(da, function(a2b) {
    peaks <- rownames(a2b$bulk.mat) %>% utilsFanc::loci.2.gr()
    lapply(names(motifs.grl), function(motif.name) {
      motif.gr <- motifs.grl[[motif.name]]
      motif.in.peaks <- subsetByOverlaps(motif.gr, peaks, ignore.strand = T)

      utilsFanc::write.zip.fanc(
        motif.in.peaks,
        out.file = paste0(out.dir, "/", root, "_", a2b$root.name, "_", motif.name, ".bed"))

    })

  }, threads = threads)

}

da.motif.pct.by.size <- function(da, motifs, motif.dir = "~/motifs/homer/sth/each_motif/",
                                 n.cut = 10, samples.bg,
                                 out.dir, root = NULL, threads = 1) {
  # divide all peaks based on size (quantile cut)
  # calculate fraction of peaks overlapping certain motifs for each quantile
  # samples.bg: the same idea as the samples.bg used to generate matched background:
  # it shouldn't be all of the samples. It could be all WT or all KO.
  if (is.null(root)) root <- basename(out.dir)
  if (is.null(names(motifs))) {
    stop("motifs must be named. set names(motifs)")
  }
  motifs.grl <- lapply(names(motifs), function(motif.name) {
    motif <- motifs[motif.name]
    motif.file <- paste0(motif.dir, "/", motif, ".bed")
    if (!file.exists(motif.file)) stop(paste0("motif file not found: ", motif.file))
    rtracklayer::import(motif.file)
  })
  names(motifs.grl) <- names(motifs)

  utilsFanc::safelapply(da, function(a2b) {
    peaks <- a2b$bulkNorm$gene %>% utilsFanc::loci.2.gr()
    peaks$size <- a2b$bulkNorm[, samples.bg, drop = F] %>% utilsFanc::pmean()
    peaks$cut <- gtools::quantcut(peaks$size, n.cut)
    peaks$cut.level <- peaks$cut %>% as.integer()

    df <- lapply(names(motifs.grl), function(motif.name) {
      motif.gr <- motifs.grl[[motif.name]]
      o <- findOverlaps(peaks, motif.gr, ignore.strand = T)
      o <- queryHits(o) %>% unique()
      tab.o <- table(peaks[o]$cut)
      tab.all <- table(peaks$cut)
      round(tab.o/tab.all, digits = 3)
    }) %>% Reduce(cbind, .)
    colnames(df) <- names(motifs.grl)
    df <- as.data.frame(df)
    df$cut <- rownames(df)
    file <- paste0(out.dir, "/", root, "_", a2b$root.name, ".tsv")
    dir.create(out.dir, showWarnings = F, recursive = T)
    write.table(df, file, quote = F, row.names = F, col.names = T, sep = "\t")
    return()
  }, threads = threads)

}

da.motif.in.dar.pct <- function(da, motifs, motif.dir = "~/motifs/homer/sth/each_motif/",
                                dar.types = c("up", "down"), summary.slot = "summary",
                                out.dir, root = NULL, threads = 1) {
  if (is.null(root)) root <- basename(out.dir)
  if (is.null(names(motifs))) {
    stop("motifs must be named. set names(motifs)")
  }
  motifs.grl <- lapply(names(motifs), function(motif.name) {
    motif <- motifs[motif.name]
    motif.file <- paste0(motif.dir, "/", motif, ".bed")
    if (!file.exists(motif.file)) stop(paste0("motif file not found: ", motif.file))
    rtracklayer::import(motif.file)
  })
  names(motifs.grl) <- names(motifs)

  stats <- utilsFanc::safelapply(da, function(a2b) {
    summ <- a2b[[summary.slot]]
    if (is.null(summ)) stop(paste0("Summary slot not fount: ", summary.slot))
    peaks <- a2b$bulkNorm$gene %>% utilsFanc::loci.2.gr()

    stats <- lapply(dar.types, function(dar.type) {
      dars <- summ[[paste0(dar.type, ".genes")]] %>% utilsFanc::loci.2.gr()
      pct.bg <- round(length(dars)/length(peaks), digits = 3)

      stats <- lapply(names(motifs.grl), function(motif.name) {
        motif.gr <- motifs.grl[[motif.name]]
        motif.in.peaks <- subsetByOverlaps(motif.gr, peaks, ignore.strand = T)
        motif.in.dars <- subsetByOverlaps(motif.in.peaks, dars, ignore.strand = T)
        stats <- data.frame(
          cluster = a2b$root.name,
          motif = motif.name,
          type = dar.type,
          n.in.peaks = length(motif.in.peaks),
          n.in.dars = length(motif.in.dars),
          pct.in.dars = round(length(motif.in.dars)/length(motif.in.peaks), digits = 3),
          bg.pct = pct.bg
        )
        stats$p.binom <- pbinom(q = stats$n.in.dars, size = stats$n.in.peaks,
                                prob = stats$bg.pct, lower.tail = F)
        return(stats)
      }) %>% do.call(rbind, .)
      return(stats)
    }) %>% do.call(rbind, .)
    stats$fdr <- p.adjust(stats$p.binom, method = "fdr")
    return(stats)
  }, threads = threads) %>% do.call(rbind, .)
  file <- paste0(out.dir, "/", root, "_motif_in_da.tsv")
  dir.create(out.dir, showWarnings = F, recursive = T)
  write.table(stats, file, quote = F, row.names = F, col.names = T, sep = "\t")

  df <- stats %>% split(., f = .$cluster) %>%
    lapply(function(stats) {
      stats %>% split(., f = .$type) %>% lapply(function(stats){
        df <- stats[, c("motif", "pct.in.dar")]
        df <- rbind(df, data.frame(motif = "Bg", pct.in.dar = stats$bg.pct[1]))
        df$cluster <- stats$cluster[1]
        df$type <- stats$type[1]
        return(df)
      }) %>% do.call(rbind, .) %>% return()
    }) %>% do.call(rbind, .) %>% return()

  # df %>% split(., f = .$type) %>% lapply(function(df){
  #   df <-
  # })

  invisible(stats)
}

de.simplify.seurat.clusters <- function(de) {
  names(de) <- names(de) %>% sub("seurat_clusters_", "sc", .)
  names(de) <- names(de) %>% sub("Clusters_C", "C", .)
  de <- lapply(de, function(s2b) {
    s2b$root.name <- s2b$root.name %>% sub("seurat_clusters_", "sc", .)
    return(s2b)
  })
}

de.filter.by.overlaps <- function(de1, de2, slot, slot.out, promoter.ext.1side = 500, invert = F) {
  # example: ~/hmtp/scAR/spbm2/DC4.2
  # returns de1, with the slot (usually summary) modified to retain only those
  # that are also in de2$s2b$slot. This will be written to slot.out
  clusters <- intersect(names(de1), names(de2))
  if (length(clusters) < 1) stop("length(clusters) < 1")

  de1 <- de1[clusters]
  de2 <- de2[clusters]
  de.slot.assert(de1, slot = slot)
  de.slot.assert(de2, slot = slot)

  de1 <- lapply(clusters, function(cluster) {
    s2b1 <- de1[[cluster]]
    s2b2 <- de2[[cluster]]
    summ.filtered <- lapply(c("up", "down"), function(direc) {
      if (de.type(de1) == "de" && de.type(de2) == "da") {
        pros <- s2b1[[slot]]$promoters[[direc]]
        if (is.null(pros)) {
          stop(paste0("promoter slot empty for ", names(des)[bG],
                      ", cluster ", cluster,", slot ", slot))
        }
        peaks <- s2b2[[slot]][[paste0(direc, ".genes")]]
        pros <- GenomicRanges::resize(pros, width = 2 * promoter.ext.1side,
                                      fix = "center")
        peaks <- peaks %>% utilsFanc::loci.2.df(loci.vec = ., remove.loci.col = F, return.gr = T)

        pros.in.peaks <- pros %>% subsetByOverlaps(peaks, invert = F, ignore.strand = T)
        genes <- pros.in.peaks$gene %>% unique()

        if (invert == T) {
          pros.not.in.peaks <- pros %>% subsetByOverlaps(peaks, invert = T, ignore.strand = T)
          # simply turing invert to TRUE is not good enough. Say one gene has 3 promoters. One of them intersects
          # peaks and the other 2 do not.
          genes.not.in.peaks <- pros.not.in.peaks$gene %>% unique()
          genes <- genes.not.in.peaks %>% .[!. %in% genes]
        }
        res <- list()
        res[[paste0("n.", direc)]] <- length(genes)
        res[[paste0(direc, ".genes")]] <- genes
        return(res)
      } else {
        stop("de1 must be de, de2 must be da. Others haven't been developed.")
      }
    }) %>% Reduce(c, .)
    summ.filtered$n.de <- summ.filtered$n.up + summ.filtered$n.down
    summ.filtered$de.genes <- c(summ.filtered$up.genes, summ.filtered$down.genes)
    s2b1[[slot.out]] <- summ.filtered
    return(s2b1)
  })
  names(de1) <- clusters
  return(de1)
}

gene.motif.mat.gen <- function(genes, gtf.gr, peakset, bPreComputedMotifOverlap = F,
                                  motifs, motif.dir = "~/motifs/homer/sth/each_motif",
                                  ext.1side = 50000, out.dir, root.name = NULL,
                               print.each.motif = F, threads.motif.overlap = 1) {
  # returns a gene x motif matrix. mat[1,1] would give the number of motif 1 found within
  # ext.1side of the tss of gene 1
  if (is.null(root.name)) root.name <- basename(out.dir)
  root.name <- paste0(root.name, "_ext", ext.1side)

  pros <- gene.promoter.map(genes = genes, gtf.gr = gtf)
  pros.bk <- pros
  pros <- pros + ext.1side
  genes.promoter.not.found <- genes[!genes %in% pros$gene]
  utilsFanc::write.zip.fanc(pros, paste0(out.dir, "/", root.name, "_pros.bed"),
                            bed.shift = T)
  n.genes.promoter.not.found <- length(genes.promoter.not.found)
  write(c(n.genes.promoter.not.found,
          round(n.genes.promoter.not.found/length(genes), digits = 3),
          genes.promoter.not.found),
        paste0(out.dir, "/", root.name, "_genes_promoter_not_found.txt"),
        sep = "\n")

  if (bPreComputedMotifOverlap) {
    utilsFanc::t.stat("Using pre-computed motif overlaps")
    if (is.character(peakset[1])) {
      peakset <- readRDS(peakset)
    }
    motifs <- colnames(mcols(peakset))
    names(motifs) <- motifs
  } else {
    utilsFanc::t.stat("Loading motif scanning results")
    if (is.null(names(motifs))) {
      # names(motifs) <- sub("\\(.+\\)$", "", motifs)
      # names(motifs) <- make.names(names(motifs))
      names(motifs) <- make.names(motifs)
      utilsFanc::check.dups(names(motifs), "names(motifs)")
    }

    motifs.grl <- utilsFanc::safelapply(names(motifs), function(motif.name) {
      motif <- motifs[motif.name]
      motif.file <- paste0(motif.dir, "/", motif, ".bed")
      if (!file.exists(motif.file)) stop(paste0("motif file not found: ", motif.file))
      rtracklayer::import(motif.file)
    }, threads = threads.motif.overlap)
    names(motifs.grl) <- names(motifs)

    utilsFanc::t.stat("Overlapping motifs with peaks")
    if (is.character(peakset)) peakset <- utilsFanc::loci.2.gr(peakset)
    utilsFanc::write.zip.fanc(peakset, paste0(out.dir, "/", root.name, "_peakset.bed"),
                              bed.shift = T)
    mcols(peakset) <- NULL
    strand(peakset) <- "*"

    mcols <- utilsFanc::safelapply(names(motifs.grl), function(motif.name) {
      motif.gr <- motifs.grl[[motif.name]]
      o <- findOverlaps(peakset, motif.gr, ignore.strand = T)
      bOverlap <- (1:length(peakset)) %in% unique(queryHits(o))
      motifs.in.peakset <- motif.gr[sort(unique(subjectHits(o)))]
      if (print.each.motif)
        utilsFanc::write.zip.fanc(motifs.in.peakset, paste0(out.dir, "/", root.name, "_motifs_in_peakset_", motif.name,".bed"))
      return(as.integer(bOverlap))
    }, threads = threads.motif.overlap)
    names(mcols) <- names(motifs.grl)
    mcols(peakset) <- mcols %>% as.data.frame()
    saveRDS(peakset, paste0(out.dir, "/", root.name, "_peakset_w_motifs.Rds"))
  }

  utilsFanc::t.stat("Summarizing by gene")
  peakset$peak <- utilsFanc::gr.get.loci(peakset)
  pros <- plyranges::join_overlap_left(pros, peakset)

  df <- mcols(pros) %>% as.data.frame()
  df <- df[, c("gene", "peak", names(motifs))]
  df <- df[!duplicated(df[, c("gene", "peak")]),]

  # there are 2 scenarios that could result in NAs in pros$peak:
  # scenario 1: no overlap between a gene and peaks
  # scenario 2: usually each gene has multiple promoters in gtf.gr. Some of them overlap peaks
  # while others do not.
  # we only remove scenario 2
  df <- df %>% dplyr::arrange(gene, peak) # For each gene, records with peak == NA will appear at the last
  df <- df[!(is.na(df$peak) & duplicated(df$gene)),]
  # now for those that really don't intersect anything, we change NA to zero:
  df[is.na(df$peak), names(motifs)] <- 0

  # remove duplicates
  df <- df[!duplicated(paste0(df$gene, "..", df$peak)),]

  write.table(df, paste0(out.dir, "/", root.name, "_genePeakMotif.tsv"),
              row.names = F, col.names = T, sep = "\t", quote = F)

  # gm: gene.motif

  gm.df <- df
  gm.df$motif.total.n <- rowSums(gm.df[, names(motifs)])
  gm.df$bTotalOver0 <- gm.df$motif.total.n > 0

  gm.df <- gm.df %>% dplyr::group_by(gene, bTotalOver0) %>%
    dplyr::summarise(
      across(-peak, sum),
      peakWmotif = paste0(peak[bTotalOver0], collapse = ","),
      peakWOmotif = paste0(peak[!bTotalOver0], collapse = ",")) %>%
    dplyr::ungroup() %>% as.data.frame()
  # peakWmotif: it could be peaks or "", when no peaks have motifs
  # peakWOmotif: it could be peaks, "", or NA. "": all peaks have motifs. NA: no peaks around promoter.
  gm.df$peakWOmotif[gm.df$peakWOmotif == "NA"] <- NA
  gm.df$peakWmotif[gm.df$peakWmotif == "" & is.na(gm.df$peakWOmotif)] <- NA

  gm.df <- gm.df %>% dplyr::group_by(gene) %>%
    dplyr::mutate(motif.total.n = NULL, bTotalOver0 = NULL) %>%
    dplyr::summarise(peakWmotif = paste0(peakWmotif, collapse = ""),
                     peakWOmotif = paste0(peakWOmotif, collapse = ""),
                     across(-c(peakWmotif, peakWOmotif), sum)) %>%
    dplyr::ungroup() %>% as.data.frame()

  write.table(gm.df, paste0(out.dir, "/", root.name, "_geneMotif.tsv"),
              row.names = F, col.names = T, sep = "\t", quote = F)

  gm.mat <- gm.df %>% dplyr::mutate(peak = NULL) %>% as.matrix()

  gm.df.melt <- df %>% dplyr::filter(!is.na(peak)) %>%
    reshape2::melt(id.vars = c("gene", "peak"), measure.vars = names(motifs),
                   value.name = "n", variable.name = "motif") %>%
    dplyr::filter(n > 0) %>%
    dplyr::group_by(gene, motif) %>%
    dplyr::summarise(peak = paste0(peak, collapse = ","), n = sum(n)) %>%
    dplyr::ungroup() %>% as.data.frame()
  gm.df.melt$motif <- as.character(gm.df.melt$motif)
  # maintains compatibility with the input requirement of de.motif.n()

  write.table(gm.df.melt, paste0(out.dir, "/", root.name, "_geneMotif_melt.tsv"),
              row.names = F, col.names = T, sep = "\t", quote = F)

  res <- list(gene.peak.motif.df = df, gene.motif.df = gm.df, gene.motif.df.melt = gm.df.melt,
              gene.motif.mat = gm.mat,
              n.genes.promoter.not.found = length(genes.promoter.not.found),
              genes.promoter.not.found = genes.promoter.not.found,
              promoters = pros.bk,
              motif.names = names(motifs),
              motifs.ori = motifs,
              params = list(ext.1side = ext.1side))
  saveRDS(res, paste0(out.dir, "/", root.name, "_result.Rds"))
  invisible(res)
}

de.annotate.genes <- function(de, slot = "annot", over.ride = T, gmt.files, gene.set.names) {
  # for example, annotate surface bound proteins.
  # each gmt contain multiple gene sets. Only those corresponding to gene.set.names will be returned.
  # gene.set.names can be a named vector to convert the gene set names in gmt files into a simpler name for the de object
  gmt <- lapply(gmt.files, function(gmt) {
    clusterProfiler::read.gmt(gmtfile = gmt)
      #        ont     gene
      # 1 chr10p11 SVIL-AS1
      # 2 chr10p11 PTCHD3P1
      # 3 chr10p11     SVIL
      # 4 chr10p11   MIR604
      # 5 chr10p11   MIR938
      # 6 chr10p11  CKS1BP2
  }) %>% do.call(rbind, .)
  rownames(gmt) <- NULL

  utilsFanc::check.intersect(gene.set.names, "gene.set.names", gmt$ont, "gmts")

  gmt <- gmt %>% dplyr::filter(ont %in% gene.set.names) %>%
    dplyr::mutate(gene = toupper(gene)) %>%
    dplyr::rename(gene.uc = gene)
  # gene.uc: gene upper case
  gmt <- unique(gmt)
  # gmt <- gmt %>% split(f = gmt$ont)

  if (is.null(names(gene.set.names)))
    names(gene.set.names) <- make.names(gene.set.names)
  name.map <- data.frame(ont = gene.set.names, gs = names(gene.set.names))

  de <- lapply(de, function(s2b) {
    genes <- rownames(s2b$bulk.mat)
    if (slot %in% names(s2b)) {
      if (over.ride) {
        print("slot already exists, overriding")
        s2b[[slot]] <- NULL
      } else {
        print("slot already exists, merging")
        if (!identical(sort(genes), sort(s2b[[slot]]$gene)))
          stop("!identical(sort(genes), sort(s2b[[slot]]$gene))")
        if (any(name.map$gs %in% colnames(s2b[[slot]]))) {
          stop("any(name.map$gs %in% colnames(s2b[[slot]]))")
        }
      }
    }

    annot <- data.frame(gene = genes, gene.uc = toupper(genes))
    annot <- dplyr::left_join(annot, gmt, by = "gene.uc")
    annot <- dplyr::left_join(annot, name.map, by = "ont")
    annot <- annot[, c("gene", "gs")]
    names(annot) <- c("gene", "annot")

    annot$value <- T
    annot <- unique(annot)
    annot <- reshape2::dcast(annot, formula = gene ~ annot, value.var = "value")

    annot <- annot[, colnames(annot) != "NA"]
    annot[is.na(annot)] <- F

    if (slot %in% names(s2b)) {
      annot <- dplyr::left_join(s2b[[slot]], annot, by = "gene")
    }
    s2b[[slot]] <- annot
    return(s2b)
  })
  return(de)
}

de.summary.filter.by.annot <- function(de, summary.slot = "summary", annot.slot = "annot", annots = NULL) {
  de <- lapply(de, function(s2b) {
    utilsFanc::check.intersect(c(annot.slot, summary.slot), "required slots",
                               names(s2b), "names(s2b)")
    if (is.null(annots)) {
      annots <- colnames(s2b[[annot.slot]]) %>% .[!.%in% "gene"]
    }

    utilsFanc::check.intersect(annots, "annots",
                               colnames(s2b[[annot.slot]]), "colnames(s2b[[annot.slot]])")
    genes.pass <- s2b[[annot.slot]]$gene[rowSums(s2b[[annot.slot]][, annots, drop = F]) > 0]
    summ <- s2b[[summary.slot]]
    trends <- c("up", "down", "de")
    for (i in trends) {
      summ[[paste0(i, ".genes")]] <- summ[[paste0(i, ".genes")]] %>% .[.%in% genes.pass]
      summ[[paste0("n.", i)]] <- length(summ[[paste0(i, ".genes")]])
    }
    s2b[[summary.slot]] <- summ
    return(s2b)
  })
  return(de)
}

de.print.summary <- function(de, summ.slot = "summary") {
  slots <- paste0("n.", c("up", "down", "de"))
  lapply(de, function(s2b) {
    print(paste0("cluster: ", s2b$root.name))
    lapply(slots, function(slot) {
      print(paste0(slot, ":"))
      print(s2b[[summ.slot]][[slot]])
    })
  })
  invisible(NULL)
}

de.deg.plot.bar <- function(de, slot = "summary", palette.fc = "R4.fc2",
                            down.name = "down", up.name = "up",
                            width = 3, height = 3, out.file) {
  de.slot.assert(de = de, slot = slot)
  df <- lapply(de, function(s2b) {
    summ <- s2b[[slot]]
    data.frame(cluster = s2b$root.name, n = c(summ$n.up,  -1 * summ$n.down),
               type = c(up.name, down.name))
  }) %>% do.call(rbind, .)

  dir.create(dirname(out.file), showWarnings = F, recursive = T)
  write.table(df, paste0(tools::file_path_sans_ext(out.file), ".tsv"),
              sep = "\t", col.names = T, row.names = F, quote = F)

  df$type <- factor(df$type, levels = c(up.name, down.name))
  df$cluster <- factor(df$cluster, levels = rev(names(de)))

  p <- ggplot(df, aes(x = n, y = cluster)) +
    geom_bar(aes(fill = type), stat = "identity")


  color.map <- utilsFanc::color.hue.fc(n = 2, palette = palette.fc)
  names(color.map) <- c(up.name, down.name)

  p <- p + scale_fill_manual(values = color.map)

  p <- p %>% utilsFanc::theme.fc.1(italic.x = F, rotate.x.45 = F)

  p <- p + theme(legend.direction = "vertical",
                 plot.margin = unit(c(0.04, 0.1, 0.01, 0), units = "in"))

  ggsave(out.file, p, width = width, height = height, device = cairo_pdf)
  invisible(p)
}

de.expr.corr <- function(de, comp.df,
                         color.by = NULL, shape.by = NULL,
                         plot.out) {
  # comp.df <- data.frame(x = "Etv2", y = c("Fhl2", "Col3a1"))
  stop("development aborted")
  pl <- lapply(de, function(s2b) {
    utilsFanc::check.intersect(unlist(comp.df), "unlist(comp.df)", s2b$bulkNorm$gene, "s2b$bulkNorm$gene")

    lapply(1:nrow(comp.df), function(i) {
      comp <- comp.df[i, ]
      df <- s2b$bulkNorm$gene %>% filter(gene %in% unlist(comp))
      df <- reshape2::melt(df, id.vars = "gene", variable.name = "sample", value.name = "expr")
    })

  })
  invisible(pl)
}
