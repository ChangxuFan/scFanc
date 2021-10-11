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

FindMarkers.fanc <- function (object, ident.1 = NULL, ident.2 = NULL, group.by = NULL, active.ident = "seurat_clusters",
                              subset.ident = NULL, assay = NULL, slot = "data", reduction = NULL, 
                              features = NULL, logfc.threshold = 0.25, test.use = "wilcox", 
                              min.pct = 0.1, min.diff.pct = -Inf, verbose = TRUE, only.pos = FALSE, 
                              max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL, 
                              min.cells.feature = 3, min.cells.group = 3, pseudocount.use = 1,
                              add.cdr = F, equal.cells = F,
                              ...) {
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
                               marker.list = NULL , n.col = NULL,  
                               genes.label = NULL,
                               sub.idents=NULL, max.quantile = NULL,
                               s2b.list = NULL, p.adj.cutoff, master.dir, root.name = NULL,
                               threads = 1,
                               use.plotly = F) {
  # different from grid.diagnostics: it tries to put all clusters into 1 file
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
        # if (sub.idents == "1")
        #   browser()
        de.genes <- grids[[norm]][[test]]$grid.de %>% filter(master.ident == sub.ident) %>% 
          .[.[, paste0("p.adj_", comp)] < p.adj.cutoff, "gene"] %>% .[!is.na(.)]
        gene.list <- list(# all.genes = NULL, 
          # marker.genes = marker.genes, 
          # Rpl.genes = Rpl.genes,
          de.genes = de.genes
          # de.markers = de.markers, 
          # de.Rpl = de.Rpl
        )
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
      trash <- wrap.plots.fanc(p.list, tooltip = "gene", plot.out = plot.out, n.col = n.col)
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

deseq2.xyplot <- function(pbl, pbl.soup = NULL, de.grid = NULL, de.grid.master.ident,
                          no.highlight = F, density.filter = NULL, device = "png", plotly.label.all = F,
                          add.label = F,
                          plot.dir, comp.vec, root.name = NULL, transformation = NULL, 
                          quantile.limit = NULL, padj.cutoff = 0.05, use.p = F, threads = 8, 
                          ...) {
  # pbl: pseudobulk list object generated by bulk.list()
  # de.grid: extend support to objects generated by de.grid()
  if (add.label == T) {
    warning("currently only support adding labels to all points. use with caution")
    label.var <- "gene"
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
  
  lapply(comp.vec, function(comp) {
    sample.y <- comp %>% sub(":.+$", "",.)
    sample.x <- comp %>% sub("^.+:", "",.)
    pl <- utilsFanc::safelapply(s2b.names, function(s2b.name) {
      # s2b <- pbl[[s2b.name]]
      if (no.highlight == F) {
        if (!is.null(de.grid)) {
          pbl[[s2b.name]]$res.exp <- de.grid$grid.sum %>% 
            filter(master.ident == sub(paste0(de.grid.master.ident, "_"), "", s2b.name)) %>% 
            rename(log2FoldChange = logFC, padj = p.adj, pvalue = p) %>% 
            mutate(log2FoldChange = log2(2.71828^log2FoldChange))
        }
        
        if (use.p == F) {
          de.df <- pbl[[s2b.name]]$res.exp %>% filter(padj < padj.cutoff) %>% select(gene, log2FoldChange)
        } else {
          de.df <- pbl[[s2b.name]]$res.exp %>% filter(pvalue < padj.cutoff) %>% select(gene, log2FoldChange)
        }
        
        de.df$log2FoldChange[de.df$log2FoldChange > 0] <- "up"
        de.df$log2FoldChange[de.df$log2FoldChange < 0] <- "down"
        
        highlight.var <- "gene"
      } else {
        highlight.var <- NULL
      }

      pbl.2 <- c(pbl[s2b.name], pbl.soup)
        
      pl <- lapply(names(pbl.2), function(s2b.name.2) {
        s2b <- pbl.2[[s2b.name.2]]
        df <- s2b$bulkNorm 
        colnames(df) <- sub("^bulkNorm_", "", colnames(df))
        if (no.highlight == F)
          df <- left_join(df, de.df)
        p <- xy.plot(df = df, x = sample.x, y = sample.y, transformation = transformation,
                     quantile.limit = quantile.limit, density.filter = density.filter,
                     highlight.var = highlight.var, highlight.values = de.df$gene, 
                     label.var = label.var, use.geom_text = T, 
                     plotly.var = "gene", highlight.color.var = "log2FoldChange", 
                     plotly.label.all = plotly.label.all,
                     show.highlight.color.var = F) +
          ggtitle(s2b.name.2)

        return(p)
      })
      return(pl)
    }, threads = threads) %>% Reduce(c, .)
    system(paste0("mkdir -p ", plot.dir))
    p.flag <- "padj"
    if (use.p == T)
      p.flag <- "p"
    plot.out <- paste0(c(root.name, paste0(comp, "_", p.flag, "_", padj.cutoff), device), collapse = ".") %>% 
      paste0(plot.dir, "/", .)
    wrap.plots.fanc(plot.list = pl, plot.out = plot.out, n.split = n.split, tooltip = "gene")
    saveRDS(pl, paste0(tools::file_path_sans_ext(plot.out), ".Rds"))
    return()
  })
  return()
}
deseq2.de.from.density <- function(pbl, comp.vec, transformation = NULL, 
                                   density.filter = 0.01, 
                                   plot.dir = NULL, plot.root.name = NULL,
                                   threads = 8) {
  utilsFanc::safelapply(pbl, function(s2b) {
    if (is.null(s2b$de.dens))
      s2b$de.dens <- list()
    for (comp in comp.vec) {
      sample.y <- comp %>% sub(":.+$", "",.)
      sample.x <- comp %>% sub("^.+:", "",.)
      # dens.df <- s2b$bulkNorm
      dens.df <- xy.plot(df = s2b$bulkNorm, x = sample.x, y = sample.y,
                         transformation = transformation, density.filter = density.filter, 
                         return.df = T)
      if (!is.null(plot.dir)) {
        if (!is.null(plot.root.name))
          root.name <- paste0(plot.dir, "/", plot.root.name,"..",s2b$root.name)
        else
          root.name <- paste0(plot.dir, "/", s2b$root.name)
        
        plot.out <- paste0(root.name, "_comp_", density.filter, ".png")
        trash <- xy.plot(df = s2b$bulkNorm, x = sample.x, y = sample.y,
                         transformation = transformation, outfile = plot.out,
                         highlight.var = "gene", highlight.values = dens.df$gene,
                         return.df = F)
        
      } 
      s2b$de.dens[[comp]] <- list()
      s2b$de.dens[[comp]]$de <- dens.df$gene
      s2b$de.dens[[comp]]$filter <- density.filter
    }
    
    return(s2b)
  })
}
deseq2.summary <- function(pbl, padj.cutoff = 0.05, force = F, stats.file = NULL, gene.out.dir = NULL) {
  pbl <- pbl[sapply(pbl, is.list)]
  s2b.names <- names(pbl)
  pbl <- lapply(s2b.names, function(s2b.name) {
    s2b <- pbl[[s2b.name]]
    if (is.null(s2b$summary) || force == T) {
      de.df <- s2b$res.exp %>% filter(padj < padj.cutoff) %>% arrange(padj) %>% select(gene, log2FoldChange)
      s2b$summary <- list(padj.cutoff = padj.cutoff, 
                          n.de = nrow(de.df),
                          n.up = nrow(de.df %>% filter(log2FoldChange > 0)),
                          n.down = nrow(de.df %>% filter(log2FoldChange < 0)),
                          de.genes = de.df$gene,
                          up.genes = de.df %>% filter(log2FoldChange > 0) %>% pull(gene),
                          down.genes = de.df %>% filter(log2FoldChange < 0) %>% pull(gene))
      if (!is.null(gene.out.dir)) {
        system(paste0("mkdir -p ", gene.out.dir))
        lapply(c("de.genes", "up.genes", "down.genes"), function(type) {
          write(s2b$summary[[type]], paste0(gene.out.dir, "/", s2b.name, "_", type, ".txt"), sep = "\n")
          return()
        })
      }
    }
    return(s2b)
  })
  names(pbl) <- s2b.names
  
  if (!is.null(stats.file)) {
    stats <- lapply(names(pbl), function(name) {
      stats <- pbl[[name]]$summary[c("n.de", "n.up", "n.down")] %>% as.data.frame()
      stats <- utilsFanc::add.column.fanc(df1 = stats, df2 = data.frame(name = name), pos = 1)
      return(stats)
    }) %>% Reduce(rbind, .)
    system(paste0("mkdir -p ", dirname(stats.file)))
    write.table(stats, stats.file, row.names = F, col.names = T, sep = "\t", quote = F)
    
  }
  return(pbl)
}

deseq2.plot.panel <- function(pbl, so, padj.cutoff, max.plots = 50, plot.type = c("de.genes", "up.genes", "down.genes"),
                              order = F, assay, reduction = "umap", cluster.ident =NULL,
                              sample.order = NULL,
                              plot.dir, root.name = NULL, threads = 6,
                              violin = F, violin.clusters = NULL) {
  if (!is.null(sample.order)) {
    so$sample <- factor(so$sample, levels = sample.order)
  }
  n.split <- so$sample %>% unique() %>% length()
  trash <- utilsFanc::safelapply(names(pbl), function(name) {
    s2b <- pbl[[name]]
    df <- s2b$res.exp %>% filter(padj < padj.cutoff) %>% 
      select(gene, padj, log2FoldChange) %>% arrange(padj)
    gene.list <- list()
    gene.list$de.genes <- df$gene
    gene.list$up.genes <- df %>% filter(log2FoldChange > 0) %>% pull(gene)
    gene.list$down.genes <- df %>% filter(log2FoldChange < 0) %>% pull(gene)
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
                                 plot.out = plot.out, raster = T,
                                 pt.size = 1.2, auto.adjust.raster.pt.size = F, page.limit = 20, max.quantile = 0.99,
                                 n.split = n.split, reduction = reduction, ident = cluster.ident)
        
      } else {
        print(genes)
        if (!is.null(root.name)) {
          name <- paste0(root.name, "..", name)
        }
        plot.out <- paste0(plot.dir, "/", name, "..", assay, "..", type, "..", order, ".vln.pdf")
        trash <- plot.panel.list(panel.list = genes, obj = so, order = order, assay = assay, 
                                 split.by = "sample", split.order = sample.order,
                                 ident = cluster.ident, cluster = violin.clusters,
                                 n.col = 6,
                                 plot.out = plot.out, 
                                 raster = T,
                                 pt.size = 0.5, auto.adjust.raster.pt.size = F, 
                                 page.limit = 6, max.quantile = 0.99, violin = T,add.median = F,
                                 sub.width = 15, sub.height = 3.5)
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

deseq2.write.results <- function(pb.m, work.dir) {
  system(paste0("mkdir -p ", work.dir, "/table/"))
  trash <- lapply(names(pb.m), function(name) {
    s2b <- pb.m[[name]]
    write.table(s2b$res.exp, paste0(work.dir, "/table/", name, ".tsv"), quote = F, row.names = F,
                col.names = T, sep = "\t")
  })
}

deseq2.pipe.ez <- function(soi, pb.m = NULL, samples, pbl.soup, assay, reduction = "UAP", 
                           cluster.ident, clusters, violin.clusters = NULL,
                           work.dir, plot.only = F, plot.xy = T, plot.panel = T, plot.violin = T,
                           padj.cutoff = 0.05, threads = 1) {
  if (plot.only == F) {
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
      deseq2.xyplot(pbl = pb.m, pbl.soup = pbl.soup, 
                    plot.dir = paste0(work.dir, "/plots_xy_w_soup/"), 
                    comp.vec = c("tumor_rep1:ctrl_rep1", "tumor_rep3:ctrl_rep3"), 
                    root.name = paste0("linear_0.99"), quantile.limit = 0.99, 
                    padj.cutoff = padj.cutoff, threads = threads)
    })
    
    try({
      deseq2.xyplot(pbl = pb.m, pbl.soup = pbl.soup, 
                    plot.dir = paste0(work.dir, "/plots_xy_w_soup/") , 
                    comp.vec = c("tumor_rep1:ctrl_rep1", "tumor_rep3:ctrl_rep3"), 
                    root.name = paste0("log2"), quantile.limit = NULL, 
                    transformation = log2,
                    padj.cutoff = padj.cutoff, threads = threads)
    })
  }
  
  if (plot.panel == T) {
    try({
      deseq2.plot.panel(pbl = pb.m, so = soi, padj.cutoff = padj.cutoff, plot.type = c("up.genes", "down.genes"), 
                        sample.order = samples, cluster.ident = cluster.ident,
                        max.plots = 50, order = F, assay = assay, reduction = reduction,
                        threads = 4, plot.dir = paste0(work.dir, "/panel/"))
    })
  }
  if (plot.violin == T) {
    try({
      deseq2.plot.panel(pbl = pb.m, so = soi, padj.cutoff = padj.cutoff, plot.type = c("up.genes", "down.genes"), 
                        sample.order = samples,
                        max.plots = 10, order = F, assay = assay, reduction = reduction,
                        plot.dir = paste0(work.dir, "/violin/"), violin = T, 
                        cluster.ident = cluster.ident, violin.clusters = violin.clusters, 
                        threads = threads.plotting)
    })
  }
  return()
}

deseq2.assess.core <- function(dds, plot.out, pca.ntop = 10000, 
                               pca.groupings = c("type", "rep"),
                               return.plot.list = F, blind = T,
                               add.hm = F, add.3d = F) {
  dir.create(dirname(plot.out), showWarnings = F, recursive = T)
  if (is.null(pca.ntop))
    pca.ntop <- 10000
  pl <- list()
  ntd <- normTransform(dds)
  rld <- rlog(dds,blind = blind, fitType = "local")
  pl$ntd <- mean.sd.plot(mat = assay(ntd))
  pl$rld <- mean.sd.plot(mat = assay(rld))
  
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
  if (flag == T)
    pc.list <- list(c(1,2))
  
  for (i in pc.list) {
    comp <- paste0("PC", i)
    p.name <- paste0(comp[1], ".", comp[2])
    if (length(pca.groupings) >= 2) {
      pl[[p.name]] <- 
        ggplot(pcaData, aes_string(comp[1], comp[2], color=pca.groupings[1],
                                                shape=pca.groupings[2])) +
        geom_point(size=3) +
        xlab(paste0(comp[1], ": ",percentVar[i[1]],"% variance")) +
        ylab(paste0(comp[2], ": ",percentVar[i[2]],"% variance")) + 
        theme_bw() +
        theme(aspect.ratio = 1) +
        ggtitle("ntop: " %>% paste0(pca.ntop))
    } 
    if (!is.na(pca.groupings[3])) {
      pl[[p.name]] <-  pl[[p.name]] + 
        geom_text(aes_string(label = pca.groupings[3]), color = "black")
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
  trash <- wrap.plots.fanc(pl, plot.out = plot.out)
  
  if (return.plot.list == T)
    return(pl)
  else
    return()
}

deseq2.homer <- function(dds, fg.vec, n.bg.each = 25, no.replace = F,
                         samples.for.bg, scatter.xy.df,
                      work.dir, write.log = T,
                      genome, size = "given", thread, denovo = F, other.params = "",
                      find.motifs.path = "/bar/cfan/software/homer/bin/findMotifsGenome.pl",
                      homer.path = "/bar/cfan/software/homer/bin/", run = T) {
  norm.mat <- DESeq2::counts(dds, normalized = T)
  bg.dir <- paste0(work.dir, "/bg/")
  prep.dir <- paste0(work.dir, "/prep/")
  homer.dir <- paste0(work.dir, "/homer/")
  system(paste0("mkdir -p ", bg.dir, " ", prep.dir, " ", homer.dir))
  
  bg.vec <- bg.gen.2(mat = norm.mat, fg.vec = fg.vec, n.bg.each = n.bg.each, samples.use = samples.for.bg,
                     no.replace = no.replace)
  write(bg.vec, paste0(bg.dir, "/bg.txt"), sep = "\n")
  bg.assess(mat = norm.mat, fg.vec = fg.vec, bg.vec = bg.vec, scatter.xy.df = scatter.xy.df,
            plot.out = paste0(bg.dir, "/bg_assess.png" ))
  bg.assess(mat = norm.mat, fg.vec = fg.vec, bg.vec = bg.vec, scatter.xy.df = scatter.xy.df,
            transformation = function(x) return(log2(x + 1)), 
            plot.out = paste0(bg.dir, "/bg_assess_log2.png" ))
  
  # prepare for homer:
  fg.bed <- paste0(prep.dir, "/fg.bed")
  trash <- utilsFanc::write.zip.fanc(df = utilsFanc::loci.2.df(loci.vec = fg.vec)[, 1:3], out.file = fg.bed)
  bg.bed <- paste0(prep.dir, "/bg.bed")
  trash <- utilsFanc::write.zip.fanc(df = utilsFanc::loci.2.df(loci.vec = bg.vec)[, 1:3], out.file = bg.bed)
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

a2bl.homer.pipe <- function(a2bl, work.dir,
                            slot, rank.key, desc = F, top.n.vec = c(100, 200, 300, 500),  # foreground choice
                            quantile.interval = NULL, quantile.samples,
                            n.bg.each = 25, no.replace = F, samples.for.bg, scatter.xy.df, # background choice
                            threads.a2bl = 1, threads.titrate = 4,  threads.homer = 1,
                            genome, size = "given", denovo = F, other.params = "", # homer params
                            find.motifs.path = "/bar/cfan/software/homer/bin/findMotifsGenome.pl", # homer params
                            homer.path = "/bar/cfan/software/homer/bin/", homer.run = T)  {# homer params
  # slot could be either res.exp or res.shrink.exp
  trash <- utilsFanc::safelapply(names(a2bl), function(a2b.name) {
    a2b <- a2bl[[a2b.name]]
    utilsFanc::safelapply(top.n.vec, function(top.n) {
      lapply(c("up", "down"), function(trend) {
        if (trend == "up")
          trend.f <- `>`
        else
          trend.f <- `<`
        
        work.dir <- paste0(work.dir, "/", a2b.name, "/", slot, "_",rank.key, "_top_", top.n,"_", trend, 
                           paste0(quantile.interval, collapse = ""),"/")
        prep.dir <- paste0(work.dir, "/prep/")
        system(paste0("mkdir -p ", work.dir, " ", prep.dir))
        fg.df <- a2b[[slot]]
        if (!is.null(quantile.interval)) {
          filter.vec <- fg.df[, quantile.samples] %>% rowMeans()
          limits <- quantile(filter.vec, quantile.interval)
          fg.df <- fg.df[filter.vec > limits[1] & filter.vec < limits[2] , ]
        }
        fg.df <- fg.df %>% filter(trend.f(log2FoldChange, 0)) %>% .[order(.[, rank.key]),]
        if (desc == T)
          fg.df <- fg.df[nrow(fg.df):1,]
        fg.df <- fg.df[1:top.n, ]
        write.table(fg.df, paste0(prep.dir, "/fg.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
        fg.vec <- fg.df$gene
        deseq2.homer(dds = a2b$dds, fg.vec = fg.vec, n.bg.each = n.bg.each, no.replace = no.replace,
                     samples.for.bg = samples.for.bg, 
                     scatter.xy.df = scatter.xy.df, work.dir = work.dir, write.log = T, genome = genome, 
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
