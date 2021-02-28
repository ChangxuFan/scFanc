de.grid <- function(so, master.ident="seurat_clusters", assay, slot, master.ident.groups = NULL, group.by, 
                    ident.df, outlist=NULL, thread = 12, test.use = "wilcox", latent.vars = NULL, equal.cells,
                    logfc.threshold = 0.25, min.pct = 0.1) {
  # ident.df: 2 columns, must be ident.1 and ident.2
  Idents(so) <- master.ident
  if (is.null(master.ident.groups))
    master.ident.groups <- so@active.ident %>% as.character() %>% unique()
  grid.de.and.sum <- mclapply(master.ident.groups, function(i) {
    de.and.sum <- ident.df %>% split(., f= factor(.$ident.1, levels = .$ident.1)) %>% 
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
  }, mc.cores = thread)
  
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
  if (is.null(sub.idents))
    sub.idents <- so@meta.data %>% factor2character.fanc() %>% pull(!!master.ident) %>% unique()
  if (is.null(threads))
    threads <- length(sub.idents)
  mclapply(sub.idents, function(sub.ident) {
    
    # marker.genes <- marker.list[[paste0(master.ident, "_", sub.ident)]]
    
    p.list <- lapply(norms, function(norm) {
      # first we prepare a dataframe for scatter plots:
      mean.df <- so.2.bulk(so = so, assay = norm, slot = "data", sub.idents = sub.ident, 
                           ident = master.ident, group.by = "sample", take.mean = T, coldata.columns = NULL)$bulk.mat %>% as.data.frame()
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
            
            gene.list <- list(all.genes = NULL, 
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
  }, mc.cores = threads)

  return()
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


