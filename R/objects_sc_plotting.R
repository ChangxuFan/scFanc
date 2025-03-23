feature.plot.comp <- function(marker.list, obj.list, plot.dir=NULL, n.marker, must.present.in=NULL) {

  p.list <- split(marker.list, f = marker.list$cluster) %>%
    lapply(function(x) {
      print(paste0("cluster: ",x$cluster[1]))
      x <- arrange(x,p_val_adj)
      if (!is.null(must.present.in)) {
        for (i in must.present.in)
          x <- x %>% filter(gene %in% rownames(obj.list[[i]]))
      }
      y <- x[1:n.marker, "gene"]
      y <- y[!is.na(y)]
      p <- lapply(seq_along(obj.list), function(i) {
        print(paste0("plotting from obj", i))
        pi <- FeaturePlot(obj.list[[i]], features = y, reduction = "umap", label=T, combine = F) %>%
          cowplot::plot_grid(plotlist = ., ncol = 1)
        return(pi)
      }) %>% cowplot::plot_grid(plotlist = ., ncol =length(obj.list))
      if (!is.null(plot.dir)) {
        system(paste0("mkdir -p ", plot.dir))
        ggsave(paste0(plot.dir, "/comp_cluster_", x$cluster[1], ".png"), p, units = "in", dpi = 200,
               width = length(obj.list)*5, height = n.marker*5, device = "png", limitsize = F)
      }

      return(p)
    })
}


wrap.plots.fanc <- function(plot.list, n.split=1, plot.out=NULL, n.col=NULL, page.limit = 100,
                            col.row.ratio=2, sub.width = 4, sub.height = 4, dpi = 300,
                            tooltip = NULL, pdf.by.row = T, pdf.use.cairo = T) {
  # col.row.ratio: in the final "big plot", # columns/# rows.
  # new behaviour: if there are too many plots (>page.limit), separate into multiple plots

  if (page.limit < n.split)
    stop("page.limit must be equal to or larger than n.split")
  n.plots <- length(plot.list)

  if (!is.null(plot.out) && grepl("(.pdf$)|(.html$)", plot.out)) {
    # if .pdf file is supplied, the multipage mode will be triggered
    system(paste0("mkdir -p ", dirname(plot.out)))
    page.limit <- floor(page.limit/n.split) * n.split
    n.sub <- min(page.limit, n.plots)
    if (is.null(n.col)) {
      n.col <- floor(col.row.ratio*(n.sub/col.row.ratio)^(0.5))
      n.col <- n.col - (n.col %% n.split)
      if (n.col == 0)
        n.col <- n.split
    }

    if (grepl("html$", plot.out)) {
      tag.list <- lapply(plot.list, function(p) {
        p.plotly <- suppressWarnings(plotly::ggplotly(p = p, tooltip = tooltip,
                                                      height = sub.height * 100,
                                                      width = sub.width * 100) %>%
                                       plotly::layout(showlegend = FALSE))
        widget <- plotly::as_widget(p.plotly)

        tag <- htmltools::div(widget, style = paste0("float:left;width:", floor(100/n.col), "%;"))
        return(tag)
      })
      doc <- do.call(what = htmltools::tagList, args = tag.list)
      new.dir <- sub(".html$", "", plot.out)
      dir.create(new.dir, showWarnings = F, recursive = T)
      plot.out <- paste0(new.dir, "/", basename(plot.out))
      htmltools::save_html(html = doc, file = plot.out)
    } else {
      n.row <- ceiling(n.sub/n.col)
      if (!pdf.use.cairo) {
        suppressMessages(suppressWarnings(extrafont::loadfonts()))
      }
      if (pdf.by.row) {
        layout_matrix <- suppressWarnings(matrix(1:length(plot.list), nrow = n.row, byrow = TRUE))
        layout_matrix[duplicated(as.vector(layout_matrix))] <- NA
        # Somehow this part will error out unless you have plotted something during this session
        p <- gridExtra::marrangeGrob(plot.list, layout_matrix = layout_matrix, 
                                     nrow = n.row, ncol = n.col,
                                     top = NULL)
        ggsave(plot.out, p, units = "in", device = ifelse(pdf.use.cairo, cairo_pdf, pdf),
               dpi = dpi, width = n.col*sub.width,
               height = n.row*sub.height, limitsize = F)
        
      } else {
        p <- gridExtra::marrangeGrob(plot.list, nrow = n.col, ncol = n.row, top = NULL)
        ggsave(plot.out, p, units = "in", device = ifelse(pdf.use.cairo, cairo_pdf, pdf),
               dpi = dpi, width = n.row*sub.width,
               height = n.col*sub.height, limitsize = F)
      }

      return(NULL)
    }


  }  else  {
    # first figure our how many pages are necessary:
    page.limit <- floor(page.limit/n.split) * n.split
    n.pages <- (n.plots/page.limit) %>% ceiling()
    split.vec <- ceiling((1:n.plots)/page.limit)

    page.list <- plot.list %>% split(f  = split.vec) %>%
      mapply(., unique(split.vec), SIMPLIFY = F, FUN = function(sub.list, id) {

        n.sub <- length(sub.list)
        if (is.null(n.col)) {
          n.col <- floor(col.row.ratio*(n.sub/col.row.ratio)^(0.5))
          n.col <- n.col - (n.col %% n.split)
          if (n.col == 0)
            n.col <- n.split
        }
        n.row <- ceiling(n.sub/n.col)
        p <- cowplot::plot_grid(plotlist = sub.list, ncol = n.col)
        if (!is.null(plot.out)) {
          system(paste0("mkdir -p ", dirname(plot.out)))
          if (n.plots > page.limit || id > 1) {
            plot.out <- sub("\\.([^\\.]+)$", paste0("_", id,".\\1"), plot.out)
          }
          ggsave(plot.out, p, units = "in", device = "png", dpi = dpi, width = n.col*sub.width,
                 height=n.row*sub.height, limitsize = F)
        }
        return(p)
      }
      )
    if (length(page.list) < 2)
      page.list <- page.list[[1]]
    return(page.list)
  }
}


plot.panel.list.bk <- function(panel.list, obj, cluster =NULL, order = T, assay,split.by=NULL,label = T, b2m = T,
                            n.split=1, ident = "seurat_clusters", plot.out = NULL, violin = F,
                            to.human = F, limits = NULL, return.list = F, pt.size = 0.05, raster = F,
                            threads = 1) {
  cells <- NULL
  if (!is.null(cluster)) {
    cells <- colnames(obj)[obj@meta.data$seurat_clusters %in% cluster]
  }
  # stop("currently not working for splits.")
  # n.split is only used for featurePlot.
  Idents(obj) <- ident
  DefaultAssay(obj) <- assay
  if (violin == T)
    n.split <- 1
  all.markers <- panel.list %>% unlist()

  avail.genes <- rownames(obj)
  avail.genes <- rbind(avail.genes, colnames(obj@meta.data))
  if (b2m == T)
    all.markers[!all.markers %in% avail.genes] <- "B2m"

  if (to.human == T) {
    all.markers <- toupper(all.markers)
  }
  n.plots <- length(all.markers)
  n.plots <- n.plots * n.split
  n.col <- floor(2*(n.plots/2)^(0.5))
  n.col <- n.col - (n.col %% n.split)
  if (n.col == 0)
    n.col <- n.split
  n.row <- ceiling(n.plots/n.col)
  in.per.plot <- 4
  if (is.null(split.by)) {
    if (violin == T) {
      p <- lapply(all.markers, function(x) {
        p.sub <- VlnPlot(object = obj, features = x, assay = assay) +
          stat_summary(fun=median, geom="point", size=3, color="red")
      })
    } else {
      p <- FeaturePlot(obj, all.markers, cells = cells, label = label,
                       order = order, pt.size = pt.size, combine = F, split.by = split.by, raster = raster)
      if (!is.null(limits)) {
        p <- lapply(p, function(p) {
          p <- p + scale_color_gradientn(colors = c("lightgrey", "blue"), limits = limits)
          return(p)
          })

      }

      p <- lapply(p, function(x) return(list(x)))
    }
  }

  else {
    p <- mclapply(all.markers, function(m) {

      if (violin == T) {
        p.sub <- VlnPlot(obj, m, split.by = split.by,
                         # cols = c("red", "orange", "blue", "green", "gray", "yellow"),
                         split.plot = F,)
        p.sub <- list(p.sub)
      } else {
        p.sub <- FeaturePlot(obj, m, label = label, cells = cells, order = order,
                             combine = F, split.by = split.by, pt.size = pt.size, raster = raster)
        if (is.null(limits)) {
          scale.max <- lapply(p.sub, function(psi) {
            g <- ggplot_build(psi)
            return(g$plot$scales$scales[[3]]$range$range[[2]])
          }) %>% unlist() %>% max()

          limits <- c(0, scale.max)
        }

        p.sub <- lapply(p.sub, function(psi) {
          psi <- suppressMessages(psi + ggtitle(m) +
            scale_color_gradientn(colors = c("lightgrey", "blue"), limits = limits) +
            theme(axis.title=element_blank(),
                  axis.text=element_blank(),
                  axis.ticks=element_blank(),
                  title = element_text(size = 8),
                  axis.title.y.right = element_text(color = "blue", size = 8),
                  legend.key.size = unit(0.3, "cm"),
                  legend.text = element_text(size = 5)))
          return(psi)
        })
      }
      return(p.sub)
    }, mc.cores = threads)
  }

 # print("miao")
  names(p) <- all.markers

  if (return.list == F) {

    p <- p %>% Reduce(c, .)
    p <- cowplot::plot_grid(plotlist = p, ncol = n.col)
    if (!is.null(plot.out)) {
      plot.out <- sub("\\.png", paste0("_", paste0(cluster, collapse = "_"),"_",order, ".png"), plot.out)
      system(paste0("mkdir -p ", dirname(plot.out)))
      ggsave(plot.out, p, units = "in", device = "png", dpi = 300, width = n.col*in.per.plot,
             height=n.row*in.per.plot, limitsize = F)
    }
  }
  return(p)
}


plot.panel.list <- function(panel.list, obj, cluster =NULL, sample = NULL, order = T, assay, 
                            split.by=NULL, split.order = NULL, use.split.as.title = F,
                            highlight.list = NULL, is.motif = F, motif.species, motif.map = ARCHETYPE.MAP,
                            label = T, label.size=4, b2m = T, italic.title = T,
                            hide.legend = F,
                            n.split=1, ident = "seurat_clusters", 
                            violin = F, violin.adjust = 1, violin.no.noise = T, violin.remove.x = F, violin.color.map = NULL,
                            ymax = NULL, ymin = 0,
                            plot.out = NULL, root.name = NULL,
                            to.human = F, limits = NULL, return.list = F, pt.size = 0.05, page.limit = 200, n.col = NULL,
                            raster = T, auto.adjust.raster.pt.size = F,
                            binarize.panel = F, binarize.items = NULL,
                            reduction = "umap", cells = NULL, max.quantile = NULL,add.median = T,
                            subset.ident = NULL, subset.idents,
                            polygon.df = NULL,
                            ao.embedding = "UMAP",
                            invisible = F,
                            publication = F, 
                            threads = 1, ...) {

  version <- obj@version %>% stringr::str_extract("\\d") %>% as.numeric()
  if (is.motif == T) {
    panel.list <- motif.2.genes(panel.list %>% unlist(), species = motif.species, map = motif.map) %>% unlist()
  }
  if (is.null(plot.out) && !is.null(root.name)) {
    plot.out <- paste0(root.name, "..", assay, "..", order, ".pdf")
  }
  if ("ArchRProject" %in% class(obj)) {
    obj <- fake.so.gen(ao = obj, ao.embedding = ao.embedding)
  }
  if (!is.null(highlight.list)) {
    panel.list <- names(highlight.list)
    obj <- add.highlight(so = obj, add.cells.list = highlight.list, is.Rds = F)
  }
  
  # filter functions:
  if (!is.null(cluster)) {
    cells <- colnames(obj)[obj@meta.data[, ident] %in% cluster]
  }
  
  if (!is.null(sample)) {
    sample.col <- "sample"
    if (is.null (obj@meta.data[, sample.col])) {
      sample.col <- "Sample"
      if (is.null (obj@meta.data[, sample.col])) {
        stop("is.null (obj@meta.data[, sample.col])")
      }
    }
    cells <- colnames(obj)[obj@meta.data[, sample.col] %in% sample]
  }
  
  if (!is.null(subset.ident)) {
    cells <- colnames(obj)[obj@meta.data[, subset.ident] %in% subset.idents]
  }
  

  if (!is.null(split.by)) {
    df <- obj@meta.data
    if (!is.null(sample)) {
      df <- df[df[, sample.col] %in% sample,]
    }
    n.split.internal <- (df[, split.by] %>% as.character() %>% unique() %>% length())
    if ( n.split.internal < 2)
      split.by <- NULL
    else if (n.split == 1) # the default n.split is 1. If the user specified a number other than 1, it will be respected regardless.
      n.split <- n.split.internal
  }
  if (!is.null(split.by) && !is.null(split.order)) {
    avail <- obj@meta.data[, split.by] %>% unique()
    if (any(! avail %in% split.order)) {
      obj <- obj[, obj@meta.data[, split.by] %in% split.order]
    }
    obj@meta.data[, split.by] <- obj@meta.data[, split.by] %>% factor(., split.order)
  }
  if (raster == T && auto.adjust.raster.pt.size == T) {
    pt.size = 1.2
  }
  if (!is.null(plot.out))
    dir.create(dirname(plot.out), showWarnings = F, recursive = T)
  if (is.null(max.quantile)) {
    max.quantile <-  1
  }
  # cells <- NULL
  # stop("currently not working for splits.")
  # n.split is only used for featurePlot.
  DefaultAssay(obj) <- assay
  Idents(obj) <- ident
  if (violin == T)
    n.split <- 1
  all.markers <- panel.list %>% unlist()
  if (binarize.panel == T) {
    obj@meta.data <- obj@meta.data %>% binarize.columns(cols = all.markers, include.items = binarize.items)
    all.markers.pattern <- paste0(paste0(all.markers, "\\.\\."), collapse = "|")
    all.markers <- colnames(obj@meta.data) %>% .[grepl(all.markers.pattern,.)] %>%
      .[gtools::mixedorder(sub(all.markers.pattern, "", .))]
  }
  avail.genes <- rownames(obj)
  avail.genes <- c(avail.genes, colnames(obj@meta.data))
  if (b2m == T) {
    B2m <- ifelse(to.human, "B2M", "B2m")
    Actb <- ifelse(to.human, "ACTB", "Actb")
    if (B2m %in% avail.genes) {
      all.markers[!all.markers %in% avail.genes] <- B2m
    } else if (Actb %in% avail.genes) {
      all.markers[!all.markers %in% avail.genes] <- Actb
    } else {
      if (any(!all.markers %in% avail.genes))
        stop("B2m and Actb are both missing from the matrix.")
    }
  }
    

  if (to.human == T) {
    all.markers <- toupper(all.markers)
  }

  #>>>>>>>>>>>>>>> plotting starts
  p <- mclapply(all.markers, function(m) {
    
    if (m %in% colnames(obj@meta.data)) {
      vec <- obj@meta.data[, m]
    }
    else {
      if (version >= 5) {
        vec <- obj@assays[[assay]]$data[m, ]
      } else {
        vec <- obj@assays[[assay]]@data[m, ]
      }
      # mat <- Seurat::GetAssayData(object = obj, slot = "data", assay = assay)
      # vec <- mat[m, ] %>% as.vector()
    }
    
    if (!(m %in% colnames(obj@meta.data) && !is.numeric(obj@meta.data[, m]))) {
      data.min <- min(vec)
      data.max <- max(vec)
    }
    
    if (is.null(limits)) {
      if (is.numeric(vec)) {
        if (is.null(ymax)) {
          ymax <- quantile(vec, max.quantile)
        }
        limits <- c(ymin, ymax)
      } else
        limits <- NULL
    }

    if (violin == T) {
      if (!is.null(cells))
        obj <- obj[, cells]
      p.sub <- VlnPlot.fanc(obj, m, split.by = split.by, adjust = violin.adjust,
                       # cols = c("red", "orange", "blue", "green", "gray", "yellow"),
                       split.plot = F, y.max = ymax, pt.size = pt.size, no.noise = violin.no.noise,
                       cols = violin.color.map)
      if (add.median == T)
        p.sub <- p.sub + stat_summary(fun=median, geom="point", size=3, color="red")
      
      if (publication) {
        p.sub <- p.sub + 
          theme(plot.title = element_text(size=6, margin = margin(b = -0.02, unit = "in"), 
                                          face = ifelse(italic.title, 'italic', 'plain'), family = "Arial"),
                text = element_text(size = 5,  color = "black", family = "Arial"), 
                axis.title= element_blank(),
                axis.line = element_line(size = 0.2),
                axis.ticks = element_line(size = 0.2),
                axis.text = element_text(size = 5,  color = "black", family = "Arial"),
                plot.margin = margin(t = 0.01, r = 0.12, b = 0, l =0, unit = "in"),
                legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"), 
                legend.box.margin = margin(t = -0.12, r = -0.1, b = 0, 
                                           l = ifelse(!is.null(split.by) && !use.split.as.title, -0.3, -0.15),
                                           unit = "in"), 
                legend.background = element_blank(), 
                legend.spacing = unit(0.02, "in"), legend.key.size = unit(0.05, "in"), 
                legend.box = "vertical", legend.title = element_text(size = 5, family = "Arial"),
                legend.text = element_text(size = 5, family = "Arial")
          )
        if (violin.remove.x) {
          p.sub <- p.sub + theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank())
        }
        if (hide.legend) {
          p.sub <- p.sub + theme(legend.position = "none")
        }
      }
      p.sub <- list(p.sub)
    } else {
      print(m)
      if (m %in% colnames(obj@meta.data) && !is.numeric(obj@meta.data[, m])) {
        Idents(obj) <- m
        if (!is.null(split.by))
          shuffle <- F
        else
          shuffle <- T
        # bug in seurat. If shuffle = TRUE, split would be incorrect
        p.sub <- DimPlot(object = obj, pt.size = pt.size, reduction = reduction,
                         split.by = split.by, label = T,#  label.size = 10,
                         combine = F, raster = raster, shuffle = shuffle, cells = cells)
        p.sub <- lapply(p.sub, function(pi) {
          pi <- pi + theme(text = element_text(size=8)) +
            theme(legend.key.size = unit(0.1, 'cm')) +
            theme(aspect.ratio = 1)
          return(pi)
        })

      } else {
        p.sub <- FeaturePlot(obj, m, label = label, label.size = label.size ,cells = cells, order = order, reduction = reduction,
                             combine = F, split.by = split.by, pt.size = pt.size, raster = F)

        # if (is.null(limits) && length(p.sub) > 1 && is.null(max.quantile)) {
        #   scale.max <- lapply(p.sub, function(psi) {
        #     g <- ggplot_build(psi)
        #     return(g$plot$scales$scales[[3]]$range$range[[2]])
        #   }) %>% unlist() %>% max()
        #
        #   limits <- c(0, scale.max)
        # }

        color.map <- c(data.min, limits, data.max)
        color.map <- c(data.min, max(data.min, min(limits)), max(limits), data.max)
        range <- data.max - data.min
        color.map[2] <- color.map[2] + 0.001 * range
        color.map[3] <- color.map[3] - 0.001 * range
        # There can't be ties. ggplot2 will generate a bug if there are ties.
        # Doing the above to get rid of ties

        p.sub <- lapply(p.sub, function(psi) {
          if (use.split.as.title) {
            tmp <- ggplot2::ggplot_build(psi)
            title <- tmp$layout$panel_params[[1]]$y.sec$name
            psi <- suppressWarnings({
              psi + scale_y_continuous()
            })
          } else {
            title <- m
          }
          psi <- suppressMessages({
            psi + ggtitle(title) +
              scale_color_gradientn(colors = c("lightgrey", "lightgrey", "blue", "blue"),
                                    values = scales::rescale(color.map),
                                    limits = c(data.min, data.max))
          })

          if (publication == T) {
            # if (publication.w.coordinates == T) {
            #   psi <- suppressMessages(psi + theme(
            #     axis.title.x = element_text(size = 8, face = "bold"),
            #     axis.title.y = element_text(size = 8, face = "bold"),
            #     axis.ticks=element_blank(),
            #     axis.text = element_blank(),
            #     axis.line = element_line(size = 0.1, colour = "grey50"),
            #     aspect.ratio = 1
            #   ))
            # } else {
            #
            # }
            # psi <- psi + theme_void()
            psi <- suppressMessages(psi + theme(
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks=element_blank(),
              axis.text = element_blank(),
              axis.line = element_blank(),
              aspect.ratio = 1
            ))
            psi <- psi +
              theme(plot.title = element_text(size=6, margin = margin(b = -0.02, unit = "in"),
                                              face = ifelse(italic.title, 'italic', 'plain'), family = "Arial"),
                    text = element_text(size = 5,  color = "black", family = "Arial"),
                    plot.margin = margin(t = 0, r = 0.12, b = 0, l =0, unit = "in"),
                    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"),
                    legend.box.margin = margin(t = -0.12, r = -0.1, b = 0,
                                               l = ifelse(!is.null(split.by) && !use.split.as.title, -0.3, -0.15),
                                               unit = "in"),
                    legend.background = element_blank(),
                    legend.spacing = unit(0.02, "in"), legend.key.size = unit(0.05, "in"),
                    legend.box = "vertical", legend.title = element_text(size = 5, family = "Arial"),
                    legend.text = element_text(size = 5, family = "Arial"),
                    axis.title.y.right = element_text(size = 6, family = "Arial", vjust = -1),
                    panel.border = element_blank()
              ) +
              guides(color = guide_colourbar(barwidth = 0.3, barheight = 2))
            if (hide.legend) {
              psi <- psi + theme(legend.position = "none")
            }

          } else {
            psi <- suppressMessages(psi +
                                      theme(axis.title=element_blank(),
                                            axis.text=element_blank(),
                                            axis.ticks=element_blank(),
                                            title = element_text(size = 8),
                                            axis.title.y.right = element_text(color = "blue", size = 8),
                                            legend.key.size = unit(0.3, "cm"),
                                            legend.text = element_text(size = 5),
                                            aspect.ratio = 1))

          }
          if (!is.null(polygon.df))
            psi <- psi + geom_polygon(data = polygon.df, mapping = aes(x = x, y = y), fill = NA, color = "red")
          if (raster) {
            psi <- ggrastr::rasterize(psi, dpi=300)
          }
          return(psi)
        })
      }

    }
    return(p.sub)
  }, mc.cores = threads)

  # print("miao")
  names(p) <- all.markers

  if (return.list == F) {
    
    p <- p %>% Reduce(c, .)
    p <- wrap.plots.fanc(plot.list = p, n.split = n.split, page.limit = page.limit, plot.out = plot.out,
                         n.col = n.col, ...)
  }
  if (invisible == T)
    invisible(p)
  else
    return(p)
}


plot.panel.list.2 <- function(panel.list, obj, cluster=NULL, assay,split.by=NULL,
 n.split=1, ident = "seurat_clusters", plot.out = NULL, violin = F, limits = NULL, ...) {
  lapply(c(T, F), function(i) {
    # plot.out.sub <- sub("\\.png", paste0("_", i, ".png"), plot.out)
    plot.out <- paste0(tools::file_path_sans_ext(plot.out), "_", i, ".", tools::file_ext(plot.out))
    trash <- plot.panel.list(panel.list = panel.list, obj = obj, order = i, cluster = cluster, assay = assay, split.by = split.by,
                             n.split=n.split, ident = ident, plot.out = plot.out, violin = violin, limits = limits, ...)
  })
  return()
}

plot.panel.list.m <- function(panel.list, obj, cluster=NULL, order, assay, split.by=NULL, ident="seurat_clusters", rearrange=T, plot.out = NULL, raster = F,
                              to.human = F, limits = NULL, return.list = F, return.list.melt = F, return.list.melt.melt=F, pt.size = 0.05,
                              threads.master = 1, threads.sub = 1, page.limit = 20) {
  # the violin part of the function is not tested at all
  # to make things uniform, all parameters should be passed as either a single object (which will be replicated), or a list
  ##never pass a list of objects in the form of vectors
  # the returned list will be a list with 2 inner levels: first layer on genes, second layer on objects.
  panel.list <- unlist(panel.list)
  if (!is.list(obj)) {
    obj <- list(obj)
  }
  n.obj <- length(obj)

  # first make sure all arguments are in the correct format
  args.check <- c("cluster", "order", "assay", "split.by", "ident", "to.human", "limits", "pt.size", "raster")


  params <- lapply(args.check, function(arg) {
    value <- get(arg)
    if (!is.list(value))
      value <- rep(list(value), n.obj)
    if (length(value) != n.obj) {
      stop(paste0("valueument ", arg, " has the length: ", length(value),
                  ", and does not agree with the number of objects passed, which is ", n.obj))
    }
    return(value)
  })

  names(params) <- args.check

  p.list <- mclapply(1:n.obj, function(i) {
    params.sub <- lapply(params, function(x) return(x[[i]]))
    params.sub <- c(obj = obj[[i]], panel.list = list(panel.list),
                    plot.out = NULL, violin = F, threads = threads.sub, return.list = T,
                    params.sub )
    p.list <- do.call(what = plot.panel.list, args = params.sub)
  }, mc.cores = threads.master, mc.cleanup = T)
  # rearrange p.list
  if (rearrange == T) {
    n.genes <- sapply(p.list, length)
    gene.names <- names(p.list[[1]])

    if (length(unique(n.genes)) != 1)
      stop("each so must return the same number of genes")
    n.genes <- unique(n.genes)
    p.list <- lapply(1:n.genes, function(i) {
      p.list.sub <- lapply(1:n.obj, function(j) {
        return(p.list[[j]][[i]])
      })
      return(p.list.sub)
    })

    names(p.list) <- gene.names

  }

  p.list.melt <- lapply(p.list, function(x) return(unlist(x, recursive = F)))
  n.split <- sapply(p.list.melt, length)
  if (length(unique(n.split)) != 1)
    stop("error when melting the list: every gene doesn't have the same number of plots associated with it")
  n.split <- unique(n.split)
  p.list.melt.melt <- p.list.melt %>% unlist(recursive = F)

  if (!is.null(plot.out)) {
    p <- wrap.plots.fanc(plot.list = p.list.melt.melt, n.split = n.split,
                             plot.out = plot.out, page.limit = page.limit)
  }

  if (return.list == T) {
    return(p.list)
  }

  if (return.list.melt == T) {
    return(p.list.melt)
  }
  if (return.list.melt.melt == T) {
    return(p.list.melt.melt)
  }
  return(p)
}



dimplot.3d <- function(so) {
  yourseuratobject <- so
  umap_1 <- yourseuratobject[["umap"]]@cell.embeddings[,1]
  umap_2 <- yourseuratobject[["umap"]]@cell.embeddings[,2]
  umap_3 <- yourseuratobject[["umap"]]@cell.embeddings[,3]

  plot.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))

  # Make a column of row name identities (these will be your cell/barcode names)
  plot.data$label <- paste(rownames(plot.data))

  # Plot your data, in this example my Seurat object had 21 clusters (0-20)
  plot_ly(data = plot.data,
          x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
          color = ~seurat_clusters,
          colors = c("lightseagreen",
                     "gray50",
                     "darkgreen",
                     "red4",
                     "red",
                     "turquoise4",
                     "black",
                     "yellow4",
                     "royalblue1",
                     "lightcyan3",
                     "peachpuff3",
                     "khaki3",
                     "gray20",
                     "orange2",
                     "royalblue4",
                     "yellow3",
                     "gray80",
                     "darkorchid1",
                     "lawngreen",
                     "plum2",
                     "darkmagenta"),
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 5, width=2), # controls size of points
          text=~label, #This is that extra column we made earlier for which we will use for cell ID
          hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names
}


pc.vln <- function(so, plot.dir, root.name = NULL, dims = NULL, assay, ident="sample", ...) {
  # dims should be a vector such as 1:30 or c(1,3,4)
  # ... passed to wrap.plots.fanc()
  Idents(so) <- ident

  if (is.null(root.name))
    root.name <- so@project.name
  dims.avail <- so@reductions$pca %>% colnames() %>% sub("PC_", "", .) %>% as.numeric()
  if (!is.null(dims))
    dims <- dims[dims %in% dims.avail]
  else
    dims <- dims.avail
  pl <- lapply(dims, function(x) {
    p <- VlnPlot(object = so, features = paste0("PC_", x), pt.size = 0.5, assay = assay) +
    stat_summary(fun=median, geom="point", size=3, color="red")
    return(p)
  })

  ps <- wrap.plots.fanc(plot.list = pl, col.row.ratio = 1, sub.width = 5, sub.height = 4, dpi = 150,
                        plot.out = paste0(plot.dir, "/", root.name, "_pc_vln_grid_",assay,".png"))
  return()
}

slice.top.features <- function(x, start, end) {
  if (is.list(x)) {
    x[[1]] <- x[[1]][start:end]
    x[[2]] <- x[[2]][start:end]
    return(x)
  } else stop("has to be a list")
}

hmtp.common.markers <- function(obj, plot.dir, root = "", threads = 1, to.human = F, assay = "RNA",
                                panels.to.plot = NULL, ...) {
  if (is.null(threads))
    threads <- 1
  panel.list <- readRDS("~/R_packages/scFanc/hmtp_panel.list.Rds")
  if (to.human == T) {
    t.f.upper <- function(x) {
      lapply(x, function(x) {
        if (!is.list(x))
          x <- toupper(x)
        else {
          x <- t.f.upper(x)
        }
      }) %>% return()
    }
    panel.list <- t.f.upper(panel.list)
  }
  utilsFanc::safelapply(seq_along(panel.list), function(i) {
    panel.name <- names(panel.list)[[i]]
    if (!is.null(panels.to.plot) && !panel.name %in% panels.to.plot)
      return()
    name <- paste0(names(panel.list)[[i]], "_panel.png")
    try(plot.panel.list.2(panel.list[[i]], obj = obj, assay = assay, to.human = to.human,
                          plot.out = paste0(plot.dir, "/", root, name), ...))
    return()
  }, threads = threads)

  # try(plot.panel.list.2(panel.list[["general"]], obj = obj, assay = "RNA",
  #                   plot.out = paste0(plot.dir, "/panel.png")))
  # try(plot.panel.list.2(panel.list[["cycle"]], obj, assay = "RNA",
  #                   plot.out = paste0(plot.dir, "/cycle_panel.png")))
  # try(plot.panel.list.2(panel.list[["lineage"]], obj, assay = "RNA",
  #                   plot.out = paste0(plot.dir, "/lineage_panel.png")))
  #
  # try(plot.panel.list.2(panel.list[["lsk"]], obj, assay = "RNA",
  #                   plot.out = paste0(plot.dir, "/lsk_panel.png")))
  # try(plot.panel.list.2(panel.list[["clp"]], obj, assay = "RNA",
  #                   plot.out = paste0(plot.dir, "/clp_panel.png")))
  #
  # try(plot.panel.list.2(panel.list[["GMP"]], obj, assay = "RNA",
  #                   plot.out = paste0(plot.dir, "/GMP_panel.png")))
  #
  # try(plot.panel.list.2(panel.list[["mono"]], obj, assay = "RNA",
  #                   plot.out = paste0(plot.dir, "/mono_panel.png")))
  #
  # try(plot.panel.list.2(panel.list[["dc"]], obj, assay = "RNA",
  #                   plot.out = paste0(plot.dir, "/dc_panel.png")))
  #
  # try(plot.panel.list.2(panel.list[["ery"]], obj, assay = "RNA",
  #                   plot.out = paste0(plot.dir, "/ery_panel.png")))
  # try(plot.panel.list.2(panel.list[["mega"]], obj, assay = "RNA",
  #                   plot.out = paste0(plot.dir, "/mega_panel.png")))
  # try(plot.panel.list.2(panel.list[["me"]], obj, assay = "RNA",
  #                            plot.out = paste0(plot.dir, "/me_panel.png")))
  return()
}

venn.core <- function(x.list, name.vec, fill=NULL) {
  if (is.null(fill))
    fill <- utilsFanc::gg_color_hue(length(name.vec))
  venn <- VennDiagram::venn.diagram(
    x = x.list,
    category.names = name.vec,
    filename = NULL,
    lwd = 1.5,
    # lty = 'blank',
    fill = fill,

    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",

    # Set names
    cat.cex = 2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
  )
  return(venn)
}

multi.venn <- function(x.list.list, name.vec.list=NULL, height.sub = 500, plot.titles.vec= NULL,
                       width.sub = 500, dpi = 150, use.upset = F,debug = F,
                       venn.fill = NULL, outlist = NULL, out.file, mainbar.y.max = NULL,
                       set_size.scale_max = NULL, order.by = c("freq", "degree")) {

  n.plots <- length(x.list.list)
  n.col <- sqrt(n.plots) %>% ceiling()
  n.row <- (n.plots/n.col) %>% ceiling()
  if (use.upset == T) {
    venns <- list()
    for (i in 1:length(x.list.list)) {
      if (!is.null(name.vec.list))
        names(x.list.list[[i]]) <- name.vec.list[[i]]
      x.list.list[[i]] <- lapply(x.list.list[[i]], function(x) {
        if (length(x) == 0)
          x <- "pseudo"
        return(x)
      })

      UpSetR::upset(UpSetR::fromList(x.list.list[[i]]), set_size.show = T ,keep.order = T, sets = names(x.list.list[[i]]),
                    mainbar.y.max = mainbar.y.max,
                    set_size.scale_max = set_size.scale_max, order.by = order.by) %>% print()
      grid::grid.edit('arrange', name=paste0("list ", i))
      grid::grid.text(ifelse(is.null(plot.titles.vec), i, plot.titles.vec[i]),
                      x = 0.1, y=0.95, gp=gpar(fontsize=10))
      vp <- grid::grid.grab()
      venns[[i]] <- vp
    }
  } else {
    venns <- lapply(seq_along(x.list.list),
                    function(i) venn.core(x.list = x.list.list[i], name.vec = name.vec.list[i], fill = venn.fill))

  }

  if (!is.null(outlist)) {
    out.file <- utilsFanc::plot.name.construct(outlist, sub.name = ifelse(use.upset == T,"_mupset.png","_mvenn.png"))
  }
  png(out.file,
      width = n.col*width.sub,
      height = n.row*height.sub, res = dpi)

  try(gridExtra::grid.arrange(grobs = venns,ncol = n.col))
  dev.off()
  return()
}

vln.core.fanc <- function(so, features,group.by = NULL, group.order=NULL, ident, ident.order = NULL,
                          pt.size = 1, assay = NULL, split.by = NULL, adjust = 1, y.max = NULL, same.y.lims = F,
                          slot = "data", split.plot = F, combine = F, outlist=NULL, ...) {

  Idents(so) <- ident

  if (!is.null(ident.order)) {
#    so <- so[, so@meta.data[,ident] %in% ident.order]
    so@meta.data[,ident] <- so@meta.data[,ident] %>% factor(levels = ident.order)
  }

  if (!is.null(group.order)) {
    so@meta.data[,group.by] <- so@meta.data[,group.by] %>% factor(levels = group.order)
  }

  p <- VlnPlot(object = so, pt.size = pt.size, assay = assay, split.by = split.by,
               group.by = group.by, features = features,
               idents = ident.order,
               adjust = adjust, y.max = y.max, same.y.lims = same.y.lims,
               slot = slot, split.plot = split.plot, combine = combine, ...)
  p <- lapply(p, function(x) return(
    x +  stat_summary(fun=median, geom="point", size=3, color="red")
  ))

  if (length(features) >= 3) {
    features <- features[1:2] %>% c("etc")
  }
  root.name.internal <- paste0("vln_", paste0(features, collapse = ":"))
  sub.name <- paste0("_",ident, "_", paste0(ident.order, collapse = ":"), "_", assay, "_", slot, "_",
                               group.by, "_", split.by, ".png")


  if (!is.null(outlist))
    wrap.plots.fanc(p, plot.out = utilsFanc::plot.name.construct(outlist = outlist, root.name.internal = root.name.internal,
                                                                         sub.name = sub.name))
  return(p)
}


rank.plot <- function(df, vars, rank.method = NULL, transformation = NULL,
                      x.limit = c(0, NA), y.limit = c(0, NA), quantile.limit.y = NULL,
                      outfile = NULL, return.df = F,
                      label.var = NULL, labels = NULL, label.size = 4,
                      pt.size = 0.5,
                      add.h.line = NULL, add.v.line = NULL, title = NULL,
                      height = 4, width = 4,
                      publication = F) {
  if (!is.data.frame(df))
    df <- as.data.frame(df)
  df.core <- df[, ! colnames(df) %in% vars, drop = F]
  df.melt <- lapply(vars, function(var) {
    df.melt <- df.core
    df.melt$value <- df[, var]
    if (!is.null(transformation))
      df.melt$value <- transformation(df.melt$value)
    if (!is.null(rank.method))
      df.melt$rank <- rank(df.melt$value , ties.method = rank.method)
    else
      df.melt$rank <- rank(df.melt$value)

    df.melt$type <- var
    return(df.melt)
  }) %>% Reduce(rbind, .)
  if (!is.null(quantile.limit.y)) {
    df.melt <- df.melt %>% filter(value < quantile(df.melt$value, quantile.limit.y))
  }
  p <- ggplot(df.melt, aes(x = rank, y = value, color = type)) +
    geom_point(size = pt.size, alpha = 0.3) +
    theme_classic()
  if (!is.null(x.limit)) 
    p <- p + xlim(x.limit)
  if (!is.null(y.limit))
    p <- p + ylim(y.limit)
  
  if (!is.null(label.var)) {
    label.data <- df.melt
    if (!is.null(labels))
      label.data <- label.data %>% .[.[,label.var] %in% labels,]

    p <- p + ggrepel::geom_text_repel(data = label.data, mapping = aes_string(x = "rank", y = "value",
                                                               color = "type", label = label.var),
                                      size = label.size)
  }

  if (!is.null(add.h.line)) {
    if (!is.null(transformation))
      add.h.line <- transformation(add.h.line)
    p <- p + geom_hline(yintercept = add.h.line)
  }

  if (!is.null(add.v.line)) {
    p <- p + geom_vline(xintercept = add.v.line)
  }
  if (!is.null(title))
    p <- p + ggtitle(title)

  p <- p + theme(aspect.ratio = 1)
  
  if (publication) {
    p <- ggrastr::rasterize(p, dpi=300)
    p <- p %>% utilsFanc::theme.fc.1(italic.x = F)
    p <- p + theme(
      legend.position = "none"
    )
  }
  
  if (!is.null(outfile)) {
    trash <- wrap.plots.fanc(plot.list = list(p), plot.out = outfile, sub.height = height, sub.width = width)
  }
  if (return.df == T)
    return(df.melt)
  return(p)
}


density.plot <- function(df, vars, transformation = NULL, x.limit =NULL, y.limit = NULL,
                         outfile = NULL, return.df = F) {
  df.melt <- reshape2::melt(df, measure.vars = vars)

  if (!is.null(transformation))
    df.melt$value <- transformation(df.melt$value)

  p <- ggplot(df.melt, aes(x = value, fill = variable)) +
    geom_density(alpha = 0.5) +
    theme_classic()

  if (!is.null(x.limit))
    p <- p + xlim(x.limit)
  if (!is.null(y.limit))
    p <- p + ylim(y.limit)

  if (!is.null(outfile)) {
    trash <- wrap.plots.fanc(plot.list = list(p), plot.out = outfile)
  }
  if (return.df == T)
    return(df.melt)
  return(p)
}

xy.plot <- function(df, x, y, is.regex = F, collapse.fun = utilsFanc::pmean, x.lab = NULL, y.lab = NULL,
                    transformation = NULL, transformation.x = NULL, transformation.y = NULL,
                    x.limit = NULL, y.limit = NULL, quantile.limit = NULL,
                    theme.text.size = NULL, pt.size = 0.05, pt.shape = 19, color.var = NULL, pt.color = "grey50",
                    density.filter = NULL, density.filter.2sided = F,
                    color.density = F, show.color.var = F,
                    highlight.var=NULL, highlight.values=NULL, highlight.color.var = NULL, color.map = NULL,
                    show.highlight.color.var = F, highlight.ptsize = NULL,
                    highlight.only = F,
                    label.var = NULL, label.values = NULL, italic.label = F, 
                    use.repel = F, hjust = 0, vjust = 0, repel.direction = "both",
                    nudge_x = 0, nudge_y = 0, text.size = 5, text.color = "midnightblue", use.geom_text = F,
                    plotly.var = NULL, plotly.label.all = F,
                    add.corr=F, add.highlight.corr = F,
                    add.smooth = F, smooth.method = NULL, add.highlight.smooth = F,
                    add.abline = T, raster = F, raster.dpi = 200,
                    add.v.line = NULL, add.h.line = NULL,
                    title = NULL,
                    outfile=NULL, return.df =F ) {
  if ("GRanges" %in% class(df) )
    df <- df %>% `names<-`(NULL) %>% as.data.frame()
  if ("RangedSummarizedExperiment" %in% class(df))
    df <- utilsFanc::se.simplify(se = df)
  if (!is.data.frame(df)) {
    df <- as.data.frame(df)
  }
  if (!suppressWarnings(is.na(as.numeric(x)))) {
    df <- utilsFanc::change.name.fanc(df = df, cols.from = x, cols.to = paste0("X", x))
    x <- paste0("X", x)
  }

  if (!suppressWarnings(is.na(as.numeric(y)))) {
    df <- utilsFanc::change.name.fanc(df = df, cols.from = y, cols.to = paste0("X", y))
    y <- paste0("X", y)
  }

  if (is.null(transformation.x))
    transformation.x <- transformation
  if (is.null(transformation.y))
    transformation.y <- transformation
  x.ori <- x
  y.ori <- y
  if (is.regex == T) {
    x <- colnames(df) %>% .[grepl(x,.)]
    y <- colnames(df) %>% .[grepl(y,.)]
  }

  if (length(x) > 1) {
    if (!is.null(x.lab))
      new.x <- x.lab
    else
      new.x <- "x.collpase"
    df[, new.x] <- df[, x]  %>% collapse.fun()
    x <- new.x
  }
  if (length(y) > 1) {
    if (!is.null(y.lab))
      new.y <- y.lab
    else
      new.y <- "y.collpase"
    df[, new.y] <- df[, y]  %>% collapse.fun()
    y <- new.y
  }


  if (!is.null(transformation.x))
    df[, x ] <- transformation.x(df[, x])
  if (!is.null(transformation.y))
    df[, y ] <- transformation.y(df[, y])

  if (!is.null(density.filter)) {
    df$density <- utilsFanc::getDensity(x = df[, x], y = df[, y])$density
    if (density.filter.2sided == T) {
      df <- df %>%  mutate(..sign = if_else((df[, y] - df[, x]) > 0, "+", "-")) %>%
        group_by(..sign) %>% filter(density < quantile(density, density.filter)) %>%
        ungroup() %>% mutate(..sign = NULL) %>% as.data.frame()
    } else {
      df <- df[df$density < quantile(df$density, density.filter),]
    }
  }

  if (is.null(color.var)) {
    if (color.density == T) {
      color.var <- "density"
      df$density <- utilsFanc::getDensity(x = df[, x], y = df[, y])$density
    } else {
      df$.color <- pt.color
      color.var <- ".color"
    }
  }

  if (!is.null(quantile.limit)) {
    x.limit <- c(0, quantile(df[, x], quantile.limit))
    y.limit <- c(0, quantile(df[, y], quantile.limit))
  }

  if (highlight.only == T) {
    df <- df %>% .[.[, highlight.var] %in% highlight.values, ]
  }

  if (!is.null(plotly.var) && plotly.label.all == T)
    p <- ggplot(df, aes_string(x = x, y = y, color = color.var, text = plotly.var))
  else
    p <- ggplot(df, aes_string(x = x, y = y, color = color.var))

  # p <- ggplot(df, aes_string(x = x, y = y, color = color.var))

  if (raster == T) {
    p <- p + ggrastr::rasterise(geom_point(
        size = pt.size, alpha = 0.5, show.legend=show.color.var, shape = pt.shape,
        ), 
      dpi = raster.dpi)
  } else {
    p <- p + geom_point(size = pt.size, alpha = 0.5, show.legend=show.color.var, shape = pt.shape)
  }

  if (color.density == T) {
    p <- p + scale_color_gradient(low = "grey70", high = "skyblue1")
  } else if(color.var == ".color") {
    p <- p + scale_color_manual(values = pt.color)
  }


  if (add.abline ==T ) {
    p <- p + geom_abline(slope = 1, intercept = 0, color = "blue", size = 0.1)
    # updated 2023-03-17: size used to be 2*pt.size
  }
  p <- p +
    theme_classic()



  if (!is.null(x.limit))
    p <- p + xlim(x.limit)
  if (!is.null(y.limit))
    p <- p + ylim(y.limit)

  if (add.corr == T) {
    p <- p + ggpubr::stat_cor(show.legend=show.color.var)
  }

  if (add.smooth == T) {
    p <- p + geom_smooth(method = smooth.method,show.legend=show.color.var)
  }
  if (!is.null(highlight.var)) {
    if (is.null(highlight.ptsize))
      highlight.ptsize <- pt.size
    df.h <- df
    if (!is.null(highlight.values))
      df.h <- df.h %>% .[.[, highlight.var] %in% highlight.values, ]
    else
      df.h <- df.h %>% .[!is.na(.[, highlight.var]),]
    if (is.null(highlight.color.var)) {
      df.h$hl.color <- "red"
      highlight.color.var <- "hl.color"
    }
    p <- p + geom_point(data = df.h, aes_string(x = x, y = y, fill = highlight.color.var,text = plotly.var),
                        size = 10 * highlight.ptsize, stroke = 0, shape = 21,
                        show.legend=show.highlight.color.var)
    if (!is.null(color.map)) {
      p <- p + scale_fill_manual(values = color.map)
    }
    
    if (raster == T) {
      p <- ggrastr::rasterize(p, dpi=300)
    }
    # fill is used to control color instead of "color".
    # source: https://stackoverflow.com/questions/34398418/geom-point-borders-in-ggplot2-2-0-0
    if (add.highlight.corr == T) {
      p <- p + ggpubr::stat_cor(show.legend=F, data =  df.h)
    }
    if (add.highlight.smooth == T) {
      p <- p + geom_smooth(data = df.h, method = smooth.method,show.legend=F, color = "green")
    }
  }
  if (!is.null(label.var)) {
    label.df <- df
    if (!is.null(label.values)) {
      label.df <- label.df %>% .[.[, label.var] %in% label.values, ]
    } else {
      label.df <- label.df %>% .[!is.na(.[, label.var]), ]
    }

    if (italic.label == T) {
      label.df[, label.var] <- label.df[, label.var] %>% paste0("italic('", ., "')")
    }
    label.df <- label.df %>% dplyr::mutate(hjust = abs(hjust), vjust = abs(vjust))
    label.df$hjust[label.df[, y] < label.df[, x]] <- -1 * label.df$hjust[label.df[, y] < label.df[, x]]
    label.df$vjust[label.df[, y] > label.df[, x]] <- -1 * label.df$vjust[label.df[, y] > label.df[, x]]
    if (use.repel == T) {
      p <- p + ggrepel::geom_text_repel(data = label.df,
                                        mapping = aes_string(x = x, y = y, label = label.var,
                                                             hjust = "hjust", vjust = "vjust"),
                                        size = text.size, color = text.color, parse = italic.label,
                                        segment.size = 0.2, direction = repel.direction,
                                        max.overlaps = 100,
                                        min.segment.length = 0)
    } else if (use.geom_text == F) {
      p <- p + geom_label(data = label.df, mapping = aes_string(x = x, y=y, label = label.var,
                                                                hjust = "hjust", vjust = "vjust"),
                          nudge_x = nudge_x, nudge_y = nudge_y, size = text.size, color = text.color,
                          fill = "white", alpha = 0.5, fontface = "bold", label.size = NA,  parse = italic.label)
    }  else {
      p <- p + geom_text(data = label.df, mapping = aes_string(x = x, y=y, label = label.var,
                                                               hjust = "hjust", vjust = "vjust"),
                         nudge_x = nudge_x, nudge_y = nudge_y, size = text.size, color = text.color,
                         fill = "white", alpha = 0.5, fontface = "bold", label.size = NA,  parse = italic.label)
    }

  }
  p <- p + theme(aspect.ratio = 1)

  if (show.highlight.color.var == F)
    p <- p + theme(legend.title = element_blank())

  if (!is.null(x.lab))
    p <- p + xlab(x.lab)
  if (!is.null(y.lab))
    p <- p + ylab(y.lab)

  if (!is.null(add.v.line)) {
    p <- p + geom_vline(xintercept = add.v.line)
  }
  if (!is.null(add.h.line)) {
    p <- p + geom_hline(yintercept = add.h.line)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (!is.null(theme.text.size)) {
    p <- p + theme(text = element_text(size = theme.text.size))
  }

  if (!is.null(outfile)) {
    trash <- wrap.plots.fanc(plot.list = list(p), plot.out = outfile, tooltip = plotly.var)
  }
  if (return.df == T)
    return(df)
  return(p)
}

archr.quick.plot <- function(ao, meta, embedding = "UMAP", plot.out = NULL, max.quantile = NULL,
                             pt.size = 0.05, simple.coloring = T, brewer.pallete = "Spectral", ...) {
  if (is.null(max.quantile))
    max.quantile <- 1
  umap <- ao@embeddings[[embedding]]$df %>% `colnames<-`(c("UMAP1", "UMAP2")) %>%
    mutate(., cell = rownames(.))
  meta.df <- getCellColData(ao, select = meta) %>% as.data.frame() %>%
    mutate(., cell = rownames(.))
  plot.df <- left_join(umap, meta.df)
  p <- ggplot(plot.df, aes_string(x = "UMAP1", y = "UMAP2", color = meta)) +
    geom_point(size = pt.size) +
    ggtitle(meta) +
    theme_classic() +
    theme(aspect.ratio = 1)
  if (simple.coloring == F) {
    if (is.numeric(meta.df[, meta]))
      p <- p + scale_color_gradientn(colors = rev(brewer.pal(5,"Spectral")),na.value = "grey70")
    else {
      p <- p + scale_color_brewer(palette = brewer.pallete)
    }
  }
  p <- wrap.plots.fanc(plot.list = list(p), plot.out = plot.out, ...)
  return(p)
}

# cell.number.distro <- function(df, cluster.ident, clusters = NULL,
#                                group.ident, groups = NULL,
#                                out.file = NULL, out.plot = NULL) {
#   if (!is.null(clusters))
#     df <- df[df[, cluster.ident] %in% clusters,]
#   if (!is.null(groups))
#     df <- df[df[, group.ident] %in% groups, ]
#   df %>% group_by(!!as.name(cluster.ident), !!as.name()) %>% summarise()
# }

# plot.cell.number.distro <- function(count.o, normalize.by)

# confusion.plot <- function(so, )
highlight.cells <- function(so, hl.list, is.Rds, plot.out, order, ...) {
  so <- add.highlight(so = so, add.cells.list = hl.list, is.Rds = is.Rds)
  p <- plot.panel.list(names(hl.list), obj = so, plot.out = plot.out, assay = "RNA",
                       order = order, ...)
  return(p)
}

mean.sd.plot <- function(mat, use.mean.rank = T, use.diff = F, use.sd.over.mean = F,
                         bLog2p1 = F, color.density = T,
                         plot.out = NULL, ...) {
  # columns are samples rows are genes
  mat <- as.matrix(mat)
  df <- data.frame(mean = rowMeans(mat))
  if (use.diff && ncol(mat) == 2) {
    df$diff <- abs(mat[, 1] - mat[, 2])
    measure <- "diff"
  } else {
    df$sd <- rowSds(mat)
    measure <- "sd"
  }

  if (use.sd.over.mean) {
    measure.old <- measure
    measure <- paste0(measure.old, "_over_mean")
    df[, measure] <- (1 + df[, measure.old])/(1 + df$mean)
    if (any(is.na(df[, measure])))
      stop("miao")
  }

  if (bLog2p1) {
    df$mean <- log2(df$mean + 1)
    if (!use.sd.over.mean) {
      df[, measure] <- log2(df[, measure] + 1)
    }
  }

  df <- df %>% mutate(mean.rank = rank(mean, ties.method = "random"))
  x.use <- ifelse(use.mean.rank, "mean.rank", "mean")
  p <- xy.plot(df = df, x = x.use, y = measure,
               add.abline = !use.mean.rank, color.density = color.density, ...)
  if (!use.mean.rank && !use.sd.over.mean) {
    p <- p + coord_fixed()
  }
  if (!is.null(plot.out))
    trash <- wrap.plots.fanc(plot.list = list(p),plot.out = plot.out)
  return(p)
}

chromVAR.ridge <- function(obj, motif.names = NULL, motif.regex = NULL,
                           group.by, groups = NULL, assay = "z",
                           use.ridge = T, plot.out, threads = 1) {
  if (is.null(motif.names)) {
    motif.names <- rownames(obj) %>% .[grepl(motif.regex,.)]
  }
  motif.names <- motif.names %>% .[.%in% rownames(obj)]

  if (is.null(groups)) {
    groups <- obj[[group.by]] %>% unique()
  }
  if (length(motif.names) < 1)
    stop("length(motif.names) < 1")
  pl <- utilsFanc::safelapply(motif.names, function(motif) {
    df <- data.frame(dev = assays(obj)[[assay]][motif,],
                     cells = colnames(obj),
                     group = obj[[group.by]] %>% as.character())
    df <- df %>% filter(group %in% groups) %>%
      mutate(group = factor(group, levels = groups))
    if (assay == "deviations")
      assay <- "raw"

    if (use.ridge == T) {
      p <- ggplot(df, aes(x = dev, fill = group, y = group)) +
        ggridges::geom_density_ridges() +
        theme_ridges() +
        theme(legend.position = "none")
      sub.width <- 4
      sub.height <- length(groups)
    } else {
      p <- ggplot(df, aes(x = dev, fill = group)) +
        geom_density(alpha = 0.5) + theme_classic() +theme(aspect.ratio = 1)
      sub.width <- 6
      sub.height <- 5
    }
    p <- p +
      ggtitle(paste0(motif, " ", assay))
    return(p)
  }, threads = threads)
  if (use.ridge == T) {
    sub.width <- 4
    sub.height <- length(groups)
  } else {
    sub.width <- 6
    sub.height <- 5
  }

  trash <- wrap.plots.fanc(plot.list = pl, plot.out = plot.out, sub.width = sub.width,
                  sub.height = sub.height)

}

cluster.composition.bar <- function(obj, cluster.ident, clusters = NULL, split.by, split.order = NULL,
                                    plot.by = NULL, plot.order = NULL,
                                    standardize = T, plot.out) {
  # split.by determines what you have each column. This will be used as x axis
  # plot.by determines what you have each plot.
  df <- get.metadata.df(obj = obj)
  if (!is.null(clusters))
    df <- df[df[, cluster.ident] %in% clusters,]
  if (!is.null(split.order)) {
    df <- df[df[, split.by] %in% split.order,]
    df[, split.by] <-  df[, split.by] %>% as.character() %>% factor(., levels = split.order)
  }


  if (is.null(plot.by)) {
    plot.by <- "composition"
    df$composition <- "composition"
  }
  if (!is.null(plot.order)) {
    df <- df[df[, plot.by] %in% plot.order,]
    df[, plot.by] <-  df[, plot.by] %>% as.character() %>% factor(., levels = plot.order)
  }

  df <- df[!is.na(df[, split.by]),]
  df <- df[!is.na(df[, cluster.ident]),]
  df <- df[!is.na(df[, plot.by]),]

  # browser()
  pl <- df %>% split(., f = .[, plot.by]) %>%
    lapply(function(df2) {

      df3 <- df2 %>% group_by(!!as.name(split.by), !!as.name(cluster.ident)) %>%
        summarise(n = n()) %>% ungroup()
      if (standardize == T) {
        df3 <- df3 %>% group_by(!!as.name(split.by)) %>%
          summarise(n = n/sum(n)) %>% ungroup()
      }
      p <- ggplot(df3, aes_string(x = split.by, y = "n", fill = cluster.ident)) +
        geom_bar(stat = "identity") +
        theme_classic()
      return(p)
    })
  trash <- wrap.plots.fanc(plot.list = pl, plot.out = plot.out)
  return()
}

chromVAR.umap <- function(ao, motif.mat.name, impute.weights = getImputeWeights(ao),
                          embedding = "UMAP", motifs.regex = "z:",
                          motifs.use = NULL, sort.by.motifs.use = T,
                          plot.out, n = NULL, ...) {
  motifs <- getFeatures(ao, useMatrix = motif.mat.name)
  if (!is.null(motifs.regex))
    motifs <- motifs %>% .[grepl(motifs.regex, .)]

  if (!is.null(motifs.use)) {
    motifs <- motifs %>% .[.%in% motifs.use]
    if (sort.by.motifs.use == T) {
      motifs <- utilsFanc::sort.by(motifs, motifs.use)
    }
  }


  if (!is.null(n)) {
    if (length(motifs) > n)
      motifs <- motifs[1:n]
  }

  p <- plotEmbedding(
    ArchRProj = ao,
    colorBy = motif.mat.name,
    name = motifs,
    embedding = "UMAP",
    continuousSet = "whiteBlue",
    imputeWeights = impute.weights, ...
  )
  p <- lapply(p, function(x) {
    x  +
      theme_ArchR(baseSize = 6.5) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
      ) +
      theme(legend.key.height = unit(0.2, 'cm'))
  })

  trash <- wrap.plots.fanc(plot.list = p, plot.out = plot.out, page.limit = 16, n.col = 4)
  return()
}

save.base.plot <- function(p = NULL, exp, file,
                           width = 700, height = 700, res = 100) {
  dir.create(dirname(file), showWarnings = F, recursive = T)
  if (grepl("png$", file)) {
    png(file, width = width, height = height, res = res)
  } else if (grepl("pdf$", file)) {
    cairo_pdf(file = file, width = width/res, height = height/res)
  }

  try({
    if (is.null(p))
      print(eval(substitute(exp)))
    else
      print(p)
  })

  dev.off()
}


plot.motif.logo <- function(motif.names = NULL, arche.names, map = ARCHETYPE.MAP,
                            logo.db = LOGO.DB, seed = 42,
                            max.col = 10, enforce.same.row = T,
                            threads = 1,
                            plot.out = NULL, return.plot.list = F) {
  if (is.null(motif.names)) {
    motif.names <- arche.2.motifs(arche.names, map)
  }
  if (!is.list(motif.names)) {
    motif.names <- list(default = motifs.names)
  }
  if (is.character(logo.db))
    logo.db <- readRDS(logo.db)
  pl <- utilsFanc::safelapply(names(motif.names), function(type.name) {
    motif.name <- motif.names[[type.name]]
    # note motif.name should still be a vector with length >= 1.
    # I messed up the naming and used motif.names earlier
    pl <- lapply(motif.name, function(x) {
      if (enforce.same.row == T)
        title <- paste0(type.name, "::", x)
      else
        title <- x
      p <- universalmotif::view_motifs(motifs = logo.db[[x]]) +
        ggtitle(title) +
        theme(text = element_text(size=20))
      return(p)
    })
    return(pl)
  }, threads = threads)

  if (return.plot.list) {
    return(pl)
  }

  if (enforce.same.row == T) {
    pl <- lapply(pl, function(x) {
      if (length(x) > max.col) {
        set.seed(seed = seed)
        x <- x[sample(1:length(x), size = max.col, replace = F)]
      }
      return(x)
    })
    max.col <- min(max.col, max(sapply(pl, length)))
    pl <- lapply(pl, function(x) {
      if (length(x) < max.col) {
        x <- c(x, list(NULL) %>% rep(max.col - length(x)))
      }
      return(x)
    })
    n.col <- max.col
  } else {
    n.col <- NULL
  }

  pl <- unlist(pl, recursive = F)

  p <- wrap.plots.fanc(pl, n.col = n.col, page.limit = 100, sub.height = 2, sub.width = 8,
                       plot.out = plot.out)
  invisible(p)
}

deseq2.scatter <- function(dds, contrast.df = NULL, contrasts, genes.include = NULL,
                           transformation = NULL,
                           highlight.gene = F, highlight.values = NULL,
                           n.split = 2, plot.out,
                           ...) {
  # contrast.df <- data.frame(x = "MPP..ctrl_rep1", y = "MPP..tumor_rep1")
  # contrasts: c("cluster", "rep", "type")
  if (missing(plot.out))
    stop("must supply plot.out")
  highlight.var <- "gene"
  if (highlight.gene == F)
    highlight.var <- NULL
  if (is.null(contrast.df)) {
    if (length(contrasts) != 3) {
      stop("only developed for 3 term contrast")
    }
    coldata <- dds@colData %>% as.data.frame()
    coldata$name <- rownames(coldata)
    contrast.df <- coldata %>% split(., .[, contrasts[1]]) %>%
      lapply(function(df.1) {
        df.1 %>% split(., .[, contrasts[2]]) %>%
          lapply(function(df.2) {
            if (nrow(df.2) != 2)
              stop("nrow(df.2) != 2")
            df.2 <- df.2 %>% arrange(!!as.name(contrasts[3]))
            df.out <- data.frame(x = df.2$name[1], y = df.2$name[2])
            return(df.out)
          }) %>% Reduce(rbind, .) %>% return()
      }) %>% Reduce(rbind,.)
  }
  norm.mat <- counts(dds, normalized = T)
  if (!is.null(genes.include))
    norm.mat <- norm.mat[genes.include,]
  pl <- contrast.df %>% split(., f = 1:nrow(.)) %>%
    lapply(function(df) {
      norm.df <- norm.mat[, c(df$x, df$y)] %>% as.data.frame()
      norm.df$gene <- rownames(norm.df)
      p <- xy.plot(df = norm.df, x = df$x, y = df$y, transformation = transformation,
                   highlight.var = highlight.var, highlight.values = highlight.values,
                   show.highlight.color.var = F,
                   plotly.var = NULL, color.density = T, ...)
      return(p)
    })
  trash <- wrap.plots.fanc(plot.list = pl, n.split = n.split, plot.out = plot.out)
  return()
}

dim.plot.simple <- function(so, plot.out = NULL, group.by = "seurat_clusters", label = T) {
  p <- DimPlot(object = so, label = label, label.size = 6, pt.size = 0.1, shuffle = T,
               group.by = group.by) +
    theme( aspect.ratio = 1,
           axis.title.x = element_text(size = 10, face = "bold"),
           axis.title.y = element_text(size = 10, face = "bold"),
           axis.ticks=element_blank(),
           axis.text = element_blank(),
           axis.line = element_line(size = 0.1, colour = "grey50"))
  if (!is.null(plot.out)) {
    dir.create(dirname(plot.out), showWarnings = F, recursive = T)
    ggsave(filename = plot.out,plot = p, width = 5, height = 5, units = "in", dpi = 100)
  }

  return(p)
}

get.tss.region <- function(genes, genome, buffer = NULL, buffer.up, buffer.down,
                           out.file = NULL) {
  if (!is.null(buffer)) {
    buffer.up <- buffer
    buffer.down <- buffer
  }
  tss.df <- utilsFanc::import.bed.fanc(paste0("~/genomes/", genome, "/", genome, "_TSS.bed"))
  loci <- mapply(function(gene, up, down) {
    if (!grepl(":", gene)) {
      gene <- tss.df %>% dplyr::filter(forth == gene) %>% utilsFanc::gr.get.loci() %>% .[1]
    }
    loci <- utilsFanc::loci.2.df(loci.vec = gene, remove.loci.col = T)
    loci <- loci %>% dplyr::mutate(start = start - up, end = end + down)
    return(loci)
  }, genes, buffer.up, buffer.down, SIMPLIFY = F) %>% Reduce(rbind, .)
  if (!is.null(out.file)) {
    utilsFanc::write.zip.fanc(loci, out.file = out.file, bed.shift = T)
  }
  return(loci)
}

ranked.bar <- function(mat, rows, rank.by = NULL, bDescending = T,
                       width = NULL, width.each.bar = 0.3, height = 5,
                       hide.x.axis = F,
                       plot.out = NULL) {
  # take in a matrix (eg: gene x cell). rows: eg: c(Cebpb, Gata1)
  # each bar will be a cell, with Cebpb and Gata1 values one on top of the other (stack, not dodge)
  # bars will be sorted by the value of rank.by.
  rows <- intersect(rows, rownames(mat))
  if (length(rows) < 1) {
    stop("none of the rows specified could be found in the rownames of mat")
  }
  if (is.null(rank.by)) {
    rank.by <- rows
  }
  rank.by <- intersect(rank.by, rows)
  if (length(rank.by) < 1) {
    stop("none of rank.by was found in rows")
  }
  mat <- mat[rows, ]
  df <- reshape2::melt(mat, varnames = c("row", "col"), value.name = "value")
  pl <- lapply(rank.by, function(rank.by.sub) {
    col.levels <- colnames(mat)[order(mat[rank.by.sub, ])]
    if (bDescending) {
      col.levels <- rev(col.levels)
    }
    row.levels <- c(rank.by.sub, rows) %>% unique()
    df$col <- factor(df$col, levels = col.levels)
    df$color <- df$row
    df$row <- factor(df$row, levels = row.levels)
    p <- ggplot(df, aes(x = col, y = value, fill = color)) +
      geom_bar(stat = "identity") +
      xlab("") + theme_classic()
    if (hide.x.axis)
      p <- p + theme(axis.ticks.x = element_blank(),
                     axis.text.x = element_blank())
    return(p)
  })
  if (is.null(width)) {
    width <- ncol(mat) * width.each.bar
  }
  p <- wrap.plots.fanc(plot.list = pl, sub.width = width,
                       sub.height = height, plot.out = plot.out, n.col = 1)
  invisible(p)
}

rank.plot.synced <- function(df, rank.by = NULL, bDescending = T,
                             width = 15, height = 5, n.col = 1,
                             plot.out = NULL) {
  # note: all columns of df will be plotted. be careful what you give to it...
  if (!is.data.frame(df)) {
    df <- as.data.frame(df)
  }
  if (is.null(rank.by)) {
    rank.by <- colnames(df)
  }
  rank.by <- intersect(rank.by, colnames(df))
  if (length(rank.by) < 1) {
    stop("rank.by not found in colnames(df)")
  }
  pl <- lapply(rank.by, function(rank.by.sub) {
    order <- order(df[, rank.by.sub])
    if (bDescending)
        order <- rev(order)
    df <- df[order, ]
    df$rank <- 1:nrow(df)
    df <- reshape2::melt(df, id.vars = "rank", variable.name = "var", value.name = "value")
    p <- ggplot(df, aes(x = rank, y = value, group = var)) +
      geom_line(aes(color = var)) +
      xlab("") + theme_classic()
    return(p)
  })
  p <- wrap.plots.fanc(plot.list = pl, sub.width = width,
                       sub.height = height, plot.out = plot.out, n.col = n.col)
  invisible(p)
}


umap.axis.schema <- function(p, x.slot = "UMAP_1", y.slot = "UMAP_2",
                             x.shift.pct = 0.08, y.shift.pct = 0.08, 
                             x.length.pct = 0.3, y.length.pct = 0.3,
                             line.size = 0.2, 
                             # x.text.nudge = 0, y.text.nudge = 0.2,
                             text.size = 2) {
  x.range <- c(min(p$data[[x.slot]]), max(p$data[[x.slot]]))
  x.len <- x.range[2] - x.range[1]
  y.range <- c(min(p$data[[y.slot]]), max(p$data[[y.slot]]))
  y.len <- y.range[2] - y.range[1]
  
  d <- data.frame(x = c(x.range[1], x.range[1] + x.shift.pct * x.len), 
                  y = c(y.range[1] + y.shift.pct * y.len, y.range[1]),
                  xend = c(x.range[1] + x.length.pct * x.len, x.range[1] + x.shift.pct * x.len), 
                  yend = c(y.range[1] + y.shift.pct * y.len, y.range[1] + y.length.pct * y.len))
  
  p <- p + geom_segment(data = d, aes(x = x, y = y, xend = xend, yend = yend),
                        inherit.aes = F,
                        size = line.size)
  dt <- data.frame(x = c(x.range[1] + (x.shift.pct + 0.05) * x.len, x.range[1]), 
                  y = c(y.range[1], y.range[1] + (y.shift.pct + 0.05) * y.len),
                  # hjust = c(-1 * x.text.nudge, x.shift.pct * 0.5),
                  # vjust = c(y.shift.pct * 0.5, y.text.nudge),
                  hjust = c(0, 0), 
                  vjust = c(0, 1),
                  angle = c(0, 90),
                  label = c(x.slot, y.slot))
  # hjust and vjust: for UMAP_1, position determined by lower left corner
  # for UMAP_2, position determined by upper left corner
  # then UMAP_2 gets rotated 90 degrees. Note we first determine its position, before rotating!
  # it is rotated using the upper left corner as the pivot.
  p <- p + geom_text(data = dt, 
                     aes(x = x, y = y, label = label, angle = angle,
                         hjust = hjust, vjust = vjust),
                     inherit.aes = F,
                     size = 2)
  
  return(p)
}


legend.plot <- function(df, cluster.col = "cell_type", color.col = "color",
                        text.size = 5, pt.size = 1, 
                        # text size is the normal size. the code will do 0.36 * text.size
                        # for geom_text()
                        width = 1, height = 2,
                        out.file) {
  # df: 2 columns, cluster and color.
  if (is.character(df)) {
    df <- read.table(df, header = T, comment.char = "", sep = "\t")
  }
  
  utilsFanc::check.intersect(c(cluster.col, color.col), "required columns", 
                             colnames(df), "colnames(df)")
  
  df <- df[, c(cluster.col, color.col)]
  colnames(df) <- c("cluster", "color")
  df <- df[nrow(df):1,]
  
  df$y <- 1:nrow(df)
  df$color[df$color == ""] <- "white"
  
  pdf <- df[, c("y", "color")]
  pdf$x <- 1
  
  cdf <- df[, c("y", "cluster")]
  cdf$x <- 1.25
  cdf$x[pdf$color == "white"] <- 1
  
  color.vec <- df$color %>% unique()
  names(color.vec) <- color.vec
  
  p <- ggplot(pdf, aes(x = x, y = y)) +
    geom_point(aes(color= color), size = pt.size) +
    geom_text(data = cdf, aes(x = x, y = y, label = cluster, hjust = 0), 
              size = 0.36 * text.size, family = "Arial",
              inherit.aes = F) +
    scale_color_manual(values = color.vec) + 
    xlim(c(1, 2.5)) +
    scale_y_continuous(expand = expansion(add = 1)) +
    theme_void() +
    theme(legend.position = "none")
  dir.create(dirname(out.file), showWarnings = F, recursive = T)
  ggsave(out.file, p, device = cairo_pdf, width = width, height = height,
         dpi = 300)
  invisible(p)
}


umap.fanc <- function(obj, 
                      # metadata mode:
                      group.by, groups = NULL, 
                      label.groups = F, label.size = 2,
                      plot.title = NULL, title.size = 6,
                      highlight.mode = F, 
                      # common
                      remove.outlier = c(0, 0, 0, 0),
                      reverse.point.order = F,
                      split.by = NULL, splits = NULL,
                      polygon.df = NULL, # polygon untested
                      cols = NULL, pt.size = 0.05, pt.shape = 18,
                      show.legends = F, 
                      width = 2, height = 2, dpi = 300, 
                      axis.type = "schema",
                      plot.out = NULL, ...) {
  # remove.outlier: c(top, right, bottom, left). if top is 1, then remove the top most point.
  if ("Seurat" %in% class(obj)) {
    df <- obj@reductions$umap@cell.embeddings %>% as.data.frame()
    
    utilsFanc::check.intersect(group.by, "group.by", colnames(obj@meta.data), "colnames(obj@meta.data)")
    df$group <- obj@meta.data[, group.by]
    
    if (!is.null(split.by)) {
      utilsFanc::check.intersect(split.by, "split.by", colnames(obj@meta.data), "colnames(obj@meta.data)")
      df$split <- obj@meta.data[, split.by]
    }
    
  } else if ("ArchRProject" %in% class(obj)) {
    df <- obj@embeddings$UMAP$df
    colnames(df) <- c("UMAP_1", "UMAP_2")
    df$group <- ArchR::getCellColData(obj, select = c(group.by), drop = T)
    if (!is.null(split.by)) {
      df$split <- ArchR::getCellColData(obj, select = c(split.by), drop = T)
    }
  }
  
  if (is.null(split.by)) {
    df$split <- ""
  }
  
  if (!is.null(groups)) {
    if (highlight.mode == T) {
      df$group <- as.character(df$group)
      groups <- as.character(groups)
      df$group[df$group %in% groups] <- "inHl"
      df$group[!df$group %in% "inHl"] <- "outHl"
      
      show.legends <- F
      
      if (is.null(cols)) {
        cols <- c("blue","grey80")
        names(cols) <- c("inHl", "outHl")
      } else if (!identical(sort(names(cols)), c("inHl", "outHl"))) {
        stop("highlight mode. cols must be named inHl and outHl")
      }
    } else {
      df <- df %>% dplyr::filter(group %in% groups)
    }
  }
  
  if (!is.null(splits)) {
    df <- df %>% dplyr::filter(split %in% splits)
    splits <- splits %>% .[.%in% df$split]
    df$split <- factor(df$split, levels = splits)
  } else {
    splits = ""
  }
  
  if (sum(remove.outlier) > 0) {
    df <- df %>% dplyr::mutate(
      top = max(rank(UMAP_2)) - rank(UMAP_2) + 1,
      right = max(rank(UMAP_1)) - rank(UMAP_1) + 1,
      bottom = rank(UMAP_2),
      left = rank(UMAP_1)
      )
    names(remove.outlier) <- c("top", "right", "bottom", "left")
    for (i in c("top", "right", "bottom", "left")) {
      df <- df[df[, i] > remove.outlier[i],]
    }
  }
  
  if (nrow(df) < 1) {
    stop("nrow(df) < 1")
  }
  
  xmax <- max(df$UMAP_1)
  xmin <- min(df$UMAP_1)
  ymax <- max(df$UMAP_2)
  ymin <- min(df$UMAP_2)
  
  n.splits <- df$split %>% unique() %>% length()
  pl <- df %>% split(f = df$split) %>% 
    lapply(function(df) {
      if (reverse.point.order) 
        df <- df[nrow(df):1,]
      if (label.groups) {
        if (!is.factor(df$group)) stop("groups must be factorized for label.groups to work  ")
      }
      p <- ggplot(df, aes_string(x = "UMAP_1", y = "UMAP_2")) +
        ggrastr::rasterise(geom_point(aes_string(color = "group", fill = "group"), 
                   size = pt.size, shape = pt.shape, stroke = 0.05), dpi = 500)
      if (label.groups) {
        p <- Seurat::LabelClusters(p, "group", size = label.size, repel = F, family = "Arial")
      }
      if (!is.null(cols)) {
        p <- p + scale_color_manual(values = cols)
      }
      
      p <- p + 
        xlim(xmin, xmax) + 
        ylim(ymin, ymax) 
      
      if (!is.null(polygon.df))
        p <- p + geom_polygon(data = polygon.df, mapping = aes(x = x, y = y),
                              fill = NA, color = "red")
      
      if (axis.type %in% c("schema", "nothing")) {
        p <- p + theme_void()
        if (axis.type == "schema")
          p <- umap.axis.schema(p)
      } else if (axis.type == "number") {
        p <- p + theme_bw()
      } else {
        stop("axis.type must be schema, number, or nothing")
      }
      
      if (n.splits > 1 ||  !is.null(plot.title)) {
        if (n.splits > 1) {
          p <- p + ggtitle(paste0(df$split[1]))
        }
        if (!is.null(plot.title)) {
          p <- p + ggtitle(paste0(plot.title))
        }
        p <- p + theme(plot.title = element_text(
          size = title.size, family = "Arial", color = "black", hjust = 0.5))
      }
      
      if (!show.legends) {
        p <- p + theme(legend.position = "none")
      }
      p <- p + theme(aspect.ratio = 1)
      return(p)
    })
  
  if (length(pl) > 1) {
    p <- wrap.plots.fanc(plot.list = pl, n.split = length(splits), plot.out = plot.out, ...)
  } else {
    p <- pl[[1]]
    if (!is.null(plot.out)) {
      dir.create(dirname(plot.out), showWarnings = F, recursive = T)
      ggsave(plot.out, p, width = width, height = height, dpi = dpi, 
             device = ifelse(grepl("pdf$", plot.out), cairo_pdf, tools::file_ext(plot.out)))
    }
  }

  invisible(p)
}

cluster.freq.bar <- function(soi, x, xs = NULL, group.by, groups = NULL, 
                             average.across, average.groups = NULL, average.method = "pct",
                             color.map = NULL, publication = T,
                             out.file = NULL) {
  # it plots the type of bar graph where each column sums to 1
  # x: what you want the x axis to be
  # xs: the x's that you want to plot. For example, you might only want to plot 2 out of 6 samples
  # group.by: e.g. seurat_clusters or cell_type
  # average.across: for example, "rep", which means you average across each replicate
  # average.method: pct or number. If pct, we first compute the percentage of each cluster 
  # for each group, and then average across, say, replicates. if number, we first merge all replicates
  # and then calculate percentage of each cluster
  # color.map: a vector of colors. names(color.map) are group names
  utilsFanc::check.intersect(c(x, group.by, average.across), "required columns",
                             colnames(soi@meta.data), "colnames(soi@meta.data)")
  meta <- soi@meta.data[, c(x, group.by, average.across)]
  colnames(meta) <- c("x", "group.by", "average.across")
  
  if (!is.null(xs)) {
    meta <- meta %>% dplyr::filter(x %in% xs)
  }
  if (!is.null(groups)) {
    meta <- meta %>% dplyr::filter(group.by %in% groups)
    meta$group.by <- factor(meta$group.by, levels = unique(groups))
  }
  if (!is.null(average.groups))
    meta <- meta %>% dplyr::filter(average.across %in% average.groups)
  
  if (nrow(meta) < 1) {
    stop("after filtering no cells were left")
  }
  
  plot.df <- meta %>% split(., f = .$x) %>% lapply(function(xdf) {
    tmp <- xdf$x[1]
    if (average.method == "pct") {
      xdf <- xdf %>% dplyr::group_by(group.by, average.across) %>% dplyr::summarise(n = n()) %>% 
        dplyr::ungroup() %>% dplyr::group_by(average.across) %>% dplyr::mutate(pct = n/sum(n)) %>% 
        dplyr::ungroup() %>% dplyr::group_by(group.by) %>% 
        dplyr::summarise(pct = mean(pct))
    } else if (average.method == "number") {
      xdf <- xdf %>% dplyr::group_by(group.by) %>% dplyr::summarise(n = n()) %>% 
        dplyr::ungroup() %>% dplyr::mutate(pct = n/sum(n))
    } else {
      stop("average.method must be pct or number")
    }
    xdf <- dplyr::ungroup(xdf) %>% as.data.frame() 
    xdf$x <- tmp
    return(xdf)
  }) %>% do.call(rbind, .)

  p <- ggplot(plot.df, aes(x = x, y = pct)) +
    geom_bar(aes(fill = group.by), stat = "identity", show.legend = !publication)
  
  if (!is.null(color.map)) {
    utilsFanc::check.intersect(unique(plot.df$group.by), "groups",
                               names(color.map), "names(color.map)")
    p <- p + scale_fill_manual(values = color.map)
  }
  if (publication && !is.null(out.file)) {
    p <- p %>% utilsFanc::theme.fc.1(rotate.x.45 = F, italic.x = F, rm.x.ticks = T)
    p <- p + coord_flip()
    p <- p + scale_y_continuous(breaks = c(0, 1)) + theme(axis.ticks.y = element_blank())
    dir.create(dirname(out.file), showWarnings = F, recursive = T)
    ggsave(out.file, p, device = cairo_pdf, width = 1.2, height = 1)
    
    # plot legend:
    color.map <- color.map %>% .[names(.) %in% unique(plot.df$group.by)]
    color.df <- data.frame(cell_type = names(color.map), color = color.map)
    legend.plot(df = color.df, pt.size = 1.5,
                out.file = paste0(tools::file_path_sans_ext(out.file), "_legend.pdf" ),
                width = 1, height = 0.1 * length(color.map))
  } else {
    trash <- wrap.plots.fanc(list(p), plot.out = out.file)
  }
  
  invisible(p)
}

plot.boxplot <- function(df, x, x.include = NULL, y, group.by, groups = NULL,
                         normalize.to.x = NULL,
                         out.dir = NULL, root.name) {
  # example: ~/hmtp/scAR/spbm2/DC/DC6.1_DEG_as_module_2024-04-02.R
  utilsFanc::check.intersect(c(x, y, group.by), "required columns",
                             colnames(df), "colnames(df)")
  df <- df[, c(x, y, group.by)]
  colnames(df) <- c("x", "y", "group.by")
  if (!is.null(x.include)) {
    x.include <- intersect(x.include, unique(df$x))
    df <- df[df$x %in% x.include,]
    df$x <- factor(as.character(df$x), levels = x.include)
  }
    
  if (!is.null(groups)) {
    groups <- intersect(groups, unique(df$group.by))
    df <- df[df$group.by %in% groups,]
    df$group.by <- factor(as.character(df$group.by), levels = groups)
  }
  if (nrow(df) < 1) stop("nrow(df) < 1")
  
  if (!is.null(normalize.to.x)) {
    utilsFanc::check.intersect(normalize.to.x, "normalize.to.x",
                               df$x, "df[, x]")
    df <- df %>% split(f = df$group.by) %>% lapply(function(df) {
      norm.factor <- mean(df$y[df$x == normalize.to.x])
      df$y <- df$y/norm.factor
      return(df)
    }) %>% do.call(rbind, .)
  }
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_boxplot(aes(color = group.by)) +
    xlab(x) + ylab(y)
  if (!is.null(out.dir)) {
    dir.create(out.dir, recursive = T, showWarnings = F)
    n.x <- df$x %>% unique() %>% length()
    ggsave(paste0(out.dir, "/", root.name, ".pdf"), p, device = cairo_pdf, width = n.x, height = 2)
  }
  invisible(p)
}



plot.line <- function(df, x, x.include = NULL, y, group.by, groups = NULL,
                         normalize.to.x = NULL,
                         out.dir = NULL, root.name) {
  # example: ~/hmtp/scAR/spbm2/DC/DC6.1_DEG_as_module_2024-04-02.R
  utilsFanc::check.intersect(c(x, y, group.by), "required columns",
                             colnames(df), "colnames(df)")
  df <- df[, c(x, y, group.by)]
  colnames(df) <- c("x", "y", "group.by")
  if (!is.null(x.include)) {
    x.include <- intersect(x.include, unique(df$x))
    df <- df[df$x %in% x.include,]
    df$x <- factor(as.character(df$x), levels = x.include)
  }
  
  if (!is.null(groups)) {
    groups <- intersect(groups, unique(df$group.by))
    df <- df[df$group.by %in% groups,]
    df$group.by <- factor(as.character(df$group.by), levels = groups)
  }
  if (nrow(df) < 1) stop("nrow(df) < 1")
  
  if (!is.null(normalize.to.x)) {
    utilsFanc::check.intersect(normalize.to.x, "normalize.to.x",
                               df$x, "df[, x]")
    df <- df %>% split(f = df$group.by) %>% lapply(function(df) {
      norm.factor <- mean(df$y[df$x == normalize.to.x])
      df$y <- df$y/norm.factor
      return(df)
    }) %>% do.call(rbind, .)
  }
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_boxplot(aes(color = group.by)) +
    xlab(x) + ylab(y)
  if (!is.null(out.dir)) {
    dir.create(out.dir, recursive = T, showWarnings = F)
    n.x <- df$x %>% unique() %>% length()
    ggsave(paste0(out.dir, "/", root.name, ".pdf"), p, device = cairo_pdf, width = n.x, height = 2)
  }
  invisible(p)
}

plot.confusion.mat <- function(mat, show_column_dend = F, 
                               no.col.cluster = T, no.row.rank = T,
                               plot.out,
                               width = 2, height = 2, ...) {
  rows <- rownames(mat)
  mat <- diag(1/rowMax(mat)) %*% mat
  rownames(mat) <- rows
  plot.mat.rank.row(mat = mat, no.ranking = no.row.rank, no.col.cluster = no.col.cluster,
                    hm.colors = c("midnightblue", "yellow"),
                    hm.values = c(0, 1), show_column_names = T, show_row_names = T, 
                    show_column_dend = show_column_dend, plot.out = plot.out,
                    width = width, height = height, row_names_side = "left", ...)
  return()
}

plot.single.cell.exp.mat <- function(so, assay, layer, binarize = T, cells, metas = NULL, gene.list, 
                                     plot.out, width = 8, height = 8,
                                     show.row.names = T, show.column.names = F, annotate.rows = F, ...) {
  # gene.list = list(myeloid= c(gene1, gene2), lymphoid = c(gene3, gene4))
  mat <- Seurat::GetAssayData(so, assay = assay, layer = layer)
  utilsFanc::check.intersect(cells, "cells", colnames(mat), paste0("assay ", assay, " and layer ", layer))
  utilsFanc::check.intersect(unlist(gene.list), "genes", rownames(mat), paste0("assay ", assay, " and layer ", layer))
  utilsFanc::check.intersect(metas, "metas", colnames(so@meta.data), "colnames(meta.data)")
  mat <- mat[unlist(gene.list), cells]
  
  if (binarize) {
    mat[mat > 0] <- 1
    mat[mat <= 0] <- 0
  } else {
    stop("only binarize == T has been developed")
  }
  
  # suppressMessages(extrafont::loadfonts())
  # ht_opt$HEATMAP_LEGEND_PADDING <- unit(0.1, "in")
  # ht_opt$DENDROGRAM_PADDING <- unit(0, "in")
  
  
  col_fun <- circlize::colorRamp2(c(0, 1), c("white", "black"))
  hm.param <- list(matrix = as.matrix(mat), col = col_fun,
                   show_row_names = show.row.names, show_column_names = show.column.names,
                   # row_names_gp = gpar(fontsize = 6), 
                   # heatmap_legend_param = list(
                   #   title = "", labels_gp = gpar(fontsize = 5), title_gp = gpar(fontsize = 6),
                   #   legend_height = unit(0.3, "in"), grid_width = unit(0.05, "in"), gap = unit(2, "in")
                   # ),
                   column_order = cells, ...)
  
  if (!is.null(metas)) {
    
    anno.df <- so@meta.data[cells, metas, drop = F]
    colors <- lapply(anno.df, function(meta) {
      meta <- unique(meta)
      color.map <- utilsFanc::gg_color_hue(length(meta))
      names(color.map) <- meta %>% as.character()
      return(color.map)
    })
    legends <- lapply(anno.df, function(x) return(list(title = "")))
    sa <- ComplexHeatmap::HeatmapAnnotation(df = anno.df, col = colors, 
                                            annotation_legend_param = legends)
    
    hm.param[["top_annotation"]] <- sa
  }
  
  hm.param[["row_order"]] <- unlist(gene.list)
  
  if (annotate.rows) {
    stop("annotate.rows have not been well developed")
    row.anno <- rep(names(gene.list), sapply(gene.list, length))
    row.anno <- data.frame(gene.type = row.anno)
    
    ha <- ComplexHeatmap::HeatmapAnnotation(
      df = row.anno, annotation_name_side = "top", which = "row", 
      annotation_name_rot = 0, show_annotation_name = F)
    hm.param[["left_annotation"]] <- ha
    
  }

  hm <- do.call(ComplexHeatmap::Heatmap, hm.param)
  if (!grepl(".pdf$", plot.out)) {
    stop("plot.out must be pdf")
  }
  system(paste0("mkdir -p ", dirname(plot.out)))
  
  cairo_pdf(filename = plot.out, width = width, height = height)
  try({print(hm)})
  dev.off()
  
  ht_opt(RESET = T)
  
  invisible(mat)
}

