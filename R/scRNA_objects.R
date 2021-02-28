# #
# path <- .libPaths()
# .libPaths(c("/bar/cfan/R/weird_packages/3.5/", path))
# options(stringsAsFactors = F)

# find.elbow <- function(df, first.deriv.filter=NULL) {
#   # first column is x, second column is y
#   if (ncol(df) != 2)
#     stop("find.elbow: number of columns has to be exactly 2. The first column is x and the second is y")
#   colnames(df) <- c("x", "y")
#   df <- df %>% arrange(x)
#   f <- splinefun(x= df$x, y= df$y)

#   if (is.numeric(first.deriv.filter)) {
#     first.deriv <- f(df$x, deriv = 1)
#     elbow.x <- df$x[which(first.deriv >=first.deriv.filter)[1]]
#   } else {
#     second.deriv <- f(df$x, deriv = 2)
#     elbow.x <- df$x[which(second.deriv==max(second.deriv))[1]]
#   }

#   elbow.y <- df %>% filter(x==elbow.x) %>% pull(y)
#   if (length(elbow.y) !=1)
#     stop(paste0("somehow elbow.x: ",elbow.x," corresponds to ",  length(elbow.y)," y values"))
#   return(list(x=elbow.x, y = elbow.y))
# }

# # find.elbow(t %>% filter (percent.mt > 2, percent.mt < 12) %>% mutate(x=rank.Feature, y=percent.mt) %>% select(x, y), 1.5)

# split.by.metadata <- function(so, sample.info,reduction, outdir=NULL, metas.include=NULL, ...) {
#   if (is.character(sample.info))
#     sample.info <- read.table(sample.info, as.is =T, header = T)
#   metas <- colnames(sample.info)[colnames(sample.info) != "dir"]
#   if (!is.null(metas.include))
#     metas <- metas[metas %in% metas.include]
#   n.plots <- length(metas)
#   if (n.plots < 1)
#     stop("nothing to plot")
#   plot.list <- lapply(metas, function (meta) {
#     p.s <- DimPlot(so, reduction = reduction, label = T, group.by = meta, split.by = meta, ...)
    
#     n.sub.plots <- sample.info[,meta] %>% unique() %>% length()
#     print(n.sub.plots)
#     p.g <- lapply(sample.info[,meta] %>% unique(), function(x) {
#       DimPlot(so, reduction = reduction, label = T, group.by = meta, order = x,repel = T, label.size = 5,...)
#     })
#     p.g <- cowplot::plot_grid(plotlist = p.g, nrow = 1)
#     if (!is.null(outdir)) {
#       ggsave(paste0(outdir, "/", reduction, "_split_", meta, ".png"), p.s, device = "png", width = 5*(n.sub.plots), height = 5,
#              dpi = 100, units = "in", limitsize = F)
#       ggsave(paste0(outdir, "/", reduction, "_group_", meta, ".png"), p.g, device = "png", width = 7*n.sub.plots, height = 5,
#              dpi = 100, units = "in", limitsize = F)
#     }
#     return(list(ps=p.s, pg=p.g))
#   })
  
#   names(plot.list) <- metas
#   return(plot.list)
# }

# feature.plot.comp <- function(marker.list, obj.list, plot.dir=NULL, n.marker, must.present.in=NULL) {

#   p.list <- split(marker.list, f = marker.list$cluster) %>%
#     lapply(function(x) {
#       print(paste0("cluster: ",x$cluster[1]))
#       x <- arrange(x,p_val_adj)
#       if (!is.null(must.present.in)) {
#         for (i in must.present.in)
#           x <- x %>% filter(gene %in% rownames(obj.list[[i]]))
#       }
#       y <- x[1:n.marker, "gene"]
#       y <- y[!is.na(y)]
#       p <- lapply(seq_along(obj.list), function(i) {
#         print(paste0("plotting from obj", i))
#         pi <- FeaturePlot(obj.list[[i]], features = y, reduction = "umap", label=T, combine = F) %>%
#           cowplot::plot_grid(plotlist = ., ncol = 1)
#         return(pi)
#       }) %>% cowplot::plot_grid(plotlist = ., ncol =length(obj.list))
#       if (!is.null(plot.dir)) {
#         system(paste0("mkdir -p ", plot.dir))
#         ggsave(paste0(plot.dir, "/comp_cluster_", x$cluster[1], ".png"), p, units = "in", dpi = 200,
#                width = length(obj.list)*5, height = n.marker*5, device = "png", limitsize = F)
#       }

#       return(p)
#     })
# }

# plot.panel.list <- function(panel.list, obj, split.by=NULL, n.split=1, ident = "seurat_clusters", plot.out = NULL, violin = F) {
#   # stop("currently not working for splits.")
#   # n.split is only used for featurePlot. 
#   Idents(obj) <- ident
  
#   if (violin == T)
#     n.split <- 1
#   all.markers <- panel.list %>% unlist()
#   n.plots <- length(all.markers)
#   n.plots <- n.plots * n.split
#   n.col <- floor(2*(n.plots/2)^(0.5))
#   n.col <- n.col - (n.col %% n.split)
#   if (n.col == 0)
#     n.col <- n.split
#   n.row <- ceiling(n.plots/n.col)
#   in.per.plot <- 4
#   if (is.null(split.by)) 
#     p <- FeaturePlot(obj, all.markers, label = T, order = T, pt.size = 0.05, combine = F, split.by = split.by)
#   else {
#     p <- lapply(all.markers, function(m) {
#       if (violin == T) {
#         p.sub <- VlnPlot(obj, m, split.by = split.by, split.plot = T, 
#                          cols = c("red", "orange", "blue", "green", "gray", "yellow")) 
#         p.sub <- list(p.sub)
#       } else {
#         p.sub <- FeaturePlot(obj, m, label = T, order = T, pt.size = 0.05, combine = F, split.by = split.by)
#         p.sub <- lapply(p.sub, function(psi) {
#           psi <- psi + ggtitle(m) +   
#             theme(axis.title=element_blank(),
#                   axis.text=element_blank(),
#                   axis.ticks=element_blank())
#           return(psi)
#         })
#       }
#       return(p.sub)
#     }) %>% Reduce(c, .)
#   }
  
#   p <- cowplot::plot_grid(plotlist = p, ncol = n.col)
#   if (!is.null(plot.out)) {
#     system(paste0("mkdir -p ", dirname(plot.out)))
#     ggsave(plot.out, p, units = "in", device = "png", dpi = 300, width = n.col*in.per.plot,
#            height=n.row*in.per.plot, limitsize = F)
#   }
#   return(p)
# }



