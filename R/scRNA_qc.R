# sc.multi.qc <- function(sample.info, work.dir, plot.dir, project.name, mt.pattern,
#                         construct.FUN = NULL, construct.arg.list=NULL, elbow.mt) {
#   system(paste0("mkdir -p ", work.dir))
#   system(paste0("mkdir -p ", plot.dir))
#   if (is.character(sample.info))
#     sample.info <- read.table(sample.info, header = T, as.is = T, sep = "\t")



#   sol <- split(sample.info, f=sample.info$sample %>% factor(.,levels=.)) %>%
#     lapply(function(x) {
#       if (!is.null(construct.FUN))
#         so <- do.call(construct.FUN, c(list(x$dir), construct.arg.list))
#       else
#         so <- Read10X(data.dir = x$dir) %>%
#           CreateSeuratObject(project = project.name)
#       so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = mt.pattern)
#       metadata <- colnames(x)[!colnames(x) %in% c("dir")]
#       for (m in metadata) {
#         so[[m]] <- x[1,m]
#       }
#       return(so)
#     })
#   n.gene <- lapply(seq_along(sol), function(i) {
#     p <- sol[[i]]@meta.data %>% mutate(rank.Feature = rank(nFeature_RNA)) %>%
#       ggplot(aes(x=rank.Feature, y=(nFeature_RNA))) +
#       geom_point(size = 0.1) +
#       ylim(0, 10000) +
#       ggtitle(sample.info$sample[i])
#     return(p)
#   })

#   cowplot::plot_grid(plotlist = n.gene) +
#     ggsave(paste0(plot.dir,"/",project.name,  "_gene_number_distro.png"), units = "in",
#            width = 10, height = 10, device = "png", dpi = 100)

#   pct.mito.plots <- list()
#   for (i in seq_along(sol) ) {
#     # sol[[i]]@misc[["elbow.mt"]] <- sol[[i]]@meta.data %>%
#     #   mutate(x = rank(percent.mt), y=percent.mt) %>% select(x,y) %>% arrange(x) %>%
#     #   find.elbow()
#     sol[[i]]@misc[["elbow.mt"]] <- list(y=elbow.mt)
#     p <- sol[[i]]@meta.data %>% mutate(rank.Feature = rank(percent.mt)) %>%
#       ggplot(aes(x=rank.Feature, y=(percent.mt))) +
#       geom_point(size = 0.1) +
#       ylim(0, 100) +
#       geom_hline(yintercept = sol[[i]]@misc[["elbow.mt"]]$y, color = "red") +
#       ggtitle(sample.info$sample[i])
#     pct.mito.plots[[i]] <- p
#   }

#   # plot mito
#   cowplot::plot_grid(plotlist = pct.mito.plots) +
#     ggsave(paste0(plot.dir,"/",project.name,  "_pct_mito.png"), units = "in",
#            width = 10, height = 10, device = "png", dpi = 100)

#   saveRDS(sol, paste0(work.dir, "/", project.name, "_sol.Rds"))
#   return(sol)

# }

# sc.qc.construct.so <- function(sample.info, project.name, mt.pattern,
#                                construct.FUN = NULL, construct.arg.list=NULL) {
#   if (is.character(sample.info))
#     sample.info <- read.table(sample.info, header = T, as.is = T, sep = "\t")

#   sol <- split(sample.info, f=sample.info$sample %>% factor(.,levels=.)) %>%
#     lapply(function(x) {
#       if (!is.null(construct.FUN))
#         so <- do.call(construct.FUN, c(list(x$dir), construct.arg.list))
#       else
#         so <- Read10X(data.dir = x$dir) %>%
#           CreateSeuratObject(project = project.name)
#       so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = mt.pattern)
#       metadata <- colnames(x)[!colnames(x) %in% c("dir")]
#       for (m in metadata) {
#         so[[m]] <- x[1,m]
#       }
#       return(so)
#     })
#   return(sol)
# }

# sc.qc.ngene <- function(sol, plot.dir=NULL, project.name) {
#   n.gene.df <- lapply(seq_along(sol), function(i) {
#     sol[[i]]@meta.data %>% mutate(rank.Feature = rank(nFeature_RNA)) %>% select(rank.Feature, nFeature_RNA) %>%
#       arrange(rank.Feature) %>%
#     return()
#   })

#   n.gene.plots <- lapply(n.gene.df, function(x) {
#     p <- x %>%
#       ggplot(aes(x=rank.Feature, y=(nFeature_RNA))) +
#       geom_point(size = 0.1) +
#       ylim(0, 10000) +
#       ggtitle(sample.info$sample[i])
#     return(p)
#   })

#   p.all <- cowplot::plot_grid(plotlist = n.gene.plots)
#   if (!is.null(plot.dir))
#     p.all + ggsave(paste0(plot.dir,"/",project.name,  "_gene_number_distro.png"), units = "in",
#            width = 10, height = 10, device = "png", dpi = 100)
#   return(list(sol, n.gene.df, p.all))
# }

# sc.qc.mito <- function(sol, plot.dir=NULL, project.name, elbow.mt.list = NULL) {
#   pct.mito.plots <- list()
#   pct.mito.df <- list()
#   for (i in seq_along(sol) ) {
#     pct.mito.df[[i]] <- sol[[i]]@meta.data %>% mutate(rank.Feature = rank(percent.mt))
#     if (is.null(elbow.mt.list[[sol[[i]]@meta.data$sample[1]]])) {
#       sol[[i]]@misc[["elbow.mt"]] <- sol[[i]]@meta.data %>%
#         mutate(x = rank(percent.mt), y=percent.mt) %>% select(x,y) %>% arrange(x) %>%
#         find.elbow()
#     } else {
#       sol[[i]]@misc[["elbow.mt"]] <- list(y=elbow.mt.list[[sol[[i]]@meta.data$sample[1]]])
#     }

#     p <- pct.mito.df[[i]] %>%
#       ggplot(aes(x=rank.Feature, y=(percent.mt))) +
#       geom_point(size = 0.1) +
#       ylim(0, 100) +
#       geom_hline(yintercept = sol[[i]]@misc[["elbow.mt"]]$y, color = "red") +
#       ggtitle(sample.info$sample[i])
#     pct.mito.plots[[i]] <- p
#   }

#   # plot mito
#   p.all <- cowplot::plot_grid(plotlist = pct.mito.plots)
#   if (!is.null(plot.dir))
#     p + ggsave(paste0(plot.dir,"/",project.name,  "_pct_mito.png"), units = "in",
#            width = 10, height = 10, device = "png", dpi = 100)
# }

# sc.qc.metadata.find.elbow.core <- function(sol, elbow.list=NULL, meta.term) {
#   meta.sum.df <- list()
#   meta.rank <- paste0("rank_", meta.term)
#   meta.elbow <- paste0("elbow_", meta.term)
#   for (i in seq_along(sol) ) {
#     # meta.sum.df[[i]] <- sol[[i]]@meta.data %>% mutate(rank.Feature = rank(percent.mt))

#     meta.sum.df[[i]] <- sol[[i]]@meta.data
#     meta.sum.df[[i]][, meta.rank] <- rank(meta.sum.df[[i]][, meta.term])
#     meta.sum.df[[i]] <- meta.sum.df[[i]][, c(meta.rank, meta.term)]
#     meta.sum.df[[i]] <- meta.sum.df[[i]][order(meta.sum.df[[i]][, meta.rank]),]
#     colnames(meta.sum.df[[i]]) <- c("x", "y")

#     if (is.null(elbow.list[[sol[[i]]@meta.data$sample[1]]])) {
#       sol[[i]]@misc[[meta.elbow]] <- meta.sum.df[[i]] %>%
#         find.elbow()
#     } else {
#       sol[[i]]@misc[[meta.elbow]] <- list(y=elbow.list[[sol[[i]]@meta.data$sample[1]]])
#     }
#   }
#   return(list(sol,meta.sum.df))

# }

# # t <- sol

# # tt <- sc.qc.metadata.find.elbow.core(t[1], meta.term = "percent.mt")
# # debugonce(find.elbow)
# # work.dir <- "int_4_27_20"
# # plot.dir <- "int_4_27_20/plots/"
# # project.name <- "test_int"
# # sample.tsv <- "int_4_27_20/kc_sample_info.tsv"
# #
# # system(paste0("mkdir -p ", plot.dir))
# #
# # sample.info <- read.table(sample.tsv, header = T, as.is = T, sep = "\t")

# # dir.list <- c("~/hmtp/kc/lodge/derek/Seq_200113/RAW/filtered_WT",
# #               "~/hmtp/kc/lodge/derek/Seq_200113/RAW/filtered_pymt+",
# #               "~/hmtp/kc/lodge/derek/Seq_200122/RAW/8wk_mtx/filtered/",
# #               "~/hmtp/kc/lodge/derek/Seq_200122/RAW/6mo_mtx/filtered/",
# #               "~/hmtp/kc/lodge/derek/Seq_200129/RAW/9wk_pyMT-/filtered/",
# #               "~/hmtp/kc/lodge/derek/Seq_200129/RAW/9wk_pyMT+/filtered/")
# #
# # sample.names <- c("WT_26w", "PyMT_24w", "PyMT_8w", "PyMT_25w", "WT_9w", "PyMT_9w")

# # so: seurat object; l: list
# # sol <- split(sample.info, f=sample.info$sample %>% factor(.,levels=.)) %>%
# #   lapply(function(x) {
# #     so <- Read10X(data.dir = x$dir) %>%
# #       CreateSeuratObject(project = "test_int")
# #     so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt")
# #     metadata <- colnames(x)[!colnames(x) %in% c("dir")]
# #     for (m in metadata) {
# #       so[[m]] <- x[1,m]
# #     }
# #     return(so)
# #   })
# #
# # # plot number of features
# # n.gene <- lapply(seq_along(sol), function(i) {
# #   p <- sol[[i]]@meta.data %>% mutate(rank.Feature = rank(nFeature_RNA)) %>%
# #     ggplot(aes(x=rank.Feature, y=(nFeature_RNA))) +
# #     geom_point(size = 0.1) +
# #     ylim(0, 10000) +
# #     ggtitle(sample.info$sample[i])
# #   return(p)
# # })
# #
# # cowplot::plot_grid(plotlist = n.gene) +
# #   ggsave(paste0(plot.dir,"/",project.name,  "_gene_number_distro.png"), units = "in",
# #          width = 10, height = 10, device = "png", dpi = 100)
# #
# # pct.mito.plots <- list()
# # for (i in seq_along(sol) ) {
# #   # sol[[i]]@misc[["elbow.mt"]] <- sol[[i]]@meta.data %>%
# #   #   mutate(x = rank(percent.mt), y=percent.mt) %>% select(x,y) %>%
# #   #   filter (y > 5, y < 12) %>% find.elbow()
# #   sol[[i]]@misc[["elbow.mt"]] <- list(y=8)
# #   p <- sol[[i]]@meta.data %>% mutate(rank.Feature = rank(percent.mt)) %>%
# #     ggplot(aes(x=rank.Feature, y=(percent.mt))) +
# #     geom_point(size = 0.1) +
# #     ylim(0, 100) +
# #     geom_hline(yintercept = sol[[i]]@misc[["elbow.mt"]]$y, color = "red") +
# #     ggtitle(sample.info$sample[i])
# #   pct.mito.plots[[i]] <- p
# # }
# #
# # # plot mito
# # cowplot::plot_grid(plotlist = pct.mito.plots) +
# #   ggsave(paste0(plot.dir,"/",project.name,  "_pct_mito.png"), units = "in",
# #          width = 10, height = 10, device = "png", dpi = 100)
# #
# # saveRDS(sol, paste0(work.dir, "/sol.Rds")

# # deriv.df$y <- tt[[2]][[1]]$y
# # which(deriv.df$y >= 4)[1] #1749

# # deriv.df %>% ggplot(aes(x=x, y = second.deriv)) +
# #   geom_point() +
# #   xlim(1700, 1800) +
# #   geom_vline(xintercept = 1749)

# # deriv.df %>% ggplot(aes(x=x, y = first.deriv)) +
# #   geom_point() +
# #   xlim(1700, 1800) +
# #   geom_vline(xintercept = 1749)

# # deriv.df %>% ggplot(aes(x=x, y = y)) +
# #   geom_point() +
# #   xlim(1700, 1800) +
# #   geom_vline(xintercept = 1749)
# # f.s <- smooth.spline(x = deriv.df$x, y = deriv.df$y)
# # deriv.df$y.s <- f.s(deriv.df$x)


# # data.frame(x=f.s$x, y = f.s$y) %>% ggplot(aes(x=x, y = y)) +
# #   geom_point() +
# #   # xlim(1700, 1800) +
# #   geom_vline(xintercept = 1749)

# # smooth.df <- data.frame(x=f.s$x, y = f.s$y)
# # smooth.df.80 <- smooth.df %>% filter(x > 0.8*max(x))

# # find.elbow(smooth.df)


# # smooth.deriv.df %>% ggplot(aes(x=x, y = second.deriv)) +
# #   geom_point() +
# # #  xlim(1700, 1800) +
# #   geom_vline(xintercept = 1749)

# # smooth.deriv.df %>% ggplot(aes(x=x, y = predict.2)) +
# #   geom_point() +
# #   #  xlim(1700, 1800) +
# #   geom_vline(xintercept = 1749)
# # smooth.deriv.df %>% ggplot(aes(x=x, y = first.deriv)) +
# #   geom_point() +
# # #  xlim(1700, 1800) +
# #   geom_vline(xintercept = 1749)

# # smooth.deriv.df %>% ggplot(aes(x=x, y = y)) +
# #   geom_point() +
# #   xlim(1700, 1800) +
# #   geom_vline(xintercept = 1749)

# # f.s$w %>% str()

# # t.predict <-  predict(f.s$fit, 1:1642, deriv = 2)
# # data.frame(x=t.predict$x, y = t.predict$y) %>% ggplot(aes(x=x, y = y)) +
# #   geom_point() +
# #   xlim(1700, 1800) +
# #   geom_vline(xintercept = 1749)
