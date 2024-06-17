# #>>>>>>>> Old code
# lin.marker.gen <- function(so, markers=NULL, cluster.ident, lin.df, fg.pct.cutoff,  each.bg.pct, each.bg.logFC,
#                            bg.pct.cutoff, bg.logFC.cutoff, out.Rds = NULL, threads = 1, out.gmt = NULL) {
#   # the format of lin.df: name, fg, optional: fg.pct.cutoff, bg.pct, bg.pct.cutoff; bg.logFC, bg.logFC.cutoff.
#   # at least one of bg.pct and bg.logFC should be specified.
#   # in bg and fg, multiple clusters can be given, separated by ","
#   if (is.character(so))
#     so <- readRDS(so)
#   if (is.character(markers))
#     markers <- readRDS(markers)
#   markers$cluster <- as.character(markers$cluster)
#   fields <- c("fg.pct.cutoff","bg.pct", "bg.logFC", "bg.pct.cutoff", "bg.logFC.cutoff")
#   for (f in fields) {
#     if (! f %in% colnames(lin.df)) {
#       lin.df[, f] <- NA
#     }
#   }
#   # browser()
#   res <- lin.df %>% split(., f= factor(.$name, levels = unique(.$name))) %>% 
#     mclapply(function(ldf) {
#       df.list <- ldf %>% split(., f = factor(.$fg, levels = unique(.$fg))) %>% 
#         lapply(function(cdf) {
#           # browser()
#           fg <- cdf$fg %>% strsplit(",") %>% unlist()
#           # bg <- cdf$bg %>% strsplit(",") %>% unlist()
#           # if (is.null(markers)) {
#           #   stop("calling findmarkers on the fly is not yet implemented")
#           # }
#           if (is.na(cdf$bg.pct) && is.na(cdf$bg.logFC) ) {
#             stop("at least one of bg.pct and bg.logFC should be specified")
#           }
#           
#           if (!(is.na(cdf$bg.pct.cutoff) )) {
#             bg.pct.cutoff <- cdf$bg.pct.cutoff
#           }
#           if (!(is.na(cdf$fg.pct.cutoff) )) {
#             fg.pct.cutoff <- cdf$fg.pct.cutoff
#           }
#           
#           if (!(is.na(cdf$bg.logFC.cutoff) )) {
#             bg.logFC.cutoff.cutoff <- cdf$bg.logFC.cutoff
#           }
#           
#           m.cl <- markers %>% filter(cluster %in% fg, avg_logFC > 0)
#           
#           
#           m.cl$pass <- T # just to make the "pass" column show up earlier
#           
#           pass.vec <- rep(T, nrow(m.cl))
#           if (!is.na(cdf$bg.pct)) {
#             m.cl$bg.pct.clusters <- cdf$bg.pct
#             bg <- cdf$bg.pct %>% strsplit(",") %>% unlist()
#             pct <- pct.detect(so = so, cluster.ident = cluster.ident, clusters = bg, genes = m.cl$gene)
#             m.cl[, "bg.pct"] <-  pct
#             pass.vec <- pass.vec & (pct <= bg.pct.cutoff)
#             
#             if (each.bg.pct == T) {
#               for (bg.sub in bg) {
#                 pct <- pct.detect(so = so, cluster.ident = cluster.ident, clusters = bg.sub,  genes = m.cl$gene)
#                 m.cl[, paste0("bg_", bg.sub, ".pct")] <-  pct
#                 pass.vec <- pass.vec & (pct <= bg.pct.cutoff)
#               }
#             } 
#           }
#           
#           
#           if (!is.na(cdf$bg.logFC)) {
#             # browser()
#             m.cl$bg.logFC.clusters <- cdf$bg.logFC
#             bg <- cdf$bg.logFC %>% strsplit(",") %>% unlist()
#             
#             logFC <- lapply(list(fg,bg), function(clusters) {
#               mean <- so.2.bulk(so = so, assay = "SCT", slot = "data", genes = m.cl$gene,
#                                 ident = cluster.ident, sub.idents = clusters, 
#                                 take.mean = T, coldata.columns = NULL)$bulk.mat[, 1]
#               # print(head(mean))
#               return(mean)
#             }) %>% Reduce(`/`, .) %>% log()
#             m.cl$bg.logFC <- logFC
#             pass.vec <- pass.vec & (logFC >= bg.logFC.cutoff)
#             if (each.bg.logFC == T) {
#               for (bg.sub in bg) {
#                 logFC <- lapply(list(fg,bg.sub), function(clusters) {
#                   mean <- so.2.bulk(so = so, assay = "SCT", slot = "data", genes = m.cl$gene,
#                                     ident = cluster.ident, sub.idents = clusters, 
#                                     take.mean = T, coldata.columns = NULL)$bulk.mat[, 1]
#                   # print(head(mean))
#                   return(mean)
#                 }) %>% Reduce(`/`, .) %>% log()
#                 m.cl[, paste0("bg_", bg.sub, ".logFC")] <- logFC
#                 pass.vec <- pass.vec & (logFC >= bg.logFC.cutoff)
#               }
#             }
#           }
# 
#           m.cl$pass <- pass.vec & (m.cl$pct.1 > fg.pct.cutoff)
#           return(m.cl)
#         }) 
#       genes.list <- lapply(df.list, function(df) return(df[df$pass == T,"gene"] %>% unique()))
#       names(genes.list) <- names(df.list)
#       genes <- unlist(genes.list) %>% unique()
#       return(list(df.list = df.list, genes.list = genes.list, genes = genes))
#     }, mc.cores = threads, mc.cleanup = T)
#   if (!is.null(out.Rds))
#     saveRDS(res, out.Rds)
#   if (!is.null(out.gmt)) {
#     marker.2.gmt(lin.marker.o = res, out.gmt = out.gmt)
#   }
#   return(res)
# }
# 
# 
# lin.marker.diag.featureplot <- function(so, lin.marker.o, out.dir, threads.master = 1, pass.only = T, use.pdf = F,
#                                         threads.sub = 10, out.root.name = NULL, split.by = NULL, n.split = NULL, ...) {
#   system(paste0("mkdir -p ", out.dir))
#   if (is.null(out.root.name))
#     out.root.name <- "/"
#   else
#     out.root.name <- paste0(out.root.name, "..")
#   mclapply(names(lin.marker.o), function(lin) {
#     df <- lin.marker.o[[lin]]$df.list %>% Reduce(rbind, .) 
#     df <- df[!duplicated(df$gene), ]
#     df$status <- NA
#     df$status[df$pass == T] <- "pass"
#     df$status[df$pass == F] <- "fail"
#     # browser()
#     df %>% split(f= factor(df$status, levels = c("pass", "fail"))) %>% lapply(function(df.sub) {
#       status <- df.sub$status[1]
#       if (pass.only == T && status == "fail") {
#         return()
#       }
#       p.list <- plot.panel.list(panel.list = df.sub$gene, obj = so, order = F, raster = T, return.list = T,
#                                 assay = "SCT", threads = threads.sub, split.by = split.by, ...)
#       p.list <- lapply(1:nrow(df.sub), function(i) {
#         p.sub.list <- lapply(p.list[[1]], function(p) {
#           title <- paste0(df.sub[i, "cluster"], "  ", df.sub[i, "gene"], "  ", df.sub[i, "pass"], "  ", format(df.sub[i, "pct.1"], digits = 2))
#           if (!is.null(df.sub$bg.pct))
#             title <- paste0(title, "  bg.pct: ", format(df.sub[i, "bg.pct"], digits = 2))
#           if (!is.null(df.sub$bg.logFC))
#             title <- paste0(title, "  bg.logFC: ", format(df.sub[i, "bg.logFC"], digits = 2))
#           p <- p + ggtitle(title)
#           return(p)
#         })
#         return(p.sub.list)
#       }) %>% Reduce(c, .)
#       # browser()
#       format <- "png"
#       if (use.pdf == T) {
#         format <- "pdf"
#       }
#       trash <- wrap.plots.fanc(plot.list = p.list, plot.out = paste0(out.dir, "/", out.root.name, lin, "..",status,".", format),
#                                n.split = n.split)
#     })
# 
#     
#     return(NULL)
#   }, mc.cores = threads.master, mc.cleanup = T) 
# }
# 
# lin.marker.diag.hm <- function(lin.marker, cluster.ident, bulk.list, sample, out.file) {
#   genes <- lapply(lin.marker, function(x) return(x$genes)) %>% unlist()
#   bulk.mat <- lapply(bulk.list, function(x) {
#     exp <- x$bulkNorm %>% filter(gene %in% genes) %>% pull(!!as.name(paste0("bulkNorm_", sample)))
#     return(exp)
#   } ) %>% Reduce(cbind, )
#   colnames(bulk.mat) <- names(bulk.list) %>% sub(cluster.ident, "",.) %>% sub("^_", "",.)
#   rownames(bulk.mat) <- genes
#   png(out.file)
#   try(print(complexHeatmap::Heatmap(bulk.mat)))
#   dev.off()
#   return()
# }
# 
# marker.2.gmt <- function(marker.list, lin.marker.o=NULL, out.gmt, gmt.name) {
#   if (!is.null(lin.marker.o)) {
#     marker.list <- lapply(lin.marker.o, function(x) return(x$genes))
#     names(marker.list) <- lin.marker.o %>% names()
#   }
#     
#   outs <- lapply(names(marker.list), function(marker.name) {
#     markers <- marker.list[[marker.name]] %>% toupper()
#     out <- c(marker.name, utilsFanc::bash2ftp(out.gmt), markers)
#     out <- paste0(out, collapse = "\t")
#     return(out)
#   }) %>% unlist()
#   system(paste0("mkdir -p ", dirname(out.gmt)))
#   write(outs, out.gmt, sep = "\n")
#   names(out.gmt) <- gmt.name
#   return(out.gmt)
# }
# 
# lin.marker.gen.ez <- function(so, idents.list, group.by, assay, ident, subset.ident = NULL,
#                               p.adj.cutoff, logfc.cutoff, ...) {
#   # example of idents.list: list(granu = list(ident.1 = "0", ident.2 = "15"))
#   DefaultAssay(so) <- assay
#   Idents(so) <- ident
#   markers.list <- lapply(idents.list, function(x) {
#     m <- FindMarkers(so, ident.1 = x$ident.1, ident.2 = x$ident.2,
#                      group.by = group.by,
#                      subset.ident = subset.ident, ...)
#     m <- utilsFanc::add.column.fanc(df1 = m, df2 = data.frame(gene = rownames(m)), pos = 1)
#     m <- m %>% filter(p_val_adj < p.adj.cutoff, abs(avg_logFC) > logfc.cutoff)
#     res <- list()
#     res$up <- m %>% filter(avg_logFC > 0)
#     res$down <- m %>% filter(avg_logFC < 0)
#     return(res)
#   })
#   return(markers.list)
# }
# #<<<<<<<< END: Old code

