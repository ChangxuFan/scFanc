ao.so.UMAP.cpm <- function(ao, so, ao.embedding = "UMAP", 
                           cell.list = NULL, plot.each.cell = F,
                           group.by, groups = NULL, split.by = NULL, splits = NULL, n.each = 20,
                           plot.out = NULL) {
  bIden <- identical(sort(ao$cellNames), sort(colnames(so)))
  if (!bIden) {
    stop("cells contained in ao and so are not identical")
  }
  if (is.null(cell.list)){
    cell.list <- get.cell.list(obj = so, is.ao = F, style = "ArchR",
                               n.cells.each = n.each, split.by = split.by, splits = splits,
                               group.by = group.by, groups = groups, return.named.vec = F)
  }
  aof <- fake.so.gen(ao = ao, ao.embedding = ao.embedding)
  pl <- lapply(names(cell.list), function(name) {
    if (plot.each.cell) {
      cells <- cell.list[[name]]
    } else {
      cells <- cell.list[name]
    }
    
    lapply(cells, function(cell) {
      lapply(c("so", "aof"), function(type) {
        if (type == "so") {
          obj <- so
        } else {
          obj <- aof
        }
        hl.list <- list(plotting = cell)
        p <- highlight.cells(so = obj, hl.list = hl.list,
                             plot.out = NULL, is.Rds = F, order = T, pt.size = 0.5,
                             return.list = T)[[1]][[1]] + 
          ggtitle(paste0(name, "\n", ifelse(plot.each.cell, paste0(cell, "\n"), ""), type))
        return(p)
      }) %>% return()
    }) %>% Reduce(c, .) %>% return()
  }) %>% Reduce(c, .) %>% return()
  
  p <- wrap.plots.fanc(plot.list = pl, n.split = 2, plot.out = plot.out)
  invisible(p)
}

ao.so.cluster.prime <- function(ao, so, rna.cluster.ident, atac.cluster.ident,
                               direction.df, out.file, plot.prime = F) {
  # df should be something like soi@meta.data
  # direction.df data.frame(direction = "ATAC->RNA", source = "C6|C8",
  # target = "9|1", primetarget = "0|2|7|8|11|13|14|17|21")
  # here we are trying to assess if there are cells that are HSPCs by ATAC identity
  # but erythroid by RNA identity
  so <- seurat.add.archr.meta(so = so, ao = ao, metas = atac.cluster.ident, overwrite = T)
  df <- so@meta.data
  if (!"cells" %in% colnames(df)) {
    df$cells <- rownames(df)
  }
  df <- utilsFanc::factor2character.fanc(df)
  fields.required <- c("sample", "cells", rna.cluster.ident, atac.cluster.ident)
  utilsFanc::check.intersect(fields.required, "requied fields", colnames(df), "colnames(df)")
  res <- direction.df %>% split(., f=  1:nrow(.)) %>% lapply(function(direc) {
    source.ident <- ifelse(direc$direction == "ATAC->RNA", atac.cluster.ident, rna.cluster.ident)
    source.clusters <- direc$source %>% strsplit("\\|") %>% unlist()
    target.ident <- ifelse(direc$direction == "ATAC->RNA", rna.cluster.ident, atac.cluster.ident)
    target.clusters <- direc$target %>% strsplit("\\|") %>% unlist()
    primetarget.clusters <- direc$primetarget %>% strsplit("\\|") %>% unlist()
    df <- df[df[, source.ident] %in% source.clusters, ]
    df$type <- NA
    df$type[df[, target.ident] %in% target.clusters] <- "target"
    df$type[df[, target.ident] %in% primetarget.clusters] <- "primetarget"
    
    res <- df %>% dplyr::group_by(sample) %>% 
      dplyr::summarise(n.source = n(), n.target = sum(type %in% "target"),
                       n.primetarget = sum(type %in% "primetarget"),
                       prime.cells = paste0(cells[type %in% "primetarget"], collapse = ",")) %>% 
      dplyr::ungroup() %>% as.data.frame()
    res <- cbind(direc, res)
    return(res)
  }) %>% do.call(rbind, .)
  dir.create(dirname(out.file), showWarnings = F, recursive = T)
  write.table(res, out.file, sep = "\t", col.names = T, row.names = F, quote = F)
  if (plot.prime) {
    cell.list <- res$prime.cells %>% strsplit(",")
    names(cell.list) <- paste0(res$sample, " ", res$source, "---", res$primetarget)
  
    plot.out <- paste0(tools::file_path_sans_ext(out.file), "_prime.png")
    ao.so.UMAP.cpm(ao = ao, so = so, cell.list = cell.list, plot.each.cell = F, 
                   plot.out = plot.out)
  }
  return(res)
}