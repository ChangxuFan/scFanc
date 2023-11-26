cc.state.mat <- function(soi, cc.state.var, cc.state.order = NULL,
                         cluster.by, clusters = NULL, group.by = "sample", groups = NULL,
                         out.dir, root.name = NULL) {
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  cols <- c(cc.state.var, cluster.by, group.by)
  utilsFanc::check.intersect(x = cols, x.name = "required columns", 
                             y = colnames(soi@meta.data), y.name = "colnames(soi@meta.data)",
                             n.examples = 5, warning.only = F)
  df <- soi@meta.data[, cols]
  if (!is.null(clusters)) {
    df <- df[df[, cluster.by] %in% clusters,]
  }
  if (!is.null(groups)) {
    df <- df[df[, group.by] %in% groups,]
  }
  
  if (nrow(df) < 1) {
    stop("nrow(df) < 1")
  }
  
  res <- df %>% split(., f = factor(.[, cluster.by], unique(.[, cluster.by]))) %>% 
    lapply(function(df) {
      tabmat <- table(df[, c(group.by, cc.state.var)]) %>% as.matrix()
      pctmat <- round(diag(1/rowSums(tabmat)) %*% tabmat, digits = 3)
      rownames(pctmat) <- rownames(tabmat)
      col_fun <- circlize::colorRamp2(c(0, 1), c("white", "red"))
      if (is.null(cc.state.order)) {
        cc.state.order <- sort(unique(colnames(pctmat)))
      } else {
        utilsFanc::check.intersect(colnames(pctmat), "cc.states", 
                                   cc.state.order, "cc.state.order")
        cc.state.order <- cc.state.order %>% .[. %in% colnames(pctmat)]
      }
      p <- ComplexHeatmap::Heatmap(pctmat, col = col_fun, column_order = cc.state.order)
      save.base.plot(p = p, file = paste0(out.dir, "/", root.name, "_", df[1, cluster.by], ".png"), 
                     width = 500, height = 500)
      tab <- as.data.frame(pctmat)
      tab <- cbind(data.frame(V1 = rownames(tab)), tab)
      colnames(tab)[1] <- group.by
      rownames(tab) <- NULL
      return(tab)
    }) %>% do.call(rbind, .)
  write.table(res, paste0(out.dir, "/", root.name, ".tsv"))
  invisible(res)
}





cc.state <- function(soi, cluster.by, clusters = NULL, cc.state.var,
                     group.by = "sample", groups = NULL,
                     pct.counter.rule.list = NULL,
                     out.dir, root.name = NULL) {
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  dir.create(out.dir,recursive = T, showWarnings = F)
  df <- soi@meta.data[, c(group.by, cluster.by, cc.state.var)]
  colnames(df) <- c("group", "cluster", "state")
  if (!is.null(clusters)) {
    df <- df %>% filter(cluster %in% clusters)
  }
  if (!is.null(groups)) {
    df <- df %>% filter(group %in% groups)
  }
  
  df <- df %>% group_by(cluster,group, state) %>% 
    summarise(n = n()) %>% as.data.frame()
  
  write.table(df, paste0(out.dir, "/", root.name, "_raw_uncast.tsv"), 
              sep = "\t", quote = F, row.names = F, col.names = T)
  
  raw.df <- reshape2::dcast(data = df, formula = group + cluster ~ state, value.var = "n")
  raw.df <- raw.df[utilsFanc::sort.by(raw.df$group, groups, return.order = T),]
  rownames(raw.df) <- NULL
  write.table(raw.df, paste0(out.dir, "/", root.name, "_raw.tsv"), 
              sep = "\t", quote = F, row.names = F, col.names = T)
  res <- list(df = df, raw.df = raw.df)
  if (!is.null(pct.counter.rule.list)) {
    rule.dfs <- lapply(names(pct.counter.rule.list), function(rule.name) {
      rule <- pct.counter.rule.list[[rule.name]]
      df <- df %>% filter(cluster %in% rule)
      rule.df <- df %>% group_by(group) %>% 
        dplyr::mutate(pct = n/sum(n)) %>% ungroup() %>% as.data.frame() %>% 
        reshape2::dcast(formula = group + cluster ~ state, value.var = "pct")
      rule.df.summary <- rule.df %>% group_by(group) %>% 
        summarise_if(is.numeric, sum)
      rule.df.summary$cluster <- "total"
      rule.df <- rbind(rule.df, rule.df.summary) %>% 
        .[utilsFanc::sort.by(.$group, groups, return.order = T),]
      write.table(rule.df, paste0(out.dir, "/", root.name, "_pct_rule_", rule.name, ".tsv"),
                  sep = "\t", quote = F, row.names = F, col.names = T)
      return(rule.df)
    })
    names(rule.dfs) <- names(pct.counter.rule.list)
    res <- c(res, rule.dfs)
  }
  
  return(res)
}