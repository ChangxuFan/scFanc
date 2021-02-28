peak.corr.plots <- function(df, cmp.list, id, min=NULL, max=NULL, n = NULL, sync.axis=F,
                            x.lim = c(0, NA), y.lim = c(0, NA), plot.out, add.corr = F) {
  if (!is.null(min))
    df <- df[df[,id] > min,]
  if (!is.null(max))
    df <- df[df[,id] < max,]
  if (!is.null(n) && n < nrow(df)) {
    set.seed(seed = 42)
    df <- df[sample(1:nrow(df), n, replace = F),]
  }
  
  p.list <- lapply(seq_along(cmp.list), function(i) {
    cmp <- names(cmp.list)[i]
    x <- cmp.list[[i]]$x
    y <- cmp.list[[i]]$y
    p <- ggpubr::ggscatter(data = df, x = x, y = y, size = 0.3, ellipse.alpha = 0.5, color = "blue") 
    if (add.corr == T) {
      p <- p+ ggpubr::stat_cor(method = "pearson")
    }  
    
    if (sync.axis == T) {
      axis.max <- max(df[,x], df[,y])
      p <- p + xlim(0, axis.max) + ylim(0, axis.max)
    } else {
      p <- p+xlim(x.lim) +ylim(y.lim)
    }
    return(p)
  })
  
  p <- scFanc::wrap.plots.fanc(plot.list = p.list, plot.out = plot.out)
  return(p)
}

peak.corr.bin.plot <- function(df, cmp.list, id, min=NULL, max=NULL, window.size, n = NULL, plot.out) {
  if (!is.null(min))
    df <- df[df[,id] > min,]
  
  if (!is.null(max))
    df <- df[df[,id] < max,]
  if (!is.null(n) && n < nrow(df)) {
    set.seed(seed = 42)
    df <- df[sample(1:nrow(df), n, replace = F),]
  }
  
  # min <- min(df[, id])
  # max <- max(df[, id])
  
  size.corr <- lapply(seq_along(cmp.list), function(i) {
    cmp <- names(cmp.list)[i]
    x <- cmp.list[[i]]$x
    y <- cmp.list[[i]]$y
    split.vec <- split.by.window(df[, id], window.size = window.size)
    # the last 2 groups should be merged together:
    max.level <- max(split.vec)
    split.vec[split.vec == max.level] <- max.level - 1
    
    subset.to <- min(table(split.vec))

    size.corr <- df %>% split(., f = split.vec) %>% 
      lapply(function(dff) {
        set.seed(seed = 21)
        dff <- dff[sample(1:nrow(dff), subset.to, replace = F), ]
        size = max(dff[, id])
        corr = cor(x = dff[, x], y = dff[, y])
        dff2 <- dff[, c(x,y)]
        colnames(dff2) <- c("x", "y")
        # p <- ggplot(dff2, aes(x=x, y = y)) +
        #   geom_point() +
        #   xlab(x) + ylab(y) +
        #   ggtitle(size)
        # print(p)
        return(data.frame(size = size, corr = corr, cmp = cmp, n = nrow(dff)))
      }) %>% `names<-`(NULL) %>% Reduce(rbind, .)
    
    
    # p <- ggpubr::ggscatter(data = df, x = x, y = y, size = 0.3, ellipse.alpha = 0.5, color = "blue") +
    #   ggpubr::stat_cor(method = "pearson")
    return(size.corr)
  }) %>% Reduce(rbind, .) %>% na.omit()
  
  p <- ggline(data = size.corr, x = "size", y = "corr", color = "cmp")
  trash <- scFanc::wrap.plots.fanc(plot.list = list(p), plot.out = plot.out)
  print(p)
  return(size.corr)
  
}

split.by.window <- function(vec, window.size) {
  vec.shift <- vec - min(vec) + 1
  f <- ceiling(vec.shift/window.size)
  return(f)
}
