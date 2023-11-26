adt.norm <- function(so, assay, slot, norm.method, plot.out = NULL) {
  so <- NormalizeData(object = so, assay = assay, normalization.method = norm.method)
  if (slot == "scale.data") {
    so <- Seurat::ScaleData(so, assay = assay)
  }
  features <- rownames(GetAssayData(so, slot = slot, assay = assay))
  if (!is.null(plot.out)) {
    pl <- lapply(features, function(feature) {
      df <- data.frame(value = GetAssayData(object = so, slot = slot, assay = assay)[feature, ])
      p <- ggplot(df, aes(x = value)) +
        geom_density() +
        ggtitle(paste0(feature, " ", assay, " ", slot))
      return(p)
    })
    wrap.plots.fanc(plot.list = pl, plot.out = plot.out)
  }
  
  return(so)
}

feature.scatter <- function(so, x, y, 
                            x.assay, x.slot, y.assay, y.slot,
                            gate.file = NULL, x.gate = NULL, y.gate = NULL,
                            split.by = NULL, splits = NULL, 
                            color.by = NULL, colors = NULL,
                            pt.size = 0.5, gate.line.size = 0.1,
                            plot.out) {
  # gate.file: tsv with columns: sample, x, y. Here x and y must be the actual names
  # x.gate and y.gate: either 1 single scalar value or a named vector for each split
  # note sample is usually the split.by that I use. However, even if split.by changes, 
  # the name of the first column in gate.file should still be "sample"
  
  meta.add <- c("nFeature_RNA", "percent.mt") # not yet useful
  
  if (missing(plot.out)) {
    stop("plot.out must be specified!")
  }
  root <- tools::file_path_sans_ext(plot.out)
  t.dir <- root
  so@meta.data$cells <- rownames(so@meta.data)
  df <- so@meta.data[, c("cells", meta.add, split.by, color.by), drop = F]
  if (!is.null(splits)) {
    splits <- splits %>% .[. %in% df[, split.by]] %>% unique()
    if (length(splits) < 1) stop("none of the splits is found")
    df <- df[df[, split.by] %in% splits,]
  } else {
    split.by <- "pd.split"
    splits <- "no.split"
    df$pd.split <- "no.split"
  }
  if (!is.null(colors)) {
    colors <- colors %>% .[. %in% df[, color.by]] %>% unique()
    if (length(colors) < 1) stop("none of the colors is found")
    df <- df[df[, color.by] %in% colors,]
  }
  
  if (!is.null(gate.file)) {
    gates <- read.table(gate.file, header = T)
    X <- x %>% gsub("\\-", ".", x)
    Y <- x %>% gsub("\\-", ".", y)
    cols <- c("sample", X, Y)
    n.f <- cols %>% .[!.%in% colnames(gates)]
    if (length(n.f) > 0) {
      stop(paste0("required columns not found in gate.file: ",
                  paste0(n.f, collapse = ", ")))
    }
    n.f <- splits %>% .[!. %in% gates$sample]
    if (length(n.f) > 0) {
      stop(paste0("some splits not found in the 'sample' column of gate.file: ",
                  paste0(n.f, collapse = ", ")))
    }
    dups <- gates$sample %>% .[duplicated(.)]
    if (length(dups) > 0) {
      stop(paste0("some elements in the 'sample' column of gate.file are duplicated: ",
                  paste0(dups, collapse = ", ")))
    }
    x.gate <- gates[, X]
    names(x.gate) <- gates[, "sample"]
    y.gate <- gates[, Y]
    names(y.gate) <- gates[, "sample"]
    
  }
  
  if (length(x.gate) == 1) {
    x.gate <- rep(x.gate, length(splits))
    names(x.gate) <- splits
  }
  
  if (length(y.gate) == 1) {
    y.gate <- rep(y.gate, length(splits))
    names(y.gate) <- splits
  }
  
  
  
  df$x <- GetAssayData(so, assay = x.assay, slot = x.slot)[x, df$cells]
  df$y <- GetAssayData(so, assay = y.assay, slot = y.slot)[y, df$cells]
  
  pl <- df %>% split(., f = factor(.[, split.by], levels = splits)) %>% 
    lapply(function(df) {
      split.name <- df[, split.by][1]
      p <- xy.plot(df = df, x = "x", y = "y", 
                   x.lab = paste0(x, " ", x.assay, " ", x.slot),
                   y.lab = paste0(y, " ", y.assay, " ", y.slot),
                   add.abline = F, pt.size = pt.size, color.var = color.by, 
                   show.color.var = T,
                   title = paste0(split.name))
      if (!is.null(x.gate))
        p <- p + geom_vline(xintercept = x.gate[split.name], linetype = "dashed", size = gate.line.size)
      if (!is.null(y.gate))
        p <- p + geom_hline(yintercept = y.gate[split.name], linetype = "dashed", size = gate.line.size)
      
      res <- list(p)
      if (!is.null(x.gate) && !is.null(y.gate)) {
        x.gate <- x.gate[split.name]
        y.gate <- y.gate[split.name]
        df$x.hi <- df$x > x.gate
        df$y.hi <- df$y > y.gate
        df2 <- df %>% dplyr::group_by(x.hi, y.hi, !!as.name(split.by)) %>% 
          dplyr::summarise(n = n()) %>% dplyr::ungroup() %>% 
          as.data.frame() %>% dplyr::mutate(pct = round(n/sum(n), digits = 3))
        
        df2[, paste0("x.", x, ".gate")] <- x.gate
        df2[, paste0("y.", y, ".gate")] <- y.gate
        
        t.df <- df2 %>% filter(x.hi == T)
        df2$y.hi.in.x.hi <- round(t.df[t.df$y.hi, "n"]/sum(t.df$n), digits = 3)
        t.df <- df2 %>% filter(y.hi == T)
        df2$x.hi.in.y.hi <- round(t.df[t.df$x.hi, "n"]/sum(t.df$n), digits = 3)
        
        if (nrow(df2) != 4) {
          stop("nrow(df2) != 4")
        }
        t.tsv <- paste0(t.dir, "/", split.name, "_gates.tsv")
        dir.create(t.dir, showWarnings = F, recursive = T)
        write.table(df2, t.tsv, sep = "\t", col.names = T, row.names = F, quote = F)
        
        p.hm <- ggplot(df2, aes(x = x.hi, y = y.hi, fill = pct)) + 
          geom_tile() + scale_fill_gradient(low = "white", high = "red") +
          geom_text(aes(label = pct)) + 
          theme(aspect.ratio = 1)
        res <- c(res, list(p.hm))
      }
      return(res)
    }) %>% Reduce(c, .)
  
  trash <- wrap.plots.fanc(plot.list = pl, n.split = 2, plot.out = plot.out)
  tsv.out <- paste0(root, "_gate.tsv")
  gate.df <- lapply(splits, function(x) {
    read.table(paste0(t.dir, "/", x, "_gates.tsv"), header = T)
  }) %>% do.call(rbind, .)
  system(paste0("rm -rf ", t.dir))
  write.table(gate.df, tsv.out, sep = "\t", col.names = T, 
              row.names = F, quote = F)
  return()
}
