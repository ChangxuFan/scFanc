VlnPlot.fanc <- function (object, features, cols = NULL, pt.size = 1, idents = NULL, 
                          sort = FALSE, assay = NULL, group.by = NULL, split.by = NULL, 
                          adjust = 1, y.max = NULL, same.y.lims = FALSE, log = FALSE, 
                          ncol = NULL, slot = "data", split.plot = FALSE, stack = FALSE, 
                          combine = TRUE, fill.by = "feature", flip = FALSE,
                          violin.scale = "width", no.noise = T) {
  # copied from Seurat
  if (!is.null(x = split.by) & getOption(x = "Seurat.warn.vlnplot.split", 
                                         default = TRUE)) {
    message("The default behaviour of split.by has changed.\n", 
            "Separate violin plots are now plotted side-by-side.\n", 
            "To restore the old behaviour of a single split violin,\n", 
            "set split.plot = TRUE.\n      \nThis message will be shown once per session.")
    options(Seurat.warn.vlnplot.split = FALSE)
  }
  return(ExIPlot.fanc(object = object, type = ifelse(test = split.plot, 
                                                     yes = "splitViolin", no = "violin"), features = features, 
                      idents = idents, ncol = ncol, sort = sort, assay = assay, 
                      y.max = y.max, same.y.lims = same.y.lims, adjust = adjust, 
                      pt.size = pt.size, cols = cols, group.by = group.by, 
                      split.by = split.by, log = log, slot = slot, stack = stack, 
                      combine = combine, fill.by = fill.by, flip = flip, 
                      violin.scale = violin.scale, no.noise = no.noise))
}


ExIPlot.fanc <- function (object, features, type = "violin", idents = NULL, ncol = NULL, 
                          sort = FALSE, assay = NULL, y.max = NULL, same.y.lims = FALSE, 
                          adjust = 1, cols = NULL, pt.size = 0, group.by = NULL, split.by = NULL, 
                          log = FALSE, slot = "data", stack = FALSE, combine = TRUE, 
                          fill.by = NULL, flip = FALSE, violin.scale = "width", no.noise = T) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  if (isTRUE(x = stack)) {
    if (!is.null(x = ncol)) {
      warning("'ncol' is ignored with 'stack' is TRUE", 
              call. = FALSE, immediate. = TRUE)
    }
    if (!is.null(x = y.max)) {
      warning("'y.max' is ignored when 'stack' is TRUE", 
              call. = FALSE, immediate. = TRUE)
    }
  }
  else {
    ncol <- ncol %||% ifelse(test = length(x = features) > 
                               9, yes = 4, no = min(length(x = features), 3))
  }
  data <- FetchData(object = object, vars = features, slot = slot)
  features <- colnames(x = data)
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  }
  else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% 
                                                 idents])
  }
  data <- data[cells, , drop = FALSE]
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  }
  else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  if (is.null(x = split.by)) {
    split <- NULL
  }
  else {
    split <- object[[split.by, drop = TRUE]][cells]
    if (!is.factor(x = split)) {
      split <- factor(x = split)
    }
    if (is.null(x = cols)) {
      cols <- scales::hue_pal()(length(x = levels(x = idents)))
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    else if (length(x = cols) == 1 && cols == "interaction") {
      split <- interaction(idents, split)
      cols <- scales::hue_pal()(length(x = levels(x = idents)))
    }
    else {
      cols <- Col2Hex(cols)
    }
    if (length(x = cols) < length(x = levels(x = split))) {
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
    names(x = cols) <- levels(x = split)
    if ((length(x = cols) > 2) & (type == "splitViolin")) {
      warning("Split violin is only supported for <3 groups, using multi-violin.")
      type <- "violin"
    }
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  if (isTRUE(x = stack)) {
    return(MultiExIPlot(type = type, data = data, idents = idents, 
                        split = split, sort = sort, same.y.lims = same.y.lims, 
                        adjust = adjust, cols = cols, pt.size = pt.size, shape = pt.shape, 
                        log = log, fill.by = fill.by, flip = flip))
  }
  plots <- lapply(X = features, FUN = function(x) {
    return(SingleExIPlot.fanc(type = type, data = data[, x, drop = FALSE], 
                              idents = idents, split = split, sort = sort, y.max = y.max, 
                              adjust = adjust, cols = cols, pt.size = pt.size, 
                              log = log, violin.scale = violin.scale, no.noise = no.noise))
  })
  label.fxn <- switch(EXPR = type, violin = if (stack) {
    xlab
  } else {
    ylab
  }, splitViolin = if (stack) {
    xlab
  } else {
    ylab
  }, ridge = xlab, stop("Unknown ExIPlot type ", type, call. = FALSE))
  for (i in 1:length(x = plots)) {
    key <- paste0(unlist(x = strsplit(x = features[i], split = "_"))[1], 
                  "_")
    obj <- names(x = which(x = Key(object = object) == key))
    if (length(x = obj) == 1) {
      if (inherits(x = object[[obj]], what = "DimReduc")) {
        plots[[i]] <- plots[[i]] + label.fxn(label = "Embeddings Value")
      }
      else if (inherits(x = object[[obj]], what = "Assay")) {
        next
      }
      else {
        warning("Unknown object type ", class(x = object), 
                immediate. = TRUE, call. = FALSE)
        plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
      }
    }
    else if (!features[i] %in% rownames(x = object)) {
      plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
    }
  }
  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = ncol)
    if (length(x = features) > 1) {
      plots <- plots & NoLegend()
    }
  }
  return(plots)
}


SingleExIPlot.fanc <- function (data, idents, split = NULL, type = "violin", sort = FALSE, 
                                y.max = NULL, adjust = 1, pt.size = 0, cols = NULL, seed.use = 42, 
                                violin.scale = "width", no.noise = T,
                                log = FALSE) {
  # copied from Seurat
  if (pt.size <- 0.1) {
    pt.shape <- 18
    line.size <- 0.1
    print("changing pt.shape to 18 for finer points")
  } else {
    pt.shape <- NULL
    line.size <- NULL
  }
  print("using SingleExIPlot.fanc")
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!is.data.frame(x = data) || ncol(x = data) != 1) {
    stop("'SingleExIPlot requires a data frame with 1 column")
  }
  feature <- colnames(x = data)
  data$ident <- idents
  if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
    data$ident <- factor(x = data$ident, levels = names(x = rev(x = sort(x = tapply(X = data[, feature], INDEX = data$ident, FUN = mean), decreasing = grepl(pattern = paste0("^", tolower(x = sort)), x = "decreasing")))))
  }
  if (log) {
    noise <- rnorm(n = length(x = data[, feature]))/200
    data[, feature] <- data[, feature] + 1
  }
  else {
    noise <- rnorm(n = length(x = data[, feature]))/100000
  }
  if (all(data[, feature] == data[, feature][1])) {
    warning(paste0("All cells have the same value of ", feature, 
                   "."))
  }
  else {
    if (!no.noise)
      data[, feature] <- data[, feature] + noise
  }
  axis.label <- "Expression Level"
  y.max <- y.max %||% max(data[, feature][is.finite(x = data[, 
                                                             feature])])
  if (type == "violin" && !is.null(x = split)) {
    data$split <- split
    vln.geom <- geom_violin
    fill <- "split"
  }
  else if (type == "splitViolin" && !is.null(x = split)) {
    data$split <- split
    vln.geom <- geom_split_violin
    fill <- "split"
    type <- "violin"
  }
  else {
    vln.geom <- geom_violin
    fill <- "ident"
  }
  switch(EXPR = type, violin = {
    x <- "ident"
    y <- paste0("`", feature, "`")
    xlab <- "Identity"
    ylab <- axis.label
    geom <- list(vln.geom(scale = violin.scale, adjust = adjust, size = line.size,
                          trim = TRUE), 
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    if (is.null(x = split)) {
      jitter <- geom_jitter(height = 0, size = pt.size, shape = pt.shape)
    } else {
      jitter <- geom_jitter(position = position_jitterdodge(jitter.width = 0.4, 
                                                            dodge.width = 0.9), size = pt.size, shape = pt.shape)
    }
    log.scale <- scale_y_log10()
    axis.scale <- ylim
  }, ridge = {
    x <- paste0("`", feature, "`")
    y <- "ident"
    xlab <- axis.label
    ylab <- "Identity"
    geom <- list(geom_density_ridges(scale = 4), theme_ridges(), 
                 scale_y_discrete(expand = c(0.01, 0)), scale_x_continuous(expand = c(0, 
                                                                                      0)))
    jitter <- geom_jitter(width = 0, size = pt.size, shape = pt.shape)
    log.scale <- scale_x_log10()
    axis.scale <- function(...) {
      invisible(x = NULL)
    }
  }, stop("Unknown plot type: ", type))

  plot <- ggplot(data = data, mapping = aes_string(x = x, y = y, fill = fill)[c(2, 3, 1)]) + 
    labs(x = xlab, y = ylab, 
         title = feature, fill = NULL) + 
    cowplot::theme_cowplot() + 
    theme(plot.title = element_text(hjust = 0.5))
  plot <- do.call(what = "+", args = list(plot, geom))
  plot <- plot + if (log) {
    log.scale
  }
  else {
    axis.scale(min(data[, feature]), y.max)
  }
  if (pt.size > 0) {
    plot <- plot + jitter
  }
  if (!is.null(x = cols)) {
    if (!is.null(x = split)) {
      idents <- unique(x = as.vector(x = data$ident))
      splits <- unique(x = as.vector(x = data$split))
      labels <- if (length(x = splits) == 2) {
        splits
      }
      else {
        unlist(x = lapply(X = idents, FUN = function(pattern, 
                                                     x) {
          x.mod <- gsub(pattern = paste0(pattern, "."), 
                        replacement = paste0(pattern, ": "), x = x, 
                        fixed = TRUE)
          x.keep <- grep(pattern = ": ", x = x.mod, fixed = TRUE)
          x.return <- x.mod[x.keep]
          names(x = x.return) <- x[x.keep]
          return(x.return)
        }, x = unique(x = as.vector(x = data$split))))
      }
      if (is.null(x = names(x = labels))) {
        names(x = labels) <- labels
      }
    }
    else {
      labels <- levels(x = droplevels(data$ident))
    }
    plot <- plot + scale_fill_manual(values = cols, labels = labels)
  }
  return(plot)
}