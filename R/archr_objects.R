seurat.name.2.archr <- function(so, cells) {
  meta.df <- so@meta.data %>% factor2character.fanc()
  meta.df <- meta.df[cells, ]
  paste0(meta.df$sample, "#", sub("_.+$", "", rownames(meta.df))) %>% return()
}

chr.sort2mixsort <- function(x) {
  lengths <- x %>% .[which(.==1)-1] 
  lengths <- c(lengths, x[length(x)])
  chr.names <- paste0("chr", c(1, 10:19,2:9, 20:length(lengths)))
  chr.names[length(chr.names)] <- "chrX"
  chr <- mapply(function(name, rep) {
    return(rep(name, rep))
  }, chr.names, lengths, SIMPLIFY = F) %>% unlist() %>% `names<-`(NULL)
  
  order <- gtools::mixedorder(chr)
  return(order)
  
}
archr.get.frag <- function(obj, cellNames, return.df = F, out.bed = NULL, return.grange = F, ...) {
  rhdf5::h5disableFileLocking()
  arrow.file <- ArchR::getArrowFiles(obj)
  
  df <- mclapply(seq_along(arrow.file), function(i) {
    sample <- names(arrow.file)[i]
    file <- arrow.file[i]
    if (sum(grepl(paste0(sample, "#"), cellNames)) > 0) {
      cells <- cellNames[grepl(paste0(sample, "#"), cellNames)]
      frags <- getFragmentsFromArrow(ArrowFile = file, 
                                      cellNames = cells, 
                                      logFile = NULL, ...)
      df <- as.data.frame(frags)
      return(df)
    }
    return(NULL)
  }, mc.cores = length(arrow.file)) %>% Reduce(rbind, .)
  
  
  # frags <- getFragmentsFromArrow(ArrowFile = arrow.file, 
  #                                cellNames = cellNames, 
  #                                logFile = NULL, ...) 
  # df <- as.data.frame(frags)
  
  # this step is basically to use the 2 ends of a fragment as 2 separate rows (target pipe)
  # according to 10x, the original fragment file is already adjusted for cut site.
  # importantly, archr adjusted the ranges to convert the 0 based 10x BED-like fragment file to 
  ##grange. I adjusted it back. 
  # getFragmentsFromArrow(ArrowFile = "~/hmtp/scAR/test_run/sth/archr_int/WT.arrow", 
  #                       cellNames = "WT#GTGCACGGTGAAACAA-1" ,
  #                       logFile = NULL, ...)
  # (jupyter) cfan@stout:~/hmtp/scAR/test_run/sth/count/WT/outs$ zcat atac_fragments.tsv.gz | grep "GTGCACGGTGAAACAA-1" | less
  
  df.insert <- df %>% mutate(pseudo.umi = 1:nrow(df), start = start-1, end = end + 1) %>% 
    reshape2::melt(id.vars = c("seqnames", "RG", "pseudo.umi"),
                   measure.vars = c("start", "end"),
                   variable.name = ("start.end"),
                   value.name = "left") %>% 
    mutate(right = left + 1, strand = "+") %>% 
    dplyr::select(seqnames, left, right, RG, pseudo.umi, strand) %>% 
    arrange(RG, seqnames, left)
  if (!is.null(out.bed)) {
    system(paste0("mkdir -p ", dirname(out.bed)))
    write.table(df.insert, out.bed, sep = "\t", quote = F, col.names = F, row.names = F)
  }
  # if (return.grange==T)
  #   return(frags)
  if (return.df == T)
    return(list(df=df, df.insert = df.insert))
  return(NULL)
}


serial.macs2 <- function(obj, cellNames.list, sub.names= NULL,master.dir, root.name  = NULL, genome,
                         get.bed = T, run.macs2 = T,
                         thread.master, thread.sub=2, spmr = F, ...) {
  rhdf5::h5disableFileLocking()
  if (is.null(root.name)) 
    root.name <- names(cellNames.list)
  
  if (length(root.name) == 1)
    root.name <- rep(root.name, length(cellNames.list))
  
  if (grepl("mm", genome))
    genome <- "mm"
  if (grepl("hg", genome))
    genome <- "hs"
  # first get all the bed 
  if (!is.null(sub.names)) {
    if (length(cellNames.list) != length(sub.names))
      stop("sub.names, if provided, must be of the same length as cellNames.list")
    
  } else {
    lengths <- sapply(cellNames.list, length)
    sub.n <- rep(1, length(cellNames.list))
    
    if (length(sub.n) > 1) {
      for (i in 2:length(sub.n)) {
        if (lengths[i] == lengths[i-1])
          sub.n[i] <-  sub.n[i-1]+1
      }
    }

    
    sub.names <- paste0(lengths, "_rep_", sub.n)
  }
  system(paste0("mkdir -p ", master.dir))
  
  
  jsongen <- mclapply(seq_along(cellNames.list), function(i) {
    x <- cellNames.list[[i]]
    sub.name <- sub.names[i]
    sub.root.name <- paste0(root.name[i], "_", sub.name)
    sub.dir <- paste0(master.dir, "/", sub.root.name, "/")
    raw.bed <- paste0(sub.dir, "/", sub.root.name, ".bed")
    if (get.bed == T)
      archr.get.frag(obj = obj, cellNames = x, out.bed = raw.bed)
    if (run.macs2 == T)
      rnaSeqFanc::macs2.atac.callpeak(infile = raw.bed, root.name = sub.root.name, 
                                      outdir = sub.dir, format = "BED", genome = genome, 
                                      q.cutoff = 0.1, shift = -75, ext = 150, spmr = spmr, subsummit = F, 
                                      bdgcmp = T, bdgcmp.method = "ppois", zip = T, 
                                      write.log = T,
                                      thread = thread.sub, run=T)
    jsongen <- data.frame(name = paste0(sub.root.name, c("_treat_pileup.bdg.gz", "_control_lambda.bdg.gz",
                                                         "_ppois.bdg.gz", "_peaks.narrowPeak.gz")),
                          type = c("bedGraph", "bedGraph", "bedGraph", "bed"))
    jsongen$url <- jsongen$name %>% paste0(sub.root.name, "/", .)
    return(jsongen)
  }, mc.cores = thread.master) %>% Reduce(rbind, .)
  json <- jsongen %>% jsonlite::toJSON() %>% jsonlite::prettify()
  write(json, paste0(master.dir, "/macs_tracks.json"))
  
  narrowpeaks <- lapply(jsongen$url, function(x) {
    read.table(x %>% sub(".gz", "", .) %>% paste0(master.dir, "/",.), as.is=T) %>% 
      return()
  })
  
  names(narrowpeaks) <- jsongen$name
  return(narrowpeaks)
}

seurat.plot.archr <- function(so, ao, ao.embedding, ao.colorBy, ao.name, plot.dir, 
                              width = 10, height = 10, ...) {
  
  embedding <- ArchR::getEmbedding(ArchRProj = ao, embedding = ao.embedding, returnDF = T)
  # rownames(embedding) <- sub(".+#", "", rownames(embedding))
  cells <- rownames(embedding)
  
  look.up.table <- so@meta.data %>% mutate(cell = rownames(so@meta.data)) %>% select(cell, sample) %>% 
    mutate (cell.ori = sub("_.+$", "", cell)) %>% 
    mutate(new.name = paste0(sample, "#", cell.ori))
  rownames(look.up.table) <- look.up.table$new.name
  embedding <- embedding[look.up.table$new.name,]
  
  colnames(embedding) <- paste0(ao.embedding, "_", c(1,2))

  so@reductions <- lapply(so@reductions, function(x) {
    x@key <- "miao"
    return(x)
  })
  
  rownames(embedding) <- colnames(so)
  
  so[["pseudo"]] <- CreateDimReducObject(embeddings = as.matrix(embedding),
                                         key = paste0(ao.embedding, "_"),
                                         assay = DefaultAssay(so))
  
  so[["to_plot"]] <- ArchR::getCellColData(ArchRProj = ao, select = ao.name, drop = T)
  
  DimPlot(so, reduction = "pseudo", pt.size = 0.5, label = T, label.size = 10, ...) + 
    ggsave(paste0(plot.dir, "/ao_", ao.embedding, "_", ao.name, ".png"), 
           device = "png", width = width, height = height, dpi = 200, units = "in")
} 


seurat.plot.archr.2 <- function(so, ao, ao.embedding, ao.plot.list, ao.name, plot.dir, 
                              width = 10, height = 10, ...) {
  
 so <- seurat.add.archr.embed(so = so, ao = ao, ao.embedding = ao.embedding)
 
 
} 
archr.get.peak.mat <- function(ao, mat.name = "PeakMatrix") {
  #stop("untested")
  mat <- ArchR::getMatrixFromProject(ao, useMatrix = mat.name)
  # stop("need to add sort2mixedsort")
  
  order <- chr.sort2mixsort(x = rowData(mat)$idx)
  mat <- mat[order, ]
  mat <- assay(mat)
  return(mat)
}

# t <- ArchR::getCellColData(ao.int.prog,drop = F)
archr.get.cells.grid <- function(ao, cluster.ident, clusters=NULL, group.ident, groups = NULL, melt =T) {
  df <- ArchR::getCellColData(ao,drop = F, select = c(cluster.ident, group.ident)) %>% 
    as.data.frame() %>%  factor2character.fanc() 
  if (!is.null(clusters))
    df <- df[df[, cluster.ident] %in% clusters,]
  if (!is.null(groups))
    df <- df[df[, group.ident] %in% groups,]
  
  cell.list <- df %>% split(., f = factor(.[,cluster.ident], levels = gtools::mixedsort(unique(.[,cluster.ident])))) %>% 
    lapply(function(cl.df) {
      cl <- cl.df[, cluster.ident][1]
      cell.list.cl <- cl.df %>% split(., f = factor(.[,group.ident], levels = gtools::mixedsort(unique(.[,group.ident])))) %>% 
        lapply(function(group.df) {
          cells <- rownames(group.df)
          return(cells)
        })
      names(cell.list.cl) <- paste0(cluster.ident, "_",cl, "..", names(cell.list.cl))
      return(cell.list.cl)
    })
  
  if (melt == T) {
    cell.list <- Reduce(c, cell.list)
  }
  return(cell.list)
}

vectorization.core <- function(ao, cell.list=NULL,cluster.ident, clusters=NULL, group.ident, groups = NULL, mat= NULL,
                               mat.name, deseq2.norm = F) {
  if (is.null(mat)) {
    mat <- ArchR::getMatrixFromProject(ao, useMatrix = mat.name)
    # stop("need to add sort2mixedsort")
    
    order <- chr.sort2mixsort(x = rowData(mat)$idx)
    mat <- mat[order, ]
    mat <- assay(mat)
  }
  
  if (mat.name == "PeakMatrix") {
    gr <- ArchR::getPeakSet(ao)
  } else {
    stop("the rest has not been developed")
  }
  if (is.null(cell.list))
    cell.list <- archr.get.cells.grid(ao = ao, cluster.ident = cluster.ident, clusters = clusters,
                                      group.ident = group.ident, groups = groups, melt = T)
  
  
  for (i in 1:length(cell.list)) {
    cells <- cell.list[[i]]
    name <- names(cell.list)[i]
    vec <- mat[,cells] %>% Matrix::rowSums()
    gr <- gr %>% plyranges::mutate(!!name := vec)
  }
  # names(gr) <- NULL
  # df <- as.data.frame(gr)
  # gr.bk <- gr
  if (deseq2.norm == T) {
    exp.df <- gr %>%  `names<-`(NULL) %>% as.data.frame() %>%
      mutate(peak = paste0(seqnames, ":", start, "-", end)) %>% 
      .[, c("peak", names(cell.list))]
    rownames(exp.df) <- exp.df$peak
    exp.df$peak <- NULL
    coldata <- data.frame(sample = colnames(exp.df))
    rownames(coldata) <- coldata$sample

    exp.mat <- as.matrix(exp.df)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = exp.mat, colData = coldata, design = ~sample)
    dds <- DESeq2::estimateSizeFactors(dds, type = "ratio")
    
    counts <- SummarizedExperiment::assay(dds) %>% as.matrix()
    norm.counts <- counts %*% diag((1/dds@colData$sizeFactor))
    norm.df <- as.data.frame(norm.counts)
    colnames(norm.df) <- as.character(dds@colData$sample) 
    
    for (i in colnames(norm.df)) {
      gr <- gr %>% plyranges::mutate(!!paste0(i, "_raw") := !!as.name(i)) %>% 
        plyranges::mutate(!!i := norm.df[, i])
    }
  }
  
  return(gr)
  
}

vectorization.m <- function(ao,cluster.ident, clusters=NULL, group.ident, groups = NULL, mat= NULL,
                            mat.name, deseq2.norm = F) {
  if (is.null(mat)) {
    mat <- ArchR::getMatrixFromProject(ao, useMatrix = mat.name)
    # stop("need to add sort2mixedsort")
    
    order <- chr.sort2mixsort(x = rowData(mat)$idx)
    mat <- mat[order, ]
    mat <- assay(mat)
  }
  

  cell.list <- archr.get.cells.grid(ao = ao, cluster.ident = cluster.ident,
                                    group.ident = group.ident, melt = F)
  gr.list <- lapply(seq_along(cell.list), function(i) {
    # cluster <- names(cell.list)[i]
    cells <- cell.list[[i]]
    gr <- vectorization.core(ao = ao, cell.list = cells, mat = mat, mat.name = mat.name,
                             deseq2.norm = deseq2.norm)
    return(gr)
  })
  clusters <- names(cell.list)
  names(gr.list) <- paste0(cluster.ident, "_", clusters)
  return(gr.list)
}



per.cluster.macs2 <- function(so, ao, meta, sort=F, subset = NULL, group.by=NULL, groups = NULL, spmr,get.bed = T,
                              master.dir, run = T, genome, thread.master, thread.sub) {
  df <- so@meta.data[, c(meta, "sample", group.by) %>% unique(), drop = F] %>% factor2character.fanc() 
  
  if (!is.null(subset))
    df <- df %>% .[.[,meta] %in% subset,]
  df$cells <- rownames(df) %>% sub("_.*$", "", .)
  df <- df %>% mutate(cells = paste0(sample, "#", cells))
  if (meta == "seurat_clusters")
    df$seurat_clusters <- paste0("cluster_", df$seurat_clusters)
  if (sort == T) {
    df <- df[gtools::mixedorder(df[, meta]), ]
  }
  
  cell.list <- df %>% split(., f=factor(.[,meta], levels = unique(.[,meta]))) %>% 
    lapply(function(x) { 
      meta.value <- x[, meta][1]

      if (!is.null(group.by)) {
        if (!is.null(groups)) {
          x <- x[x[, group.by] %in% groups,]
        }
        cells.by.group <- x %>% split(., f = factor(.[,group.by], levels = unique(sort(.[,group.by])))) %>% 
          lapply(function(y) {
            return(y[,"cells"])
          })
        names(cells.by.group) <- paste0(meta.value, "_",names(cells.by.group))
      } else {
        cells.by.group <- list(x[, "cells"])
        names(cells.by.group) <- meta.value
      }

      
      return(cells.by.group)
    })
  cell.list <- Reduce(c, cell.list)
  
  trash <- serial.macs2(obj = ao, cellNames.list = cell.list, genome = genome,
                        master.dir = master.dir, root.name = names(cell.list), 
                        get.bed = get.bed, thread.master = thread.master, thread.sub = thread.sub, spmr = spmr, run.macs2 = run)
  
  return(cell.list)
}

archr.add.seurat <- function(ao, so, meta = "seurat_clusters") {
  ao <- addCellColData(ArchRProj = ao, data = so@meta.data[, meta] %>% as.character(),
                       cells = paste0(so$sample, "#", sub("_.+$","",colnames(so))),
                       name = meta, force = T)
  return(ao)
}


so.2.se <- function(so, assay, archr.names) {
  se <- Seurat::as.SingleCellExperiment(x = so, assay = assay)
  if (archr.names == T) {
    colnames(se) <- get.cell.names.seurat(so = so, cells = colnames(se), style = "ArchR")
  }
  
  return(se)
}

so.2.se.archr <- function(so, h5 = H5.VEC, assay, slot ) {
  
  
  spMat <- Seurat::GetAssayData(object = so, slot = slot, assay = assay)
  colnames(spMat) <- get.cell.names.seurat(so = so, cells = colnames(spMat), style = "ArchR")
  
  sub.f.feature.gen <- function(h5) {
    features <- h5read(h5, "/matrix/features")
    features$`_all_tag_keys` <- NULL
    features$name[duplicated(features$name)] <- features$name[duplicated(features$name)] %>% paste0(".1")
    features <- lapply(seq_along(features), function(x) {
      df <- DataFrame(x = features[[x]])
      colnames(df) <- names(features)[x]
      df
    })
    features <- Reduce("cbind", features[!unlist(lapply(features, 
                                                        is.null))])
    features <- features %>% .[.$feature_type == "Gene Expression",]
    
    return(features)
  }
  
  ft <- lapply(h5, sub.f.feature.gen)
  ft <- Reduce(rbind, ft) %>% unique()
  rownames(ft) <- ft$name
  ft <- ft[rownames(spMat),]
  features <- ft
  
  rownames(features) <- NULL
  rownames(spMat) <- NULL
  # Browse[3]> rownames(spMat) %>% .[!.%in% as.character(ft$name)]
  # [1] "Ptp4a1.1"        "Gm16701.1"       "Pakap.1"         "Fam220a.1"       "Aldoa.1"         "Dpep2.1"         "Gm16364.1"      
  # [8] "4930594M22Rik.1" "Snhg4.1"         "Gm35438.1"       "Gm35558.1"
  se <- SummarizedExperiment(assays = SimpleList(counts = spMat), 
                             rowData = features)
  rownames(se) <- features$name
  
  if ("interval" %in% colnames(rowData(se))) {
    idxNA <- which(rowData(se)$interval == "NA")
    if (length(idxNA) > 0) {
      se <- se[-idxNA, ]
    }
    rr <- GRanges(paste0(rowData(se)$interval))
    mcols(rr) <- rowData(se)
    se <- SummarizedExperiment(assays = SimpleList(counts = assay(se)), 
                               rowRanges = rr)
  }
  idxDup <- which(rownames(se) %in% rownames(se[duplicated(rownames(se))]))
  names(idxDup) <- rownames(se)[idxDup]
  if (length(idxDup) > 0) {
    dupOrder <- idxDup[order(Matrix::rowSums(assay(se[idxDup])), 
                             decreasing = TRUE)]
    dupOrder <- dupOrder[!duplicated(names(dupOrder))]
    se <- se[-as.vector(idxDup[idxDup %ni% dupOrder])]
  }
  gc()
  se
}


da.list.gen <- function(ao, cluster.ident, clusters = NULL, group.ident, groups = NULL, mat = NULL,  n.reps = 5) {
  cell.list <- archr.get.cells.grid(ao = ao, cluster.ident = cluster.ident, clusters = clusters, group.ident = group.ident,
                                    groups = groups, melt = F)
  
  gr.list <- lapply(seq_along(cell.list), function(i) {
    cluster <- names(cell.list)[i]
    sub.list <- cell.list[[i]]
    return.list <- lapply(sub.list %>% seq_along(), function(j) {
      group.name <- names(sub.list)[j]
      cells <- sub.list[[j]]
      lapply(1:n.reps, function(r) {
        set.seed(seed = 100 + r)
        sub.1 <- sample(x = cells, size = floor(length(cells)/2), replace = F)
        # sub.2 <- cells[-sub.1]
        sub.2 <- cells[!cells %in% sub.1]
        return.list <- list(sub_1 = sub.1, sub_2 = sub.2)
        names(return.list) <- paste0(group.name, "..rep_", r,"..", names(return.list))
        return(return.list)
      }) %>% Reduce(c, .)
    }) %>% Reduce(c, .)
    
    print("miao")
    
    gr <- vectorization.core(ao = ao, cell.list = return.list, mat = mat, mat.name = "PeakMatrix", deseq2.norm = T)
  })
  
  names(gr.list) <- names(cell.list)
  return(gr.list)
}

da.plot.gen <- function(da.list, plot.dir, comp, threads = NULL, ... ) {
  if (is.null(threads)) {
    threads <- length(da.list)
  }
  
  if (threads > 12)
    threads <- 12
  
  mclapply(da.list, function(gr) {
    df <- gr %>% `names<-`(NULL) %>% as.data.frame()
    # get the number of replicates
    comp.cols <- colnames(df) %>% .[grepl("\\.\\.rep_\\d+\\.\\.", .)]
    n.reps <- comp.cols %>% 
      sub("^.+rep_", "", .) %>% sub("\\.\\..+$", "", .) %>% as.numeric() %>% max()
    p.list <- lapply(1:n.reps, function(i) {
      sample.1 <- sub(":.+$", "", comp)
      sample.2 <- sub("^.+:", "", comp)
      comp.list <- list(list(sample.1, sample.2), 
                       list(c(sample.1, "sub_1"), c(sample.1, "sub_2")),
                       list(c(sample.2, "sub_1"), c(sample.2, "sub_2")),
                       list(c(sample.1, "sub_1"), c(sample.2, "sub_1")),
                       list(c(sample.1, "sub_2"), c(sample.2, "sub_2")))
      p.list.comp <- lapply(comp.list, function(comp.sub) {
        comp.df <- lapply(comp.sub, function(comp.axis) {
          axis <- comp.cols %>% .[!grepl("_raw$", .)] 
          for (m in comp.axis) {
            axis <- axis %>% .[grepl(m, .)]
          }
          axis <- df[, axis, drop = T]
          if (is.data.frame(axis)) {
            axis <- Reduce(`+`, as.list(axis))
          }
          return(axis)
        }) %>% as.data.frame()
        
        axis.names <- sapply(comp.sub, function(x) return(paste0(x, collapse = "_")))
        colnames(comp.df) <- axis.names
        
        p <- xy.plot(df = comp.df, x = axis.names[2], y = axis.names[1], add.corr = T, ...)
        return(p)
        
      })
      return(p.list.comp)
      
    }) %>% Reduce(c, .)
    
    cluster <- comp.cols[1] %>% sub("\\.\\..+$", "",.)
    trash <- wrap.plots.fanc(plot.list = p.list, n.col = n.reps, plot.out = paste0(plot.dir, "/", cluster, ".png"), 
                             page.limit = 100)
    return()
     
  }, mc.cores = threads)
  
  return()
  
}

seurat.add.archr.meta <- function(so, ao, metas = NULL) {
  so.meta <- so@meta.data
  so.meta[,"cells"] <- paste0(so.meta$sample, "#", sub("_.+$", "", colnames(so)))
  ao.meta <- getCellColData(ArchRProj = ao) %>% as.data.frame()
  ao.meta <- ao.meta %>% utilsFanc::change.name.fanc(cols.from = "Sample", cols.to = "sample")
  metas.avail <- colnames(ao.meta)
  
  if (is.null(metas)) {
    metas <- metas.avail
    metas <- metas[!metas %in% c("Sample", "seurat_clusters")]
  }
  if (sum(!metas %in% metas.avail) > 0) {
    not.in <- metas[!metas %in% metas.avail]
    warning(paste0("The following metas are not available from the ao provided: ",
                   paste0(not.in, collapse = "; ")))
    metas <- metas[metas %in% metas.avail]
  }
  if (length(metas) < 1) {
    stop("at least one metadata term needs to be provided")
  }
  
  ao.meta <- ao.meta[, metas, drop = F]
  ao.meta$cells <- rownames(ao.meta)
  meta.new <- left_join(so.meta, ao.meta)
  rownames(meta.new) <- rownames(so.meta)
  so@meta.data <- meta.new
  
  return(so)
}

seurat.add.archr.embed <- function(so, ao, ao.embedding, embedding.name = NULL) {
  embedding <- ArchR::getEmbedding(ArchRProj = ao, embedding = ao.embedding, returnDF = T)
  # rownames(embedding) <- sub(".+#", "", rownames(embedding))
  cells <- rownames(embedding)
  
  look.up.table <- so@meta.data %>% mutate(cell = rownames(so@meta.data)) %>% select(cell, sample) %>% 
    mutate (cell.ori = sub("_.+$", "", cell)) %>% 
    mutate(new.name = paste0(sample, "#", cell.ori))
  rownames(look.up.table) <- look.up.table$new.name
  embedding <- embedding[look.up.table$new.name,]
  
  colnames(embedding) <- paste0(ao.embedding, "_", c(1,2))
  
  # so@reductions <- lapply(so@reductions, function(x) {
  #   x@key <- "miao"
  #   return(x)
  # })
  # 
  rownames(embedding) <- colnames(so)
  if (is.null(embedding.name)) {
    embedding.name <- paste0("ArchR", ao.embedding)
  }
  
  so[[embedding.name]] <- CreateDimReducObject(embeddings = as.matrix(embedding),
                                         key = paste0(embedding.name , "_"),
                                         assay = DefaultAssay(so))
  return(so)
}

.saveUWOT.fanc <- function(model, file){
  file <- file.path(normalizePath(dirname(file)), basename(file))
  wd <- getwd()
  mod_dir <- tempfile(pattern = "dir")
  dir.create(mod_dir)
  uwot_dir <- file.path(mod_dir, "uwot")
  dir.create(uwot_dir)
  model_tmpfname <- file.path(uwot_dir, "model")
  .safeSaveRDS(model, file = model_tmpfname)
  metrics <- names(model$metric)
  n_metrics <- length(metrics)
  for (i in seq_len(n_metrics)) {
    nn_tmpfname <- file.path(uwot_dir, paste0("nn", i))
    if (n_metrics == 1) {
      if (!is.null(model$nn_index$ann)) {
        model$nn_index$ann$save(nn_tmpfname)
        model$nn_index$ann$unload()
        model$nn_index$ann$load(nn_tmpfname)
      } else {
        model$nn_index$save(nn_tmpfname)
        model$nn_index$unload()
        model$nn_index$load(nn_tmpfname)
      }

    }
    else {
      if (!is.null(model$nn_index[[i]]$ann)) {
        model$nn_index[[i]]$ann$save(nn_tmpfname)
        model$nn_index[[i]]$ann$unload()
        model$nn_index[[i]]$ann$load(nn_tmpfname)
      } else {
        model$nn_index[[i]]$save(nn_tmpfname)
        model$nn_index[[i]]$unload()
        model$nn_index[[i]]$load(nn_tmpfname)
      }

    }
  }
  setwd(mod_dir)
  system2("tar", "-cvf uwot.tar uwot", stdout = NULL, stderr = NULL)
  o <- .fileRename("uwot.tar", file)
  setwd(wd)
  if (file.exists(mod_dir)) {
    unlink(mod_dir, recursive = TRUE)
  }
  return(o)
}

seurat.sync.archr <- function(so, ao) {
  df <- so@meta.data
  df$so.cells <- rownames(df)
  df$cells <- paste0(df$sample, "#", sub("_.+$", "", df$so.cells))
  ao.cells <- rownames(ao)
  df.intsct <- df %>% filter(cells %in% ao.cells)
  return(list(so = so[,df.intsct$so.cells], ao = ao[df.intsct$cells,]))
}