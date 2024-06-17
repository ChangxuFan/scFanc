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
archr.get.frag.2 <- function(obj, cell.list, bCount.only = F, threads = 1) {
  # written on 2022-08-17
  # supposed to be a simplified + paralleled version of archr.get.frag()
  # didn't finish writing...
  
}
archr.get.frag <- function(obj, cellNames, n.frag = NULL, seed = 42,  return.df = F, out.bed = NULL,
                           out.bdg = NULL, ext.left = 75, ext.right = 75,
                           return.grange = F, ...) {
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
                   value.name = "start") %>% 
    mutate(end = start + 1, strand = "+") %>% 
    dplyr::select(seqnames, start, end, RG, pseudo.umi, strand) %>% 
    arrange(RG, seqnames, start)
  if (!is.null(n.frag)) {
    if (n.frag >= nrow(df.insert)) {
      n.frag <- nrow(df.insert)
      warning("requesting more fragments than available")
    }
    set.seed(seed = seed)
    df.insert <- df.insert[sample(1:nrow(df.insert), n.frag, replace = F), ]
  }
  if (!is.null(out.bed)) {
    system(paste0("mkdir -p ", dirname(out.bed)))
    write.table(df.insert, out.bed, sep = "\t", quote = F, col.names = F, row.names = F)
  }
  if (!is.null(out.bdg)) {
    system(paste0("mkdir -p ", dirname(out.bdg)))
    # browser()
    gr.insert.ext <- df.insert %>% mutate(start = pmax(1, start - ext.left), end = end + ext.right) %>% 
      makeGRangesFromDataFrame()
    cov <- GenomicRanges::coverage(gr.insert.ext)
    cov <- cov[sort(names(cov))]
    bdg <- lapply(names(cov), function(chr) {
      rle <- cov[[chr]]
      bdg <- data.frame(chr = chr, start = start(rle) - 1, end = end(rle), score = runValue(rle)) %>% 
        filter(score > 0)
      return(bdg)
    }) # %>% Reduce(rbind, .)
    # bdg <- Reduce(rbind, bdg)
    
    for (i in length(bdg)) {
      if (i == 1)
        append <- F
      else
        append <- T
      readr::write_tsv(bdg[[i]], out.bdg, col_names = F, append = append)
    }
    
    
    
  }
  # if (return.grange==T)
  #   return(frags)
  if (return.df == T)
    return(list(df=df, df.insert = df.insert))
  return(NULL)
}




insert.2.bdg <- function(in.bed, genome, left.ext = 75, right.ext = 75, 
                         ext.bed = NULL, out.bdg = NULL, 
                         bedtools = "/bar/cfan/anaconda2/envs/jupyter/bin/bedtools") {
  stop("this function is too slow")
  genome <- paste0("/bar/cfan/genomes/", genome, "/", genome, ".chrom.sizes")
  if (is.null(ext.bed))
    ext.bed <- sub(".bed$", ".ext.bed", in.bed)
  if (is.null(out.bdg))
    out.bdg <- sub(".bed$", ".ext.bedgraph", in.bed)
  
  if (is.character(in.bed))
    in.bed <- read.table(in.bed, as.is = T, header = F, quote = "")
  
  ext.df <- in.bed
  ext.df[, 2] <- pmax(0, ext.df[, 2] - left.ext)
  ext.df[, 3] <- ext.df[, 3] + right.ext
  
  utilsFanc::write.zip.fanc(ext.df, ext.bed, bed.shift = F, zip = T)
  cmd <- paste0(bedtools, " genomecov -bg -i ",ext.bed, "-g ", genome, " > ", out.bdg)
  system(cmd)
  return(out.bdg)
}

serial.macs2 <- function(obj, cellNames.list, n.frag.vec = NULL, sub.names= NULL,master.dir, root.name  = NULL, genome,
                         get.bed = T, run.macs2 = T, get.bdg = T,
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
  # browser()
  
  jsongen <- mclapply(seq_along(cellNames.list), function(i) {
    n.frag <- n.frag.vec[i]
    x <- cellNames.list[[i]]
    sub.name <- sub.names[i]
    sub.root.name <- paste0(root.name[i], "_", sub.name)
    sub.dir <- paste0(master.dir, "/", sub.root.name, "/")
    raw.bed <- paste0(sub.dir, "/", sub.root.name, ".bed")
    bdg <- sub(".bed$", ".ext.bedgraph", raw.bed)
    if (get.bed == T)
      archr.get.frag(obj = obj, cellNames = x, out.bed = raw.bed, n.frag = n.frag.vec)
    else if (get.bdg == T) {
      archr.get.frag(obj = obj, cellNames = x, out.bdg = bdg)
    }
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



serial.macs2.2 <- function(obj, cellNames.list, sub.names= NULL,master.dir, root.name  = NULL, genome,
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
  
  
  browser()
  
  archr.get.frag.m(obj = obj, cell.list = cellNames.list, 
                   out.bdg.files = paste0(master.dir, "/", root.name, "_", sub.names, 
                                          "/", root.name, "_", sub.names, ".bedgraph"))
  
  jsongen <- mclapply(seq_along(cellNames.list), function(i) {
    x <- cellNames.list[[i]]
    sub.name <- sub.names[i]
    sub.root.name <- paste0(root.name[i], "_", sub.name)
    sub.dir <- paste0(master.dir, "/", sub.root.name, "/")
    raw.bed <- paste0(sub.dir, "/", sub.root.name, ".bed")
    if (get.bed == T) {
      archr.get.frag(obj = obj, cellNames = x, out.bed = raw.bed)
      if (run.macs2 == T) {
        rnaSeqFancLite::macs2.atac.callpeak(infile = raw.bed, root.name = sub.root.name, 
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
      }
      
    }
  }, mc.cores = thread.master) %>% Reduce(rbind, .)
  json <- jsongen %>% jsonlite::toJSON() %>% jsonlite::prettify()
  write(json, paste0(master.dir, "/macs_tracks.json"))
  
}

fake.so.gen <- function(ao, ao.embedding, n.gene = 5) {
  embedding <- ArchR::getEmbedding(ArchRProj = ao, embedding = ao.embedding, returnDF = T)
  # rownames(embedding) <- sub(".+#", "", rownames(embedding))
  cells <- rownames(embedding)
  
  # mat <- diag(length(cells))
  # rownames(mat) <- paste0("fakeGene-", 1:nrow(mat))
  # colnames(mat) <- cells
  
  n.cells <- length(cells)
  mat <- matrix(data = sample(c(0,1), x = n.gene * n.cells, replace = T),
                nrow = n.gene, ncol = n.cells)
  rownames(mat) <- paste0("fakeGene-", 1:nrow(mat))
  colnames(mat) <- cells
  
  fake.so <- CreateSeuratObject(counts = mat, project = "fake", min.cells = 0, min.features = 0)

  if (!identical(colnames(fake.so), cells)) {
    stop ("!identical(colnames(fake.so), cells)")
  }
  fake.so[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding),
                                              key = paste0(ao.embedding, "_"),
                                              assay = DefaultAssay(fake.so))
  fake.so@meta.data <- ArchR::getCellColData(ArchRProj = ao, drop = F) %>% as.data.frame()
  fake.so$Clusters <- factor(fake.so$Clusters, levels = gtools::mixedsort(unique(fake.so$Clusters)))
  fake.so[["sample"]] <- fake.so$Sample
  return(fake.so)
}  

fake.so.gen.2 <- function(cells = NULL, n.gene = 5, meta.data = NULL) {
  # note: meta.data must be a seurat object like meta.data, with rownames indicating cell names
  if (is.null(cells))
    cells <- rownames(meta.data)
  n.cells <- length(cells)
  mat <- matrix(data = sample(c(0,1), x = n.gene * n.cells), nrow = n.gene, ncol = n.cells)
  rownames(mat) <- paste0("fakeGene-", 1:nrow(mat))
  colnames(mat) <- cells
  fake.so <- CreateSeuratObject(counts = mat, project = "fake", min.cells = 0, min.features = 0)
  if (!is.null(meta.data)) {
    if (!identical(cells, rownames(meta.data)))
      stop("!identical(cells, rownames(meta.data))")
    fake.so@meta.data <- meta.data
  }
  
  return(fake.so)
  
}


seurat.plot.archr <- function(so, ao, ao.embedding, 
                              ao.name, use.levels = NULL, 
                              plot.dir = NULL, root.name = "ao",
                              label = T, pt.size = 0.5, cols = NULL,
                              width = 10, height = 10, ...) {
  
  embedding <- ArchR::getEmbedding(ArchRProj = ao, embedding = ao.embedding, returnDF = T)
  # rownames(embedding) <- sub(".+#", "", rownames(embedding))
  cells <- rownames(embedding)
  
  mat <- diag(length(cells))
  rownames(mat) <- paste0("fakeGene-", 1:nrow(mat))
  colnames(mat) <- cells
  
  fake.so <- CreateSeuratObject(counts = mat, project = "fake", min.cells = 0, min.features = 0)
  if (!identical(colnames(fake.so), cells)) {
    stop ("!identical(colnames(fake.so), cells)")
  }
  
  
  fake.so[["pseudo"]] <- CreateDimReducObject(embeddings = as.matrix(embedding),
                                         key = paste0(ao.embedding, "_"),
                                         assay = DefaultAssay(fake.so))
  
  fake.so[[ao.name]] <- ArchR::getCellColData(ArchRProj = ao, select = ao.name, drop = T)
  if (!is.null(use.levels)) {
    fake.so@meta.data[[ao.name]] <- factor(fake.so@meta.data[[ao.name]], levels = use.levels)
    cols <- utilsFanc::gg.color.keeplevel(fac.vec = fake.so@meta.data[[ao.name]])
  }
  p <- DimPlot(fake.so, reduction = "pseudo",group.by = ao.name,
          pt.size = pt.size, label = label, cols = cols,
          ...) 
  if (!is.null(plot.dir))  {
    dir.create(plot.dir, showWarnings = F,recursive = T)
    ggsave(paste0(plot.dir, "/",root.name,"_", ao.embedding, "_", ao.name, ".png"), p, 
           device = "png", width = width, height = height, dpi = 200, units = "in")
  }
  
  return(p)
} 


seurat.plot.archr.bk <- function(so, ao, ao.embedding, ao.colorBy, ao.name, plot.dir, 
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


archr.get.peak.mat <- function(ao, mat.name = "PeakMatrix", mat = NULL, return.full = F) {
  #stop("untested")
  warning("please don't run these function more than once.")
  if (is.null(mat))
    mat <- ArchR::getMatrixFromProject(ao, useMatrix = mat.name)
  # stop("need to add sort2mixedsort")
  
  order <- chr.sort2mixsort(x = rowData(mat)$idx)
  mat <- mat[order, ]
  if (return.full == F)
    mat <- assay(mat)
  return(mat)
}

# t <- ArchR::getCellColData(ao.int.prog,drop = F)
archr.get.cells.grid <- function(ao, cluster.ident, clusters=NULL, group.ident, groups = NULL,
                                 match.bg.df = NULL, match.bg.outdir, bias = c("TSSEnrichment", "log10(nFrags)"), root.name = "match_bg",
                                 filter.by = NULL, filter.limits, bQuantileFilter = c(F, F), filter.plot.dir = NULL,
                                 melt =T, pseudo.gen = F, seed = NULL, nprep = 2) {
  # the format of match.bg.df: fg, bg. If multiple are specified for bg, use "," to separate.
  # fg can only be one sample. This is determined by ArchR
  # filter.limits: length of 2. lower and upper bound. bQuantileFilter indicate if they are
  #quantiles
  print("using archr.get.cells.grid. A similar function is get.cell.list()")
  df <- ArchR::getCellColData(ao,drop = F, select = c(cluster.ident, group.ident, filter.by)) %>% 
    as.data.frame() %>%  factor2character.fanc() 
  df <- df %>% na.omit()
  
  if (length(group.ident) > 1) {
    df$mult <- Reduce(function(x,y) paste(x, y, sep = ".."),
                      df[, group.ident] %>% as.list())
    if (!is.null(groups)) {
      if (!is.list(groups) || length(groups) != length(group.ident)) {
        stop("when length(group.ident) > 1 and !is.null(groups), groups must be a list, each element corresponding to group.ident")
      }
      groups <- lapply(seq_along(groups), function(i) {
        x <- groups[[i]]
        if (is.null(x)) {
          x <- df[, group.ident[i]] %>% unique() %>% .[!is.na(.)]
        }
        return(x)
      })
      groups <- Reduce(function(x, y) outer(x, y, FUN = function(x, y) paste(x,y, sep = "..")),
                       groups) %>% as.character()
    }
    group.ident <- "mult"
  }
  if (!is.null(clusters))
    df <- df[df[, cluster.ident] %in% clusters,]
  if (!is.null(groups))
    df <- df[df[, group.ident] %in% groups,]
  
  cell.list <- df %>% split(., f = factor(.[,cluster.ident], levels = gtools::mixedsort(unique(.[,cluster.ident])))) %>% 
    lapply(function(cl.df) {
      cl <- cl.df[, cluster.ident][1]
      cell.list.cl <- cl.df %>% split(., f = factor(.[,group.ident], levels = gtools::mixedsort(unique(.[,group.ident])))) %>% 
        lapply(function(group.df) {
          group <- group.df[, group.ident][1]
          if (!is.null(filter.by)) {
            if (bQuantileFilter[1] == T) {
              filter.limits[1] <- quantile(group.df[, filter.by], filter.limits[1])
            }
            if (bQuantileFilter[2] == T) {
              filter.limits[2] <- quantile(group.df[, filter.by], filter.limits[2])
            }
            if (!is.null(filter.plot.dir)) {
              filter.plot <- paste0(filter.plot.dir, "/", filter.by, "..", cl, "..", group, "..", 
                                    paste0(filter.limits, collapse = "_"), ".png")
              trash <- rank.plot(df = group.df, vars = filter.by,
                                 add.h.line = filter.limits)  %>% list() %>% 
                wrap.plots.fanc(plot.out = filter.plot)
            }
            group.df <- group.df %>% .[.[, filter.by] > filter.limits[1] & 
                                         .[, filter.by] < filter.limits[2], ]
          }
          cells <- rownames(group.df)  
          return(cells)
        })
      # we could just use this cell.list.cl. However, if you want to match bg by archr,
      #we redefine cell.list.cl. 
      
      cl.df <- cl.df[unlist(cell.list.cl),]
      # just to take advantage of filtering
      
      if (!is.null(match.bg.df)) {
        ao.sub <- ao[ao$cellNames %in% rownames(cl.df),]
        cell.list.cl <- match.bg.df %>% split(., f= 1:nrow(.)) %>% 
          lapply(function(comp) {
            colDat <- getCellColData(ao.sub)
            fg.groups <- comp$fg
            bg.groups <- comp$bg %>% strsplit(split = ",") %>% unlist()
            n.fg <- sum(cl.df[, group.ident] %in% fg.groups)
            n.bg <- sum(cl.df[, group.ident] %in% bg.groups)
            if (n.fg > n.bg) {
              swap <- comp$fg
              comp$fg <- comp$bg
              comp$bg <- swap
              fg.groups <- comp$fg
              bg.groups <- comp$bg %>% strsplit(split = ",") %>% unlist()
            }
            matchObj <- .matchBiasCellGroups(input = colDat, groups = getCellColData(ao.sub, group.ident, drop = T), 
                                             useGroups = fg.groups, bgdGroups = bg.groups, bias = bias, 
                                             k = 100, n = 10000, bufferRatio = 0.8, logFile = NULL)
            res <- list(ao.sub$cellNames[matchObj$matchbgd[[1]]$cells %>% unique()],
                        ao.sub$cellNames[matchObj$matchbgd[[1]]$bgd %>% unique()])
            # I have taken the cell list and manually compute nFrags and TSS stats to make sure the 
            #cell names I get is actually the ones indicated in matchObj
            dir.create(match.bg.outdir, showWarnings = F, recursive = T)
            saveRDS(matchObj, paste0(match.bg.outdir, "/", root.name,"_", comp$fg, "..", comp$bg, ".Rds"))
            names(res) <- c(comp$fg, comp$bg)
            return(res)
          }) %>% Reduce(c, .)
      } 
    
      # browser()
      if (pseudo.gen == T) {
        cell.list.cl <- lapply(seq_along(cell.list.cl), function(i) {
          cells <- cell.list.cl[[i]]
          name <- names(cell.list.cl)[i]
          set.seed(seed = seed)
          split.id <- sample(1:nprep, size = length(cells), replace = T)
          cells.prep <- cells %>% split(f = split.id)
          names(cells.prep) <- paste0(name, "..prep_", 1:length(cells.prep))
          return(cells.prep)
        }) %>% Reduce(c, .)
      }
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

per.cluster.macs2.2 <- function(so=NULL, ao, meta, meta.is.ao = F, sort=F, subset = NULL, 
                                group.by=NULL, groups = NULL, equal.depth = F, spmr,get.bed = T,
                              master.dir, run = T, genome, thread.master, thread.sub) {
  
  if (is.null(so) || meta.is.ao == T) {
    df <- ArchR::getCellColData(ao, select = c(meta, group.by)) %>% as.data.frame() %>% 
      factor2character.fanc() %>% mutate(., cells = rownames(.))
  } else {
    df <- so@meta.data[, c(meta, "sample", group.by) %>% unique(), drop = F] %>% factor2character.fanc() 
    df$cells <- rownames(df) %>% sub("_.*$", "", .)
    df <- df %>% mutate(cells = paste0(sample, "#", cells))
  }
  
  if (!is.null(subset))
    df <- df %>% .[.[,meta] %in% subset,]
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
        if (equal.depth == T) {
          depths <- sapply(cells.by.group, function(x) {
            d <- getCellColData(ao, select = c("nFrags")) %>% as.data.frame() %>% .[x,"nFrags"] %>% sum()
            return(d)
          })
          depth <- min(depths)
          names(cells.by.group) <- paste0("depth", depth, "@", names(cells.by.group))
        }
      } else {
        cells.by.group <- list(x[, "cells"])
        names(cells.by.group) <- meta.value
      }
      
      
      return(cells.by.group)
    })
  cell.list <- Reduce(c, cell.list)
  
  if (equal.depth == T) {
    depth.vec <- names(cell.list) %>% stringr::str_extract(pattern = "depth\\d+@") %>% 
      stringr::str_extract("\\d+") %>% as.numeric()
  } else {
    depth.vec <- NULL
  }
  
  trash <- serial.macs2(obj = ao, cellNames.list = cell.list, genome = genome, n.frag.vec = depth.vec,
                        master.dir = master.dir, root.name = names(cell.list), 
                        get.bed = get.bed, thread.master = thread.master, thread.sub = thread.sub, spmr = spmr, run.macs2 = run)
  
  return(cell.list)
}

archr.add.seurat <- function(ao, so, meta = "seurat_clusters", as.factor = F) {
  ao.cells <- getCellNames(ao)
  so.df <- data.frame(cells = get.cell.names.seurat(so = so, style = "ArchR"),
                      meta = so@meta.data[, meta] %>% as.character())
  so.df <- so.df %>% dplyr::filter(cells %in% ao.cells)
  if (length(ao.cells %>% .[!.%in% so.df$cells]) > 0) {
    add.df <- data.frame(cells = ao.cells %>% .[!.%in% so.df$cells], 
                         meta = NA)
    so.df <- rbind(so.df, add.df)
  }

  ao <- addCellColData(ArchRProj = ao, data = so.df$meta,
                       cells = so.df$cells,
                       name = meta, force = T)
  if (as.factor == T) {
    ao@cellColData[[meta]] <- factor(ao@cellColData[[meta]],
                                     levels = gtools::mixedsort(unique(ao@cellColData[[meta]])))
  }
  return(ao)
}
archr.add.seurat.m <- function(ao, so, metas = "seurat_clusters", as.factor = F) {
  for (meta in metas) {
    ao <- archr.add.seurat(ao = ao, so = so, meta = meta, as.factor = as.factor)
  }
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
  se <- se[order(rowRanges(se)),]
  # this line is important in ensuring that rows won't be permuted when 
  # geneExpressionMatrix is added to / retrieved from arrow files.
  # it seems that when things are added to arrow files they are automatically sorted
  # by ranges. But what lives within the archr object is not automatically sorted
  # for debug info: RNA_vs_promoter-2021-08-25.R in rep2_3
  return(se)
}
t.f.compare.se <- function(se.a, se.f, a.is.mat = F, f.is.mat = F, 
                           gene.a = "Hlf", gene.f = "Hlf",
                       transformation.a = NULL,
                       cor.method = "spearman",
                       scatter = F) {
  if (is.null(rownames(se.a)))
    rownames(se.a) <- rowData(se.a)$name
  
  shared <- intersect(colnames(se.f), colnames(se.a))
  se.a <- se.a[, shared]
  se.f <- se.f[, shared]
  if (a.is.mat == F)
    a <- assay(se.a)[gene.a,]
  else
    a <- se.a[gene.a, ]
  if (f.is.mat == F)
    f <- assay(se.f)[gene.f,]
  else
    f <- se.f[gene.f, ]
  
  if (!is.null(transformation.a)) {
    a.ori <- a
    a <- transformation.a(a)
  }
    
  identical(a, f) %>% print()
  cor(a, f, method = cor.method) %>% print()
  df <- data.frame(a = a, f = f, a.ori = a.ori)
  if (scatter == T)
    xy.plot(df = df, x = "f", y = "a") %>% print()
  return(df)
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

seurat.add.archr.meta <- function(so, ao, metas = NULL, overwrite = F) {
  so.meta <- so@meta.data
  so.meta[,"cells"] <- get.cell.names.seurat(so = so, style = "ArchR")
  ao.meta <- getCellColData(ArchRProj = ao) %>% as.data.frame()
  ao.meta <- ao.meta # %>% utilsFanc::change.name.fanc(cols.from = "Sample", cols.to = "sample")
  metas.avail <- colnames(ao.meta)
  
  if (is.null(metas)) {
    metas <- metas.avail
    # metas <- metas[!metas %in% c("Sample", "seurat_clusters")]
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
  if (overwrite == T) {
    so.meta <- so.meta[, c("cells", colnames(so.meta) %>% .[! .%in% colnames(ao.meta)])]
  } else {
    ao.meta <- ao.meta[, c("cells", colnames(ao.meta) %>% .[! .%in% colnames(so.meta)])]
  }
  
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
    # mutate (cell.ori = sub("_.+$", "", cell)) %>% 
    mutate(new.name = get.cell.names.seurat(so, style = "ArchR"))
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

seurat.add.archr.matrix <- function(so, mat.se, assays.to.add = NULL, 
                                    assay.root.name = NULL,
                                    overwrite = F) {
  if (!identical(sort(colnames(so)), sort(colnames(mat.se)))) {
    stop("!is.null(identical(colnames(so), colnames(mat.se)))")
  }
  if (is.null(assays.to.add)) {
    assays.to.add <- names(assays(mat.se))
  }

  utilsFanc::check.intersect(assays.to.add, "assays.to.add", names(assays(mat.se)), "names(assays(mat.se))")

  so.assays <- names(so@assays)
  new.assays <- assays.to.add
  if (!is.null(assay.root.name)) {
    new.assays <- paste0(assay.root.name, "_", assays.to.add)
  }
  bExist <- new.assays %in% so.assays
  
  if (sum(bExist) > 0) {
    if (!overwrite) stop(paste0("Some of the assays already exists, try overwrite = T"))
    for (a in new.assays[bExist]) {
      so[[a]] <- NULL
    }
  }
  
  for (a in assays.to.add) {
    mat <- assays(mat.se)[[a]]
    a.name <- a
    if (!is.null(assay.root.name)) a.name <- paste0(assay.root.name, "_", a)
    so[[a.name]] <- CreateAssayObject(data = mat)
  }
  
  return(so)
}

seurat.add.embed <- function(so, embed.df, embedding.name) {
  # embed.df: rownames of embed.df is taken as the cell names.
  # the first 2 columns are taken as the UMAP coordinates.
  # a matrix is also okay
  warning("this is a fast function that might not be the most reliable")
  embed.df <- embed.df[rownames(embed.df) %in% colnames(so), 1:2]
  so[[embedding.name]] <- CreateDimReducObject(embeddings = as.matrix(embed.df),
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
  df$cells <- get.cell.names.seurat(so = so, cells = df$so.cells, style = "ArchR")
  ao.cells <- rownames(ao)
  df.intsct <- df %>% filter(cells %in% ao.cells)
  return(list(so = so[,df.intsct$so.cells], ao = ao[df.intsct$cells,]))
}

gr.sync <- function(gr1, gr2, minoverlap, out.dir) {
  o <- findOverlaps(gr1, gr2, minoverlap = minoverlap)
  stats <- data.frame(total = c(length(gr1), length(gr2)),
                      hit = c(sum(!duplicated(queryHits(o))), sum(!duplicated(subjectHits(o)))),
                      dup = c(sum(duplicated(queryHits(o))), sum(duplicated(subjectHits(o)))))
  rownames(stats) <- c("gr1", "gr2")
  gr1.filter <- gr1[queryHits(o) %>% unique(),]
  gr2.filter <- gr2[subjectHits(o) %>% unique(),]
  gr.merge <- c(gr1.filter, gr2.filter) %>% GenomicRanges::reduce()
  res.list <- list(gr1.filter = gr1.filter, gr2.filter = gr2.filter,
                   stats = stats, gr.merge = gr.merge, gr1 = gr1, gr2 = gr2)
  trash <- mclapply(c("gr1.filter", "gr2.filter", "gr.merge", "gr1", "gr2"), function(x) {
    gr <- res.list[[x]]
    utilsFanc::write.zip.fanc(df = gr, out.file = paste0(out.dir, "/",x,".bed"), bed.shift = T)
    return()
  }, mc.cores = 10)
  
  res.list$gr1 <- NULL
  res.list$gr2 <- NULL
  trash <- utilsFanc::write.zip.fanc(df = res.list$stats, 
                                     out.file = paste0(out.dir, "/stats.tsv"), zip = F,
                                     row.names = T, col.names = T)
  return(res.list)
}

archr.split.pileup <- function(aoi, split.by, splits = NULL) {
  if (is.null(splits))
    splits <- aoi@cellColData[[split.by]] %>% as.character() %>% unique()
}


getGroupBW.fanc <- function (ArchRProj = NULL, cell.list = NULL, groupBy, groups = NULL, split.by, splits = NULL,
                             normMethod = "ReadsInTSS", out.dir,  
                             tileSize = 100, maxCells = 1000, ceiling = 4, verbose = TRUE, 
                             threads = getArchRThreads(), logFile = createLogFile("getGroupBW")) {
  cellGroups <- cell.list
  dir.create(out.dir, showWarnings = F, recursive = T)
  # .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  # .validInput(input = groupBy, name = "useMatrix", valid = c("character"))
  # .validInput(input = normMethod, name = "groupBy", valid = c("character"))
  # .validInput(input = tileSize, name = "divideN", valid = c("integer"))
  # .validInput(input = maxCells, name = "scaleTo", valid = c("integer", 
  #                                                           "null"))
  # .validInput(input = ceiling, name = "ceiling", valid = c("integer", 
  #                                                          "null"))
  # .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  # .validInput(input = threads, name = "threads", valid = c("integer"))
  # .validInput(input = logFile, name = "logFile", valid = c("character"))
  tstart <- Sys.time()
  normMethods <- c("None", "ReadsInTSS", "nCells", "ReadsInPromoter", 
                   "nFrags")
  if (tolower(normMethod) %ni% tolower(normMethods)) {
    stop(paste0("normMethod (", normMethod, ") not in supported normMethods : ", 
                paste0(normMethods, collapse = ", ")))
  }
  # .startLogging(logFile = logFile)
  # .logThis(normMethod, "normMethod", logFile = logFile)
  ArrowFiles <- getArrowFiles(ArchRProj)
  if (tolower(normMethod) %in% tolower(c("ReadsInTSS", "ReadsInPromoter", 
                                         "nFrags"))) {
    normBy <- getCellColData(ArchRProj = ArchRProj, select = normMethod)
  }
  else {
    normBy <- NULL
  }
  if (is.null(cellGroups)) {
    cellGroups <- get.cell.list(obj = ArchRProj, is.ao = T, split.by = split.by,
                                splits = splits, group.by = groupBy, groups = groups)
    names(cellGroups) <- paste0(names(cellGroups), "..", sapply(cellGroups, length))
  } 
  
  if (!is.null(maxCells)) {
    gnames <- names(cellGroups)
    cellGroups <- lapply(seq_along(cellGroups), function(x) {
      if (length(cellGroups[[x]]) > maxCells) {
        sample(cellGroups[[x]], maxCells)
      }
      else {
        cellGroups[[x]]
      }
    })
    names(cellGroups) <- gnames
  }
  cellsInArrow <- split(rownames(getCellColData(ArchRProj)), 
                        getCellColData(ArchRProj)$Sample)
  availableChr <- .availableSeqnames(head(getArrowFiles(ArchRProj)))
  chromLengths <- getChromLengths(ArchRProj)
  chromSizes <- getChromSizes(ArchRProj)
  tiles <- unlist(slidingWindows(chromSizes, width = tileSize, 
                                 step = tileSize))
  if (threads > 1) {
    h5disableFileLocking()
  }
  covFiles <- c()
  for (i in seq_along(cellGroups)) {
    # browser()
    o <- .createGroupBW(i = i, cellGroups = cellGroups, 
                        ArrowFiles = getArrowFiles(ArchRProj), cellsInArrow = cellsInArrow, 
                        availableChr = availableChr, chromLengths = chromLengths, 
                        normMethod = normMethod, normBy = normBy, ceiling = ceiling, 
                        tiles = tiles, tileSize = tileSize, bwDir = out.dir, 
                        threads = threads, verbose = verbose, tstart = tstart)
    covFiles <- c(covFiles, o)
  }
  if (threads > 1) {
    h5enableFileLocking()
  } 
  
  return()
}


archr.diffbind.gen <- function(ao = NULL, cluster.ident, clusters, 
                               group.ident, groups, 
                               macs2.dir, file.out, 
                               threads.clusters = 1, threads.groups = 4, threads.bam.sort = 1,
                               genome, force.bam.gen = F, 
                               bedToBam = "/bar/cfan/anaconda2/envs/jupyter/bin/bedToBam", 
                               samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
  # groups: basically samples. the information of replicates are parsed from the sample names
  if (!is.null(ao))
    stop("this part of the function has not been developed")
  sample.sheet <- utilsFanc::safelapply(clusters, function(cluster) {
    utilsFanc::safelapply(groups, function(group) {
      dir.name <- paste0(macs2.dir, "/", cluster, "_", group, "_*") %>% Sys.glob()
      root.name <- basename(dir.name)
      factor <- sub("_rep.+", "", group)
      if (length(dir.name) != 1)
        stop("length(dir.name) != 1")
      rep <- stringr::str_extract(group, "rep\\d+") %>% sub("rep", "", .) %>% as.numeric()
      peak.file <- paste0(dir.name, "/", root.name, "_peaks.narrowPeak")
      if (!file.exists(peak.file))
        stop(paste0(peak.file, " doesn't exist"))
      bam.file <- paste0(dir.name, "/", root.name, ".bam")
      if (!file.exists(bam.file) || force.bam.gen == T) {
        bed.file <- paste0(dir.name, "/", root.name, ".bed")
        cmd <- paste0(bedToBam, " -i ", bed.file, " -g ", "~/genomes/", genome, "/", genome, ".chrom.sizes",
                      " -mapq 42 | ", samtools, " sort - -m 4G -o ",
                      bam.file, " -@ ", threads.bam.sort)
        print(cmd); system(cmd)
        
        cmd <- paste0(samtools, " index ",bam.file )
        print(cmd); system(cmd)
        
        if (!file.exists(bam.file))
          stop(paste0(bam.file, " was not successfully created"))
      }
      bam.index <- paste0(bam.file, ".bai")
      if (!file.exists(bam.index)) {
        cmd <- paste0(samtools, " index ", bam.file)
        print(cmd); system(cmd)
      } 
      
      df <- data.frame(SampleID = group, Factor = factor, Replicate = rep, bamReads = bam.file,
                       Peaks = peak.file, PeakCaller = "narrow")
      return(df)
    }, threads = threads.groups) %>% Reduce(rbind, .) %>% return()
  }, threads = threads.clusters) %>% Reduce(rbind, .)
  
  dir.create(dirname(file.out), showWarnings = F, recursive = T)
  write.table(sample.sheet, file.out, quote = F, sep = "\t", row.names = F, col.names = T)
  return(sample.sheet)
}

archr.enrich.2.df <- function(archr.se, assays = NULL) {
  if (is.null(assays)) {
    assays <- names(assays(archr.se))
  }
  assays <- assays %>% .[.%in%names(assays(archr.se))]
  df <- lapply(assays, function(assay) {
    return(assays(archr.se)[[assay]][ ,1])
  }) %>% `names<-` (assays) %>% as.data.frame()
  # browser()
  # df <- utilsFanc::add.column.fanc(df1 = df, df2 = data.frame(feature = rownames(archr.se)),
  #                                  pos = 1)
  return(df)
}

save.archr.fanc <- function(ao, new.work.dir, copy.arrow = T, pseudo.copy = F, over.write = F) {
  new.work.dir <- normalizePath(new.work.dir, mustWork = F)
  if (dir.exists(new.work.dir) && over.write == F)
    stop(paste0(new.work.dir, " already exists"))
  dir.create(new.work.dir, showWarnings = F, recursive = T)
  ao@projectMetadata$outputDirectory <- new.work.dir
  samples <- ao$Sample %>% unique()
  ao@sampleColData <- ao@sampleColData[rownames(ao@sampleColData) %in% samples,, drop = F]
  if (copy.arrow == T) {
    new.arrow.dir <- paste0(new.work.dir, "/ArrowFiles")
    dir.create(new.arrow.dir)
    old.arrow <- ao@sampleColData$ArrowFiles
    new.arrow <- paste0(new.arrow.dir, "/", basename(old.arrow))
    if (pseudo.copy == F) {
      utilsFanc::safelapply(1:length(old.arrow), function(i) {
        system(paste0("rsync --partial ", old.arrow, " ", new.arrow))
        return()
      }, threads = length(old.arrow))
    }
    ao@sampleColData$ArrowFiles <- new.arrow
  }
  saveRDS(ao, paste0(new.work.dir, "/ao.Rds"))
  return(ao)
}

get.named.coldata <- function(ao, col, na.rm = F) {
  if (length(col) != 1)
    stop("length(col) != 1")
  res <- getCellColData(ArchRProj = ao, select = col, drop = T)
  names(res) <- ao$cellNames
  if (na.rm == T)
    res <- res %>% .[!is.na(.)]
  return(res)
}

archr.find.closest.cells <- function(cell.list, ao,
                                     n.bg.each,
                                     per.sample = T,
                                     plot.dir = NULL, root.name) {
  if (!is.list(cell.list)) {
    vec.flag <- T
    cell.list <- list(cell.list)
  } else {
    vec.flag <- F
  }
  if (is.null(names(cell.list)))
    names(cell.list) <- paste0("cells_", 1:length(cell.list))
  res.list <- lapply(names(cell.list), function(name) {
    cells <- cell.list[[name]]
    utilsFanc::check.intersect(x = cells, x.name = paste0("cells in ", name),
                               y = ao$cellNames, y.name = "cells in ao")
    df <- ao@embeddings$UMAP$df
    if (per.sample == T) {
      cells.by.sample <- cells %>% split(., f = sub("#.+$", "", .))
      df.by.sample <- df %>% split(., f = sub("#.+$", "", rownames(.)))
      res <- lapply(names(cells.by.sample), function(sample) {
        cells <- cells.by.sample[[sample]]
        df <- df.by.sample[[sample]]
        res <- bg.gen.2(mat = df, fg.vec = cells, n.bg.each = n.bg.each, no.replace = T, 
                        method = "euclidean", scale = "none")
        return(res)
      }) %>% unlist()
    } else {
      stop("only support per.sample matching so far")
    }
    ao$tmp <- NA
    ao$tmp[ao$cellNames %in% cells] <- "fg"
    ao$tmp[ao$cellNames %in% res] <- "bg"
    if (!is.null(plot.dir)) {
      plot.panel.list(panel.list = "tmp", ao, order = T, assay = "RNA", binarize.panel = T,
                      split.by = "sample", invisible = T,
                      plot.out = paste0(plot.dir, "/", root.name, "_", name, ".png"))
    }
    return(res)
  })
  names(res.list) <- names(cell.list)
  if (length(res.list) == 1) {
    res.list <- unlist(res.list)
  } 
  if (!is.null(plot.dir)) {
    saveRDS(res.list, paste0(plot.dir, "/", root.name, ".Rds"))
  }
  return(res.list)
}


footprint.motif.dar <- function(motifPos, da, slot = "summary") {
  res <- lapply(da, function(a2b) {
    res <- lapply(c("up", "down"), function(type) {
      dars <- a2b[[slot]][[paste0(type, ".genes")]] %>% 
        utilsFanc::loci.2.gr()
      pos <- lapply(motifPos, function(motifs.all) {
        strand(motifs.all) <- "*"
        motifs.sub <- subsetByOverlaps(motifs.all, dars)
        return(motifs.sub)
      })
      return(pos)
    })
    names(res) <- c("up", "down")
    return(res)
  })
  return(res)
}

archr.peak.motif.df.gen <- function(motif.gr.list, peaks = NULL, aoi,
                                    threads = 12, out.rds) {
  if (missing(out.rds)) {
    stop("missing(out.rds)")
  }
  if (!is.null(peaks)) {
    if (is.character(peaks)) {
      peaks <- utilsFanc::loci.2.gr(peaks)
    }
  } else {
    peaks <- getPeakSet(aoi)
  }
  mcols(peaks) <- NULL
  names(peaks) <- NULL
  strand(peaks) <- "*"
  df <- utilsFanc::safelapply(names(motif.gr.list), function(motif.name) {
    motif.gr <- motif.gr.list[[motif.name]]
    if (is.null(motif.gr$score)) {
      stop("the 'score' column must be present in motif.gr")
    }
    
    gr <- peaks %>% plyranges::join_overlap_inner(motif.gr)
    df <- data.frame(gene = utilsFanc::gr.get.loci(gr = gr),
                     motif = motif.name, score = gr$score)
    return(df)
  }, threads = threads) %>% do.call(rbind, .)
  df <- df %>% dplyr::group_by(gene, motif) %>% 
    dplyr::summarise(score = max(score), n.motif = n()) %>% 
    dplyr::ungroup() %>% 
    as.data.frame()
  rownames(df) <- NULL
  dir.create(dirname(out.rds), showWarnings = F, recursive = T)
  saveRDS(df, out.rds)
  invisible(df)
}

archr.gene.motif.df.gen <- function(promoters.gr,
                                    peak.motif.df, motifs.use = NULL, motif.min.score  = NULL,
                                    method = "distance", distance.1side = 50000) {
  # promoter.gr: requires the "gene" column for gene names
  utilsFanc::check.intersect(c("gene"), "required fields",
                             colnames(mcols(promoters.gr)), "colnames(mcols(promoters.gr))")
  utilsFanc::check.intersect(c("gene", "motif", "score"), "required fields",
                             colnames(peak.motif.df), "colnames(peak.motif.df)")
  
  if (!is.null(motifs.use)) {
    peak.motif.df <- peak.motif.df %>% dplyr::filter(motif %in% motifs.use)
  }
  if (!is.null(motif.min.score)) {
    peak.motif.df <- peak.motif.df %>% dplyr::filter(score >= motif.min.score)
  }
  motifs <- peak.motif.df$motif %>% unique()
  if (nrow(peak.motif.df) < 1) {
    stop("nrow(peak.motif.df) < 1")
  }
  colnames(peak.motif.df)[colnames(peak.motif.df) == "gene"] <- "peak"
  
  peak.motif.gr <- utilsFanc::loci.2.df(
    df = peak.motif.df, loci.col.name = "peak", remove.loci.col = F, return.gr = T)
  
  if (method == "distance") {
    promoters.gr <- resize(promoters.gr, width = 2 * distance.1side, fix = "center")
    df <- plyranges::join_overlap_left(promoters.gr, peak.motif.gr) %>% 
      as.data.frame() %>% dplyr::select(gene, motif, peak) %>% unique() %>% na.omit()
    df <- df %>% dplyr::group_by(gene, motif) %>% 
      dplyr::summarise(n = n(), peak = paste0(peak, collapse = ",")) %>% 
      as.data.frame()
    attr(df, "distance.1side") <- 50000
  }
  attr(df, "method") <- method
  attr(df, "motifs.use") <- motifs
  attr(df, "motif.min.score") <- ifelse(is.null(motif.min.score), "not specified", motif.min.score)
  return(df)
}

footprint.decomposition <- function(seFoot, flank = 250, flankNorm = 50, # norm.factors = NULL,
                                    normMethod = "none",
                                    footprint.smoothWindow = 20, pileup.smoothWindow = 100,
                                    pal = NULL, baseSize = 5, ymax.fix = NULL,
                                    out.dir, root = NULL, title.root = "", remove.legend = F,
                                    debug = F, debug.xlim = 60) {
  # Motif composition doesn't really make sense
  # Eventually this function was only useful for the debug mode where I assessed what the 
  # background bias looked like.
  
  
  # seFoot <- getFootprints(ArchRProj = aoi, positions = motifPositions[greped[3]],
  #                         groupBy = "cluster_type", useGroups = groups[[3]])
  if (is.null(root)) root <- basename(out.dir)
  motifs <- names(assays(seFoot))
  df <- lapply(motifs, function(name) {
    smoothWindow <- footprint.smoothWindow
    # copying from .ggFootprint from ArchR.
    errorList <- list()
    rowDF <- SummarizedExperiment::rowData(seFoot)
    footMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[, 2] == "footprint"), ], name)
    biasMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[, 2] == "bias"), ], name)
    footDF <- rowDF[BiocGenerics::which(rowDF[, 2] == "footprint"), ]
    # FANC added:
    printDF <- footDF
    pileupDF <- footDF
    #
    biasDF <- rowDF[BiocGenerics::which(rowDF[, 2] == "bias"), ]
    
    errorList$footMat <- footMat
    errorList$biasMat <- biasMat
    errorList$footDF <- footDF
    errorList$biasDF <- biasDF
    
    if (!is.null(smoothWindow)) {
      footMat <- apply(footMat, 2, function(x) .centerRollMean(x, smoothWindow))
      biasMat <- apply(biasMat, 2, function(x) .centerRollMean(x, smoothWindow))
    }
    
    # if (!is.null(norm.factors)) {
    #   
    # }
    
    if (!is.null(flankNorm)) {
      idx <- which(abs(footDF$x) >= flank - flankNorm)
      footMat <- t(t(footMat)/colMeans(footMat[idx, , drop = FALSE]))
      biasMat <- t(t(biasMat)/colMeans(biasMat[idx, , drop = FALSE]))
      errorList$footMatNorm <- footMat
      errorList$biasMatNorm <- footMat
    }
    
    if (tolower(normMethod) == "none") {
      title <- ""
    }
    else if (tolower(normMethod) == "subtract") {
      title <- "Tn5 Bias Subtracted\n"
      footMat <- footMat - biasMat
    }
    else if (tolower(normMethod) == "divide") {
      title <- "Tn5 Bias Divided\n"
      footMat <- footMat/biasMat
    }
    else {
      stop("normMethod not recognized!")
    }
    
    # FANC: add pileupMat. This gives the general trend line
    pileupMat <- apply(footMat, 2, function(x) .centerRollMean(x, pileup.smoothWindow))
    printMat <- footMat - pileupMat
    # 
    
    footMatMean <- .groupMeans(footMat, SummarizedExperiment::colData(seFoot)$Group)
    footMatSd <- .groupSds(footMat, SummarizedExperiment::colData(seFoot)$Group)
    biasMatMean <- .groupMeans(biasMat, SummarizedExperiment::colData(seFoot)$Group)
    biasMatSd <- .groupSds(biasMat, SummarizedExperiment::colData(seFoot)$Group)
    smoothFoot <- rowMaxs(apply(footMatMean, 2, function(x) .centerRollMean(x, 11))) # FANC: only used for quantile computation later to set y axis limits
    
    # FANC: added:
    printMatMean <- .groupMeans(printMat, SummarizedExperiment::colData(seFoot)$Group)
    printMatSd <- .groupSds(printMat, SummarizedExperiment::colData(seFoot)$Group)
    pileupMatMean <- .groupMeans(pileupMat, SummarizedExperiment::colData(seFoot)$Group)
    pileupMatSd <- .groupSds(pileupMat, SummarizedExperiment::colData(seFoot)$Group)
    #
    
    errorList$footMatMean <- footMatMean
    errorList$footMatSd <- footMatSd
    errorList$biasMatMean <- biasMatMean
    errorList$biasMatSd <- biasMatSd
    errorList$smoothFoot <- smoothFoot
    
    plotIdx <- seq_len(nrow(footMatMean))
    
    plotFootDF <- lapply(seq_len(ncol(footMatMean)), function(x) {
      data.frame(x = footDF$x, 
                 mean = footMatMean[, x], 
                 sd = footMatSd[, x], 
                 group = colnames(footMatMean)[x]
                 )[plotIdx, , drop = FALSE]
    }) %>% Reduce("rbind", .)
    plotFootDF$group <- factor(paste0(plotFootDF$group), levels = unique(gtools::mixedsort(paste0(plotFootDF$group))))
    
    
# >     plotFootDF %>% head()
#      x        mean          sd group
# 1 -250 -0.03867912 0.002682626 3..KO
# 2 -249 -0.03867912 0.002682626 3..KO
# 3 -248 -0.03867912 0.002682626 3..KO
# 4 -247 -0.03867912 0.002682626 3..KO
# 5 -246 -0.03867912 0.002682626 3..KO
# 6 -245 -0.03867912 0.002682626 3..KO
    
    plotBiasDF <- lapply(seq_len(ncol(biasMatMean)), function(x) {
      data.frame(x = biasDF$x, 
                 mean = biasMatMean[, x], 
                 sd = biasMatSd[, x], 
                 group = colnames(biasMatMean)[x]
                 )[plotIdx, , drop = FALSE]
    }) %>% Reduce("rbind", .)
    plotBiasDF$group <- factor(paste0(plotBiasDF$group), levels = unique(gtools::mixedsort(paste0(plotBiasDF$group))))
    
    plotPrintDF <- lapply(seq_len(ncol(printMatMean)), function(x) {
      data.frame(x = printDF$x, 
                 mean = printMatMean[, x], 
                 sd = printMatSd[, x], 
                 group = colnames(printMatMean)[x]
      )[plotIdx, , drop = FALSE]
    }) %>% Reduce("rbind", .)
    plotPrintDF$group <- factor(paste0(plotPrintDF$group), levels = unique(gtools::mixedsort(paste0(plotPrintDF$group))))
    
    plotPileupDF <- lapply(seq_len(ncol(pileupMatMean)), function(x) {
      data.frame(x = pileupDF$x, 
                 mean = pileupMatMean[, x], 
                 sd = pileupMatSd[, x], 
                 group = colnames(pileupMatMean)[x]
      )[plotIdx, , drop = FALSE]
    }) %>% Reduce("rbind", .)
    plotPileupDF$group <- factor(paste0(plotPileupDF$group), levels = unique(gtools::mixedsort(paste0(plotPileupDF$group))))
    
    
    errorList$plotFootDF <- plotFootDF
    errorList$plotBiasDF <- plotBiasDF
    errorList$plotPrintDF <- plotPrintDF
    
    if (is.null(pal)) {
      pal <- paletteDiscrete(values = gtools::mixedsort(SummarizedExperiment::colData(seFoot)$Group))
    }
#     plotMax <- plotFootDF[order(plotFootDF$mean, decreasing = TRUE), 
#                           ]
#     plotMax <- plotMax[abs(plotMax$x) > 20 & abs(plotMax$x) < 
#                          50, ]
#     plotMax <- plotMax[!duplicated(plotMax$group), ]
#     plotMax <- plotMax[seq_len(ceiling(nrow(plotMax)/4)), 
#                        ]
#     plotMax$x <- 25 # remember that each footprint plot generated in ArchR, only the highest curve is labeled?
#     # This is how!
# # > plotMax %>% head()
# #      x     mean         sd group
# # 725 25 2.036531 0.05803274 3..WT
    
    # FANC: refuse to write each plot individually as in ArchR. Put into a list!
    if (debug) {
      df.list <- list(foot = plotFootDF, bias = plotBiasDF)
    } else {
      df.list <- list(foot = plotFootDF, pileup = plotPileupDF,
                      print = plotPrintDF, printZoom = plotPrintDF %>% filter(abs(x) <= 50))
    }
    
    
    pl <- lapply(names(df.list), function(type) {
      ymin <- quantile(df.list[[type]]$mean, 0)
      ymax <- 1.15 * quantile(smoothFoot, 0.999)
      xlims <- c(min(plotFootDF$x), max(plotFootDF$x))
      
      if (!is.null(ymax.fix)) ymax <- ymax.fix
      
      if (type %in% c("print", "printZoom")) {
        ymax = 1.2 * max(df.list[[type]]$mean)
      }
      
      if (type == "printZoom") {
        xlims <- c(min(df.list$printZoom$x), max(df.list$printZoom$x))
      }
      
      p <- ggplot(df.list[[type]], aes(x = x, y = mean, color = group)) + 
        geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), 
                    alpha = 0.4, show.legend = ifelse(type == "foot",T,F)) + 
        geom_line(size = 0.15) + 
        scale_color_manual(values = pal) + 
        scale_fill_manual(values = pal) + 
        xlab("Distance to motif center (bp)") + 
        coord_cartesian(expand = FALSE, 
                        ylim = c(ymin, ymax), 
                        xlim = xlims) + 
        # scale_x_continuous(n.breaks = 3) +
        guides(fill = FALSE) + guides(color = FALSE)
      
      #>>>>>>>> testing zone
      if (debug) {
        p <- ggplot(df.list[[type]], aes(x = x, y = mean, color = group)) + 
          geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), 
                      alpha = 0.4, show.legend = ifelse(type == "foot",T,F)) + 
          geom_line(size = 0.15) + 
          # geom_vline(xintercept = seq(from = -20, to = 20, by = 2), size = 0.1, color = "grey50", linetype = "dashed") +
          geom_vline(xintercept = c(-30, -20, 20, 30), size = 0.5, color = "red", linetype = "dashed") +
          scale_color_manual(values = pal) + 
          scale_fill_manual(values = pal) + 
          xlab("Distance to motif center (bp)") + 
          coord_cartesian(expand = FALSE,
                          # ylim = c(ymin, ymax),
                          xlim = c(-debug.xlim, debug.xlim)) +
          # scale_x_continuous(n.breaks = 3) +
          guides(fill = FALSE) + guides(color = FALSE)
          
        
        # file <- paste0(out.dir, "/", root, "_", name, "_", type , "_zoomInView.pdf")
        # dir.create(out.dir, showWarnings = F, recursive = T)
        # ggsave(file, p, device = cairo_pdf, width = 4, height = 2)
        
        if (type == "foot") {
          p <- p + ggtitle(name)
        }
        return(p)
      }
      #<<<<<<<< END testing zone
      
      if (type == "foot")
        p <- p + ggtitle(paste0(title.root, sub("\\(.+$", "", name)))
      
      p <- utilsFanc::theme.fc.1(p, rm.x.ticks = F, rotate.x.45 = F, text.size = baseSize, italic.x = F)
      if (type %in% c("foot", "pileup"))
        p <- p + theme(aspect.ratio = 1)
      if (type %in% c("print", "printZoom")) {
        p <- p + scale_y_continuous(n.breaks = 2)
      }
      
      if (remove.legend)
        p <- p + theme(legend.position = "none")
      # dir.create(out.dir, showWarnings = F, recursive = T)
      # file <- paste0(out.dir, "/", root, "_", name, "_", type, ".pdf")
      return(p)
      # ggrepel::geom_label_repel(data = plotMax, aes(label = group), size = 3, xlim = c(75, NA))
    })
    names(pl) <- names(df.list)
    if (debug) {
      p <- pl$foot + pl$bias + patchwork::plot_layout(ncol = 1)
      file <- paste0(out.dir, "/", root, "_", name, "_", "_alined.pdf")
      dir.create(out.dir, showWarnings = F, recursive = T)
      ggsave(file, p, device = cairo_pdf, width = 4, height = 4)
      return()
    } 
    
    
    design <- "
#B
AB
AB
AB
AB
AB
AC
#D
"
    p <- pl$foot + pl$pileup + pl$print + pl$printZoom + patchwork::plot_layout(design = design)
    file <- paste0(out.dir, "/", root, "_", name, ".pdf")
    dir.create(out.dir, showWarnings = F, recursive = T)
    ggsave(file, p, device = cairo_pdf, width = 2, height = 2)
    rds <- paste0(out.dir, "/", root, "_", name, "_pl.Rds")
    saveRDS(pl, rds)
    # bar plots
    # center mean
    center <- ceiling(nrow(footMat)/2)
    center.range <- (center - 25):(center+25)
    print("For center mean, using range:")
    print(paste0(min(center.range), " -> ", max(center.range)))
    meanMat <- footMat[center.range, ] %>% colMeans() %>% t()
    
    df <- colData(seFoot) %>% as.data.frame()
    df[colnames(meanMat),"center.mean"] <- meanMat[1, ]
    df <- df[, c("Group", "Name", "center.mean")]
    df$motif <- name
    max.value <- df %>% group_by(Group) %>% summarise(mean = mean(center.mean)) %>% pull(mean) %>% max()
    df$center.mean <- df$center.mean/max.value
    
    return(df)
    
  }) %>% do.call(rbind, .)
  if (debug) {
    return()
  }
  df$motif <- sub("\\(.+$", "", df$motif)
  file <- paste0(out.dir, "/", root, "_stats.tsv")
  write.table(df, file, sep = "\t", row.names = F, col.names = T, quote = F)
  groups <- df$Group %>% unique()
  p <- utilsFanc::barplot.pub.3(
    df = df, x = "motif", y = "center.mean", color.by = "Group",
    spread.width = 0.25, spread.bin.size = 1, bar.line.size = 0.5,
    bar.width = 0.70, dodge.width = 0.8, pt.size = 1, palette.fc = "R4.fc1",
    add.pval = T, pval.adjust = "fdr", pval.group.1 = groups[1], 
    pval.group.2 = groups[2]) %>% 
    utilsFanc::theme.fc.1(rm.x.ticks = F, italic.x = F, rotate.x.45 = F) +
    coord_cartesian(ylim = c(0.6, 1.05)) +
    ggsave(paste0(out.dir, "/", root, "_center_mean_bar", ".pdf"),
           device = cairo_pdf, width = 1, height = 1, dpi = 300)
  return()
}
t.f.footprint.arrange <- function(in.dir, locations, motifs.by.cluster, 
                                  plot.use = "foot",
                                  out.dir, root = NULL) {
  if (is.null(root)) root <- basename(out.dir)
  
  lapply(names(motifs.by.cluster), function(cluster) {
    
    stats.files <- paste0(in.dir, "/", locations, "_", cluster, "_none_stats.tsv")
    utilsFanc::check.file.exists(stats.files, "File")
    
    stats <- lapply(1:length(stats.files), function(i) {
      df <- read.table(stats.files[i], sep = "\t", header = T)
      df$location <- locations[i]
      return(df)
    }) %>% do.call(rbind, .)
    
    motifs <- motifs.by.cluster[[cluster]]
    lapply(motifs, function(motif) {
      motif.short <- motif %>% sub("\\(.+$", "", .)
      
      utilsFanc::check.intersect(motif.short, "motif.short", stats$motif, "stats$motif")
      
      stats <- stats %>% filter(motif == motif.short)
      Groups <- stats$Group %>% unique()
      locations.short <- locations %>% sub("Promoter", "Pro.", .) %>% 
        sub("Intergenic", "Interg.", .) %>% 
        sub("Intron", "Intr.", .)
      stats$location <- stats$location %>% sub("Promoter", "Pro.", .) %>% 
        sub("Intergenic", "Interg.", .) %>% 
        sub("Intron", "Intr.", .)
      stats$location <- factor(stats$location, levels = locations.short)
      
      p <- utilsFanc::barplot.pub.3(
        df = stats, x = "location", y = "center.mean", color.by = "Group",
        spread.width = 0.25, spread.bin.size = 1, bar.line.size = 0.5,
        bar.width = 0.70, dodge.width = 0.8, pt.size = 1, palette.fc = "ArchR",
        add.pval = T, pval.adjust = "fdr", pval.group.1 = Groups[1], alternative = "greater",
        pval.group.2 = Groups[2]) %>% 
        utilsFanc::theme.fc.1(rm.x.ticks = F, italic.x = F, rotate.x.45 = F, text.size = 5) +
        theme(aspect.ratio = 1, legend.position = "none") +
        coord_cartesian(ylim = c(0.6, 1.2))
      # ggsave(paste0(out.dir, "/", root, "_center_mean_bar", ".pdf"),
      #        device = cairo_pdf, width = 1, height = 1, dpi = 300)
      
      
      pl.rds <- paste0(in.dir, "/", locations, "_", cluster, "_none_", motif, "_pl.Rds")
      utilsFanc::check.file.exists(x = pl.rds, x.name = "File")
      
      pl <- lapply(pl.rds, function(x) return(readRDS(x)[[plot.use]]))
      names(pl) <- locations
      pl <- lapply(locations, function(location) {
        p <- pl[[location]]
        p <- p + labs(title = NULL, subtitle = paste0("    ", location)) +
          theme(legend.position = "none", 
                plot.subtitle = element_text(margin = margin(b = -7)))
        return(p)
      })
      
      pl$stats <- p
      
      p <- cowplot::plot_grid(plotlist = pl, align = "hv")
      file <- paste0(out.dir, "/", root, "_", cluster, "_", motif, "_byLocation.pdf")
      dir.create(out.dir, showWarnings = F, recursive = T)
      ggsave(file, p, device = cairo_pdf, width = 1.8, height = 1.8, dpi = 300)
      return()
    })
    
  })
}

footprint.df.gen <- function(seFoots, flank = 1000, flankNorm = 100, 
                             out.dir, root = NULL) {
  if (is.null(root)) root <- basename(out.dir)
  
  if (is.null(names(seFoots))) {
    stop("seFoots must be named")
  }
  clusters <- names(seFoots)
  
  dfs <- lapply(clusters, function(cluster) {
    seFoot <- seFoots[[cluster]]
    motifs <- names(assays(seFoot))
    dfs <- lapply(motifs, function(name) {
      # smoothWindow <- footprint.smoothWindow
      # copying from .ggFootprint from ArchR.
      errorList <- list()
      rowDF <- SummarizedExperiment::rowData(seFoot)
      footMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[, 2] == "footprint"), ], name)
      biasMat <- .getAssay(seFoot[BiocGenerics::which(rowDF[, 2] == "bias"), ], name)
      footDF <- rowDF[BiocGenerics::which(rowDF[, 2] == "footprint"), ]
      # # FANC added:
      # printDF <- footDF
      # pileupDF <- footDF
      # #
      biasDF <- rowDF[BiocGenerics::which(rowDF[, 2] == "bias"), ]
      
      errorList$footMat <- footMat
      errorList$biasMat <- biasMat
      errorList$footDF <- footDF
      errorList$biasDF <- biasDF
      
      # if (!is.null(smoothWindow)) {
      #   footMat <- apply(footMat, 2, function(x) .centerRollMean(x, smoothWindow))
      #   biasMat <- apply(biasMat, 2, function(x) .centerRollMean(x, smoothWindow))
      # }
      
      # if (!is.null(norm.factors)) {
      #   
      # }
      
      if (!is.null(flankNorm)) {
        idx <- which(abs(footDF$x) >= flank - flankNorm)
        footMat <- t(t(footMat)/colMeans(footMat[idx, , drop = FALSE]))
        biasMat <- t(t(biasMat)/colMeans(biasMat[idx, , drop = FALSE]))
        errorList$footMatNorm <- footMat
        errorList$biasMatNorm <- footMat
      }
      
      # if (tolower(normMethod) == "none") {
      #   title <- ""
      # }
      # else if (tolower(normMethod) == "subtract") {
      #   title <- "Tn5 Bias Subtracted\n"
      #   footMat <- footMat - biasMat
      # }
      # else if (tolower(normMethod) == "divide") {
      #   title <- "Tn5 Bias Divided\n"
      #   footMat <- footMat/biasMat
      # }
      # else {
      #   stop("normMethod not recognized!")
      # }
      
      # FANC: add pileupMat. This gives the general trend line
      # pileupMat <- apply(footMat, 2, function(x) .centerRollMean(x, pileup.smoothWindow))
      # printMat <- footMat - pileupMat
      # 
      
      footMatMean <- .groupMeans(footMat, SummarizedExperiment::colData(seFoot)$Group)
      footMatSd <- .groupSds(footMat, SummarizedExperiment::colData(seFoot)$Group)
      biasMatMean <- .groupMeans(biasMat, SummarizedExperiment::colData(seFoot)$Group)
      biasMatSd <- .groupSds(biasMat, SummarizedExperiment::colData(seFoot)$Group)
      #    smoothFoot <- rowMaxs(apply(footMatMean, 2, function(x) .centerRollMean(x, 11))) # FANC: only used for quantile computation later to set y axis limits
      
      # # FANC: added:
      # printMatMean <- .groupMeans(printMat, SummarizedExperiment::colData(seFoot)$Group)
      # printMatSd <- .groupSds(printMat, SummarizedExperiment::colData(seFoot)$Group)
      # pileupMatMean <- .groupMeans(pileupMat, SummarizedExperiment::colData(seFoot)$Group)
      # pileupMatSd <- .groupSds(pileupMat, SummarizedExperiment::colData(seFoot)$Group)
      # #
      
      # errorList$footMatMean <- footMatMean
      # errorList$footMatSd <- footMatSd
      # errorList$biasMatMean <- biasMatMean
      # errorList$biasMatSd <- biasMatSd
      # errorList$smoothFoot <- smoothFoot
      
      plotIdx <- seq_len(nrow(footMatMean))
      
      plotFootDF <- lapply(seq_len(ncol(footMatMean)), function(x) {
        data.frame(x = footDF$x, 
                   mean = footMatMean[, x], 
                   sd = footMatSd[, x], 
                   group = colnames(footMatMean)[x]
        )[plotIdx, , drop = FALSE]
      }) %>% Reduce("rbind", .)
      plotFootDF$group <- factor(paste0(plotFootDF$group), levels = unique(gtools::mixedsort(paste0(plotFootDF$group))))
      
      
      # >     plotFootDF %>% head()
      #      x        mean          sd group
      # 1 -250 -0.03867912 0.002682626 3..KO
      # 2 -249 -0.03867912 0.002682626 3..KO
      # 3 -248 -0.03867912 0.002682626 3..KO
      # 4 -247 -0.03867912 0.002682626 3..KO
      # 5 -246 -0.03867912 0.002682626 3..KO
      # 6 -245 -0.03867912 0.002682626 3..KO
      
      plotBiasDF <- lapply(seq_len(ncol(biasMatMean)), function(x) {
        data.frame(x = biasDF$x, 
                   mean = biasMatMean[, x], 
                   sd = biasMatSd[, x], 
                   group = colnames(biasMatMean)[x]
        )[plotIdx, , drop = FALSE]
      }) %>% Reduce("rbind", .)
      plotBiasDF$group <- factor(paste0(plotBiasDF$group), levels = unique(gtools::mixedsort(paste0(plotBiasDF$group))))
      
      coldata <- colData(seFoot) %>% as.data.frame()
      res <- list(foot.mat = footMat, bias.mat = biasMat, foot.meanSd = plotFootDF, bias.meanSd = plotBiasDF,
                  coldata = coldata)  
      
      return(res)
    }) 
    names(dfs) <- motifs
    return(dfs)
  })
  
  names(dfs) <- clusters
  file <- paste0(out.dir, "/", root, "_dfs.Rds")
  dir.create(out.dir, showWarnings = F, recursive = T)
  saveRDS(dfs, file)
  return(dfs)
}

footprint.center.vs.flank <- function(fpo,
                                      flank.right, flank.left = NULL, 
                                      center) {
# fpo: list(foot.mat, bias.mat, foot.meanSd, bias.meanSd, coldata)
# you get it from footprint.df.gen
  required <- c("foot.mat", "bias.mat", "foot.meanSd", "bias.meanSd", "coldata")
  utilsFanc::check.intersect(required, "required elements", names(fpo), "names(fpo)")
# format of coldata: from ArchR:
#   Browse[1]> coldata
#                 Group            Name
# 0..KO._.KO_rep3 0..KO 0..KO._.KO_rep3
# 0..KO._.KO_rep2 0..KO 0..KO._.KO_rep2
# 0..WT._.WT_rep3 0..WT 0..WT._.WT_rep3
# 0..WT._.WT_rep2 0..WT 0..WT._.WT_rep2
#                                                                                                                                                File
# 0..KO._.KO_rep3 /scratch/fanc/others/jun/2022-01-25/sync_all/atac/merge_v1/GroupCoverages/cluster_type_9mer/X0..KO._.KO_rep3.insertions.coverage.h5
# 0..KO._.KO_rep2 /scratch/fanc/others/jun/2022-01-25/sync_all/atac/merge_v1/GroupCoverages/cluster_type_9mer/X0..KO._.KO_rep2.insertions.coverage.h5
# 0..WT._.WT_rep3 /scratch/fanc/others/jun/2022-01-25/sync_all/atac/merge_v1/GroupCoverages/cluster_type_9mer/X0..WT._.WT_rep3.insertions.coverage.h5
# 0..WT._.WT_rep2 /scratch/fanc/others/jun/2022-01-25/sync_all/atac/merge_v1/GroupCoverages/cluster_type_9mer/X0..WT._.WT_rep2.insertions.coverage.h5
#                 nCells nInsertions
# 0..KO._.KO_rep3    500    31871130
# 0..KO._.KO_rep2    362    25468138
# 0..WT._.WT_rep3    500    36840432
# 0..WT._.WT_rep2    500    41786660
  if (is.null(flank.left)) flank.left <- sort(-1 * flank.right)
  ranges <- list(center = list(center), flank = list(flank.left, flank.right))
  
  sum <- lapply(c("foot", "bias"), function(x) {
    mat <- fpo[[paste0(x, ".mat")]]
    meanSd <- fpo[[paste0(x, ".meanSd")]]
    x.range <- meanSd$x %>% unique() %>% sort()
    if (length(x.range) != nrow(mat)) stop("length(x.range) != nrow(mat)")
    rownames(mat) <- x.range
    sum <- lapply(c("center", "flank"), function(r) {
      sum <- lapply(ranges[[r]], function(r) {
        r <- sort(r)
        rows <- as.character(r[1]:r[2])
        mat <- mat[rows, ]
        sum <- colSums(mat)
        return(sum)
      }) %>% Reduce(`+`, .)
      sum <- t(sum)
      return(sum)
    }) %>% do.call(rbind, .)
    rownames(sum) <- paste0(x, ".", c("center", "flank"))
    return(sum)
  }) %>% do.call(rbind, .)
  sum <- t(sum)
  if (!identical(sort(rownames(sum)), sort(fpo$coldata$Name)))
    stop("!identical(sort(rownames(sum)), sort(fpo$coldata$Name))")
  
  sum <- sum[fpo$coldata$Name,]
  df <- fpo$coldata[, c("Name", "Group")]
  df <- cbind(df, as.data.frame(sum))
  df$foot.protect <- 1 - (df$foot.center / df$foot.flank)
  df$bias.protect <- 1 - (df$bias.center / df$bias.flank)
  df$foot.vs.bias <- round(df$foot.protect / df$bias.protect, digits = 2)
  
  return(df)
}

