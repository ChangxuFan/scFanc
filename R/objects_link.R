read10x.link <- function(bedpe, filter = F, pos.only = F) {
  if (is.character(bedpe)) {
    bedpe.ori <- bedpe
    bedpe <- read.table(bedpe, sep = "\t", quote = "", as.is = T)
    colnames(bedpe) <- c("chr1", "start1", "end1", "chr2", "start2", "end2",
                         "name", "score", "strand1", "strand2", "p", "dist", "type")
  }
  if (filter == T)
    bedpe <- bedpe %>% filter(type %in% c("peak-gene", "gene-peak"))
  if (pos.only == T)
    bedpe <- bedpe %>% filter(score > 0)
  bedpe$direction <- NA
  bedpe$direction[bedpe$score > 0] <- "pos_corr"
  bedpe$direction[bedpe$score < 0] <- "neg_corr"
  return(bedpe)
  
}
join.peaks.10x <- function(peaks.10x.bed, out.bed = NULL) {
  peaks <- lapply(peaks.10x.bed, rtracklayer::import.bed)
  peak.j <- Reduce(c, peaks)
  peak.m <- reduce(peak.j)
  if (!is.null(out.bed)) {
    system(paste0("mkdir -p ", dirname(out.bed)))
    rtracklayer::export.bed(object = peak.m, con = out.bed)
    system(paste0("/bar/cfan/scripts/bed_browser_v2.sh ", out.bed))
  }
  return(peak.m)
}
qc.10x.link <- function(bedpe, qc.dir, nCRE.quantile = 0.99, nGene.quantile = 0.99) {
  system(paste0("mkdir -p ", qc.dir))
  bedpe <- read10x.link(bedpe = bedpe, filter = T, pos.only = F)
  # distance stats
  t.f.qc <- function(df, x, title, quantile, qc.dir) {
    title <- gsub(" +", "_", title)
    max.quantile <- df[, x] %>% quantile(quantile)
    p.density <- df %>% ggplot(aes(x= !!as.name(x), fill = direction)) + 
      geom_density(alpha = 0.2) + xlim(c(0, max.quantile)) +
      ggtitle(paste0(title, "; density", ", quantile ", quantile)) +
      ggsave(paste0(qc.dir, "/", title, ".density.quantile.", quantile,".png" ))
    p.rank.pos <- df %>% filter(direction == "pos_corr") %>% 
      rank.plot(vars = x, rank.method = "random") + ylim(c(0, max.quantile)) +
      ggtitle(paste0(title, "; rank. pos", ", quantile ", quantile)) +
      ggsave(paste0(qc.dir, "/", title, ".rank.pos.quantile.", quantile,".png" ))
    p.rank.neg <- df %>% filter(direction == "neg_corr") %>% 
      rank.plot(vars = x, rank.method = "random") + ylim(c(0, max.quantile)) +
      ggtitle(paste0(title, "; rank. neg", ", quantile ", quantile)) +
      ggsave(paste0(qc.dir, "/", title, ".rank.neg.quantile.", quantile,".png" ))
    return(list(p.density, p.rank.pos, p.rank.neg))
  }
  
  
  p.list.dist <- t.f.qc(df = bedpe, x = "dist", quantile = 1, qc.dir = qc.dir, title = "distance") 
  
  dorc <- bedpe.2.dorc(bedpe = bedpe, reverse = F, gr = F)
  n.CRE <- dorc %>% group_by(gene, direction) %>% summarise(n = n()) %>% as.data.frame()
  p.list.nCRE <- lapply(c(1, nCRE.quantile), function(q) {
    t.f.qc(df = n.CRE, x = "n", title = "n.CRE", quantile = q, qc.dir = qc.dir) %>% return()
  }) %>% Reduce(c, .)
  
  dorc.r <- bedpe.2.dorc(bedpe = bedpe, reverse = T, gr = F)
  n.CRE.r <- dorc.r %>% group_by(gene, direction) %>% summarise(n = n()) %>% as.data.frame()
  p.list.nGene <- lapply(c(1, nGene.quantile), function(q) {
    t.f.qc(df = n.CRE.r, x = "n", title = "n.Gene", quantile = q, qc.dir = qc.dir) %>% return()
  }) %>% Reduce(c, .)
  
  p.p <- xy.plot(df = bedpe %>% filter(score > 0), x = "dist", y = "p") +
    ggtitle("p_vs_dist") +
    ggsave(paste0(qc.dir, "/p_vs_dist.png"))
  p.p <- list(p.p)
  
  p.list <- Reduce(c, list(p.list.dist, p.list.nCRE, p.list.nGene, p.p))
  p.master <- wrap.plots.fanc(p.list, plot.out = paste0(qc.dir, "/qc.png"), n.col = 6)
  return()
}


bedpe.2.browser <- function(bedpe, filter = F, pos.only = T, use.score = F, out.file = NULL,
                            highlight.idx = NULL, highlight.only = F, genes = NULL,
                            # the viriables below are for peak watching
                            watch.id = NULL, n, min=NULL, max=NULL, seed = 42 ) {
  if (is.character(bedpe)) {
    bedpe.ori <- bedpe
    bedpe <- read.table(bedpe, sep = "\t", quote = "", as.is = T)
    colnames(bedpe) <- c("chr1", "start1", "end1", "chr2", "start2", "end2",
                         "name", "score", "strand1", "strand2", "p", "dist", "type")
  }
  
  if (use.score == T)
    bedpe$p <- bedpe$score
  if (filter == T)
    bedpe <- bedpe %>% filter(type %in% c("peak-gene", "gene-peak"))
  if (pos.only == T)
    bedpe <- bedpe %>% filter(score > 0)
  
  if (!is.null(watch.id)) {
    highlight.idx <- peak.watch.core(df = bedpe, id = watch.id, max = max, 
                                     min = min, n = n, seed = seed, return.idx = T)
  }
  
  if (!is.null(genes)) {
    gp <- bedpe$type == "gene-peak"
    pg <- bedpe$type == "peak-gene"
    bedpe$gene <- NA
    bedpe[gp, "gene"] <- bedpe[gp, "name"] %>% sub("(<.+>)(<.+>)", "\\1",.) %>% 
      sub("<", "",.) %>% sub(">", "",.)
    bedpe[pg, "gene"] <- bedpe[pg, "name"] %>% sub("(<.+>)(<.+>)", "\\2",.) %>% 
      sub("<", "",.) %>% sub(">", "",.)
    highlight.idx <- which(bedpe$gene %in% genes)
  }
  
  # print(highlight.idx)
  if (!is.null(highlight.idx)) {
    bedpe$p <- abs(bedpe$p)
    bedpe$p[highlight.idx] <- -1 * bedpe$p[highlight.idx]
  }
  if (highlight.only == T) {
    bedpe <- bedpe[highlight.idx,]
  }

  half.1 <- bedpe %>% mutate(add = paste0(chr2, ":", start2, "-", end2, ",", p)) %>% 
    select(chr1, start1, end1, add)
  colnames(half.1) <- c("chr", "left", "right", "int")
  
  half.2 <- bedpe %>% mutate(add = paste0(chr1, ":", start1, "-", end1, ",", p)) %>% 
    select(chr2, start2, end2, add)
  colnames(half.2) <- c("chr", "left", "right", "int")
  lr <- rbind(half.1, half.2)
  if (is.null(out.file)) {
    out.file <- bedpe.ori %>% paste0(".lr")
  } 
  system(paste0("mkdir -p ", dirname(out.file)))
  write.table(lr, out.file, sep = "\t", row.names = F, col.names = F, quote = F)
  cmd <- paste0("/bar/cfan/scripts/bed_browser_v2.sh ", out.file)
  print(cmd); system(cmd)
  
  if (!is.null(highlight.idx)) {
    watch.tsv <- paste0(out.file, ".watch.tsv")
    write.table(bedpe[highlight.idx,], watch.tsv, sep = "\t", row.names = F, col.names = F, quote = F)
    # cmd <- paste0("/bar/cfan/scripts/bed_browser_v2.sh ", watch.tsv)
    # print(cmd); system(cmd)
  }
  
  return(utilsFanc::bash2ftp(paste0(out.file, ".gz")))
  
}
# parse.int.name <- function(int.names) {
#   df <- data.frame(gene.1.proto = sub("(<.+>)(<.+>)", "\\1", int.names) %>% sub("<", "",.) %>% sub(">", "", .),
#                    gene.2.proto = sub("(<.+>)(<.+>)", "\\2", int.names) %>% sub("<", "",.) %>% sub(">", "", .)) %>% 
#     mutate(gene.1.proto = sapply(gene.1.proto, function(x) {
#       
#     }))
#   
#   protos <- lapply(1:2, function(i) {
#     proto <- sub("(<.+>)(<.+>)", paste0("\\", i), int.names) %>% sub("<", "",.) %>% sub(">", "", .)
#     pro <- !grepl("_",proto)
#     proto[pro] <- paste0(proto[pro], "_TSS")
#     name <- sub("(.*)_(.*)", "\\1", proto)
#     type <- sub("(.*)_(.*)", "\\2", proto)
#     return(proto)
#   })
#   
# }

bedpe.2.dorc <- function(bedpe, out.bed = NULL, pos.only = F, reverse = F, gr = F) {
  bedpe <- read10x.link(bedpe, pos.only = pos.only) %>%
    select(chr1, start1, end1, chr2, start2, end2, name, score, direction, dist, p, type)
  if (reverse == T) {
    pg <- bedpe$type == "peak-gene"
    gp <- bedpe$type == "gene-peak"
    bedpe[pg, "type"] <- "gene-peak"
    bedpe[gp, "type"] <- "peak-gene"
  }
  bedpe.1 <- bedpe %>% filter(type == "gene-peak") %>%
    mutate(name = sub("(<.+>)(<.+>)", "\\2\\1", name))
  colnames(bedpe.1) <- c("chr2","start2", "end2", "chr1", "start1", "end1", "name", "score", "direction","p", "dist", "type")
  bedpe.2 <- bedpe %>% filter(type == "peak-gene")
  dorc <- rbind(bedpe.2, bedpe.1) %>% arrange(chr2, start2, end2) %>% 
    mutate(gene = sub("(<.+>)(<.+>)", "\\2", name) %>% sub("<", "",.) %>% sub(">", "",.)) %>% 
    select(chr1, start1, end1, chr2, start2, end2, gene, name, type, score, direction, dist,p)
  if (!is.null(out.bed)) {
    bed <- dorc %>% mutate(V4 = paste0(chr2, ":", start2, "-", end2, ";", gene,";",score, ";",p))
    write.table(bed, out.bed, col.names = F, row.names = F, quote = F, sep = "\t")
    cmd <- paste0("/bar/cfan/scripts/bed_browser_v2.sh ", out.bed)
    print(cmd); system(cmd)
  } 
  if (gr == T) {
    dorc <- dorc %>% rename(chr = chr1, start = start1, end = end1) %>% 
      mutate(prox = paste0(chr2, ":", start2, "-", end2)) %>% 
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  }
  return(dorc)
}

bedpe.2.dorc.m <- function(bedpe.vec, sample.names=NULL, pos.only,
                           merge = F, reverse = F) {
  if (is.null(sample.names))
    sample.names <- names(bedpe.vec)
  if (is.null(sample.names))
    stop("sample.names must be given or contained in the names of bedpe.vec")
  dorc.list <- lapply(seq_along(bedpe.vec), function(i) {
    bedpe <- bedpe.vec[i]
    dorc <- bedpe.2.dorc(bedpe = bedpe, reverse = reverse, gr = T, pos.only = pos.only)
    dorc$sample = sample.names[i]
    return(dorc)
  })
  names(dorc.list) <- sample.names
  return(dorc.list)
}

tss.gex.gen <- function(s2b.list, tss.df = TSS.CELLRANGER.RDS, samples, master.dir, 
                        thread = 6, zip.only = F) {
  if (is.character(s2b.list))
    s2b.list <- readRDS(s2b.list)
  if (is.character(tss.df))
    tss.df <- readRDS(tss.df)

  system(paste0("mkdir -p ", master.dir))
  # re-do 
  if (zip.only == F) {
    bdgs <- mclapply(seq_along(s2b.list), function(i) {
      s2b <- s2b.list[[i]]
      cluster <- names(s2b.list)[i]
      lapply(samples, function(sample) {
        gex.df <- data.frame(GEX = s2b$bulk.mat[, sample],
                             gene = rownames(s2b$bulk.mat), row.names = NULL)
        gex.df$gene <- rownames(s2b$bulk.mat)
        bdg <- left_join(tss.df, gex.df) %>% select("chr", "start", "end", "GEX", "gene")
        bdg.file <- paste0(master.dir, "/", sample, "..", cluster, ".bdg")
        write.table(bdg, bdg.file, sep = "\t", quote = F, row.names = F, col.names = F)
        cmd <- paste0("/bar/cfan/scripts/bed_browser_v2.sh ", bdg.file)
        print(cmd); system(cmd)
        return(NULL)
      })
      return()
    }, mc.cleanup = T, mc.cores = thread)    
  }
  cmd <- paste0("/bar/cfan/R_for_bash/json_dir2.R -d", master.dir, " -f ", "\"bdg.gz$\"")
  print(cmd); system(cmd)
  return(utilsFanc::bash2ftp(paste0(master.dir, "/files.json")))

}

da.diag.link <- function(da.gr.list, cmp, de.grid.sum, dorc, cluster.ident, clusters = NULL,
                         plot.out=NULL, threads = 6, de.p.filter=1, de.p.adj.filter,
                         return.df = F, return.p = F, return.p.list = F,
                         violin = F, max.quantile = 1, min.size = 0,
                         min.quantile = 0, violin.color.by.sample = F) {
  sample.y <- sub(":.+$", "", cmp)
  sample.x <- sub("^.+:", "", cmp)
  p.list <- mclapply(seq_along(da.gr.list), function(i) {
    da <- da.gr.list[[i]]
    cluster <- names(da.gr.list)[i] %>% sub(paste0(cluster.ident,"_"), "",.)
    if (!is.null(clusters)) {
      if (!cluster %in% clusters)
        return(NULL)
    }
    de <- de.grid.sum %>% filter(master.ident == cluster, comp == comp,
                                 p <= de.p.filter, p.adj <= de.p.adj.filter)
    y.up.genes <- de %>% filter(logFC > 0) %>% pull(gene)
    y.up.link <- dorc %>% plyranges::filter(gene %in% y.up.genes) %>% reduce()
    y.up.hits <- findOverlaps(y.up.link, da, ignore.strand = T) %>% subjectHits()
    
    x.up.genes <- de %>% filter(logFC < 0) %>% pull(gene)
    x.up.link <- dorc %>% plyranges::filter(gene %in% x.up.genes) %>% reduce()
    x.up.hits <- findOverlaps(x.up.link, da, ignore.strand = T) %>% subjectHits() %>% unique()
    hit.names <- c("y", "x")
    hits <- list(y.up.hits, x.up.hits)
    da.df <- da %>% `names<-`(NULL) %>% as.data.frame()
    if (violin == T) {
      da.df$type <- "other"
      da.df$type[y.up.hits] <- paste0(sample.y, ".up")
      da.df$type[x.up.hits] <- paste0(sample.x, ".up")
      da.df.melt <- da.df %>% mutate(peak = paste0(seqnames, ":", start, "-", end)) %>% 
        reshape2::melt(id.vars = c("peak","type"),
                       measure.vars = paste0(cluster.ident, "_", cluster, "..",c(sample.y, sample.x)),
                       variable.name = "sample", value.name = "size") %>% 
        mutate(sample = sub(paste0(cluster.ident, "_"), "", sample)) %>% 
        mutate(sample = factor(sample, levels = paste0(cluster, "..", c(sample.y, sample.x))), 
               type = factor(type, levels = c("other", paste0(c(sample.y, sample.x), ".up"))))
      max.size <- quantile(da.df.melt$size, max.quantile)
      min.size <- max(min.size, quantile(da.df.melt$size, min.quantile))
      da.df.melt <- da.df.melt %>% filter(size <= max.size, size>=min.size)
      violin.x <- "sample"
      violin.color <- "type"
      if (violin.color.by.sample == T) {
        violin.x <- "type"
        violin.color <- "sample"
      }
      p.list <- ggpubr::ggviolin(da.df.melt, x = violin.x,y = "size",
                                 color = violin.color, title = paste0(cluster.ident, "_", cluster), 
                                 #ylim = c(min.size, max.size)
                                 ) %>% list()
      return(p.list)
    } else {
      p.list <- lapply(seq_along(hits), function(i) {
        hit <- hits[[i]]
        name <- hit.names[i]
        da.df$highlight <- NA
        da.df$highlight[hit] <- "highlight"
        if (return.df == T) {
          p <- da.df
        } else {
          p <- xy.plot(df = da.df, x = paste0(cluster.ident, "_", cluster,"..", sample.x),
                       y =  paste0(cluster.ident, "_", cluster,"..", sample.y), highlight.var = "highlight", 
                       highlight.values = "highlight") +
            ggtitle(paste0(cluster.ident, "_", cluster, "; ", name, ".up"))
        }
        
        return(p)
      })
      names(p.list) <- paste0(cluster.ident, "_", cluster, "..", c("y", "x"), ".up")
    }
    return(p.list)
  }, mc.cores = threads, mc.cleanup = T) %>% Reduce(c, .)
  if (return.df == T) {
    return(p.list)
  } else {
    p <- wrap.plots.fanc(p.list, plot.out = plot.out)
    if (return.p == T)
      return(p) 
    
    if (return.p.list == T)
      return(p.list)
    else
      return()
  }

}

t.f.titrate.link <- function(da.gr.list, de.grid.sum, dorc, out.dir, threads.ld = 4, threads.lp = 4,
                             clusters, peak.size.quantiles, link.p.quantiles, link.dist.quantiles) {
  lapply(clusters, function(cluster) {
    lapply(peak.size.quantiles, function(ps.q) {
      p.list <- mclapply(link.p.quantiles, function(lp.q) {
        mclapply(link.dist.quantiles, function(ld.q) {
          # browser()
          dorc.filtered <- dorc %>% plyranges::filter(p > quantile(p, lp.q), dist < quantile(dist, ld.q))
          if (length(dorc.filtered) < 1)
            return()
          # print("miao")
          p <- da.diag.link(da.gr.list = da.gr.list, cmp = "PyMT:WT", de.grid.sum = de.grid.sum, 
                            dorc = dorc.filtered, cluster.ident = "seurat_clusters", clusters = cluster, 
                            de.p.adj.filter = 0.05, return.df = F, return.p.list = T,
                            plot.out = NULL, violin = T, threads = 6, max.quantile = 0.999,
                            min.quantile = ps.q, violin.color.by.sample = T)[[1]]
          p <- p + ggtitle(paste0("link p cutoff: ", lp.q, "\n", "link dist cutoff: ", ld.q))
          return(p)
        }, mc.cores = threads.ld, mc.cleanup = T) %>%  return() 
      }, mc.cores = threads.lp, mc.cleanup = T) %>% Reduce(c, .) %>% return() 
      # print("miao")
      trash <- wrap.plots.fanc(plot.list = p.list, n.col = length(link.p.quantiles),
                               plot.out = paste0(out.dir,"/seurat_clusters_", cluster , "/size_q..", ps.q, ".png"))
      return()
    })
  })
}

find.int.genes <- function(loci.vec, dorc) {
  loci.gr <- utilsFanc::loci.2.df(loci.vec = loci.vec, return.gr = T, remove.loci.col = T)
  join.gr <- plyranges::join_overlap_left(x = loci.gr, y = dorc)
  
  return(join.gr)
}

joined.peakwatch <- function(a2b.list, s2b.list, dorc, grid.sum,comp, bedpe.vec=NULL, 
                             bedpe.name.vec = NULL, bedpe.out.dir = NULL, highlight.only = F,
                             pw.out.dir, pw.root.name = NULL,
                             slot.name = "peakwatch", force = F,
                             id = "pvalue", sort.by = NULL,
                             cluster.ident, threads = 6, ...) {
  # ... will be passed to peak.watch.core
  dorc <- dorc[, !colnames(mcols(dorc)) %in% c("chr2", "start2", "end2")]
  dorc <- dorc %>% utilsFanc::add.column.fanc(df1 = dorc, df2 = data.frame(distal = utilsFanc::gr.get.loci(dorc)),
                                              is.gr = T, pos = 1)
  clusters <- Reduce(intersect, 
                     list(names(a2b.list), names(s2b.list))) %>% 
    sub(paste0(cluster.ident, "_"), "", .)
  a2b.list <- mclapply(clusters, function(cluster) {
    a2b <- a2b.list[[paste0(cluster.ident, "_",cluster)]]
    s2b <- s2b.list[[paste0(cluster.ident, "_",cluster)]]
    if (is.null(pw.root.name))
      pw.root.name <- slot.name
    out.bed <- paste0(pw.out.dir, "/", pw.root.name, "..", cluster.ident, "_", cluster, ".bed")
    out.pw <- paste0(pw.out.dir, "/", pw.root.name, "..", cluster.ident, "_", cluster, ".tsv")
    if (slot.name %in% names(a2b) && force != T)
      warning("indicated peakwatch slot already exist. skipping peakwatch generation")
    else
      a2b[[slot.name]] <- peak.watch.core(df = a2b$res.exp, id = id, return.df.watch = T, 
                                          out.bed = out.bed,...)
    pw <- utilsFanc::loci.2.df(df = a2b[[slot.name]], loci.col.name = "gene", 
                                             remove.loci.col = T, return.gr = T)
    pw <- plyranges::join_overlap_left(x = pw, y = dorc)
    # 
    pw <- pw %>% `names<-`(NULL) %>% as.data.frame() 
    gex.bulk <- s2b$bulkNorm
    colnames(gex.bulk) <- sub("bulkNorm", "gex", colnames(gex.bulk))
    # rownames(gex.bulk) <- gex.bulk$gene
    # gex.bulk <- gex.bulk[pw$gene,]
    pw <- left_join(pw, gex.bulk)
    gex.stats <- grid.sum %>% filter(master.ident == cluster, comp == comp) %>% 
      select(gene, p, logFC, pct.1, pct.2, p.adj)
    colnames(gex.stats)[2:ncol(gex.stats)] <- paste0("gex_", colnames(gex.stats)[2:ncol(gex.stats)])
    pw <- left_join(pw, gex.stats)
    pw <- makeGRangesFromDataFrame(pw, keep.extra.columns = T)
    a2b[[slot.name]] <- pw
    pw.df <- as.data.frame(pw %>% `names<-`(NULL))
    if (!is.null(sort.by)){
      pw.df <- pw.df %>% arrange(stratify, !!as.name(sort.by))
    }
    utilsFanc::write.zip.fanc(pw.df, out.file = out.pw, zip = F,
                              col.names = T)
    return(a2b)
  }, mc.cores = threads, mc.cleanup = T)
  
  if (!is.null(bedpe.vec)) {
    
    if (is.null(bedpe.name.vec)) {
      bedpe.name.vec <- names(bedpe.vec)
      if (is.null(bedpe.name.vec))
        stop("sample names must be provided for bedpe.vec")
    }
    pw.all <- lapply(a2b.list, function(x) return(x[[slot.name]])) %>% Reduce(c,.) %>% 
      as.data.frame() %>% select(distal, prox) %>% na.omit() %>% unique()
    trash <- lapply(seq_along(bedpe.vec), function(i) {
      bedpe.file <- bedpe.vec[i]
      bedpe.name <- names(bedpe.vec)[i]
      bedpe <- read10x.link(bedpe = bedpe.file, filter = T, pos.only = T)
      ids <- lapply(1:2, function(m) {
        i <- 3-m
        bedpe$distal <- paste0(bedpe[,paste0("chr", i)], ":",bedpe[,paste0("start", i)] ,"-",
                               bedpe[,paste0("end", i)])
        bedpe$prox <- paste0(bedpe[,paste0("chr", m)], ":",bedpe[,paste0("start", m)] ,"-",
                             bedpe[,paste0("end", m)])
        bedpe <- bedpe[, c("distal", "prox")]
        bedpe$id <- 1:nrow(bedpe)
        j <- left_join(pw.all, bedpe) %>% pull(id) %>% .[!is.na(.)] %>% unique()
        return(j)
      }) %>% Reduce(c, .)
      out.bedpe <- paste0(bedpe.out.dir, "/", bedpe.name, "_", slot.name,".lr")
      bedpe.2.browser(bedpe = bedpe.file, filter = T, pos.only = T, use.score = F, 
                      out.file = out.bedpe, highlight.only = highlight.only,
                      highlight.idx = ids)
      print(utilsFanc::bash2ftp(paste0(out.bedpe, ".gz")))
      return()
    })
  }
  return(a2b.list)
}