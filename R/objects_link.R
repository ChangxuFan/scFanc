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
                            highlight.idx = NULL, genes = NULL,
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
    select(chr1, start1, end1, chr2, start2, end2, name, score, direction, p, type)
  if (reverse == T) {
    pg <- bedpe$type == "peak-gene"
    gp <- bedpe$type == "gene-peak"
    bedpe[pg, "type"] <- "gene-peak"
    bedpe[gp, "type"] <- "peak-gene"
  }
  bedpe.1 <- bedpe %>% filter(type == "gene-peak") %>%
    mutate(name = sub("(<.+>)(<.+>)", "\\2\\1", name))
  colnames(bedpe.1) <- c("chr2","start2", "end2", "chr1", "start1", "end1", "name", "score", "direction","p", "type")
  bedpe.2 <- bedpe %>% filter(type == "peak-gene")
  dorc <- rbind(bedpe.2, bedpe.1) %>% arrange(chr2, start2, end2) %>% 
    mutate(gene = sub("(<.+>)(<.+>)", "\\2", name) %>% sub("<", "",.) %>% sub(">", "",.)) %>% 
    select(chr1, start1, end1, chr2, start2, end2, gene, name, type, score, direction, p)
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

da.diag.link <- function(da.gr.list, cmp, de.grid.sum, dorc, cluster.ident,
                         plot.out, threads = 6, de.p.filter=1, de.p.adj.filter) {
  sample.y <- sub(":.+$", "", cmp)
  sample.x <- sub("^.+:", "", cmp)
  p.list <- mclapply(seq_along(da.gr.list), function(i) {
    da <- da.gr.list[[i]]
    cluster <- names(da.gr.list)[i] %>% sub(paste0(cluster.ident,"_"), "",.)
    de <- de.grid.sum %>% filter(master.ident == cluster, comp == comp,
                                 p <= de.p.filter, p.adj <= de.p.adj.filter)
    y.up.genes <- de %>% filter(logFC > 0) %>% pull(gene)
    y.up.link <- dorc %>% plyranges::filter(gene %in% y.up.genes) %>% reduce()
    y.up.hits <- findOverlaps(y.up.link, da, ignore.strand = T) %>% subjectHits()
    
    x.up.genes <- de %>% filter(logFC < 0) %>% pull(gene)
    x.up.link <- dorc %>% plyranges::filter(gene %in% x.up.genes) %>% reduce()
    x.up.hits <- findOverlaps(x.up.link, da, ignore.strand = T) %>% subjectHits()
    hit.names <- c("y", "x")
    hits <- list(y.up.hits, x.up.hits)
    da.df <- da %>% `names<-`(NULL) %>% as.data.frame()
    p.list <- lapply(seq_along(hits), function(i) {
      hit <- hits[[i]]
      name <- hit.names[i]
      da.df$highlight <- NA
      da.df$highlight[hit] <- "highlight"
      p <- xy.plot(df = da.df, x = paste0(cluster.ident, "_", cluster,"..", sample.x),
                   y =  paste0(cluster.ident, "_", cluster,"..", sample.y), highlight.var = "highlight", 
                   highlight.values = "highlight") +
        ggtitle(paste0(cluster.ident, "_", cluster, "; ", name, ".up"))
      return(p)
    })
    return(p.list)
  }, mc.cores = threads, mc.cleanup = T) %>% Reduce(c, .)
  trash <- wrap.plots.fanc(p.list, plot.out = plot.out)
  return()
}


