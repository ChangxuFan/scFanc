de.da.intersect <- function(pbl1, pbl2, name1, name2, 
                            slot = "summary", clusters = NULL, 
                            rule = "promoter", promoter.ext.1side = 1000,
                            out.dir, root.name = NULL) {
  # pbl1 and pbl2 could be either de or da.
  # when both de: intersect by gene name
  # when de and da: convert gene name to loci, then perform loci intersection (GRanges based)
  # when both da: loci intersection (GRanges based)
  # rule is used to instruct how to tell if a DEG is intersecting a DAR.
  if (is.null(root.name)) {
    root.name <- basename(out.dir)
  }
  all.clusters <- intersect(names(pbl1), names(pbl2))
  if (length(all.clusters) < 1) {
    stop("length(all.clusters) < 1")
  }
  
  if (is.null(clusters)) {
    clusters <- all.clusters
  }
  clusters <- clusters %>% .[.%in% all.clusters]
  if (length(clusters) < 1) {
    stop("length(clusters) < 1")
  }
  pbl1 <- pbl1[clusters]
  pbl2 <- pbl2[clusters]
  
  des <- list(pbl1, pbl2)
  names(des) <- c(name1, name2)
  
  type <- sapply(des, function(de) {
    gn <- rownames(de[[1]]$bulk.mat)[1]
    return(ifelse(grepl(":", gn), "da", "de"))
  })
  bG <- which(type == "de")
  bP <- which(type == "da")
  stats <- lapply(clusters, function(cluster) {
    if (sum(type == "de") > 0 && sum(type == "da") > 0) {
      if (rule == "promoter") {
        if (length(des) != 2) stop("length(des) != 2")
        
        stats <- lapply(c("up", "down"), function(direc) {
          pros <- des[[bG]][[cluster]][[slot]]$promoters[[direc]]
          if (is.null(pros)) {
            stop(paste0("promoter slot empty for ", names(des)[bG], 
                        ", cluster ", cluster,", slot ", slot))
          }
          peaks <- des[[bP]][[cluster]][[slot]][[paste0(direc, ".genes")]]
          if (is.null(peaks)) {
            stop(paste0("peak slot empty for ", names(des)[bP], 
                        ", cluster ", cluster,", slot ", slot))
          }
          
          pros <- GenomicRanges::resize(pros, width = 2 * promoter.ext.1side, 
                                        fix = "center")
          peaks <- peaks %>% utilsFanc::loci.2.df(loci.vec = ., remove.loci.col = F, return.gr = T)
          
          o <- plyranges::join_overlap_inner(pros, peaks)
          utilsFanc::write.zip.fanc(
            o, out.file = paste0(out.dir, "/", root.name, "_", cluster, "_", direc, "_overlap.bed"), 
            bed.shift = T)
          
          stats <- data.frame(cluster = cluster, direction = direc,
                              n.de = des[[bG]][[cluster]][[slot]][[paste0("n.", direc)]],
                              n.pro = pros$gene %>% unique() %>% length(),
                              n.da = length(peaks),
                              n.de.w.da = o$gene %>% unique() %>% length()) %>% 
            dplyr::mutate(frac.de.w.da = round(n.de.w.da/n.de, digits = 3))
          return(stats)
        }) %>% do.call(rbind, .)
        return(stats)
        
      } else if (rule == "nearestGene") {
        stats <- lapply(c("up", "down"), function(direc) {
          anno <- des[[bP]][[cluster]][[slot]]$anno[[direc]]
          if (is.null(anno)) {
            stop(paste0("anno slot empty for ", names(des)[bP], 
                        ", cluster ", cluster,", slot ", slot))
          }
          deg <- des[[bG]][[cluster]][[slot]][[paste0(direc, ".genes")]]
          if (is.null(deg)) {
            stop(paste0("deg slot empty for ", names(des)[bG], 
                        ", cluster ", cluster,", slot ", slot))
          }
          
          o <- anno[anno$gene %in% deg]
          utilsFanc::write.zip.fanc(
            o, out.file = paste0(out.dir, "/", root.name, "_", rule, "_", cluster, "_", direc, "_overlap.bed"), 
            bed.shift = T)
          
          stats <- data.frame(
            cluster = cluster, direction = direc,
            n.de = des[[bG]][[cluster]][[slot]][[paste0("n.", direc)]],
            n.da = length(anno),
            n.da.w.anno = sum(!is.na(anno$gene)),
            n.da.w.de = anno$gene %>% .[!is.na(.)] %>% .[. %in% deg] %>% unique() %>% length(),
            n.de.w.da = deg %>% .[. %in% anno$gene] %>% unique() %>% length()
          ) %>% dplyr::mutate(
            frac.da.w.de = round(n.da.w.de/n.da, digits = 3),
            frac.de.w.da = round(n.de.w.da/n.de, digits = 3)
          )
        
          return(stats)
        }) %>% do.call(rbind, .)
        return(stats)
      } else {
        stop("rule can be: promoter")
      }
    } else stop("only de vs da has been developed")
  }) %>% do.call(rbind, .)
  dir.create(out.dir, showWarnings = F, recursive = T)
  write.table(stats, paste0(out.dir, "/", root.name, "_", rule, "_stats.tsv"), sep = "\t", quote = F,
              row.names = F, col.names = T)
  invisible(stats)
}

