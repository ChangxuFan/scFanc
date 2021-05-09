to.palantir <- function(mat, outdir, start.cell, extrema.df = NULL, other.params = "", run = T) {
  system(paste0("mkdir -p ", outdir))
  if (!is.null(extrema.df)) {
    write.table(extrema.df, paste0(outdir, "/terminal_states_lookup.txt") , sep = "\t", col.names = T, 
                row.names = F, quote = F)
  }
  write(start.cell, "/start_cell.txt" %>% paste0(outdir, .))
  
  write.table(x = as.data.frame(mat), file = "/mat.tsv" %>% paste0(outdir, .), sep = "\t", col.names = T, 
              row.names = T, quote = F)
  cmd <- paste0("cd ", outdir, " && /bar/cfan/anaconda2/envs/jupyter/bin/python ",
                "~/R_packages/scFanc/python/palantirPipe.py", " ", other.params)
  print(cmd)
  if (run == T)
    system(cmd)
  
  
}

load.palantir <- function(dir, extrema.lookup, mask.df, mask.threashold) {
  # extrema.lookup: dataframe, $curve and $cell
  # mask.df: NA means a cell does not belong to a certain trajectory
  # note: the order of extrema.lookup$curve must be the same as the colnames of mask.df
  
  n.curves <- nrow(extrema.lookup)
  extrema.lookup <- extrema.lookup %>% mutate(cell = gsub("[^A-Za-z0-9]", ".", cell))
  rownames(extrema.lookup) <- extrema.lookup$cell
  
  branch_probs <- read.csv(paste0(dir, "/branch_probs.csv"), header = T, row.names = 1)
  
  colnames(branch_probs) <- extrema.lookup[colnames(branch_probs), "curve"]
  branch_probs <- branch_probs[, extrema.lookup$curve]
  branch_probs$cell <- rownames(branch_probs)
  
  pseudotime <- read.csv(paste0(dir, "/pseudotime.csv"))
  
  pseudotime.df <- rep(list(pseudotime[,2]), n.curves) %>% `names<-`(extrema.lookup$curve) %>% 
    as.data.frame() %>% `rownames<-`(branch_probs$cell)
  
  mask.df[mask.df < mask.threashold] <- NA
  mask.df[!is.na(mask.df)] <- 1
  
  pseudotime.df.masked <- pseudotime.df * mask.df
  
  entropy <- read.csv(paste0(dir, "/entropy.csv")) %>% `colnames<-`(c("cell", "entropy"))
  
  waypoints <- read.csv(paste0(dir, "/waypoints.csv"))[,1] 
  
  res <- list(n.curves = n.curves, extrema.lookup = extrema.lookup, branch_probs = branch_probs, 
              pseudotime = pseudotime,
              pseudotime.df = pseudotime.df, pseudotime.df.masked = pseudotime.df.masked, 
              entropy = entropy, waypoints = waypoints)
  
  return(res)
}

palantir.rank.plot <- function(df, cluster.ident, group.by, vars, outdir) {
  ## first just do density plot
  clusters <- df[, cluster.ident] %>% unique() %>% gtools::mixedsort()
  trash <- mclapply(clusters, function(cluster) {
    df.sub <- df[df[,cluster.ident] == cluster,]
    d.p.list <- mclapply(vars, function(var) {
      p <- ggplot(df.sub, aes_string(x = var, fill = group.by)) +
        geom_density(alpha = 0.4) +
        ggtitle(paste0(cluster, "  ", var))
    })

    trash <- scFanc::wrap.plots.fanc(plot.list = d.p.list, plot.out = paste0(outdir, "/", cluster,".png"))
    return(NULL)
  }, mc.cores = 6)
  
  return()
}