addSlingShotTrajectories.fanc <- function(ArchRProj, curve.object, name = "SlingShot", force = F) {
  pt <- slingPseudotime(curve.object)
  colnames(pt) <- paste0(name, ".Curve", seq_len(ncol(pt)))
  ptn <- apply(pt, 2, ArchR:::.getQuantiles) * 100
  for (i in seq_len(ncol(ptn))) {
    ArchRProj <- addCellColData(ArchRProj = ArchRProj, data = as.vector(ptn[,i]),
                                name = colnames(ptn)[i], cells = rownames(ptn), 
                                force = force)
  }
  ArchRProj
}


