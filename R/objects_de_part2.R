de.subset.deg <- function(de, slot = "summary", 
                          method, seed = 42,
                          frac = 0.1, n.min = 1000, n.max = 3000) {
  possible.methods <- c("random", "top")
  if (!method %in% possible.methods) {
    stop("!method %in% possible.methods")
  }
  
  if (method == "random")
    set.seed(seed = seed)
  
  DARs <- lapply(de, function(s2b) {
    if (!slot %in% names(s2b)) {
      stop(paste0("slot '", slot, "' not found"))
    }
    DARs <- lapply(c("up", "down"), function(direction) {
      dar <- s2b[[slot]][[paste0(direction, ".genes")]]
      if (is.null(dar)) {
        stop("is.null(dar)")
      }
      
      len <- length(dar)
      if (len < n.min)
        return(dar)
      
      n.subset.to <- len * frac
      if (n.subset.to < n.min) {
        n.subset.to <- n.min
      }
      if (n.subset.to > n.max) {
        n.subset.to <- n.max
      }
      if (method == "random")
        id <- sample(1:len, size = n.subset.to, replace = F) %>% sort()
      else if (method == "top") {
        id <- 1:n.subset.to
      } else {
        stop("method must be random or top")
      }
      
      if(any(duplicated(id))) {
        stop("any(duplicated(id))")
      }
      
      return(dar[id])
    })
    names(DARs) <- c("up", "down")
    return(DARs)
  })
  return(DARs)
}
