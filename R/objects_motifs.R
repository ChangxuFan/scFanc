motif.2.genes <- function(motifs, species, map = ARCHETYPE.MAP) {
  # the motifs here are actually archetypes
  if (is.character(map))
    map <- read.table(map, header = T, quote = "")
  genes.list <- lapply(motifs, function(motif) {
    genes <- map %>% filter(Name == motif) %>% pull(gene) %>% 
      strsplit(split = "\\+") %>% unlist() %>% unique()
    return(genes)
  })
  
  if (species == "mouse") {
    genes.list <- lapply(genes.list, stringr::str_to_title)
  } else if (species == "human") {
    genes.list <- lapply(genes.list, toupper)
  }
  names(genes.list) <- motifs
  return(genes.list)
  
}

arche.2.motifs <- function(arche.names, map = ARCHETYPE.MAP, enforce.list = F) {
  if (is.character(map))
    map <- read.table(map, header = T, quote = "")
  
  res <- lapply(arche.names, function(arche) {
    motifs <- map %>% filter(Name == arche) %>% pull(Motif)
    return(motifs)
  })
  names(res) <- arche.names
  if (enforce.list == F && length(res) == 1)
    res <- unlist(res)
  return(res)
}