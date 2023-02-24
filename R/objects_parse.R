# mostly parse 10x statistics
cellranger.parse <- function(sample.tsv, samples.exclude = NULL, out.root.name) {
  df <- lapply(sample.tsv, function(x) return(read.table(x, header = T)[, c("dir", "sample")])) %>% 
    Reduce(rbind, .)
  if (!is.null(samples.exclude))
    df <- df %>% filter(! sample %in% samples.exclude)
  
  qc.df <- df %>% split(., f = factor(.$sample, levels = .$sample)) %>% 
    lapply(function(x) {
      x$dir <- sub("/outs.+$", "/outs/", x$dir)
      qc <- read.csv(x$dir %>% paste0("/summary.csv"), header = T)
      qc <- utilsFanc::add.column.fanc(df1 = qc, df2 = data.frame(sample = x$sample), pos = 1)
      return(qc)
    }) %>% Reduce(rbind, .)
  
  ATAC.df <- qc.df %>% .[, c("sample",colnames(.)[grepl("^ATAC" ,colnames(.))])]
  GEX.df <- qc.df %>% .[, c("sample",colnames(.)[grepl("^GEX" ,colnames(.))])]
  
  dfs <- list(all = qc.df, atac = ATAC.df, gex = GEX.df)
  lapply(names(dfs), function(x) {
    df <- dfs[[x]]
    file.out <- paste0(out.root.name, "..", x, ".xlsx")
    dir.create(dirname(file.out), showWarnings = F, recursive = T)
    xlsx::write.xlsx(df, file.out, col.names = T, row.names = F)
    return(NULL)
  })
 return(dfs)
}