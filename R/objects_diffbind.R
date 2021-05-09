quick.consensus.with.diffbind <- function(peak.files, score.col) {
  sample.df <- data.frame(SampleID = paste0("sample", 1:length(peak.files)), 
                          Peaks = peak.files)
  dba <- DiffBind::dba(minOverlap = 1, sampleSheet = sample.df, scoreCol = score.col, bRemoveM = T, 
               bRemoveRandom = T)
  conss <-  DiffBind::dba.peakset(DBA = dba, consensus = T, 
                                  bRemoveM = T, bRemoveRandom = T, minOverlap = 1,
                                  bRetrieve = T)
  mcols(conss) <- NULL
  return(conss)
}