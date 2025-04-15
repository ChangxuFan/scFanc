# deeptools.refpoint <- function(bw.vec, regions.vec, 
#                                upstream = 2000, downstream = 2000, color.map = "Reds",
#                                sort.using.samples = 1,
#                                compute.matrix = T, plot.heatmap = T,
#                                threads,
#                                out.dir, root.name = NULL,
#                                computeMatrix = "/bar/cfan/anaconda2/envs/jupyter/bin/computeMatrix",
#                                plotHeatmap = "/bar/cfan/anaconda2/envs/jupyter/bin/plotHeatmap") {
#   if (!is.null(names(bw.vec))) {
#     sample.labels <- names(bw.vec)
#   } else {
#     sample.labels <- basename(bw.vec) %>% sub("\\.bw$|\\.bigwig$|\\.bigWig$|\\.BigWig$", "", .)
#   }
#   if (length(color.map) != 1) {
#     if (length(color.map) != length(bw.vec)) {
#       stop("length(color.map) != length(bw.vec)")
#     }
#   }
#   if (is.null(root.name)) {
#     root.name <- basename(out.dir)
#   }
#   dir.create(out.dir, showWarnings = F, recursive = T)
#   prefix <- paste0(out.dir, "/", root.name)
#   mat <- paste0(prefix, ".mat.gz")
#   mat_hm <- paste0(prefix, ".tsv.gz")
#   sorted_bed <- paste0(prefix, ".sorted.bed")
#   heatmap <- paste0(prefix, ".hm.pdf")
#   
#   if (compute.matrix == T) {
#     cmd <- paste0(computeMatrix, " reference-point -S ", paste0(bw.vec, collapse = " "),
#                   " -R ", paste0(regions.vec, collapse = " "),
#                   " -o ", mat, " --missingDataAsZero ", 
#                   " -b ", downstream, " -a ", upstream, " --referencePoint center ", 
#                   " --samplesLabel ", paste0(sample.labels, collapse = " "),
#                   " -p ", threads)
#     print(cmd); system(cmd)
#   }
#   
#   if (plot.heatmap == T) {
#     cmd <- paste0(plotHeatmap, " -m ", mat, " -out ", heatmap, 
#                   " --outFileNameMatrix ", mat_hm, " --outFileSortedRegions ", sorted_bed, 
#                   " --sortUsingSamples ", sort.using.samples, " --colorMap ", color.map, 
#                   " --refPointLabel center")
#     print(cmd); system(cmd)
#   }
#   return()
#   
# }