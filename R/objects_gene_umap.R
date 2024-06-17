gene.cell.swap <- function(so, already.swapped = F,
                           pc.run = 50, pc.dim = 1:5,
                      cluster.resolution = 0.6,
                      umap.seed = 0, cluster.seed = 0,
                      work.dir, plot.dir = NULL, save.Rds = T) {
  if (is.null(plot.dir)) {
    plot.dir <- paste0(work.dir, "/plots/")
  }
  dir.create(plot.dir, showWarnings = F, recursive = T)
  dir.create(work.dir, showWarnings = F, recursive = T)
  
  if (!already.swapped) {
    mat <- so@assays$RNA@counts
    mat <- Matrix::t(mat)
    sog <- Seurat::CreateSeuratObject(counts = mat, project = "gene_cell_flip",
                                      min.cells = 0, min.features = 0)
  } else{
    sog <- so
  }
  
  assay <- "RNA"
  sog <- NormalizeData(object = sog, normalization.method = "LogNormalize", assay = assay)
  sog <- FindVariableFeatures(sog, assay = assay)
  sog <- ScaleData(object = sog, assay = assay, vars.to.regress = NULL)
  
  sog <- RunPCA(sog, npcs = pc.run, assay = assay)
  
  p <- ElbowPlot(sog)
  wrap.plots.fanc(plot.list = list(p), sub.width = 6, sub.height = 6, 
                  plot.out = paste0(plot.dir, "/elbow.png"))
  
  sog <- RunUMAP(sog, reduction = "pca", dims = pc.dim, assay = assay, seed.use = umap.seed)
  
  sog <- FindNeighbors(sog, reduction = "pca", dims = pc.dim)
  sog <- FindClusters(sog, resolution = cluster.resolution, random.seed = cluster.seed)
  
  p <- DimPlot(sog, group.by = "seurat_clusters", label = T)
  wrap.plots.fanc(plot.list = list(p), sub.width = 6, sub.height = 6, 
                  plot.out = paste0(plot.dir, "/umap.png"))
  
  sog$nCell_RNA <- sog$nFeature_RNA
  plot.panel.list(c("nCell_RNA"), sog, assay = "RNA", order = F, 
                  root.name = paste0(plot.dir, "/depth"), max.quantile = 0.99)
  
  plot.panel.list(c("nCell_RNA"), sog, assay = "RNA", order = F, 
                  violin = T,
                  root.name = paste0(plot.dir, "/depth_violin"), max.quantile = 0.99)
  
  if (save.Rds) {
    saveRDS(sog, paste0(work.dir, "/sog.Rds"))
  }
  return(sog)
}