library(patchwork)
library(Seurat)
library(harmony)

outdir = prepDir("../results/s2_harmony")

seurat = readRDS(file.path("../results/rds/", "seurat.s1.rds"))
Idents(seurat) = "sample"

seurat = seurat %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE)
seurat = seurat %>% RunPCA(pc.genes = pbmc@var.genes,
                           npcs = 40,
                           verbose = FALSE)

whichPC(seurat)

seurat  = seurat %>% RunUMAP(reduction = "pca", dims = 1:20)

seurat <-
  RunHarmony(seurat,
             "sample",
             plot_convergence = T,
             assay.use = "RNA")


DimPlot(
  object = seurat,
  reduction = "harmony",
  pt.size = .1,
  group.by = "sample"
)

seurat <- seurat %>%
  RunUMAP(reduction = "harmony",
          dims = 1:20,
          reduction.name = "harmony_umap") %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = c(0.8))

saveRDS(seurat, "../results/rds/seurat.s2.harmony.ver1.rds")

