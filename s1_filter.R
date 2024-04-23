library(tidyverse)
library(patchwork)
library(Seurat)
library(DoubletFinder)


###############################
## Filter low quality cells ##
###############################
outdir = prepDir("../results/s1_qc_harmony")
min.features = 200 # filter cell by feature
min.cells = 3 # filter feature by cell

files = list.files(path = "../data", full.names = T)
samples = gsub("../data/", "", files[1:6])
samples = gsub("_.*", "", samples)

slist = lapply(seq_along(samples), function(x) {
  id = samples[x]
  # = grep(id,files,value = T)
  hd5.f = files[x]
  print(sprintf("Reading %s", hd5.f))
  dt = Read10X_h5(hd5.f)
  s = CreateSeuratObject(
    counts = dt,
    project = id,
    min.features = min.features,
    min.cells = min.cells
  )
})

combined <- merge(slist[[1]], y = slist[-1])

combined$sample = factor(combined@active.ident, levels = samples)
Idents(combined) = "sample"

combined <-
  PercentageFeatureSet(object = combined,
                       pattern = "^mt-",
                       col.name = "percent.mt")
VlnPlot(
  combined,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.0,
  group.by = "sample",
  log = F
)

th.nFeature_RNA = 500
th.mt = 15

combined <-
  subset(combined, subset = nFeature_RNA > th.nFeature_RNA &
           percent.mt < th.mt)

seurat_list = SplitObject(combined, split.by = "sample")
table(combined$sample)
mr <- c(0.046, 0.057, 0.087, 0.061, 0.071, 0.076)


PCs = c(15, 17, 19, 13, 18, 18)
for (i in 1:6) {
  rm(seu)
  seu <- seurat_list[[i]]
  seu <- NormalizeData(seu)
  seu <-
    FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = PCs[i])
  
  seu <- RunUMAP(seu, dims = 1:PCs[i])
  seu <- FindNeighbors(seu, dims = 1:PCs[i])
  seu <- FindClusters(seu)
  
  sweep.res.list <- paramSweep_v3(seu, PCs = 1:PCs[i], sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  
  pK = bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
  pK = as.numeric(as.character(pK[[1]]))
  
  annotations <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(mr[i] * nrow(seu@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  seu <-
    doubletFinder_v3(
      seu,
      PCs = 1:PCs[i],
      pN = 0.25,
      pK = pK,
      nExp = nExp_poi,
      reuse.pANN = FALSE,
      sct = FALSE
    )
  
  write.csv(seu@meta.data,
            sprintf("%s/doublet.%s.predict.csv", outdir, names(seurat_list)[i]))
}


doublets = list.files(outdir, pattern =  "*.predict.csv", full.names = T)
doublets <- lapply(doublets, function(x)  {
  x = read.csv(x)
  x = x %>% mutate(orig.ident = X) %>% select(c("orig.ident", "sample", ncol(x)))
  colnames(x) = c("orig.ident", "sample", "doublets")
  return(x)
})
doublets = do.call("rbind", doublets)

darray = doublets$doublets
names(darray) = doublets$orig.ident
combined = AddMetaData(combined, metadata = darray, col.name = "doublets")
DefaultAssay(combined) = "RNA"

th.nFeature_RNA.doublet = 9000
th.nCount_RNA.doublet = 50000
seurat.sub = subset(combined, subset = nFeature_RNA <= th.nFeature_RNA.doublet)
seurat.sub = subset(seurat.sub, subset = nCount_RNA <= th.nCount_RNA.doublet)

seurat.sub = subset(seurat.sub, subset = doublets != "Doublet")
seurat.sub = DietSeurat(seurat.sub, scale.data = F)

saveRDS(seurat.sub, file.path("../rds", "seurat.s1.rds"))
