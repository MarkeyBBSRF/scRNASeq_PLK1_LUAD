library(patchwork)
library(Seurat)
library(harmony)

seurat = readRDS("../results/rds/seurat.s2.harmony.ver1.rds")
Idents(seurat) = "RNA_snn_res.0.8"
seurat@meta.data$main_seurat_cluster <-
  seurat@meta.data$RNA_snn_res.0.8
cluster.ids <-
  sort(as.numeric(unique(
    as.character(seurat@meta.data$RNA_snn_res.0.8)
  )))
immune_annotation <- rep("immune", length(cluster.ids))
immune_annotation[10] = "non-immune"

seurat@meta.data$immune_annotation <-
  seurat@meta.data$RNA_snn_res.0.8
seurat@meta.data$immune_annotation <-
  plyr::mapvalues(
    x = seurat@meta.data$immune_annotation,
    from = cluster.ids,
    to = immune_annotation
  )


# All immune cell types

Idents(seurat) = "immune_annotation"
immune = subset(seurat, idents = "immune")
immune = immune %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE)
immune = immune %>% RunPCA(npcs = 40, verbose = FALSE)
immune <-
  RunHarmony(immune,
             "sample",
             plot_convergence = T,
             assay.use = "RNA")

immune <- immune %>%
  RunUMAP(reduction = "harmony",
          dims = 1:15,
          reduction.name = "umap") %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = c(0.5))

markers.immune <- FindAllMarkers(object = immune, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, 
                                 do.print=T, max.cells.per.ident = 500)

Idents(immune) = "RNA_snn_res.0.5
"

IDS =  c(Tcells = c(13,14,17,18,21,25,26),
         B = c(11,27),
         Neu = c(3,10,16,23),
         DC = c(20),
         Macro_1C = c(0,1,2,4,5,6,8,9,12,19,22),
         Macro_2C = c(7,15,24)
)
old = as.numeric(IDS)
newIDS = gsub("[0-9]+$","",names(IDS))
names(newIDS) = old
immune <- RenameIdents(object = immune,newIDS)
immune$immune_subtype_annotation = Idents(immune)

# Macrophages

mf.cell.tiss <- subset(tiss_immune, idents = c("Macro_1C" ,"Macro_2C"))  

mf.cell.tiss = mf.cell.tiss %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE)
mf.cell.tiss = mf.cell.tiss %>% RunPCA(pc.genes = pbmc@var.genes, npcs = 40, verbose = FALSE) 
mf.cell.tiss <- RunHarmony(mf.cell.tiss, "sample",plot_convergence = T,assay.use = "RNA")

mf.cell.tiss <- mf.cell.tiss %>% 
  RunUMAP(reduction = "harmony", dims = 1:15,reduction.name = "umap") %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = c(0.3))
Idents(mf.cell.tiss) = "RNA_snn_res.0.3"
IDS =  c(pDC = c(10,11),
         Macro_1C = c(0,1,2,3,5,6,7),
         Macro_2C = c(4,8,9)
)
old = as.numeric(IDS)
newIDS = gsub("[0-9]+$","",names(IDS))
names(newIDS) = old
mf.cell.tiss <- RenameIdents(object = mf.cell.tiss,newIDS)
mf.cell.tiss$mf_subtype_annotation = Idents(mf.cell.tiss)


# T cells
t.cell.tiss <- subset(tiss_immune, idents = c("Tcells"))  
t.cell.tiss

t.cell.tiss = t.cell.tiss %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE)
t.cell.tiss = t.cell.tiss %>% RunPCA(pc.genes = pbmc@var.genes, npcs = 40, verbose = FALSE) 

t.cell.tiss <- RunHarmony(t.cell.tiss, "sample",plot_convergence = T,assay.use = "RNA")

t.cell.tiss <- t.cell.tiss %>% 
  RunUMAP(reduction = "harmony", dims = 1:15,reduction.name = "umap") %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = c(0.3))
Idents(t.cell.tiss) = "RNA_snn_res.1"
IDS =  c(DP_T = 14,
         CD8T = c(4,6,15),
         CD4T = c(0,2,3,5,9,12,18),
         Neu_T = c(13),
         Macro_T = c(7),
         T_unkown = c(8),
         NK = c(11),
         DN_T = c(1,16),
         Basophil = c(17),
         Empty = c(10)
         
)
old = as.numeric(IDS)
newIDS = gsub("[0-9]+$","",names(IDS))
names(newIDS) = old
t.cell.tiss <- RenameIdents(object = t.cell.tiss,newIDS)
t.cell.tiss$t_subtype_annotation = Idents(t.cell.tiss)

