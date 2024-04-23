library(tidyverse)
library(patchwork)
library(Seurat)
library(cowplot)
library(harmony)

outdir = prepDir("../results/manuscript/")


# Figure 1
fig1dir = prepDir(file.path(outdir, "fig1"))

immune$higher_annotation <-
  immune$immune_subtype_annotation
Idents(immune) = "higher_annotation"
IDS =  c(
  "T cells" = c("Tcells"),
  "B cells" = "B",
  "Neutrophils" = "Neu",
  "Dendritic cells" = "DC",
  "Macrohpages" = c("Macro_1C", "Macro_2C")
)
old = as.factor(unname(IDS))
newIDS = names(IDS)
newIDS = gsub("[0-9]+$", "", names(IDS))
names(newIDS) = old
immune <- RenameIdents(object = immune, newIDS)
immune$higher_annotation = Idents(immune)
immune$higher_annotation[colnames(immune) %in% colnames(seurat_new.clean)[seurat_new.clean$higher_annotation ==
                                                                            "pDC"]] = "Dendritic cells"
cols = scales::hue_pal()(6)

fig1b <- DimPlot(
  immune,
  label = F,
  group.by = "higher_annotation",
  reduction = "umap",
  cols = cols[1:5]
)
ggsave(file.path(fig1dir, "fig1_B.pdf"),
       fig1b,
       width = 7,
       height = 5)


#> Export table of all samples broken down by cell type

library(tidyr)
tab.S1 <-
  as.data.frame(table(immune@meta.data$sample, immune@meta.data$higher_annotation))
tab.S1 <- spread(tab.S1, Var2, Freq)
# Add Sum column
tab.S1$"Sum_Immune" <- rowSums(tab.S1[, -1])
# Add a column and annotate samples used in the fractional analysis
samples <- unique(immune$sample)
tab.S1$"Fractional" <- NA
for (i in 1:nrow(tab.S1)) {
  a <- which(samples == as.character(tab.S1$Var1[i]))
  if (length(a) != 0) {
    tab.S1$Fractional[i] <- 1
  }
}
openxlsx::write.xlsx(tab.S1,
                     file = file.path(fig1dir, "Table_of_immune_cell_types_by_sample.xlsx"))


require(ggridges)
require(ggplot2)
# # #
cell.genes <-
  openxlsx::read.xlsx("../refs/1-s2.0-S0092867420308825-mmc2.xlsx", sheet =
                        5)
cell.types  <- as.character(unique(cell.genes$Cell.Type))
coor <- Embeddings(immune, reduction = "umap")
##
data_temp <-
  as.data.frame(as.matrix(GetAssayData(object = immune, slot = "data")))

ggplot.list <- list()
# ggplot.list.2 <- list()
#
# rm(temp)
for (i in 1:length(unique(cell.types))) {
  genes <-
    as.character(cell.genes$Gene[which(cell.genes$Cell.Type == cell.types[i])])
  genes <- stringr::str_to_title(genes)
  gene.exp <-
    colMeans(as.matrix(data_temp[which(rownames(data_temp) %in% genes), ]))[row.names(coor)]
  #clusters <- immune@meta.data$RNA_snn_res.0.7
  # Make ggplot friendly
  #temp <- as.data.frame(cbind(coor, as.data.frame(gene.exp), as.data.frame(clusters)))
  temp <- as.data.frame(cbind(coor, as.data.frame(gene.exp)))
  # Plot with ggplot
  ggplot.list[[i]] <- ggplot(temp, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(colour = gene.exp)) +
    scale_colour_gradient(low = "grey95", high = "red") +
    labs(title = cell.types[i],
         subtitle = paste(genes, collapse = ", "))
  # # Boxplot per cluster
  # ggplot.list.2[[i]] <- ggplot(temp, aes(x = clusters, y = gene.exp)) +
  #                       geom_boxplot() +
  #                       ggtitle(cell.types[i]) + ylab("Average gene expression (log)")
}
# Plot all
require(grid)
require(gridExtra)
#require(gbm)
n <- length(ggplot.list)
nCol <- floor(sqrt(n))
# Expression on tSNE
pdf(
  file.path(
    fig1dir,
    "Immune_cells_UMAP_with_average_expression_of_cell_markers.pdf"
  ),
  15,
  15
)
do.call("grid.arrange", c(ggplot.list, ncol = nCol))
dev.off()
# # Expression per cluster boxplots
# pdf(file.path(dir,"plot_out/IM01/Immune_cells_boxplots_with_average_expression_of_cell_markers.pdf"),15,15)
# do.call("grid.arrange", c(ggplot.list.2, ncol=nCol))
# dev.off()



#### FigureS1

figs1dir = prepDir(file.path(outdir, "figs1"))
fig_fs1 <- DimPlot(immune,
                   label = T,
                   group.by = "RNA_snn_res.0.5",
                   reduction = "umap")

Idents(immune) = "RNA_snn_res.0.5"

markers.immune_0.1 <-
  FindAllMarkers(
    object = immune,
    only.pos = TRUE,
    min.pct = 0.3,
    thresh.use = 1,
    do.print = T,
    max.cells.per.ident = 300
  )
openxlsx::write.xlsx(
  markers.immune_0.1,
  file = file.path(figs1dir, "Immune_cells_marker_gene_1table.xlsx"),
  overwrite = T
)
markers.immune_0.1.list = split(markers.immune_0.1, f = markers.immune_0.1$cluster)

openxlsx::write.xlsx(
  markers.immune_0.1.list,
  file = file.path(figs1dir, "Immune_cells_marker_gene.xlsx"),
  overwrite = T
)


markers.small  <-
  markers.immune_0.1 %>% group_by(cluster) %>% top_n(5, avg_log2FC)
genes_to_check <- unique(markers.small)

immune <-
  ScaleData(immune, features = c(genes_to_check$gene, VariableFeatures(immune)))


pdf(file.path(figs1dir, "UMAP.pdf"), 10, 10)
DimPlot(immune, reduction = "umap", label = F)
DimPlot(immune, reduction = "umap", label = T)
DimPlot(immune,
        reduction = "umap",
        label = F,
        group.by = "sample")
dev.off()

p_heat = DoHeatmap(immune, features = genes_to_check$gene, raster = T) + NoLegend() + theme(text = element_text(size = 20))
ggsave(
  file.path(figs1dir, "heatmap.pdf"),
  p_heat + theme(text = element_text(size = 10)),
  height = 10,
  width = 10
)
tab1 <-
  cbind(
    as.data.frame(immune@meta.data$sample),
    as.data.frame(immune@meta.data$RNA_snn_res.0.5)
  )
colnames(tab1) <- c("Sample", "Immune.cluster")
ggplot(tab1) +
  aes(x = Immune.cluster, fill = factor(Sample)) +
  geom_bar(position = "fill")
dev.off()


Idents(immune) = "higher_annotation"

markers.immune_higher <-
  FindAllMarkers(
    object = immune,
    only.pos = TRUE,
    min.pct = 0.3,
    thresh.use = 1,
    do.print = T,
    max.cells.per.ident = 300
  )
openxlsx::write.xlsx(
  markers.immune_higher,
  file = file.path(
    figs1dir,
    "Immune_cells_marker_gene_1table_higherannotation.xlsx"
  ),
  overwrite = T
)

markers.immune_higher.list = split(markers.immune_higher, f = markers.immune_higher$cluster)

openxlsx::write.xlsx(
  markers.immune_higher.list,
  file = file.path(figs1dir, "Immune_cells_marker_gene_higherannotation.xlsx"),
  overwrite = T
)



# figS3

library(fgsea)
library(msigdbr)


m_df.h <- msigdbr(species = "Mus musculus", category = "H")
m_df.kegg <-
  msigdbr(species = "Mus musculus",
          category = "C2",
          subcategory = "KEGG")
gsea.sets <-
  rbind(m_df.h, m_df.kegg) %>% split(x = .$gene_symbol, f = .$gs_name)



dir.create("../results/manuscript/figS3")
ann = as.character(seurat_new.clean$higher_annotation)

immune_ann = as.character(immune$higher_annotation)
t_ann = seurat_t$t_subtype_annotation

t_immune_idx = match(names(t_ann), names(immune_ann))
t_immune = t_immune[!is.na(t_immune_idx)]
immune_ann[t_immune] = as.character(t_ann[!is.na(t_immune_idx)])

imm_all_idx = match(names(immune$higher_annotation),
                    names(seurat_new.clean$higher_annotation))
ann[imm_all_idx] = immune_ann
seurat_new.clean$tmp = ann
DimPlot(
  seurat_new.clean,
  group.by = "tmp",
  reduction = "harmony_umap",
  label = T,
  label.size = 2
)


Idents(seurat_new.clean) = "tmp"
DefaultAssay(seurat_new.clean) <- "RNA"

seurat_new.clean$bygroup <-
  paste(Idents(seurat_new.clean), seurat_new.clean$Treatment, sep = "_")


clusters = c("B cells", "CD4T", "CD8T", "Neutrophils", "Dendritic cells")

des = list()

group_id = c("KPP", "KP")
for (i in as.character(clusters)) {
  plots = list()
  Idents(seurat_new.clean) <- "bygroup"
  ids = paste(i, group_id, sep = '_')
  nCells = table(seurat_new.clean@meta.data %>% filter(bygroup %in% ids) %>% pull(sample))
  if (any(nCells < 3))
    next
  print(i)
  de = FindMarkers(
    seurat_new.clean,
    ident.1 = ids[1],
    ident.2 = ids[2],
    verbose = FALSE,
    test.use = "wilcox",
    min.pct = .1,
    logfc.threshold = 0,
    max.cells.per.ident = 500
  )
  print("done de")
  de$cluster = i
  de$genes = rownames(de)
  colnames(de)[3:4] = ids
  des[[i]] = de
  
  rank = de %>%
    mutate(order = (avg_log2FC)) %>%
    dplyr::select(genes, order)
  write.table(
    rank,
    file.path(
      "../results/manuscript/figS3/",
      sprintf("log2fc.%s.rnk", paste(ids, collapse = "_vs_"))
    ),
    row.names = F,
    col.names = F,
    sep = '\t',
    quote = F
  )
  rank <- deframe(rank)
  
  fgseaRes_rank <- fgsea(gsea.sets, stats = rank, nperm = 100000)
  fgseaRes_rank <- fgseaRes_rank %>% as_tibble() %>% arrange(padj)
  openxlsx::write.xlsx(fgseaRes_rank,
                       file.path(
                         "../results/manuscript/figS3/",
                         sprintf("fgsea.%s.xlsx", paste(ids, collapse = "_vs_"))
                       ),
                       overwrite = T)
  
  fgseaResTidy <- fgseaRes_rank %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  
  pdf(file = file.path(
    "../results/manuscript/figS3/",
    sprintf("fgsea.barplot.%s.pdf", paste(ids, collapse = "_vs_"))
  ),
  width = 12)
  p = fgseaResTidy %>% filter(padj < 0.1) %>% ggplot(aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = padj)) + scale_fill_gradient(high = "yellow",
                                                     low = "red",
                                                     na.value = NA) +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = "Pathways NES from GSEA, padj < 0.1") +
    theme_minimal()
  print(p)
  dev.off()
  
  
  filtered_pathway <- subset(fgseaResTidy, padj < 0.05)
  
  filt_p <- as.vector(filtered_pathway$pathway)
  
  pdf(
    file.path(
      "../results/manuscript/figS3/",
      sprintf("fgsea.enrichplot.%s.pdf", paste(ids, collapse = "_vs_"))
    ),
    height = 5,
    width = 7.5
  )
  for (pathw in filt_p) {
    plt <- plotEnrichment(
      pathway = gsea.sets[[pathw]],
      gseaParam = 1,
      ticksSize = 0.3,
      stats = rank
    ) +
      labs(title = pathw) + theme(plot.title = element_text(hjust = 0.5, face =
                                                              "bold"))
    print(plt)
    
  }
  dev.off()
  
  
  
  # select only de genes
  sig = de %>% filter(p_val_adj < 0.05) %>% top_n(n = 200, wt = abs(avg_log2FC))
  degenes = de %>% filter(p_val_adj < 0.05) %>% pull(genes)
  
  
  if (nrow(sig) > 1) {
    Idents(seurat_new.clean) = "tmp"
    up_genes_sig = de %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0) %>% arrange((p_val))
    down_genes_sig = de %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC < 0) %>% arrange((p_val))
    t1 = sprintf("%s %s vs %s upreg", i, ids[1], ids[2])
    t2 = sprintf("%s %s vs %s downreg", i, ids[1], ids[2])
    p1 = DoHeatmap(
      seurat_new.clean,
      features = up_genes_sig$genes,
      cells = WhichCells(seurat_new.clean, idents = i),
      group.by = "sample"
    ) + labs(title = t1) + theme(axis.text.y = element_text(size = 4))
    p2 = DoHeatmap(
      seurat_new.clean,
      features = down_genes_sig$genes,
      cells = WhichCells(seurat_new.clean, idents = i),
      group.by = "sample"
    ) + labs(title = t2) + theme(axis.text.y = element_text(size = 4))
    plots[[t1]] = p1 + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
    plots[[t2]] = p2 + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
    pdf(
      file.path(
        "../results/manuscript/figS3/",
        sprintf("DE.heatmap.%s.pdf", paste(ids, collapse = "_vs_"))
      ),
      width = 6,
      height = 15
    )
    print(plots)
    dev.off()
    
    
  }
}


openxlsx::write.xlsx(
  des,
  file.path(
    "../results/manuscript/figS3/",
    "DE_by_Treatment_collapse_replicates.xlsx"
  ),
  row.names = F
)



# macrophage
figdir = prepDir(file.path(outdir, "macrophage"))

Idents(mf.cell.tiss) = "mf_subtype_annotation"
mf.cell.tiss = subset(mf.cell.tiss, idents = "pDC", invert = T)

pdf(
  file.path(figdir, "MF_UMAP_annotation_dropDC.pdf"),
  width = 6,
  height = 5
)
DimPlot(mf.cell.tiss, reduction = "umap")
DimPlot(mf.cell.tiss, reduction = "umap", label = T)
dev.off()

markers = openxlsx::read.xlsx(file.path(figdir, "Cell_Type_Marker_tables.xlsx"))
markers = markers %>% filter(cluster != "pDC")
markers.sel = markers %>%  group_by(cluster) %>% top_n(50, wt = avg_log2FC)
mf.cell.tiss = ScaleData(mf.cell.tiss, features = c(VariableFeatures(mf.cell.tiss), markers.sel$gene))

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(256)

p3 = DoHeatmap(
  mf.cell.tiss,
  features = markers.sel$gene,
  group.by = "Treatment",
  slot = "scale.data"
) + scale_fill_gradientn(colours = rev(mapal))
ggsave(
  filename = file.path(figdir, "Top50_markers_byTreatment_v2.pdf"),
  plot = p3,
  width = 10,
  height = 10
)


FeaturePlot(mf.cell.tiss,
            features = c("Cd86", "H2-Eb1"),
            reduction = "umap") & labs(x = "UMAP1", y = "UMAP2")
ggsave(
  file.path(figdir, "Featureplot_M1_markers.pdf"),
  width = 10,
  height = 5
)

FeaturePlot(mf.cell.tiss,
            features = c("Chil3", "Pparg"),
            reduction = "umap") & labs(x = "UMAP1", y = "UMAP2")
ggsave(
  file.path(figdir, "Featureplot_M2_markers.pdf"),
  width = 10,
  height = 5
)



# T

Idents(t.cell.tiss) = "t_subtype_annotation"

t.cell.tiss.clean <-
  subset(
    t.cell.tiss,
    idents = c("doublet", "Basophil", "Dead", "Prolif_T"),
    invert = T
  )
t.cell.tiss.clean$t_subtype_annotation = droplevels(t.cell.tiss.clean$t_subtype_annotation)

t.cell.tiss.clean = FindVariableFeatures(t.cell.tiss.clean)
t.cell.tiss.clean <-
  RunHarmony(t.cell.tiss.clean,
             "sample",
             plot_convergence = T,
             assay.use = "RNA")

t.cell.tiss.clean <- t.cell.tiss.clean %>%
  RunUMAP(reduction = "harmony",
          dims = 1:20,
          reduction.name = "umap") %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = c(0.3))

DimPlot(
  t.cell.tiss.clean,
  reduction = "umap",
  label = T,
  label.size = 2,
  group.by = "t_subtype_annotation"
)

Idents(t.cell.tiss.clean) = "t_subtype_annotation"
DimPlot(t.cell.tiss.clean, reduction = "umap")
DimPlot(
  t.cell.tiss.clean,
  reduction = "umap",
  label = T,
  label.size = 2
)


library(tidyr)
tab.S1 <-
  as.data.frame(
    table(
      t.cell.tiss.clean@meta.data$sample,
      t.cell.tiss.clean@meta.data$t_subtype_annotation
    )
  )
tab.S1 <- spread(tab.S1, Var2, Freq)
tab.S1$"Sum" <- rowSums(tab.S1[, -1])
samples <- unique(t.cell.tiss.clean$sample)
tab.S1$"Fractional" <- NA
for (i in 1:nrow(tab.S1)) {
  a <- which(samples == as.character(tab.S1$Var1[i]))
  if (length(a) != 0) {
    tab.S1$Fractional[i] <- 1
  }
}
openxlsx::write.xlsx(tab.S1,
                     file = file.path(
                       "../results/manuscript/tcells/",
                       "Table_of_tcells_by_sample.xlsx"
                     ))


markers.annot <-
  FindAllMarkers(
    object = t.cell.tiss.clean,
    only.pos = TRUE,
    min.pct = 0.3,
    thresh.use = 1,
    do.print = T
  )

markers = c("Cd3d", "Cd3e", "Cd8a", "Cd4", "Nkg7", "Ncr1", "Tcrg-C1", "Foxp3")

for (i in markers) {
  p = FeaturePlot(t.cell.tiss.clean, features = i)
  print(p)
}
