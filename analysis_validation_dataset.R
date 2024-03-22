library(trqwe)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(data.table)
library(devtools)
library(ArchR)

options(future.globals.maxSize = 1e10)
options(Seurat.object.assay.version = "v5")
t2.data <- Read10X("dataset_validation/T2_RNA/outs/filtered_feature_bc_matrix/")
t2 <- CreateSeuratObject(counts = t2.data, project = "T2_ccRCC", min.cells = 3, min.features = 200)
t2[["percent.mt"]] <- PercentageFeatureSet(t2, pattern = "^MT-")

t3.data <- Read10X("dataset_validation/T3_RNA/outs/filtered_feature_bc_matrix/")
t3 <- CreateSeuratObject(counts = t2.data, project = "T3_ccRCC", min.cells = 3, min.features = 200)
t3[["percent.mt"]] <- PercentageFeatureSet(t3, pattern = "^MT-")

t4.data <- Read10X("dataset_validation/T4_RNA/outs/filtered_feature_bc_matrix/")
t4 <- CreateSeuratObject(counts = t4.data, project = "T4_ccRCC", min.cells = 3, min.features = 200)
t4[["percent.mt"]] <- PercentageFeatureSet(t4, pattern = "^MT-")

merged_obj <- merge(t2, list(t3,t4))
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
merged_obj <- FindNeighbors(merged_obj, dims = 1:30, reduction = "pca")
merged_obj <- FindClusters(merged_obj, resolution = 2, cluster.name = "unintegrated_clusters")
merged_obj <- RunUMAP(merged_obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(merged_obj, reduction = "umap.unintegrated", group.by ="orig.ident")
merged_obj <- IntegrateLayers(
  object = merged_obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
merged_obj <- IntegrateLayers(
  object = merged_obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)
merged_obj <- IntegrateLayers(
  object = merged_obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

mcsaveRDS(merged_obj,"dataset_validation/tmp_mergedT234.rds")
merged_obj <- RunUMAP(merged_obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
merged_obj <- RunUMAP(merged_obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
p_integrated <- DimPlot(merged_obj,reduction = "umap.cca" , group.by="orig.ident") + 
  DimPlot(merged_obj,reduction = "umap.rpca" , group.by="orig.ident") +
  DimPlot(merged_obj,reduction = "umap.harmony" , group.by="orig.ident")

merged_obj <- FindNeighbors(merged_obj, reduction = "harmony", dims = 1:30)

merged_obj$sample <- merged_obj$orig.ident
merged_obj <- JoinLayers(merged_obj)




cMarkers <- fread("dataset_validation/ccRCC_canonicalMarkers.csv")
cMarkers$CellType <- gsub("\\+","",cMarkers$CellType)
cMarkers$CellType <- gsub(" ","_",cMarkers$CellType)
celltypes <- unique(cMarkers$CellType)
celltypes
cMarkers[CellType== celltypes[1], Gene]
cMakrers_list <- lapply(split(cMarkers, cMarkers[[2]]), function(x) x$Gene)
cMakrers_list <- cMakrers_list[match(celltypes,names(cMakrers_list))]


merged_obj<- AddModuleScore(merged_obj,features = cMakrers_list, name= names(cMakrers_list))

j <- 1
t4 <- t2 
celltypes_meta <- celltypes[-1]

for(j in 1:length(celltypes)){
  tmp_markers <- cMarkers[CellType== celltypes_meta[j], Gene]
  
  
  merged_obj <- MetaFeature(merged_obj,features =c(tmp_markers),meta.name =paste0("score_",celltypes_meta[j]))
  
}
p_fscore_h <- FeaturePlot(merged_obj, features = paste0("score_",celltypes_meta),reduction = "umap.harmony")
p_fscore_cca <- FeaturePlot(merged_obj, features = paste0("score_",celltypes_meta),reduction = "umap.cca")
p_fscore_rpca <- FeaturePlot(merged_obj, features = paste0("score_",celltypes_meta),reduction = "umap.rpca")
dir.create("result_validation")
ggsave(filename = "result_validation/umap_Merged_RNA.png",
       plot=p_integrated, dpi=300,width=8, height = 8)
ggsave(filename = "result_validation/featureScore_Merged_RNA_harmony.png",
       plot=p_fscore_h, dpi=300,width=12, height = 12)
ggsave(filename = "result_validation/featureScore_Merged_RNA_cca.png",
       plot=p_fscore_cca, dpi=300,width=12, height = 12)
ggsave(filename = "result_validation/featureScore_Merged_RNA_rpca.png",
       plot=p_fscore_rpca, dpi=300,width=12, height = 12)
merged_obj <- FindClusters(merged_obj, resolution = 0.2, cluster.name = "harmony_clusters")
pharmony_cluster <- DimPlot(merged_obj,reduction = "umap.harmony",group.by="harmony_clusters",label=T)
pharmony_sample <- DimPlot(merged_obj,reduction = "umap.harmony",group.by="sample",label=F)

ggsave(filename = "result_validation/umap_Merged_RNA_Harmony.png",
       plot=pharmony_cluster, dpi=300,width=8, height = 8)
ggsave(filename = "result_validation/umap_Merged_RNA_Harmony_Sample.png",
       plot=pharmony_sample, dpi=300,width=8, height = 8)

cl_to_type <- data.table(c(levels(merged_obj$harmony_clusters)),
                         celltype=c("T_cell", # 0
                                    "NK_cell", #1
                                    "Macrophage", # 2
                                    "Mesangial_cell", # 3
                                    "Endothelium", # 4
                                    "ccRCC", # 5
                                    "Monocyte", # 6
                                    "CD8_T_cell", #7
                                    "Endothelium", #8
                                    "Neutrophil", # 9
                                    "B_cell",#10
                                    "B_cell",#11
                                    "Mesangial_cell",#12
                                    "Mast_cell", #13
                                    "Proliferating_CD8_T_cell", #14
                                    "Macrophage", #15
                                    "Macrophage", #16
                                    "Dendritic_cell"# 17
                         ))


new.cluster.ids <- cl_to_type$celltype
names(new.cluster.ids) <- levels(merged_obj$harmony_clusters)

Idents(merged_obj) <- "harmony_clusters"
merged_obj <- Seurat::RenameIdents(merged_obj,new.cluster.ids)




p_umap2 <- DimPlot(merged_obj, reduction = "umap.harmony",label = T,label.size = 4.5) + NoLegend()
p_umap2
p_celltype_score <- VlnPlot(merged_obj,features = paste0("score_",celltypes_meta))
ggsave(filename = "result_validation/Vln_celltype_merged_RNA.png",
       plot=p_celltype_score, dpi=500,width=16, height = 16)
ggsave(filename = "result_validation/celltype_merged_RNA.png",
       plot=p_umap2, dpi=300,width=8, height = 8)

p_lnc <- FeaturePlot(t4,features = rownames(head(deg_lnc,9)))
ggsave(filename = "analysis_lncRNA/DEG_lnc_T2_RNA.png",
       plot = p_lnc, dpi=300,width=12, height = 12)

p_all <- FeaturePlot(t4,features = rownames(head(deg_allgene,9)))
ggsave(filename = "analysis_lncRNA/DEG_all_T2_RNA.png",
       plot=p_all,dpi=300,width=12, height = 12)

trqwe::mcsaveRDS(object = merged_obj,file = "result_validation/RNA_merged_obj.rds",mc.cores = 8)
#######################################################################################################################################

setwd("result_validation/")

# ArchR::installExtraPackages()
addArchRThreads(threads = 78) 
addArchRGenome("hg38")
inputFiles <- c("../dataset_validation/T2_ATAC/outs/fragments.tsv.gz","../dataset_validation/T3_ATAC/outs/fragments.tsv.gz",
                "../dataset_validation/T4_ATAC/outs/fragments.tsv.gz")
names(inputFiles) <- c("T2_ATAC","T3_ATAC","T4_ATAC")


ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ArchR_T234_raw",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

proj2 <- filterDoublets(proj)
proj2 <- addIterativeLSI(
  ArchRProj = proj2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

proj2 <- addHarmony(
  ArchRProj = proj2,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

proj2 <- addClusters(
  input = proj2,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "H_Clusters",
  resolution = 0.8
)

proj2 <- addUMAP(
  ArchRProj = proj2, 
  reducedDims = "Harmony", 
  name = "H_UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
proj2 <- addUMAP(
  ArchRProj = proj2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
proj2 <- addClusters(
  input = proj2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)
p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "H_UMAP")
p4 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "H_Clusters", embedding = "H_UMAP")


ggAlignPlots(p1, p2,p3,p4, type = "hv")

genes <- getFeatures(proj2)
genes[grepl("^L1",genes)]
(projHeme2)

proj2 <- saveArchRProject(proj2,outputDirectory = "ArchR_T234")
proj2 <- loadArchRProject(path = "ArchR_T234")

proj_subset <- subsetArchRProject(proj2,outputDirectory = "ArchRSubset_T234",cells = proj2[proj2$predictedScore > 0.5,]$cellNames)
#######################################################################################################################################

seRNA <- mcreadRDS("result_validation/RNA_merged_obj_v3.rds")
proj <- loadArchRProject("result_validation/ArchR_T234")
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  force= TRUE,
  groupRNA = "celltype",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)
saveArchRProject(proj)
getAvailableMatrices(proj)

#######################################################################################################################################
setwd("result_validation/")

addArchRThreads(threads = 78) 
addArchRGenome("hg38")
proj_subset <- loadArchRProject( "ArchRSubset_T234")
proj_subset$predictedGroup

l1hs <- fread("../data/rmsk/L1HS_regions.txt")
tmp_regions <- GRanges(seqname= l1hs$V1, ranges = IRanges(start=l1hs$V2, end=l1hs$V3),strand = l1hs$V6)
tmp_regions
tmp_regions$gene_id <- gsub(":", "_",l1hs$V4,fixed=T)
tmp_regions$symbol <- gsub(":", "_", l1hs$V4,fixed=T)

genes_te <- do.call(c, as(list(tmp_regions),"GRangesList"))

proj_subset <- addGeneScoreMatrix(proj_subset,genes=genes_te, matrixName="TEScoreMatrix",force=T)
proj_subset <- addImputeWeights(proj_subset)


markersTE <- getMarkerFeatures(
  ArchRProj = proj_subset, 
  useMatrix = "TEScoreMatrix", 
  groupBy = "predictedGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersGS <- getMarkerFeatures(
  ArchRProj = proj_subset, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "predictedGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

TEList <- getMarkers(markersTE, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
TEList$ccRCC
TEList
markers1 <- TEList$ccRCC[abs(TEList$ccRCC$end - TEList$ccRCC$start) > 5000,"name" ]

plotEmbedding(proj_subset,embedding="H_UMAP",name = "predictedGroup")
plotEmbedding(proj_subset,embedding="UMAP",name = "predictedGroup")
plotEmbedding(proj_subset,embedding="UMAP",name = markers1, colorBy="TEScoreMatrix")

mcsaveRDS(markersTE,"./TEMarkers.rds")
mcsaveRDS(markersGS,"./GeneMarkers.rds")

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_subset, 
  useMatrix = "PeakMatrix", 
  groupBy = "predictedGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")



