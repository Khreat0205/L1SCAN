library(trqwe)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(ArchR)
library(regioneR)
library(ArchRtoSignac)
####scATAC-seq
args<-commandArgs(trailingOnly=TRUE)
sample_id <- args[1]
sprintf("Start to Create Signac object - %s ",sample_id)
### make cmd
# samples <- list.files("../dataset_discovery/",pattern="^RCC[0-9]")
# sprintf("Rscript 1.OP_eachObject.R %s &",samples)

# sample_id <- "RCC81"
set.seed(1234)
plan("multisession", workers = 5)
options(future.globals.maxSize = 50000 * 1024^2)

combined.peaks <- mcreadRDS("combined_peaks.rds")
md.RCC <- read.table(
  file = sprintf("../dataset_discovery/%s/singlecell.csv",sample_id),
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.RCC <- md.RCC[md.RCC$is__cell_barcode > 0.5, ]

frags.RCC <- CreateFragmentObject(
  path = sprintf("../dataset_discovery/%s/fragments.tsv.gz",sample_id),
  cells = rownames(md.RCC)
)

RCC.counts <- FeatureMatrix(
  fragments = frags.RCC,
  features = combined.peaks,
  cells = rownames(md.RCC)
)


metadata_RCC <- read.csv(
  file = sprintf("../dataset_discovery/%s/singlecell.csv",sample_id),
  header = TRUE,
  row.names = 1
)

RCC_assay <- CreateChromatinAssay(RCC.counts, fragments = frags.RCC)
RCC <- CreateSeuratObject(RCC_assay, assay = "ATAC", meta.data = metadata_RCC)
#################
RCC$dataset <- sample_id


annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(RCC) <- annotations

gene.activities <- GeneActivity(RCC)
RCC[['RNA']] <- CreateAssayObject(counts = gene.activities)
RCC <- NormalizeData(
  object = RCC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(RCC$nCount_RNA)
)
mcsaveRDS(RCC,file=sprintf("signac_objects/%s.rds",sample_id))
sprintf("Finished to Create Signac object - %s ",sample_id)
#######################################################################################################################################


set.seed(1234)
plan("multisession", workers = 60)
options(future.globals.maxSize = 50000 * 1024^2)

RCC81 <- mcreadRDS("signac_objects/RCC81.rds")
RCC84 <- mcreadRDS("signac_objects/RCC84.rds")
RCC86 <- mcreadRDS("signac_objects/RCC86.rds")
RCC87 <- mcreadRDS("signac_objects/RCC87.rds") 
RCC94 <- mcreadRDS("signac_objects/RCC94.rds")
RCC96 <- mcreadRDS("signac_objects/RCC96.rds") 
RCC99 <- mcreadRDS("signac_objects/RCC99.rds") 
RCC100 <- mcreadRDS("signac_objects/RCC100.rds") 
RCC101 <- mcreadRDS("signac_objects/RCC101.rds")
RCC103 <- mcreadRDS("signac_objects/RCC103.rds") 
RCC104 <- mcreadRDS("signac_objects/RCC104.rds")
RCC106 <- mcreadRDS("signac_objects/RCC106.rds")
RCC112 <- mcreadRDS("signac_objects/RCC112.rds")
RCC113 <- mcreadRDS("signac_objects/RCC113.rds")
RCC114 <- mcreadRDS("signac_objects/RCC114.rds")
RCC115 <- mcreadRDS("signac_objects/RCC115.rds")
RCC116 <- mcreadRDS("signac_objects/RCC116.rds")
RCC119 <- mcreadRDS("signac_objects/RCC119.rds")
RCC120 <- mcreadRDS("signac_objects/RCC120.rds")

###############################
combined <- merge(
  x = RCC81,
  y = list(RCC84, RCC86,RCC87, RCC94,RCC96, RCC99, RCC100, RCC101,RCC103, RCC104,RCC106,RCC112,RCC113, RCC114,RCC115,RCC116, RCC119,RCC120),
  add.cell.ids = c("RCC81", "RCC84", "RCC86", "RCC87", "RCC94", "RCC96", "RCC99", "RCC100", "RCC101", "RCC103", "RCC104", "RCC106", "RCC112", "RCC113", "RCC114", "RCC115","RCC116", "RCC119", "RCC120")
)
mcsaveRDS(combined,"./raw_merged_ATAC_hg38.rds")
print("Finished to Combine objects")

#######################################################################################################################################
set.seed(1234)
plan("multisession", workers = 60)
options(future.globals.maxSize = 50000 * 1024^2)


combined <- mcreadRDS("./raw_merged_ATAC_hg38.rds")

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(combined) <- annotations

combined <- NucleosomeSignal(object = combined)
combined <- TSSEnrichment(object = combined, fast = T)

combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments
combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')

combined <- subset(
  x = combined,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 20000 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)


combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
combined <- RunSVD(combined)

combined <- RunUMAP(object = combined, reduction = 'lsi', dims = 2:40)
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:40)
combined <- FindClusters(object = combined, resolution = 0.4, verbose = FALSE, algorithm = 3)
###여기 돌리기
mcsaveRDS(combined, file="processed_merged_ATAC_hg38.rds")
print("Finished to Process combined object")

#######################################################################################################################################
### seurat 4.9-

set.seed(1234)
plan("multisession", workers = 60)
options(future.globals.maxSize = 50000 * 1024^2)

combined <- mcreadRDS("processed_merged_ATAC_hg38.rds")
rna <- mcreadRDS("../dataset_discovery/mRNA19_ver2_rpca_annotated.RData")


transfer.anchors <- FindTransferAnchors(reference = rna, query = combined,
                                        features = VariableFeatures(object = rna),
                                        reference.assay = "RNA", query.assay = "RNA",
                                        reduction = "cca")


celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$celltype,
                                     weight.reduction = combined[["lsi"]], dims = 2:30)


combined <- AddMetaData(combined, metadata = celltype.predictions)
subset(combined, prediction.score.max > 0.5)
combined$prediction.score.max
# combined_subset <- combined[,combined$prediction.score.max > 0.5,]

DimPlot(combined,group.by="predicted.id",label=T,raster =F ) + FeaturePlot(combined,"prediction.score.max",label=T,raster =F,min.cutoff = 0.5 )
mcsaveRDS(combined, file="labeled_atac_fromRNA_hg38.rds")
fwrite(celltype.predictions,row.names = T,file = "labeled_atac_fromRNA_hg38.tsv",sep="\t")
# mcsaveRDS(combined_subset, file="subset0.5_labeled_atac_fromRNA.rds")

#######################################################################################################################################


options(future.globals.maxSize= 3000*1024^2)
# options(Seurat.object.assay.version = 'v4')
addArchRThreads(threads = 64) 
addArchRGenome("hg38")

cell_labels <- fread("labeled_atac_fromRNA_hg38.tsv")
cell_subset <- cell_labels[cell_labels$prediction.score.max > 0.5,]
cellnames_subset <- sub("_","_ver2#",cell_subset$V1)


proj <- loadArchRProject(path = "../RCCproj_Ver7")
proj2 <- loadArchRProject(path = "../RCCproj/")


length(setdiff(proj2$cellNames,cellnames_subset))
length(setdiff(cellnames_subset,proj2$cellNames))
length(setdiff(proj$cellNames,cellnames_subset))
length(setdiff(cellnames_subset,proj$cellNames))
intersections <- intersect(cellnames_subset,proj$cellNames)

proj_subset <- proj[proj$cellNames %in% intersections,]
proj_subset
seurat_atac <- mcreadRDS("subset0.5_labeled_atac_fromRNA.rds")
colnames(seurat_atac)<- sub("_","_ver2#",colnames(seurat_atac))
seurat_atac <- seurat_atac[,colnames(seurat_atac) %in% intersections]

library(SeuratObject)
seurat_atac[,seurat_atac$cellbarcode %in% intersections]
seurat_atac$cellbarcode  <- names(seurat_atac$cell_id)
seurat_atac_intersections <- subset(seurat_atac,  cellbarcode %in% intersections)
seurat_atac_intersections



gsm <- getGeneScoreMatrix(ArchRProject = proj_subset, SeuratObject = seurat_atac)
getGeneScoreMatrix
preds_data <- seurat_atac@meta.data[, grepl("predict",colnames(seurat_atac@meta.data))]
preds_data<- preds_data[rownames(preds_data) %in% intersections,]
preds_data<- preds_data[match(proj_subset$cellNames,rownames(preds_data)),]
preds_data
length(proj_subset$cellNames)

sum(proj_subset$cellNames == rownames(preds_data))
proj_subset$celltype.rna <- preds_data$predicted.id
proj_subset@projectSummary


plotEmbedding(proj_subset,name = "celltype.rna")

#######################################################################################################################################

options(future.globals.maxSize= 300*1000*1024^2)
options(Seurat.object.assay.version = 'v4')
addArchRThreads(threads = 90) 
addArchRGenome("hg38")
inputFiles <- list.files("../dataset_discovery/",full.names = T,pattern = "fragments.tsv.gz$",recursive = T)
names(inputFiles) <- basename(dirname(inputFiles))
inputFiles

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE
  # addGeneScoreMat = TRUE
)

getwd()
seurat_atac <- mcreadRDS("not_used/subset0.5_labeled_atac_fromRNA.rds")
length(colnames(seurat_atac))
colnames(seurat_atac)<- sub("_","#",colnames(seurat_atac))



proj<- ArchRProject(ArrowFiles,outputDirectory = "ArchR_RCCs/")
getAvailableMatrices(proj)
intersect(getCellNames(proj),colnames(seurat_atac))
proj_subset <- subsetArchRProject(proj,cells= intersect(getCellNames(proj),colnames(seurat_atac)), outputDirectory = "ArchRSubset")
labeled <- fread("labeled_atac_fromRNA_hg38.tsv")
labeled$V1 <- sub("_","#",labeled$V1)

labeled <- labeled[labeled$V1 %in% getCellNames(proj_subset),]
labeled <- labeled[match(getCellNames(proj_subset),labeled$V1),]

proj_subset$predicted.id <- labeled$predicted.id
proj_subset$prediction.score.max <- labeled$prediction.score.max
proj_subset$prediction.score.ccRCC <- labeled$prediction.score.ccRCC
getAvailableMatrices(proj_subset)

l1hs <- fread("../../SoloTE-main/SoloTE_1.08/solote_hg38.ucsc.L1HS_chr.txt")
tmp_regions <- GRanges(seqname= l1hs$V1, ranges = IRanges(start=l1hs$V2, end=l1hs$V3),strand = l1hs$V6)
tmp_regions
tmp_regions$gene_id <- gsub(":", "_",l1hs$V4,fixed=T)
tmp_regions$symbol <- gsub(":", "_", l1hs$V4,fixed=T)

genes_te <- do.call(c, as(list(tmp_regions),"GRangesList"))

proj_subset <- addGeneScoreMatrix(proj_subset,genes=genes_te, matrixName="TEScoreMatrix",force=T)
# proj2 <- filterDoublets(proj2)
proj_subset <- addIterativeLSI(
  ArchRProj = proj_subset,
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


proj_subset <- addHarmony(
  ArchRProj = proj_subset,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)


# rm(list=ls()[grepl("^K",ls())])
# gc()
# rm(test_obj)
proj_subset <- addClusters(
  input = proj_subset,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)
# proj2 <- addClusters(
#   input = proj2,
#   reducedDims = "IterativeLSI",
#   method = "scran",
#   name = "Clusters",
#   k= 22
# )
proj_subset <- addUMAP(
  ArchRProj = proj_subset, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
proj_subset <- addUMAP(
  ArchRProj = proj_subset, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
markersGS <- getMarkerFeatures(
  ArchRProj = proj_subset, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "predicted.id",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersTE <- getMarkerFeatures(
  ArchRProj = proj_subset, 
  useMatrix = "TEScoreMatrix", 
  groupBy = "predicted.id",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
markerList$ccRCC
TEList <- getMarkers(markersTE, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
TEList$ccRCC

proj_subset <- addImputeWeights(proj_subset)


plotEmbedding(proj_subset,embedding="UMAPHarmony")
plotEmbedding(proj_subset,embedding="UMAPHarmony",name = "predicted.id")
plotEmbedding(proj_subset,embedding="UMAP",name = "predicted.id")

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_subset, 
  useMatrix = "PeakMatrix", 
  groupBy = "predictedGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


TEList$ccRCC
markers1 <- TEList$ccRCC[abs(TEList$ccRCC$end - TEList$ccRCC$start) > 5000,"name" ]
markers1
p_tumor_te_h<- plotEmbedding(proj_subset,
                             colorBy = "TEScoreMatrix", 
                             name=markers1,
                             embedding="UMAP",
                             imputeWeights = getImputeWeights(proj)
)

# devtools::install_github("GreenleafLab/chromVARmotifs")
plotEmbedding(proj_subset,embedding="UMAP",name = "predicted.id")+ p_tumor_te_h[[1]]
plotEmbedding(proj_subset,embedding="UMAP",name = "predicted.id")+plotEmbedding(proj_subset,embedding="UMAP",name = "prediction.score.max")
proj_subset$prediction.score.max

p_tumor_te_h[[2]]
proj_subset <- addGroupCoverages(ArchRProj = proj_subset, groupBy = "predicted.id")
saveArchRProject(proj_subset)
pathToMacs2 <- findMacs2()

proj_subset <- addReproduciblePeakSet(
  ArchRProj = proj_subset, 
  groupBy = "predicted.id",
  pathToMacs2 = pathToMacs2
)

saveArchRProject(proj_subset)
proj_subset <- addPeakMatrix(proj_subset)
saveArchRProject(proj_subset)

tumorPeaks <- getMarkerFeatures(
  ArchRProj = proj_subset, 
  useMatrix = "PeakMatrix", 
  groupBy = "predicted.id",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "ccRCC"
)
getAvailableMatrix(proj_subset)

pma <- markerPlot(seMarker = tumorPeaks, name = "ccRCC", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma
proj_subset <- addMotifAnnotations(ArchRProj = proj_subset, motifSet = "cisbp", name = "Motif")

tumor_motifsUp <- peakAnnoEnrichment(
  seMarker = tumorPeaks,
  ArchRProj = proj_subset,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
tumor_motifsUp
df <- data.frame(TF = rownames(tumor_motifsUp), mlog10Padj = assay(tumor_motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
df
df[grep("ESR",df$TF),]
df[grep("FOXA1",df$TF),]
df[grep("E2F1",df$TF),]
library(ArchRtoSignac)

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

peakAnnoEnrichment(
  seMarker = tumorPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
getPeakSet(proj_subset)
getPeakSet(tumorPeaks)
genes_te@ranges@width > 5000
tmp_regions <- GRanges(seqname= l1hs$V1, ranges = IRanges(start=l1hs$V2 - 5000, end=l1hs$V3 +500),strand = l1hs$V6)
tmp_regions
tmp_regions$gene_id <- gsub(":", "_",l1hs$V4,fixed=T)
tmp_regions$symbol <- gsub(":", "_", l1hs$V4,fixed=T)
genes_te_extend <- do.call(c, as(list(tmp_regions),"GRangesList"))
de_score <- genes_te_extend[genes_te_extend$gene_id %in% markers1,]
de_score

subsetByOverlaps(de_score,getPeakSet(proj_subset))
subsetByOverlaps(getPeakSet(proj_subset),de_score)
getPeakSet(proj_subset)
overlapped <- subsetByOverlaps(getPeakSet(proj_subset),de_score)
saveArchRProject(proj_subset)


