library(Seurat)
library(ArchR)
library(regioneR)
library(trqwe)
options(future.globals.maxSize= 300*1000*1024^2)
options(Seurat.object.assay.version = 'v4')
addArchRThreads(threads = 90) 
addArchRGenome("hg38")
seurat_atac = mcreadRDS("..") # Seurat Object path
proj = loadArchRProject("..") # ArchR Project path
l1hs <- fread("./data/rmsk/L1HS_regions.txt")
l1hs_regions <- GRanges(seqname= l1hs$V1, ranges = IRanges(start=l1hs$V2, end=l1hs$V3),strand = l1hs$V6)
l1hs_regions$gene_id <- gsub(":", "_",l1hs$V4,fixed=T)
l1hs_regions$symbol <- gsub(":", "_", l1hs$V4,fixed=T)

proj_scored <- addGeneScoreMatrix(proj,genes= l1hs_regions, matrixName="L1HSScoreMatrix",force=T)
proj_scored <- addImputeWeights(proj_scored)
unique(proj_scored$predictedGroup)
markersL1HS <- getMarkerFeatures(
  ArchRProj = proj_scored, 
  useMatrix = "L1HSScoreMatrix", 
  useGroups = "ccRCC",
  groupBy = "predictedGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

list_L1HS <- getMarkers(markersL1HS, cutOff = "FDR <= 0.01 & Log2FC >= 0.55")


L1HS_full_length <- list_L1HS$ccRCC[abs(list_L1HS$ccRCC$end - list_L1HS$ccRCC$start) > 5000,"name" ]
L1HS_full_length[c(1,2,8)]

p_tumor_te_h<- plotEmbedding(proj_scored,
                             colorBy = "L1HSScoreMatrix", 
                             name=L1HS_full_length[c(1,2,8)],
                             embedding="H_UMAP",
                             quantCut = c(0.1,0.95),
                             imputeWeights = getImputeWeights(proj_scored)
)


l1hs_score_se <- ArchR::getMatrixFromProject(proj_scored,"L1HSScoreMatrix")
assays(l1hs_score_se)$L1HSScoreMatrix

getTEScoreMatrix <- function (ArchRProject, SeuratObject)
{
  print("In Progress:")
  print("Get TE Score Matrix From ArchRProject")
  TEScore_matrix <- ArchR::getMatrixFromProject(ArchRProject, 
                                                useMatrix = "TEScoreMatrix")
  gsm <- assays(TEScore_matrix)$TEScoreMatrix
  print("get TE Features From ArchRProject")
  GeneFeatures <- getFeatures(ArchRProj = ArchRProject, useMatrix = "TEScoreMatrix", 
                              select = NULL, ignoreCase = TRUE)
  colnames(gsm) <- gsub("#", "_", colnames(gsm))
  ix <- match(colnames(SeuratObject), colnames(gsm))
  gsm <- gsm[, ix]
  print("Saving TE Features From ArchRProject into TE Score Matrix")
  rownames(gsm) <- GeneFeatures
  print("Return TE Score Matrix")
  return(gsm)
}

tsm <-getTEScoreMatrix(proj_subset,seurat_atac)
gsm <- getGeneScoreMatrix(ArchRProject = proj_subset, SeuratObject = seurat_atac)



seurat_atac[['TEScore']] <- CreateAssayObject(counts = tsm)
seurat_atac[['GeneScore']] <- CreateAssayObject(counts = gsm)

# getAvailableMatrices(proj)
