# lib_dir <- c("/home/local/kyeonghunjeong_920205/R_lib")
# if(!dir.exists(lib_dir)) dir.create(lib_dir)
# .libPaths(lib_dir)
# .libPaths()

library(data.table)
ucsc <- fread("./data/rmsk/hg38.ucsc.rmsk.txt")

ucsc <- ucsc[grepl("LINE",ucsc$repClass) & grepl("L1HS",ucsc$repName),]
ucsc <- ucsc[grepl("L1HS",ucsc$repName),]
ucsc <- ucsc[!grepl("_",ucsc$genoName),]

genescore_input <- data.table(ucsc$genoName,
                       ucsc$genoStart,
                       ucsc$genoEnd,
                       paste0(ucsc$genoName,"|",
                              ucsc$genoStart,"|",
                              ucsc$genoEnd,"|",
                              ucsc$repName,":",
                              ucsc$repFamily,":",
                              ucsc$repClass,"|",
                              ucsc$strand),
                       ucsc$swScore, ucsc$strand)
fwrite(genescore_input,"./data/rmsk/L1HS_regions.txt")


