runTranslateGenes <- function () {
    idConversionTable <- read.table ( "/Users/hrdueck/Documents/Murray/Projects/20160419sampleSelection/references/mart_export.txt" , header = T , as.is = T , fill = T , sep = '\t' , row.names = 1 )
    counts <- read.table ( "/Users/hrdueck/Documents/Murray/Projects/20160422combineData/bulkRNASeq/wm9LdB_star_HTSeqCounts_full-2.tsv" , header = T , sep = '\t' , as.is = T , fill = T )
    
    # Select samples
    desSamples <- c ( "A6_bulk_NoDrug" , "C10_bulk_NoDrug" )
    keep <- counts$sampleID %in% desSamples
    counts <- counts [ keep , ]
    
    # transform table shape
    library ( reshape2 )
    countsMatrix <- acast ( counts , gene_id ~ sampleID , value.var = "counts" )
    
    # translate genes
    countsMatrix <- translateGenes ( idConversionTable , countsMatrix )
    
    # Save:
    write.table ( countsMatrix , "data/BulkRawCounts.txt" )
    
}

runNormalizeCounts <- function () {
    counts <- read.table (  "data/BulkRawCounts.txt" , as.is = T , header = T )
    
    # RPM
    rpmNorm <- apply ( counts , 2 , function ( x ) { 1000000 * x / sum ( x ) } )
    
    # GAPDH
    gapdhNorm <- apply ( counts , 2 , function ( x ) { x / x [ "GAPDH" ] } )
    
    # Save
    write.table ( rpmNorm , "data/BulkRPM.txt" )
    write.table ( gapdhNorm  , "data/BulkGAPDH.txt" )
}

runExtractFISHGenes <- function () {
    FISH <- read.table ( "~/Dropbox/paper/extractedData/data/fishSubset.txt" , header = T , as.is = T )

    rpm <- read.table ( "data/BulkRPM.txt" , header = T , as.is = T )
    gapdh <- read.table ( "data/BulkGAPDH.txt" , header = T , as.is = T )
    
    sharedGenes <- intersect ( colnames ( FISH ) , rownames ( rpm ) )
    
    write.table ( t ( rpm [ sharedGenes , ] ) , "data/BulkRPMSubset.txt" )
    write.table ( t ( gapdh [ sharedGenes , ] ) , "data/BulkGAPDHSubset.txt" )
 
 
}
### From here: /Users/hrdueck/Documents/Murray/Projects/20160419sampleSelection/getCellStats.R
# Modified to remove list of desired genes.
translateGenes <- function ( idConversionTable , counts ) {
    # Tested.
    # Start with ids present in counts
    idConversionTable <- subset ( idConversionTable, subset = rownames ( idConversionTable ) %in% rownames ( counts ) )
    # 56885 of these.
    
    # Remove duplicated symbols
    dups <- unique ( idConversionTable$Associated.Gene.Name [ duplicated ( idConversionTable$Associated.Gene.Name ) ] )
    # 668 of these
    keep <- ! ( idConversionTable$Associated.Gene.Name %in% dups )
    idConversionTable <- subset ( idConversionTable , subset = keep )
    # 52761 remaining
    
    # Then convert names
    counts <- counts [ rownames ( idConversionTable ) , ]
    rownames ( counts ) <- idConversionTable$Associated.Gene.Name
    
    # Return
    counts
}
