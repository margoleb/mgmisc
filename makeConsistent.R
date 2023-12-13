#' @export makeConsistent

makeConsistent <- function(
    asvtab = NULL,
    sampledata = NULL,
    taxonomy = NULL,
    tree = NULL,
    verbose = TRUE
){
  out <- list()
  if(is.null(asvtab)){ message('ASV (OTU) table must be present'); stop()}

  nsamples_orig <- nrow(asvtab)
  nasv_orig <- ncol(asvtab)

  if( is.null(sampledata) & is.null(taxonomy) & is.null(tree) ){message('No data sampledata, taxonomy or tree were given, nothing to do'); stop()}

  if( !is.null(sampledata)){
    asvtab <- asvtab[ rownames(asvtab) %in% rownames(sampledata), ] # only those samples, for which there are sampledata
    sampledata <- sampledata[ rownames(sampledata) %in% rownames(asvtab), ] # only those samples, for which there are ASVs
    asvtab <- asvtab[ order(rownames(asvtab)), ]
    sampledata <- sampledata[ order(rownames(sampledata)), ]
    if( !identical(rownames(asvtab), rownames(sampledata))){
      stop('Something is wrong, there was no possibility of making ASV table and sampledata rownames identical')
    }
    if( nrow(asvtab) == 0 ){
      stop('No common row names in ASV table and sampledata, check objects (and their row names)')
    }
   }
    out$sampledata <- sampledata
    nsamples <- nrow(asvtab)

    if( !is.null(taxonomy)){
      asvtab <- asvtab[ , colnames(asvtab) %in% rownames(taxonomy)] # only those ASVs for which there is taxonomy
      taxonomy <- taxonomy[ rownames(taxonomy) %in% colnames(asvtab), ]
      asvtab <- asvtab[ , order(colnames(asvtab))]
      taxonomy <- taxonomy[ order(rownames(taxonomy)), ]
      if(!identical(colnames(asvtab), rownames(taxonomy))){
        stop('Something went wrong, sets ASV names and taxa names differ')
      }
      if(ncol(asvtab) == 0){
        stop('No common names in ASVs table and tax table, check objects (OTUs vs ASVs?)')
      }
      nasv <- ncol(asvtab)
    }
    out$taxonomy <- taxonomy

    if(!is.null(tree)){
      tree <- phangorn::midpoint(tree) # midpoint rooting
      tree$edge.length[which(is.na(tree$edge.length))] <- 0 # some edge lenghts can be NA, they are coverted to 0
      asvtab <- asvtab[ , colnames(asvtab) %in% tree$tip.label ] # only ASVs that are in the tree
      tree <- ape::drop.tip(tree, tree$tip.label[!(tree$tip.label %in% colnames(asvtab))]) # only ASVs that are in the asvtab (it can be rarefied)
      if(ncol(asvtab) == 0){
        stop('No common names in ASV table and tree, check objects (ASVs vs OTUs?)')
      }
      nasv <- ncol(asvtab)
    }
    out$tree <- tree

    if(verbose){
      message(paste0("Out of ", nsamples_orig, " samples ", nsamples, " were retained.\nOut of ", nasv_orig, " ASVs (OTUs) ", nasv, " were retained."))
    }

    return(out)
}
