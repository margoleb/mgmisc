#' @export extractFeatureSet
#'

extractFeatureSet <- function(
    experiment=NULL,
    feature=NULL,
    rarefied=F,
    d=NULL # distance threshold needed if extracting OTUs
){
  outlist <- list()

  outlist$featuretab <- mgmisc1::featuretab(experiment=experiment, feature=feature, rarefied=rarefied, d=d)
  outlist$sampledata <- mgmisc1::sampledata(experiment=experiment, feature=feature, rarefied=rarefied, d=d)
  outlist$featureannot <- mgmisc1::featureannot(experiment=experiment, feature=feature, rarefied=rarefied, d=d)
  if(!(feature %in% c('original.data', 'KOs', 'pathways'))){
    outlist$tree <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, d=d, what="tree")
    outlist$tree.file <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, d=d, what="tree.file")
    outlist$fasta <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, d=d, what="fasta")
    outlist$fasta.file <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, d=d, what="fasta.file")
    outlist$count_table <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, d=d, what="count_table")
    outlist$count.file <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, d=d, what="count.file")
  }else{
    outlist$tree <- NULL
    outlist$tree.file <-""
    outlist$fasta <- NULL
    outlist$fasta.file <- ""
    outlist$count_table <- NULL
    outlist$count.file <- ""
  }

  return(outlist)
}



