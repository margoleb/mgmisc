#' @export sliceData


sliceData <- function(
  asvtab = NULL,
  design = NULL,
  tree = NULL,
  taxa = NULL,
  experiment.object = NULL,
  variable = NULL,
  levels = NULL){

  message("sliceData")
  if(!is.null(asvtab) & !is.null(experiment.object)){
    message("Either ASV table, design, taxonomy and tree or an experiment object may be processed")
    stop()
  }
  if(!is.null(asvtab) & is.null(design)){
    message("When ASV table is processed, design must be given")
  }

}
