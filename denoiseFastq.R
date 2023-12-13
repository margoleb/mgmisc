#' @importFrom dada2 learnErrors filterAndTrim derepFastq dada mergePairs makeSequenceTable removeBimeraDenovo
#' @export denoiseFastq

denoiseFastq <- function(path=NULL, novaseqLoess=NULL, save = TRUE, basename = 'bacteria', pattern="_R1.fastq.gz", truncLen=c(200,200), reads = "both"){
  path <- path

  if(!reads %in% c('both', 'r1', 'r2')){
    warning("Reads must be one of 'both', 'r1' or 'r2', unknown value encountered, assuming paired reads ('both')")
    reads = "both"
  }

    fnFs <- sort(list.files(path, pattern=pattern)) # list of R1 files
    pattern <- gsub("_R1", "_R2", pattern)
    fnRs <- sort(list.files(path, pattern=pattern)) # list of R2 files


  sample.names <- sapply(strsplit(fnFs, "_R1.fastq"), `[`, 1)
  # Specify the full path to the fnFs and fnRs
  fnFs <- file.path(path, fnFs)
  fnRs <- file.path(path, fnRs)
  filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
  filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
  out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncLen=truncLen, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
  if(!is.null(novaseqLoess)){
    errF <- dada2::learnErrors(filtFs, multithread=TRUE, nbases=1e8, errorEstimationFunction=novaseqLoess, verbose=T)
    errR <- dada2::learnErrors(filtRs, multithread=TRUE, nbases=1e8, errorEstimationFunction=novaseqLoess, verbose=T)
  }else{
    errF <- dada2::learnErrors(filtFs, multithread=TRUE, nbases=1e8, verbose=T)
    errR <- dada2::learnErrors(filtRs, multithread=TRUE, nbases=1e8, verbose=T)
  }

  derepFs <- dada2::derepFastq(filtFs, verbose=TRUE)
  derepRs <- dada2::derepFastq(filtRs, verbose=TRUE)
  # Name the derep-class objects by the sample names
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  dadaFs <- dada2::dada(derepFs, err=errF, multithread=TRUE)
  dadaRs <- dada2::dada(derepRs, err=errR, multithread=TRUE)
  mergers <- dada2::mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=10, verbose=TRUE)
  seqtab <- dada2::makeSequenceTable(mergers)
  seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  if(save){
    fname = paste0(basename, ".seqtab.nochim.rds")
    saveRDS(seqtab.nochim, fname)
  }
  dim(seqtab.nochim)
  sum(seqtab.nochim)/sum(seqtab)
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
  rownames(track) <- sample.names
  fname = paste0(basename, ".track.csv")
  write.table(track, fname, sep="\t", col.names=NA)

  return(track)
}
