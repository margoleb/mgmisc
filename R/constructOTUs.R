#' constructOTUs
#'
#' Cluster ASVs into OTUs. Opticlust algorithm (Schloss and Westcott, 2016?) is used. Tree is computed and classification of represetative sequences is performed if needed.
#'
#' @param seqtab A named integer matrix with samples in rows and ASVs in columns. Column names should be actual ASVs' sequences.
#' @param fasta.file A path to a fasta file containing ASVs' sequences. Needs to be given if there is no seqtab. Defaults to 'seqs.fasta'.
#' @param count.file A path to a count_table file with abundances of ASVs. Needs to be given if there is no seqtab. Defaults to 'seqs.count_table'.
#' @param dist.file A path to a distance matrix file in column or phylip format. Optional. Defaults to NULL.
#' @param diss A vector of numbers between 0 and 1 indicating at which dissimilarities OTUs should be constructed. Defaults to c(0.03).
#' @param aligned A logical indicating if sequences in fasta.file are aligned. Defaults to FALSE.
#' @param reference.alignment A path to a reference alignment in fasta format (e.g. SILVA or Greengenes). Defaults to NULL.
#' @param start An integer indicating beginning of a kept region of an alignment or NULL in which case it will be calculated so that min. 75\% of sequences are kept. Defaults to NULL.
#' @param stop An integer indicating end of a kept region of an alignment or NULL in which case it will be calculated so that min. 75\% of sequences are kept. Defaults to NULL.
#' @param pairwise.alignment A logical indicating if pairwise alignment should be performed. Defaults to FALSE.
#' @param classification.db A path to a classification database suitable for use with dada2's assignTaxoonmy. Defaults to NULL.
#' @param mothur A path to a MOTHUR executable. Defaults to 'mothur'.
#' @param processors An integer indicating how many processors should be used. Defaults to 1.
#'
#' @returns A list of featuresets, one for each of dissimilarities.
#'
#' @export constructOTUs


constructOTUs <- function(
    seqtab = NULL,
    fasta.file = "seqs.fasta",
    count.file = "seqs.count_table",
    dist.file = NULL, # path to a pre-calculated distance file (in column or phylip format)
    diss = c(0.03),
    aligned = F,
    reference.alignment = NULL,
    start=NULL,
    stop=NULL,
    pairwise.alignment = FALSE,
    classification.db = NULL,
    mothur = "mothur", # path to mothur executable
    processors = 1){

  if( system2(mothur, "--version") != 0 ){
    stop("constructOTUs: ", mothur, " is not a functional mothur binary")
  }

  if(!xor(is.null(fasta.file), is.null(seqtab)) ){
    message("constructOTUs: either a fasta file and a matching count file or a dada2-produced seqtab must be given")
    stop()
  }

  if(is.null(reference.alignment) & !pairwise.alignment){
    message("constructOTUs: if no reference alignment was given, pairwise.alignment must be set to TRUE")
    stop()
  }

  if(is.null(fasta.file) & !is.null(seqtab )){
    message("constructOTUs: generating fasta and count_table")
    temp <- mgmisc1::generateFasta(seqtab)
    fasta.file <- temp[1]
    count.file <- temp[3]
  }

  out <- c()
  if( is.null(dist.file) ){
    message("constructOTUs: generating distance matrix")
    res <- generateDistMat(fasta.file = fasta.file, count.file = count.file, aligned=aligned, start = start, stop = stop, format="column", pairwise.alignment = pairwise.alignment, reference.alignment = reference.alignment, mothur = mothur, processors=processors)
    dist.file <- res$distance.matrix
    count.file <- res$count.file
    alignment <- res$fasta.file
  }else{
    message("constructOTUs: constructing OTUs based on ", dist.file)
    if(grep("phylip", dist.file)){
      stop("Phylip formatted distance matrix is incompatible with the opticlust algorithm")
    }
  }

  # generating OTUs
  for(cutoff in diss){
    cutoff <- as.character(cutoff)
    message(paste0("constructOTUs: constructing OTUs for ", cutoff, " dissimilarity using ", dist.file, " and ", count.file))
    # clustering using OptiMcc
    cmd <- paste0("\"#cluster(column=", dist.file, ", count=", count.file, ", cutoff=", cutoff, ", processors=", processors, "); quit()\"")
    stdout <- system2(mothur, cmd, stdout=T)
    list.file <- stdout[ grep("\\.opti_mcc\\.list", stdout) ]
    # generating shared OTU table
    message("constructOTUs: preparing OTU table")
    cmd <- paste0("\"#make.shared(list=", list.file, ", count=", count.file, ", label=", cutoff, "); quit()\"")
    stdout <- system2(mothur, cmd, stdout=T)
    shared.file <- stdout[ grep("\\.opti_mcc\\.shared", stdout) ]
    # reading in shared file
    otutab <- as.data.frame(read.table(shared.file, header=T, sep="\t"))
    rownames(otutab) <- otutab$Group
    otutab$Group <- NULL
    otutab$numOtus <- NULL
    otutab$label <- NULL
    out[[cutoff]]$featuretab <- otutab
    # finding representative sequences
    message("constructOTUs: finding representative sequences")
    cmd <- paste0("\"#get.oturep(list=", list.file, ", column=", dist.file, ", count=", count.file, ", cutoff=", cutoff, ", fasta=", fasta.file, "); quit()\"")
    stdout <- system2(mothur, cmd, stdout=T)
    rep.fasta.file <- stdout[grep("\\.rep\\.fasta", stdout)]
    rep.count.file <- stdout[grep("\\.rep\\.count_table", stdout)]
    message(paste0("constructOTUs: representative sequences in ", rep.fasta.file, " count table: ", rep.count.file))
    # changing fasta ids from ASVs to OTUs
    f <- ShortRead::readFasta(".", rep.fasta.file)
    f <- ShortRead::ShortRead(sread=ShortRead::sread(f), id=Biostrings::BStringSet(sub("\\|.*", "", sub("ASV.*O", "O", ShortRead::id(f)))))
    ShortRead::writeFasta(f, rep.fasta.file, mode='w')
    cmd <- paste0("\"#degap.seqs(fasta=", rep.fasta.file, "); quit()\"")
    stdout <- system2(mothur, cmd, stdout=T)
    rep.fasta.ng.file <- stdout[grep("\\.rep\\.ng\\.fasta", stdout)]
    out[[cutoff]]$fasta.file <- rep.fasta.ng.file
    # preparing seqtab for classification
    message("constructOTUs: classifying representative sequences")
    f <- ShortRead::readFasta(".", rep.fasta.file)
    seqtab <- data.frame(sequence=as.vector(ShortRead::sread(f)), abundance=rep(1, length(as.vector(ShortRead::id(f)))), name=as.vector(ShortRead::id(f)))
    taxa <- dada2::assignTaxonomy(seqtab, classification.db, multithread=T)
    rownames(taxa) <- seqtab$name
    out[[cutoff]]$featureannot <- taxa
    # generating tree
    message(paste0("constructOTUs: generating tree for OTUs at ", cutoff, " dissimilarity level"))
    if(pairwise.alignment){
      message(paste0("constructOTUs: using pairwise alignment, aligning sequences in ", rep.fasta.file))
      d <- mgmisc1::generateDistMat(fasta.file=rep.fasta.file, count.file=rep.count.file, pairwise.alignment=T, aligned=F, format='phylip', processors=processors)
    }else{
      message(paste0("constructOTUs: using sequences aligned against a reference alignment in ", rep.fasta.file))
      d <- mgmisc1::generateDistMat(fasta.file=rep.fasta.file, count.file=rep.count.file, pairwise.alignment=F, aligned=T, format='phylip', processors=processors)
    }
    message(paste0("constructOTUs: distance matrix in ", d$distance.matrix))
    t <- mgmisc1::generateTree(distance.file=d$distance.matrix)
    out[[cutoff]]$tree <- t$tree
    out[[cutoff]]$tree.file <- t$tree.file
  }

  return(out)
}
