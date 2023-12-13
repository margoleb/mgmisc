#' @export constructOTUs


constructOTUs <- function(
    seqtab = NULL,
    fasta.file = "seqs.fasta",
    count.file = "seqs.count_table",
    diss = c(0.03),
    reference.alignment = NULL,
    start=NULL,
    stop=NULL,
    pairwise.alignment = FALSE,
    classification.db = NULL,
    mothur = "mothur", # path to mothur executable
    processors = 1){

  if(!is.null(seqtab)){
    bname="seqs."
  }
  bname <- gsub("fasta", "", fasta.file)
  if(!xor(is.null(fasta.file), is.null(seqtab)) ){
    message("Either fasta file or seqtab must be given")
    stop()
  }

  if(is.null(reference.alignment) & !pairwise.alignment){
    message("If no reference alignment was given, pairwise.alignment must be set to TRUE")
    stop()
  }

  if(is.null(fasta.file)){
    message("Generating fasta and count_table")
    temp <- mgmisc1::generateFasta(seqtab)
    fasta.file <- temp[1]
    count.file <- temp[3]
  }

  out <- list()

  if( pairwise.alignment){
    message("Generating OTUs using pairwise alignments")
    cmd <- paste0("\"#pairwise.seqs(fasta=", fasta.file, ", output=lt, processors=", processors, "); quit();\"")
    d <- paste0(bname, ".phylip.dist")
    list.file <- paste0(bname, ".phylip.optimcc.list")
    shared.file <- paste0(bname, ".phylip.opti_mcc.shared")
    rep.file <- paste0(bname, ".phylip.opti_mcc.")
    system2(mothur, cmd, stdout=F)
  }else{
    message("Generating OTUs using reference alignment")
    cmd <- paste0("\"#align.seqs(fasta=", fasta.file, ", reference=", reference.alignment, ", processors=", processors, "); quit();\"")
    system2(mothur, cmd, stdout=F)
    alignment <- paste0(bname, "align")
    seqs_to_remove <- paste0(bname, "flip.accnos")
    if(file.exists(seqs_to_remove)){
      cmd <- paste0(" \"#remove.seqs(fasta=", fasta.file, ", count=", count.file, ", accnos=", seqs_to_remove, "); quit();\"")
      system2(mothur, cmd, stdout=F)
      alignment <- paste0(bname, "pick.align")
      count.file <- paste0(bname, "pick.count_table")
      bname <- paste0(bname, "pick.")
    }
    # if no desired alignment regions was given (start and stop), it is calculated as 75th percentile start and stop
    if(is.null(start) | is.null(stop)){
      message("Determining start and stop positions")
      cmd <- paste0("\"#count.seqs(count=", count.file, "); quit()\"")
      count.file <- paste0(bname, "sparse.count_table");
      system2(mothur, cmd, stdout=F)
      cmd <- paste0("\"#summary.seqs(fasta=", alignment, ", count=", count.file, ", processors=", processors, "); quit();\"")
      stdout <- system2(mothur, cmd, stdout=T)
      stdout
      start <- strsplit(stdout[39], "\t")[[1]][2]
      stop <- strsplit(stdout[39], "\t")[[1]][3]
      message(paste0("Alignment trimmed to ", start, " - ", stop, " position"))
    }
    # screening for sequences covering the desired region of the alignment
    cmd <- paste0("\"#screen.seqs(fasta=", alignment, ", count=", count.file, ", start=", start, ", end=", stop, ", processors=", processors, "); quit()\"")
    system2(mothur, cmd, stdout=F)
    alignment <- paste0(bname, "good.align")
    # all-gap columns and columns containing at least one terminal gap are removed
    cmd <- paste0("\"#filter.seqs(fasta=", alignment, ", vertical=T, trump=., processors=", processors, "); quit()\"")
    system2(mothur, cmd, stdout=F)
    alignment <- paste0(bname, "good.filter.fasta")
    # distance matrix is calculated and saved in the phylip format
    cmd <- paste0("\"#dist.seqs(fasta=", alignment, ", output=lt, processors=", processors, "); quit()\"")
    system2(mothur, cmd, stdout=F)
    d <- paste0(bname, "good.filter.phylip.dist")
    list.file <- paste0(bname, ".good.filter.phylip.optimcc.list")
    shared.file <- paste0(bname, ".good.filter.phylip.opti_mcc.shared")
    rep.file <- paste0(bname, ".good.filter.phylip.opti_mcc.")
  }

  # generating OTUs
  for(cutoff in diss){
    # clustering using OptiMcc
    cmd <- paste0("\"#cluster(phylip=", d, ", cutoff=", cutoff, "); quit()\"")
    system2(mothur, cmd, stdout=F)
    # generating shared OTU table
    cmd <- paste0("\"#make.shared(list=", list.file, ", count=", count.file, "); quit()\"")
    system2(mothur, cmd, stdout=F)
    # reading in shared file
    otutab <- read.table(shared.file, header=T, sep="\t")
    rownames(otutab <- otutab$Group)
    otutab$Group <- NULL
    otutab$numOtus <- NULL
    otutab$label <- NULL
    out[paste0(cutoff)]$otutab <- otutab
    # finding representative sequences
    cmd <- paste0("\"#get.otureps(list=", list.file, ", phylip=", d, ", label=", cutoff, "fasta=", fasta.file, "); quit()\"")
    system2(mothur, cmd, stdout=F)
    rep.fasta.file <- paste0(rep.file, cutoff, ".rep.fasta")
    cmd <- paste0("\"#degap.seqs(fasta=", rep.fasta.file, "); quit()\"")
    rep.ng.fasta.file <- paste0(rep.file, cutoff, ".rep.ng.fasta")
    rep.count.file <- paste0(rep.file, cutoff, ".rep.count_table")
    system2(mothur, cmd, stdout=F)
    # preparing seqtab for classification
    f <- ShortRead::readFasta(".", rep.ng.fasta.file)
    seqtab <- data.frame(sequence=as.vector(ShortRead::sread(f)), abundance=rep(1, length(as.vector(ShortRead::id(f)))), name=as.vector(ShortRead::id(f)))
    taxa <- dada2::assignTaxonomy(seqtab, classification.db, multithread=T)
    rownames(taxa) <- seqtab$name
    out[paste0(cutoff)]$taxonomy <- taxa
    # generating tree

    tree <- mgmisc::generateTree(distance.file=d)
    out[paste0(cutoff)]$tree <- tree


  }

  return(out)
}
