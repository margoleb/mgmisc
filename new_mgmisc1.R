#' @export new_mgmisc1


new_mgmisc1 <- function(
  rand.seed = 667,
  sampledata = NULL, # path to a tsv file with rows representing samples and cols representing variables
  reads.dir = NULL, # path to a directory where fastq(.gz) files reside
  reads.used = "both", # whether reads are to be merged (both), or only r1 or r2 are to be used
  mothur.version = NULL,
  classification.db = NULL,
  seq.type = "16S"
  ){
    stopifnot(class(seqtab) == "matrix", class(design) == "data.frame", is.integer(rand.seed), rand.seed %% 2 == 1)
    date.created = date()
    l <- mgmisc1::generateFasta(seqtab)
    orig = list(
      seqtab = seqtab,
      taxa = taxa,
      tree = tree,
      fasta = l$fasta,
      fasta.filename = l$fasta.filename,
      count.filename = l$count.filename,
      samplenames = colnames(seqtab),
      asvnames = l$count$Representative_sequence
    )


  ASV <- list(
    orig = orig,
    rarefied = NULL
  )

  OTU <- list()

  value = list(
    date.created = date.created,
    date.modified = date.created,
    rand.seed = rand.seed,
    design = design,
    mothur.version = mothur.version,
    seq.type = seq.type,
    classification.db = classification.db,
    ASV = asv,
    OTU = otu
  )
  class(value) <- "mgmisc"

  return(value)
}
