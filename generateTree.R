#' @export generateTree

generateTree <- function(
    seqtab = NULL,
    fasta.file = "seqs.fasta",
    count.file = "seqs.count_table",
    aligned = FALSE,
    alignment.file = NULL,
    pairwise.alignment = FALSE,
    reference.alignment = NULL,
    start = NULL,
    stop = NULL,
    distance.file = NULL,
#    algorithm = "rNJ",
    mothur = "mothur",  # path to a mothur executable, default assumes that mothur is in $PATH
    processors = 1) {

    if(!is.null(seqtab)){
      bname="seqs."
    }
    bname <- gsub("fasta", "", fasta.file)
    if(!xor(is.null(fasta.file), is.null(seqtab)) ){
      message("Either fasta file or seqtab must be given")
      stop()
    }

    if(!aligned & is.null(reference.alignment) & !pairwise.alignment){
      message("If sequences are not aligned and no reference alignment was given, pairwise.alignment must be set to TRUE")
     stop()
    }
    tree <- NULL

    if(is.null(fasta.file)){
      message("Generating fasta and count_table")
      temp <- mgmisc1::generateFasta(seqtab)
      fasta.file <- temp[1]
      count.file <- temp[3]
    }

    if( !is.null(distance.file)){
      bname <- sub("dist", "", distance.file )
      cmd <- paste0("\"#clearcut(phylip=", distance.file, "); quit()\"")
      system2(mothur, cmd, stdout=F)
    } else {
      if( aligned ){
        message(paste0("Generating tree using aligned sequences"))
        # distance matrix is calculated and saved in the phylip format
        cmd <- paste0("\"#dist.seqs(fasta=", alignment.file, ", output=lt, processors=", processors, "); quit()\"")
        system2(mothur, cmd, stdout=F)
        bname <- sub("")
        d <- paste0(bname, ".dist")
        # tree is computed using relaxed neighbor joining algorithm (Sheneman et al.)
        cmd <- paste0("\"#clearcut(phylip=", d, "); quit();\"")
        system2(mothur, cmd, stdout=F)
      } else {
        if( pairwise.alignment){
          message("Generating tree using pairwise alignments")
          cmd <- paste0(mothur, " \"#pairwise.seqs(fasta=", fasta.file, ", output=lt, processors=", processors, "); clearcut(phylip=current ); quit();\"")
          system(cmd)
        } else {
          message("Generating tree using reference alignment")
          cmd <- paste0(mothur, " \"#align.seqs(fasta=", fasta.file, ", reference=", reference.alignment, ", processors=", processors, "); quit();\"")
          system(cmd)
          alignment <- paste0(bname, "align")
          seqs_to_remove <- paste0(bname, "flip.accnos")
          if(file.exists(seqs_to_remove)){
            cmd <- paste0(mothur, " \"#remove.seqs(fasta=", fasta.file, ", count=", count.file, ", accnos=", seqs_to_remove, "); quit();\"")
            system(cmd)
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
          bname <- paste0(bname, "good.filter.phylip.")
          # tree is computed using relaxed neighbor joining algorithm (Sheneman et al.)
          cmd <- paste0("\"#clearcut(phylip=", dist, "); quit();\"")
          system2(mothur, cmd, stdout=F)
        }
      }
    }

    tree.file <- paste0(bname, ".tre")
    tree <- ape::read.tree(treename)
    message("Tree generated")
    return(tree)
}
