#' @export generateFasta

generateFasta <- function(
  seqtab=seqtab.final,
  fasta.file="seqs.fasta",
  count.file="seqs.count_table"){

  seqs1.fasta <- data.frame(seq=colnames(seqtab), row.names=paste(">ASV", seq(1, dim(seqtab)[2], 1), "\n", sep=""))
  write.table(seqs1.fasta, fasta.file, quote=F, sep="", col.names=F)
  seqs1.count_table <- as.data.frame(t(seqtab))
  rownames(seqs1.count_table) <- paste(">ASV", seq(1, dim(seqtab)[2], 1), "\n", sep="")
  seqs1.count_table$total <- rowSums(seqs1.count_table)
  seqs1.count_table <- seqs1.count_table[, c(ncol(seqs1.count_table), 1:(ncol(seqs1.count_table)-1))]
  seqs1.count_table$Representative_Sequence <- paste("ASV", seq(1, dim(seqtab)[2], 1), sep="")
  seqs1.count_table <- seqs1.count_table[, c(ncol(seqs1.count_table), 1:(ncol(seqs1.count_table)-1))]
  write.table(seqs1.count_table, count.file, quote=F, sep="\t", row.names=F)

  out <- list(fasta.file, seqs1.fasta,  count.file, seqs1.count_table)
  return(out)
}
