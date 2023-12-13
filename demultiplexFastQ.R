#' @export demultiplexFastQ

demultiplexFastQ <- function(
    reads1, reads2,
    primerF, primerR,
    design,
    outdir,
    max_mismatch,
    ncores) {

  destinations <- gsub("-", "_", as.character(design[, 1]))

  ## get the fastq to split (raise error if fastq untrimmed/not existing)
  fastq_R1 <- reads1
  fastq_R2 <- reads2

  # register parallel backend
  param <- MulticoreParam(workers = ncores)
  if (!BiocParallel::bpisup(param)) {
    BiocParallel::bpstart(param)
    on.exit(BiocParallel::bpstop(param))
  }
  message("de-multiplexing the FASTQ file")
  ## filter and write
  info <-
    BiocParallel::bplapply(seq_along(destinations), function(i) {
      split1 <- file.path(outdir, paste0(destinations[i], "_R1.fastq.gz"))
      split2 <-
        file.path(outdir, paste0(destinations[i], "_R2.fastq.gz"))
      ## stop if the outfiles already exist (otherwise the output would be appended)
      if (file.exists(split1) | file.exists(split2)) {
        warning("Output files already exist! Overwriting...")
        if (file.exists(split1))
          file.remove(split1)
        if (file.exists(split1))
          file.remove(split2)
      }

      ## split files and save kept read number
      kept <- splitFastQNovogene(
        fastq_R1,
        fastq_R2,
        primerF = primerF,
        primerR = primerR,
        outfile_R1 = split1,
        outfile_R2 = split2,
        ind_f = design[ i, 2 ],
        ind_r = design[ i, 3 ],
        max_mismatch = max_mismatch
      )
      dfout <- data.frame(R1 = split1,
                          R2 = split2,
                          kept_reads = kept)
      return(dfout)
    }, BPPARAM = param)

}

