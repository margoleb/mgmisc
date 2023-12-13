#' @export preparePICRUSt2

preparePICRUSt2 <- function(
  asvtab = asvtab,
  sampledata = sampledata,
  variable = NULL,
  fasta.file = "seqs.fasta",
  filename.base = "bacteria",
  conda = "conda", # path to conda
  condaenv = "picrust2", # name of conda environment in which PICRUSt2 is installed
  picrust = "picrust2_pipeline.py",
  picrust.outdir = "picrust2_out",
  plot.nsti = T,
  processors = 1
){
  if(plot.nsti & is.null(variable)){
    stop("A variable must be given if an NSTI plot is to be generated")
  }
  if(dir.exists(picrust.outdir)){
    stop("PICRUSt2 output directory exists, either delete or specify other dir")
  }

  filename.biom = paste0(filename.base, ".biom")
  creation.date <- date()
  biom <- mgmisc1::prepareBIOM(asvtab=asvtab, filename=filename.biom, write=T)
  message("Biom file generated")

  if(!is.null(condaenv)){
    reticulate::use_condaenv(conda=conda, condaenv=condaenv, required=T)
    reticulate::py_config() # python initialization, setting up env vars so that contents of conda environment's bin directory could be accessed
    message(condaenv, " conda environment activated")
  }
  picrust.version = system2(picrust, "--version")
  message("Now running PICRUSt2 pipline with ", processors, " cores")
  cmd <- paste0( " -s ",  fasta.file, " -i ", filename.biom, " -o ", picrust.outdir, " -p ", processors)
  message(picrust, cmd)
  system2(picrust, cmd)

  cmd <- paste0("-i ", picrust.outdir, "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz")
  system2("add_descriptions.py", cmd)

  cmd <- paste0("-i ", picrust.outdir, "/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o pathways_out/path_abun_unstrat_descrip.tsv.gz")
  system2("add_descriptions.py", cmd)

  # KOs
  f <- paste0(picrust.outdir, "/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv")
  picrust_ko <- read.table(f, sep='\t', dec='.', header=T, row.names=1, quote='"')
  ko_descriptions <- data.frame(row.names=rownames(picrust_ko), fun=rownames(picrust_ko), desc=picrust_ko$description)
  picrust_ko$description <- NULL
  picrust_ko <- t(picrust_ko)

  # pathways
  f <- paste0(picrust.outdir, "/pathways_out/path_abun_unstrat_descrip.tsv")
  picrust_pathways <- read.table(f, sep='\t', dec='.', header=T, row.names=1, quote='"')
  pathways_descriptions <- data.frame(row.names=rownames(picrust_pathways), fun=rownames(picrust_pathways), desc=picrust_pathways$description)
  picrust_pathways$description <- NULL
  picrust_pathways <- t(picrust_pathways)

  # NSTI
  filename.nstiplot <- paste0(filename.base, "_NSTI.svg")
  nsti <- read.table("picrust2_out/KO_metagenome_out/weighted_nsti.tsv", header=T, sep="\t", dec=".")
  if(plot.nsti){
    rownames(nsti) <- nsti$sample
    nsti$variable <- variable
    svg(filename.nstiplot, width=7, height=7, pointsize=3)
    p <- ggplot( nsti, aes(x=variable, y=nsti)) + geom_boxplot() + xlab(variable) + ylab("weighted NSTI") + theme(text=element_text(size=6))
    print( p )
    dev.off()
  }

  result < list()
  result$KOs <- picrust_ko
  result$pathways <- picrust_pathways
  result$nsti <- nsti

  return(result)

}
