#' @export prepareDataFromDada2

prepareDataFromDada2 <- function(
    seqtab = seqtab.nochim,
    taxa = NULL,
    classification.db = NULL,
    design = design,
    design.rowname = NULL,
    remove.rare = 2, # removing global singletons and doubletons
    remove.taxa = c("Chloroplast", "Mitochondria"),
    keep.only = "Bacteria",
    generate.tree = TRUE,
    pairwise.alignment = FALSE,
    reference.alignment=NULL,
    rarefy = 100,
    depth = 2000,
    basename = "bacteria",
    rand.seed = 667,
    fasta.file = "seqs.fasta", # name of a fasta file with ASVs sequences
    count.file = "seqs.count_table", # name of a count file with ASVs number per sample
    OTU = FALSE,
    save = TRUE, # save output as an rda object?
    processors = 1)

{
#  out = mgmisc1::new_mgmisc1(seqtab=seqtab, sampledata=design, taxa=taxa, tree=tree, rand.seed=rand.seed, orig.fasta.file=fasta.file, orig.count.file=count.file)


  if(is.null(design.rowname)){
    message("No name of column in the design was given as samplenames (they must match rownames of the seqtab)")
    stop()
  }else if(! design.rowname %in% colnames(design)){
    message(paste0(design.rowname, "is not a name of a column in design"))
    stop()
  }else if(length(rownames(seqtab) %in% design[[design.rowname]]) == 0 ){
    message("Samplenames in the design do not match row names of the seqtab")
  }
  if( !is.numeric(remove.rare) & remove.rare != FALSE){
    message("remove.rare must be either FALSE or a non-negative integer")
    stop()
  }
  if( !is.numeric(rarefy) & rarefy != FALSE ){
    message("Number of rarefaction iterations either must be a positive integer or 'rarefy' must be FALSE")
    stop()
  } else if( !is.numeric(depth) | depth <= 1){
    message("Rarefying depth must be an integer greater than 1")
    stop()
  }
  if( is.null(taxa) & (!is.null(remove.taxa) | !is.null(keep.only)) ){
    if(!is.null(classification.db)){
      message("Classifying sequences")
      taxa <- dada2::assignTaxonomy(seqtab, classification.db, multithread=T)
      if(save){
        fname <- paste0(basename, ".taxa.rds")
        saveRDS(taxa, fname)
      }
    } else {
      message("No object with taxonomy was given, and the sequences cannot be classfied as classification.db is NULL")
      stop()
    }
  }



  set.seed(rand.seed)

  f <- mgmisc1::generateFasta(seqtab)

  asvtab <- as.data.frame(seqtab)
  colnames(asvtab) <- paste("ASV", seq(1, dim(seqtab)[2], 1), sep="")
  taxa <- as.data.frame(taxa)
  rownames(taxa) <- paste("ASV", seq(1, dim(seqtab)[2], 1), sep="")



  if(remove.rare){
    message(paste0("Removing ASVs with abundance less than ", remove.rare, " across all samples"))
    seqtab.final <- seqtab[ , colSums(seqtab) >= remove.rare ]
    asvtab <- asvtab[ , colSums(asvtab) >= remove.rare ]
    if(OTU){
     for( diss in names(outlist$otutab)){
       message(paste0("Removing ", diss, " OTUs with abundance less than ", remove.rare, "across all samples"))
       otutab.final <- outlist$otutab[[diss]][, colSums(outlist$outab[[diss]]) >= remove.rare]
       outlist$otutab[[diss]] <- otutab.final
      }
    }

  } else {
    seqtab.final <- seqtab
  }


  if( !is.null(taxa) ){
    taxa <- as.data.frame(taxa)
    rownames(taxa) <- paste("ASV", seq(1, dim(seqtab)[2], 1), sep="")
    taxa <- taxa[ rownames(taxa) %in% colnames(asvtab), ]
#    outlist$taxa <- taxa
  }

  if(!is.null(remove.taxa)){
    message("Removing unwanted taxa")
    taxa$taxonomy <- paste(taxa$Kingdom, taxa$Phylum, taxa$Class, taxa$Order, taxa$Family, taxa$Genus, sep=";")
    pat <- paste(remove.taxa, collapse="|")
    taxa <- taxa[ grep(pat, taxa$taxonomy, invert=T), ]
    message(paste0("After removing unwanted taxa ", nrow(taxa), " remained"))
    if(nrow(taxa) == 0){
      stop()
    }
    asvtab <- asvtab[ , colnames(asvtab) %in% rownames(taxa)]
  }


  if(!is.null(keep.only)){
    message(paste0("Removing taxa which are not ", paste0(keep.only)))
    taxa <- taxa[grep(keep.only, taxa$taxonomy),]
    taxa$taxonomy <- NULL
    asvtab <- asvtab[ , colnames(asvtab) %in% rownames(taxa) ]
    message(paste0("After removing ", ncol(asvtab), " left"))
  }

  if(is.logical(OTU) & OTU == TRUE){
    otus <- constructOTUs( fasta.file, count.file, 0.03 )
    outlist$otutab[["0.03"]] <- otus$otutab
    outlist$otutree[["0.03"]] <- otus$tree
  } else if(is.numeric(OTU)) {
    for( diss in OTU ){
      message(paste("Generating OTUs at", diss, "dissimilarity level"))
      otus <- constructOTUs( fasta.file, count.file, diss )
      outlist$otutab[[as.character(diss)]] <- otus$otutab
      outlist$otutree[[as.character(diss)]] <- otus$tree
      message("Done")
    }
  }

  if( generate.tree ){
    message("Generating a tree")
    tree <- mgmisc1::generateTree(fasta.file=fasta.file, count.file=count.file, reference.alignment=reference.alignment, pairwise.alignment=pairwise.alignment, processors=processors)
    tree <- phangorn::midpoint(tree)
    tree$edge.length[which(is.na(tree$edge.length))] <- 0
#    outlist$asvtree <- tree
    asvtab <- asvtab[ , colnames(asvtab) %in% tree$tip.label ] # only ASVs that are in the tree
    tree <- ape::drop.tip(tree, tree$tip.label[!(tree$tip.label %in% colnames(asvtab))]) # only ASVs that are in the seqtab (it can be rarefied)
#    outlist$asvtab <- seqtab.final
  }


  if(rarefy){
    message(paste0("Rarefying to ", depth, " reads ", rarefy, " times"))
    temp <- GUniFrac::Rarefy(otu.tab=asvtab, depth=depth)$otu.tab.rff
    for( i in 1:rarefy ){ temp <- temp + GUniFrac::Rarefy(otu.tab=asvtab, depth=depth)$otu.tab.rff }
    seqtab.final <- round( temp/100 )
    outlist$asvtab.rarefied <- seqtab.final
    message(paste0("After rarefying ", nrow(asvtab), " samples and ", ncol(asvtab), " ASVs left"))
  }



  if( !is.null(design) ){
    rownames(design) <- design[[design.rowname]]
    design <- design[ rownames(design) %in% rownames(asvtab), ]
    asvtab <- asvtab[ rownames(asvtab) %in% rownames(design), ]
    design <-design[ order(rownames(design)), ]
    asvtab <- asvtab[ order(rownames(asvtab)), ]
    stopifnot(identical(rownames(asvtab), rownames(design))) # must be TRUE!
  }




  asvtab.name <- paste0(basename, ".asvtab")
  taxa.name <- paste0(basename, ".tax")
  tree.name <- paste0(basename, ".tree")
  design.name <- paste0(basename, ".design")

  assign(asvtab.name, asvtab, envir=parent.frame())
  assign(taxa.name, taxa, envir=parent.frame())
  assign(tree.name, tree, envir=parent.frame())
  assign(design.name, design, envir=parent.frame())

  return(outlist)

}
