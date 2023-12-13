#' @export analyzeAll

analyzeAll <- function(
    asvtab = NULL,
    taxonomy = NULL,
    tree = NULL,
    design = NULL,
    picrust.path = "", # path to picrust2 scripts
    picrust.level = c('ko', 'pathway'),
    picrust.desc = NULL,
    color.var = NULL,
    shape.var = NULL,
    varpart.var = NULL,
    basename = NULL,
    rand.seed = 667,
    alpha.diversity = T,
    beta.diversity = T,
    taxonomy.plots = T,
    generate.picrust = T,
    diff.features = c("ASV", "KOs"),
    generate.OTU = F,
    OTU.levels = c('0.03') # at what dissimilarity level OTUs need to be constructed?
){
  possible_diff.features = c("ASV", "OTU", "KOs", "pathways", "taxa")
  possible_taxa_levels = c("kingdoms", "phyla", "classes", "orders", "families", "genera", "species")
  set.seed(rand.seed)
  basename <- paste0(basename, "_", shape.var);
  message(basename)
  mgmisc1::plotAlphaDiversity(asvtab, sampledata, variable=shape.var, basename=basename, comps=as.list(as.data.frame(combn(levels(as.factor(sampledata[[shape.var]])), 2))))
  mgmisc1::plotBetaDiversity(asvtab, sampledata, variable=shape.var, color.variable=color.var, varpart.variable=varpart.var, basename=basename, dists=c( 'GUniFrac'), unifrac_types=c('d_0.5'), tree=tree, statistics=T, rand.seed = rand.seed)
  if( findDRASVs ){
    mgmisc1::findDifferentiallyAbundantFeatures(as.formula(paste("~",  grouping.variable)), asvtab=asvtab, data=sampledata, basename=basename, annot=taxonomy, nproc=4, write.final=T )
  }
  if( !is.null(taxonomy) ){
    taxlist <- mgmisc1::plotTaxonomy(asvtab, data=sampledata, taxonomy=taxonomy, variable=shape.var, basename=basename )
    resdf <- data.frame(matrix(nrow=0, ncol=14))
    difftaxa_fname <- paste0(basename, ".difftaxa.csv")
    if( findDRtaxa ){
      for( tl in names(taxlist) ){
        message(paste0("Analyzing ", tl, " level"))
        atab <- taxlist[[tl]]
        nbasename <- paste0(basename, ".", tl)
        cus <- mgmisc1::findDifferentiallyAbundantFeatures(as.formula(paste("~",  shape.var)), asvtab=atab, annot=taxonomy, data=sampledata, basename=nbasename, nproc=4, write.final=F )
        row <- rep( tl, 14 )
        names(row) <- colnames(cus)
        colnames(resdf) <- colnames(cus)
        resdf <- rbind( resdf, row )
        colnames(resdf) <- colnames(cus)
        resdf <- rbind( resdf, cus )
      }
      write.table( resdf, difftaxa_fname, sep="\t", col.names=NA )
    }
  }
  if(!is.null(generate.picrust)){
    picrust=mgmisc::preparePICRUSt2()
    basename=paste0(picrust, "_", shape.var)
    message(basename)
    mgmisc::plotAlphaDiversity(picrust, sampledata, variable=shape.var, basename=basename, comps=as.list(as.data.frame(combn(levels(as.factor(sampledata[[shape.var]])), 2))))
    mgmisc::plotBetaDiversity(picrust, sampledata, variable=shape.var, color.variable=color.var, varpart.variable=varpart.var, basename=basename, dists=c( 'bray', 'horn'), unifrac_types=NULL, tree=NULL, statistics=T)
    if(findDRPicrust){
      mgmisc::findDifferentiallyAbundantFeatures(as.formula(paste("~",  grouping.variable)), asvtab=kotab, data=sampledata, basename=basename, annot=kodesc, nproc=4, write.final=T )
    }
  }
  return(taxlist)
}

