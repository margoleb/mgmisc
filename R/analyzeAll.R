#' analyzeAll
#'
#' @description
#' Analyze all features in a mg-class experiment.
#'
#' @details The function performs analyzes on various types of features resulting from one experiment. Currently the set of possible features comprises ASVs, OTUs, KO functions and pathways (from PICRUSt2 analysis). In the future other types of features will be added, such as genes or transcripts.
#'
#' @param experiment An mg-class experiment.
#' @param parameters A list or a path to a tsv file with parameters. Any other function parameter may apart from experiment may be set in this object/file. If parameters are set both in file/object, the settings can be overriden in the function call. Defaults to NULL.
#' @param dissimilarities A vector of numbers between 0 and 1 determining for which dissimilarities OTUs should be analyzed (and generated, if necessary). Defaults to c(0.03).
#' @param analyze What to analyze? A vector of strings being names of features to be analyzed. Current possibilities: 'ASVs', 'KOs', 'pathways', 'OTUs'. Defaults to c('ASVs', "OTUs', 'KOs', 'pathways').
#' @param analyzes How to analyze? A vector of strings being names of analyzes to be performed. Current possibilities: 'alpha-diversity', 'beta-diversity', 'taxonomy', 'diff-features'. Defaults to c('alpha-diversity', 'beta-diversity', 'taxonomy', 'diff-features')
#' @param processors Number of processors to use. Defaults to 1.
#' @param basename What should be a common part of produced files' names? Defaults to 'bacteria'.
#' @param outdir A path to directory where resulting plots and other files should be stored. Defaults to ".".
#' @param graphics.type Which graphics type should be used? One of 'ggplot' or 'base'. Defaults to 'ggplot'.
#' @param statistics Should statistical tests be performed? A logical. Defaults to TRUE.
#' @param rand.seed A random seed for reproducibility. Must be an odd positive integer. Defaults to 667 (a new, better beast).
#' @param significance.threshold A number between 0 and 1 determining significance threshold for statistical tests. Defaults to 0.05.
#' @param find.differing.pairs Should pairwise tests be performed to find differing group pairs if overall test is significant? A logical. Defaults to TRUE.
#' @param permu A positive integer determining the number of permutations performed during non-parametric tests such as PERMANOVA. Defaults to 999.
#' @param palette.family Which palette family should be used? Any family included in the paletteer package may be used. Defaults to 'ggsci'.
#' @param palette Which palette should be used? It should be a name of a palette from family chosen with palette.family. Any palette from the paletteer package may be used here.  Defaults to 'npg-nrc' from the 'ggsci' family.
#' @param collect A logical determining if identical legends should be collected in one place of a compound plot. Defaults to FALSE.
#' @param vertical A logical determining if compound plots should be organized vertically (i.e. their width be smaller than height). Defaults to TRUE.
#' @param width A number determining width of a single plot (in inches). Defaults to 3.5.
#' @param height A number determining height of a single plot (in inches). Defaults to 3.5.
#' @param main.variable The main variable to be analyzed. A string indicating a categorical variable from experiment's sampledata. Absolutely necessary. It will be mapped to symbol shape on beta-diversity plots and used as grouping variable in all other analyzes. Defaults to NULL.
#' @param facetting.variables Should facetted composite plots be produced? A vector of one or two variable names. Optional. The names must be chosen from colnames(sampledata) or variableNames(experiment). Defaults to NULL (no facetting).
#' @param sample.color.variable A string indicating which variable from sampledata should be mapped to color of symbols denoting samples. Optional. Must be categorical. Defaults to NULL.
#' @param sample.size.variable A string indicating which variable from sampledata should be mapped to size of symbols denoting samples. Optional. Must be continuous. Defaults to NULL.
#' @param feature.color.variable A string indicating which variable from featureannot should be mapped to color of symbols representing features. Must be categorical
#' @param chemical.variables A vector of strings indicating which variables from sampledata should be used in CCA or (db)RDA. Defaults to NULL.
#' @param varpart.variable A string indicating which variable from sampledata should be used in variance partitioning to enable assessing main variable's influence. Defaults to NULL.
#' @param alpha.div.measures Which species richness, diversity and evenness metrics should be used in alpha-diversity analysis? A vector of index names. Current possibilities: i) species richness: 'sobs' - stands for 'species observed'; 'chao' - Chao1 index (estimated total richness); 'ace' - ACE index (estimated total richness); ii) diversity: 'shannon' - Shannon's H'; 'simpson' - Simpson's D; 'invsimpson' - Simpson's 1/D; iii) evenness: 'shannoneven' - Shannon's E. Defaults to c('sobs', 'shannon', 'shannoneven').
#' @param alpha.div.points A logical determining if poits representing samples should be added to alpha diversity plots. Defaults to TRUE.
#' @param alpha.plot.style A string determining how alpha diversity will be plotted. Current possibilities: 'box', 'violin', 'dot'. Defauts to 'box'.
#' @param alpha.plot.cols A positive integer determining number of columns in compound alpha diversity plots. Defaults to 2.
#' @param beta.metrics Which beta-diversity metrics should be used? GUniFrac or any distance possible to calculate with 'vegdist' can be used. Defaults to c('bray', 'GUniFrac') - Bray-Curtis and GUniFrac.
#' @param unifracs Which variants of Generalized UniFrac should be used as beta-diversity metrics? ... Defaults to c('d_05').
#' @param beta.diversity.analyzes Which beta-diversity analyzes should be performed? Current possibilities: 'nmds', 'rda', 'dbrda'. Defaults to c('nmds', 'dbrda')
#' @param biplot A logical determining if biplots should be generated. Defaults to TRUE.
#' @param triplot A logical determining if triplots should be generated. Requires giving at least one variable in chemical.variables. Defaults to FALSE.
#' @param size.means.abundance Should size of feature symbols on a biplot convey features' abundance? Defaults to TRUE.
#' @param abundance.cutoff How many features should be displayed on a biplot? Most abundant features are displayed. An integer greater than 1. Defaults to 50.
#' @param diff.features.methods Which methods should be used to find differentially abundant features? A vector of method names. Current possibilities: 'deseq' (DESeq2 is used). Defaults to c('deseq').
#' @param log2FC.threshold A number determining a minimal absolute value of log2 fold change of a feature required to consider the feature differentially abundant. Defaults to 1 (meaning two-fold difference).
#' @param rare.threshold Threshold for identification of rare features. All features whose abundance is less than threshold will be lumped to a category named 'rare'. Must be a number between 0 and 1. Defaults to 0.01.
#'
#'
#' @returns An mg-class experiment object with 'analyzes' slots populated with results of chosen analyzes.
#'
#'
#' @export analyzeAll

analyzeAll <- function(
    # general parameters
    experiment = NULL, # an experiment object of class mg or a path to an .rds file containing one; if given, no featuretab, taxonomy, tree, design, fasta.file and count.file can be given.
    parameters = NULL,
    analyze = c('ASVs', 'OTUs', 'KOs', 'pathways'),
    dissimilarities = c(0.03),
    analyzes = c('alpha-diversity', 'beta-diversity', 'taxonomy', 'diff-features'),
    basename = 'bacteria',
    outdir = '.', # where results are to be stored
    graphics.type = 'ggplot', # if base ('base') graphics or ggplot2 ('ggplot') is to be used
    vertical = TRUE,
    collect = FALSE,
    tag_level = 'A',
    statistics = TRUE,
    rand.seed = 667,
    permu = 999,
    significance.threshold = 0.05,
    find.differing.pairs = TRUE,
    palette.family = 'ggsci', # palette family (package), must be one of those gathered in paletteer (https://github.com/EmilHvitfeldt/paletteer?tab=readme-ov-file), if present in the experiment object will be overwritten
    palette = 'nrc_npg', # palette name, must come from the package in 'palette.family', if present in the experiment object will be overwritten
    processors = 1,
    width = 3.5,
    height = 3.5,
    # variables
    main.variable = NULL,
    facetting.variables = NULL,
    varpart.variable = NULL,
    sample.color.variable = NULL,
    sample.size.variable = NULL,
    feature.color.variable = NULL,
    chemistry.variables = NULL,
    # alpha diversity
    alpha.div.measures = c('sobs', 'H', 'E'),
    alpha.div.points = T,
    ref.group = NULL,
    alpha.plot.style = 'box',
    alpha.plot.cols = 2,
    # beta diversity
    beta.diversity.metrics = c('bray', 'GUniFrac'), # one or more of metrics implemented in vegan::vegdist or GUniFrac
    unifrac.types = c('d_0.5'), # GUniFrac metrics (for ASVs and OTUs only)
    beta.diversity.analyzes = c('nmds', 'dbrda'),
    biplot=T,
    triplot=F,
    size.means.abundance=T,
    abundance.cutoff=50,
    ellipses = T,
    vector.col = 'red',
    ordistep.steps=500,
    ordistep.direction='both',
    # differentially abundant features identification parameters
    diff.features.methods=c('deseq'), # one of 'deseq', 'aldex',
    log2FC.threshold = 1,
    mergeby = 0,
    consensus = F, # include only those features which were found to be DR by all used methods in the final results
    plot.diff.features = F, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
    num.top.features = 10,
    write.single = T,
    write.final = T, # write final, combined table to disc
    mc.samples = 128, # how many Monte Carlo samples to use in ALDEx2
    # taxonomy plots parameters
    rare.threshold = 0.01,
    taxlevel.names = c('kingdoms', 'phyla', 'classes', 'orders', 'families', 'genera')
){
  argv <- as.list(match.call()[-1])
  params <- set_parameters(p=parameters, argv)

  set.seed(params$rand.seed)

  if(is.null(experiment)){
    stop("analyzeAll: an object of class 'mg' needs to be given")
  }

  if(!R.utils::isAbsolutePath(params$outdir)){
    outdir <- file.path(getwd(), params$outdir)
  }
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  setwd(outdir)

  if(!is.null(experiment)){
    if(is.character(experiment)){
      experiment <- readRDS(experiment)
    }
    if(!inherits(experiment, "mg")){
      stop("analyzeAll: an experiment to analyze must be of class 'mg'")
    }
    if(!is.null(params$palette.family)){
      experiment$metadata$palette.family <- params$palette.family
    }else if(!is.null(experiment$metadata$palette.family)){
      palette.family <- experiment$metadata$palette.family
    }else{
      palette.family <- 'ggsci'
    }
    if(!is.null(params$palette)){
      experiment$metadata$palette <- params$palette
    }else if(!is.null(experiment$metadata$palette)){
      palette <- experiment$metadata$palette
    }else{
      palette <- 'nrc_npg'
    }

    for( feature in params$analyze ){
      if( feature == 'ASVs'){
        message("analyzeAll: analyzing ASVs")
        d <- paste0("ASVs_", gsub(":", "-", gsub("\\ ", "_", as.character(date()))))
        if(!dir.exists(d)){
          dir.create(d)
        }
        bn <- paste0(params$basename, "_ASVs")
        message("analyzeAll: results will be stored in the ", d)
        setwd(d)
        tmpres <- analyzeFeature(
          experiment=experiment,
          featurename="ASVs",
          basename=bn,
          directory=".",
          analyzes=params$analyzes,
          graphics.type=params$graphics.type,
          vertical=params$vertical,
          collect=params$collect,
          palette.family=params$palette.family,
          palette=params$palette,
          tag_level=params$tag_level,
          width=params$width,
          height=params$height,
          statistics=params$statistics,
          permu=params$permu,
          significance.threshold=params$significance.threshold,
          find.differing.pairs=params$find.differing.pairs,
          rand.seed=params$rand.seed,
          variable=params$main.variable,
          facetting.variables=params$facetting.variables,
          alpha.div.measures=params$alpha.div.measures,
          alpha.plot.style=params$alpha.plot.style,
          alpha.plot.cols=params$alpha.plot.cols,
          reflevel=params$ref.group,
          alpha.subplots.annotation=params$alpha.subplots.annotation,
          points=params$alpha.div.points,
          beta.diversity.analyzes=params$beta.diversity.analyzes,
          beta.diversity.metrics=params$beta.diversity.metrics,
          unifrac.types=params$unifrac.types,
          sample.color.variable = params$sample.color.variables, #
          sample.size.variable = params$sample.size.variable, # needs to be a continuous one
          ellipses = params$ellipses,
          biplot = params$biplot,
          feature.color.variable = params$feature.color.variable,
          size.means.abundance = params$size.means.abundance,
          abundance.cutoff = params$abundance.cutoff,
          triplot = params$triplot, #' plot vectors representing chemical variables significantly influencing community structure? For this to take effect 'rda', 'dbrda' or 'cca' must be in beta.diversity.analyzes
          chemistry.variables = params$chemical.variables,
          vector.col = params$vector.col,
          ordistep.steps=params$ordistep.steps,
          ordistep.direction=params$ordistep.direction,
          varpart.variable = params$varpart.variable,
          ### differentially abundant features identification-specific parameters
          log2FC.threshold = params$log2FC.threshold,
          mergeby = params$mergeby,
          diff.abund.methods = params$diff.features.methods,
          consensus = params$consensus, # include only those features which were found to be DR by all used methods in the final results
          plot.diff.features = params$plot.diff.features, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
          num.top.features = params$num.top.features,
          write.single = params$write.single,
          write.final = params$write.final, # write final, combined table to disc
          mc.samples = params$mc.samples, # how many Monte Carlo samples to use in ALDEx2
          ### taxonomy analysis-specific parameters
          taxlevel.names = params$taxlevel.names,
          rare.threshold = params$rare.threshold
        )
        setwd(outdir)
      }else if( feature == 'OTUs' ){
        message("analyzeAll: analyzing OTUs")
        if(is.null(experiment$OTUs)){
          message("analyzeAll: generating OTUs")
          experiment$OTUs <- mgmisc::constructOTUs(
            fasta.file = experiment$ASVs$fasta.file,
            count.file = experiment$ASVs$count.file,
            diss = params$dissimilarities,
            aligned = experiment$ASVs$aligned,
            reference.alignment = experiment$metadata$reference.alignment,
            start=experiment$metadata$alignment_trimming$start,
            stop=experiment$metadata$alignment_trimming$stop,
            pairwise.alignment = experiment$metadata$pairwise.alignment,
            classification.db = experiment$metadata$classification.db,
            mothur = experiment$metadata$mothur, # path to mothur executable
            processors = params$processors
          )
          message("analyzeAll: OTUs generated")
        }

        d <- paste0("OTUs_", gsub(":", "-", gsub("\\ ", "_", as.character(date()))))
        dir.create(d)
        message("analyzeAll: results will be stored in the ", d)
        setwd(d)
        bn <- paste0(basename, "_OTUs")
        tmpres <- analyzeFeature(
          experiment=tmpres,
          featurename="OTUs",
          d=params$dissimilarities,
          basename=bn,
          directory=".",
          analyzes=params$analyzes,
          graphics.type=params$graphics.type,
          vertical=params$vertical,
          collect=params$collect,
          palette.family=params$palette.family,
          palette=params$palette,
          tag_level=params$tag_level,
          width=params$width,
          height=params$height,
          statistics=params$statistics,
          permu=params$permu,
          significance.threshold=params$significance.threshold,
          find.differing.pairs=params$find.differing.pairs,
          rand.seed=params$rand.seed,
          variable=params$main.variable,
          facetting.variables=params$facetting.variables,
          alpha.div.measures=params$alpha.div.measures,
          alpha.plot.style=params$alpha.plot.style,
          alpha.plot.cols=params$alpha.plot.cols,
          reflevel=params$ref.group,
          alpha.subplots.annotation=params$alpha.subplots.annotation,
          points=params$alpha.div.points,
          beta.diversity.analyzes=params$beta.diversity.analyzes,
          beta.diversity.metrics=params$beta.diversity.metrics,
          unifrac.types=params$unifrac.types,
          sample.color.variable = params$sample.color.variables, #
          sample.size.variable = params$sample.size.variable, # needs to be a continuous one
          ellipses = params$ellipses,
          biplot = params$biplot,
          feature.color.variable = params$feature.color.variable,
          size.means.abundance = params$size.means.abundance,
          abundance.cutoff = params$abundance.cutoff,
          triplot = params$triplot, #' plot vectors representing chemical variables significantly influencing community structure? For this to take effect 'rda', 'dbrda' or 'cca' must be in beta.diversity.analyzes
          chemistry.variables = params$chemical.variables,
          vector.col = params$vector.col,
          ordistep.steps=params$ordistep.steps,
          ordistep.direction=params$ordistep.direction,
          varpart.variable = params$varpart.variable,
          ### differentially abundant features identification-specific parameters
          log2FC.threshold = params$log2FC.threshold,
          mergeby = params$mergeby,
          diff.abund.methods = params$diff.features.methods,
          consensus = params$consensus, # include only those features which were found to be DR by all used methods in the final results
          plot.diff.features = params$plot.diff.features, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
          num.top.features = params$num.top.features,
          write.single = params$write.single,
          write.final = params$write.final, # write final, combined table to disc
          mc.samples = params$mc.samples, # how many Monte Carlo samples to use in ALDEx2
          ### taxonomy analysis-specific parameters
          taxlevel.names = params$taxlevel.names,
          rare.threshold = params$rare.threshold
        )
        setwd(outdir)
      }else if( feature == 'KOs'){
        message("analyzeAll: analyzing KOs")
        d <- paste0("KOs_", gsub(":", "-", gsub("\\ ", "_", as.character(date()))))
        dir.create(d)
        message("analyzeAll: results will be stored in the ", d)
        setwd(d)
        bn <- paste0(basename, "_KOs")
        bdm <- beta.diversity.metrics[ -grep("GUniFrac", beta.diversity.metrics) ]
        tmpres <- analyzeFeature(
          experiment=tmpres,
          featurename="KOs",
          d=NULL,
          basename=bn,
          directory=".",
          analyzes=params$analyzes,
          graphics.type=params$graphics.type,
          vertical=params$vertical,
          collect=params$collect,
          palette.family=params$palette.family,
          palette=params$palette,
          tag_level=params$tag_level,
          width=params$width,
          height=params$height,
          statistics=params$statistics,
          permu=params$permu,
          significance.threshold=params$significance.threshold,
          find.differing.pairs=params$find.differing.pairs,
          rand.seed=params$rand.seed,
          variable=params$main.variable,
          facetting.variables=params$facetting.variables,
          alpha.div.measures=params$alpha.div.measures,
          alpha.plot.style=params$alpha.plot.style,
          alpha.plot.cols=params$alpha.plot.cols,
          reflevel=params$ref.group,
          alpha.subplots.annotation=params$alpha.subplots.annotation,
          points=params$alpha.div.points,
          beta.diversity.analyzes=params$beta.diversity.analyzes,
          beta.diversity.metrics=bdm,
          unifrac.types=NULL,
          sample.color.variable = params$sample.color.variables, #
          sample.size.variable = params$sample.size.variable, # needs to be a continuous one
          ellipses = params$ellipses,
          biplot = params$biplot,
          feature.color.variable = params$feature.color.variable,
          size.means.abundance = params$size.means.abundance,
          abundance.cutoff = params$abundance.cutoff,
          triplot = params$triplot, #' plot vectors representing chemical variables significantly influencing community structure? For this to take effect 'rda', 'dbrda' or 'cca' must be in beta.diversity.analyzes
          chemistry.variables = params$chemical.variables,
          vector.col = params$vector.col,
          ordistep.steps=params$ordistep.steps,
          ordistep.direction=params$ordistep.direction,
          varpart.variable = params$varpart.variable,
          ### differentially abundant features identification-specific parameters
          log2FC.threshold = params$log2FC.threshold,
          mergeby = params$mergeby,
          diff.abund.methods = params$diff.features.methods,
          consensus = params$consensus, # include only those features which were found to be DR by all used methods in the final results
          plot.diff.features = params$plot.diff.features, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
          num.top.features = params$num.top.features,
          write.single = params$write.single,
          write.final = params$write.final, # write final, combined table to disc
          mc.samples = params$mc.samples, # how many Monte Carlo samples to use in ALDEx2
          ### taxonomy analysis-specific parameters
          taxlevel.names = params$taxlevel.names,
          rare.threshold = params$rare.threshold
        )
        setwd(outdir)
      }else if( feature == 'pathways'){
        message("analyzeAll: analyzing pathways")
        d <- paste0("pathways_", gsub(":", "-", gsub("\\ ", "_", as.character(date()))))
        dir.create(d)
        message("analyzeAll: results will be stored in ", d)
        setwd(d)
        bn <- paste0(basename, "_pathways")
        bdm <- beta.diversity.metrics[ -grep("GUniFrac", beta.diversity.metrics) ]
        tmpres <- analyzeFeature(
          experiment=tmpres,
          featurename="pathways",
          d=ds,
          basename=bn,
          directory=".",
          analyzes=analyzes,
          graphics.type=graphics.type,
          vertical=vertical,
          collect=collect,
          palette.family=palette.family,
          palette=palette,
          tag_level=tag_level,
          width=width,
          height=height,
          statistics=statistics,
          permu=permu,
          significance.threshold=significance.threshold,
          find.differing.pairs=find.differing.pairs,
          rand.seed=rand.seed,
          variable=main.variable,
          facetting.variables=facetting.variables,
          alpha.div.measures=alpha.div.measures,
          alpha.plot.style=alpha.plot.style,
          alpha.plot.cols=alpha.plot.cols,
          reflevel=ref.group,
          alpha.subplots.annotation=alpha.subplots.annotation,
          points=alpha.div.points,
          beta.diversity.analyzes=beta.diversity.analyzes,
          beta.diversity.metrics=bdm,
          unifrac.types=NULL,
          sample.color.variable = sample.color.variables, #
          sample.size.variable = sample.size.variable, # needs to be a continuous one
          ellipses = ellipses,
          biplot = biplot,
          feature.color.variable = feature.color.variable,
          size.means.abundance = size.means.abundance,
          abundance.cutoff = abundance.cutoff,
          triplot = triplot, #' plot vectors representing chemical variables significantly influencing community structure? For this to take effect 'rda', 'dbrda' or 'cca' must be in beta.diversity.analyzes
          chemistry.variables = chemical.variables,
          vector.col = vector.col,
          ordistep.steps=ordistep.steps,
          ordistep.direction=ordistep.direction,
          varpart.variable = varpart.variable,
          ### differentially abundant features identification-specific parameters
          log2FC.threshold = log2FC.threshold,
          mergeby = mergeby,
          diff.abund.methods = diff.features.methods,
          consensus = consensus, # include only those features which were found to be DR by all used methods in the final results
          plot.diff.features = plot.diff.features, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
          num.top.features = num.top.features,
          write.single = write.single,
          write.final = write.final, # write final, combined table to disc
          mc.samples = mc.samples, # how many Monte Carlo samples to use in ALDEx2
          ### taxonomy analysis-specific parameters
          taxlevel.names = taxlevel.names,
          rare.threshold = rare.threshold
        )
        setwd(outdir)
      }
      setwd(outdir)
    }
  }

  return(experiment)
}


set_parameters <- function(
    p=parameters, # a path to a tsv file with two columns: par_names and values or a named list where names are parameter names and elements are values
    a=argv){  # list of arguments of a calling function

  defaults <- list(
    analyze = c('ASVs', 'OTUs', 'KOs', 'pathways'),
    dissimilarities = c(0.03),
    analyzes = c('alpha-diversity', 'beta-diversity', 'taxonomy', 'diff-features'),
    basename = 'bacteria',
    outdir = '.', # where results are to be stored
    graphics.type = 'ggplot', # if base ('base') graphics or ggplot2 ('ggplot') is to be used
    vertical = TRUE,
    collect = FALSE,
    tag_level = 'A',
    statistics = TRUE,
    rand.seed = 667,
    permu = 999,
    significance.threshold = 0.05,
    find.differing.pairs = TRUE,
    palette.family = 'ggsci', # palette family (package), must be one of those gathered in paletteer (https://github.com/EmilHvitfeldt/paletteer?tab=readme-ov-file), if present in the experiment object will be overwritten
    palette = 'nrc_npg', # palette name, must come from the package in 'palette.family', if present in the experiment object will be overwritten
    processors = 1,
    width = 3.5,
    height = 3.5,
    # variables
    main.variable = NULL,
    facetting.variables = NULL,
    varpart.variable = NULL,
    sample.color.variable = NULL,
    sample.size.variable = NULL,
    feature.color.variable = NULL,
    chemistry.variables = NULL,
    # alpha diversity
    alpha.div.measures = c('sobs', 'H', 'E'),
    alpha.div.points = T,
    ref.group = NULL,
    alpha.plot.style = 'box',
    alpha.plot.cols = 2,
    # beta diversity
    beta.diversity.metrics = c('bray', 'GUniFrac'), # one or more of metrics implemented in vegan::vegdist or GUniFrac
    unifrac.types = c('d_0.5'), # GUniFrac metrics (for ASVs and OTUs only)
    beta.diversity.analyzes = c('nmds', 'dbrda'),
    biplot=T,
    triplot=F,
    size.means.abundance=T,
    abundance.cutoff=50,
    ellipses = T,
    vector.col = 'red',
    ordistep.steps=500,
    ordistep.direction='both',
    # differentially abundant features identification parameters
    diff.features.methods=c('deseq'), # one of 'deseq', 'aldex',
    log2FC.threshold = 1,
    mergeby = 0,
    consensus = F, # include only those features which were found to be DR by all used methods in the final results
    plot.diff.features = F, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
    num.top.features = 10,
    write.single = T,
    write.final = T, # write final, combined table to disc
    mc.samples = 128, # how many Monte Carlo samples to use in ALDEx2
    # taxonomy plots parameters
    rare.threshold = 0.01,
    taxlevel.names = c('kingdoms', 'phyla', 'classes', 'orders', 'families', 'genera')
  )
  if(is.character(p)){
    pdf <- read.table(p, header=T, sep="\t", dec=".")
    if(!identical(colnames(pdf), c("par_names", "values"))){
      stop(paste0("set_parameters: Incorrectly formatted paramters file ", p, ". Columns should be named 'par_names' and 'values'\n"))
    }
    p <- list()
    for(n in 1:nrow(pdf)){
      p[[n]] <- pdf[n,2]
    }
    names(p) <- pdf[,1]
  }
  for(n in names(p)){
    if(n %in% names(defaults)){
      if(!is.identical(defaults[[n]], p[[n]])){
        defaults[[n]] <- p[[n]]
      }
    }else{
      warning("set_parameters: unknown parameter", n, " in parameters file/object\n")
    }
  }
  a <- a[-grep("experiment", names(a))]
  a <- a[-grep("parameters", names(a))]
  for(n in names(a)){
    defaults[[n]] <- a[[n]]
  }
  return(defaults)
}
