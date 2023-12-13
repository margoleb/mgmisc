#' @export plotBetaDiversity
plotBetaDiversity <- function(
    asvtab,
    data,
    variable,
    color.variable=NULL,
    varpart.variable=NULL,
    basename,
    analyses=c('nmds', 'dbrda'),
    dists=c('bray', 'horn'),
    tree=NULL,
    unifrac_types=c('d_0.5'),
    statistics=T,
    rand.seed=667,
    palette='default',
    width=3.5,
    height=3.5,
    ggplot=F ) {
	stopifnot( identical(rownames(asvtab), rownames(data)), length(levels(as.factor(data[[variable]]))) >= 2, analyses %in% c('nmds', 'rda', 'dbrda') );
	unifracs <- NULL
	disttabs <- list()

	message("Beta-diversity analysis")
  set.seed(rand.seed)
	for( d in dists ) {
		if( d %in% c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower","morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", "robust.aitchison")) {
			disttab <- vegdist( asvtab, method=d )
			disttabs[[d]] <- disttab
		} else if( grepl("UniFrac", d) ) {
			if( is.null(unifracs) ) {
				alpha <- as.numeric(sub("d_", "", unifrac_types))
				alpha <- alpha[ !is.na(alpha) ]
				unifracs <- GUniFrac(asvtab, tree, alpha=alpha)$unifracs
			}
			for( a in  unifrac_types) {
				disttab <- as.dist(unifracs[ , , a ])
				name <- paste0(a, "GUniFrac")
				disttabs[[name]] <- disttab
			}
		} else {
			stop("Unknown distance")
		}
	}

	message("Distance matrices generated")

	for( analysis in analyses ) {
			if( analysis == 'nmds' ) {
				for( d in names(disttabs) ){
				  message(paste0(d, " NMDS analysis"))
					bname <- paste0(basename, ".", d, "_", variable )
					plotNMDS( disttabs[[d]], data, variable=variable, color.variable=color.variable, varpart.variable=varpart.variable, basename=bname, statistics=statistics, palette=palette, width=width, height=height )
				}
			} else if( analysis == 'dbrda' ) {
				for( d in names(disttabs) ){
				  message(paste0(d, " dbRDA analysis"))
					bname <- paste0(basename, ".", d, "_", variable )
					plotdbRDA( disttabs[[d]], data, variable, color.variable, varpart.variable, basename=bname, statistics=statistics, palette=palette, width=width, height=height )
				}
			} else if( analysis == 'rda' ){
			  message("RDA analysis")
				bname <- paste0(basename, "_", variable )
				plotRDA( asvtab, data, variable, color.variable, varpart.variable, basename=bname, statistics=statistics, palette=palette, width=width, height=height )
			}
		}
}


