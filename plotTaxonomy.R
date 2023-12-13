#' @export plotTaxonomy

plotTaxonomy <- function(
    asvtab,
    data,  # sample data table for each of the rows in asvtab
    taxonomy=NULL, # taxonomy table for each of the columns in asvtab
    variable, # must be one of 'data' columns, if not a factor, will be converted to one
    reflevel=NULL, #
    rare.threshold=0.01, # minimal abundance (across all data) of a taxon not falling into the 'rare' category
    basename="bacteria", # a string, base for forming file names, taxlevel and variable names will be added after underscore, e.g. bacteria_phyla_variable.svg
    taxlevels=c('phyla', 'classes', 'orders', 'families', 'genera'), # vector of taxonomic level names in plural, if not in a classical set, final 's' will be stripped to obtain singular form
    palette='default',
    width=3.5,
    height=3.5,
    ggplot=F, # whether use ggplot2 or base graphic (default)
    pointsize=3) {
stopifnot( variable %in% colnames(data), identical(rownames(asvtab), rownames(data)), levels(data[[variable]]) > 1 )
tax <- taxonomy[ rownames(taxonomy) %in% colnames(asvtab), ]

message("Plotting taxonomy")

outlist <- list()
	for( taxlevelpl in taxlevels ){
		if(taxlevelpl == 'kingdoms'){
		  taxlevelsngl == "Kingdom"
		  } else if( taxlevelpl == "phyla" ){
			taxlevelsngl <- "Phylum";
			} else if( taxlevelpl == "classes" ){
			taxlevelsngl <- "Class";
			} else if( taxlevelpl == "orders" ){
			taxlevelsngl <- "Order";
			} else if( taxlevelpl == "families" ){
			taxlevelsngl <- "Family";
			} else if( taxlevelpl == 'genera'){
			taxlevelsngl <- "Genus";
			} else if( taxlevelpl == 'species'){
			  taxlevelsngl <- "Species"
			} else {
			  taxlevelsngl <- sub("s$", "", taxlevelpl)
			}
		taxtmp <- as.data.frame(t(asvtab))
		taxtmp$taxlvl <- tax[[taxlevelsngl]];
		taxtmp <-  aggregate(. ~ taxlvl, data=taxtmp, FUN='sum');
		rownames(taxtmp) <- taxtmp$taxlvl;
		taxtmp$taxlvl <- NULL;
		taxtmp <- as.data.frame(t(taxtmp));
		outlist[[taxlevelpl]] <- taxtmp; # count matrix is returned (taxa as columns)
		message(paste0(taxlevelpl))
		taxtmp <- taxtmp/rowSums(taxtmp) # normalization
		taxtmp <- aggregate( . ~ data[[variable]], data=taxtmp, FUN='mean'); # mean over levels of the variable
		rownames(taxtmp) <- paste0(taxtmp$"data[[variable]]");
		taxtmp$"data[[variable]]" <- NULL;
		taxtmp.t <- as.data.frame(t(taxtmp));
		taxtmp.t <- taxtmp.t[ order(rowSums(taxtmp.t), decreasing=T), ];
		taxtmp <- as.data.frame(t(taxtmp.t));
		taxtmp_unclassified <- as.data.frame(taxtmp[ , grep("unknown|unclassified|_ge|_fa|_or|_cl|_ph", colnames(taxtmp)) ]);
		taxtmp_classified <- as.data.frame(taxtmp[ , grep("unknown|unclassified|_ge|_fa|_or|_cl|_ph", colnames(taxtmp), invert=T) ]);
		taxtmp_classified.nonrare <- as.data.frame(taxtmp_classified[, colSums(taxtmp_classified)/sum(taxtmp) > rare.threshold, drop=F ]);
		taxtmp_classified.rare <- as.data.frame(taxtmp_classified[, colSums(taxtmp_classified)/sum(taxtmp) <= rare.threshold]);
		taxtmp_classified.nonrare$rare <- rowSums(taxtmp_classified.rare);
		taxtmp_classified.nonrare$unclassifed <- rowSums(taxtmp_unclassified);
		taxtmp <- taxtmp_classified.nonrare;
		taxtmp_percent <- 100 * taxtmp
    if(!ggplot){
	    file <- paste0( basename, ".", taxlevelpl, "_", variable, ".svg" );
		  palette( RColorBrewer::brewer.pal(10, "Set3") );
		  svg( file, width=width, height=height, pointsize=pointsize );
		  par( mar=c(6,4,2,12)+0.1, mfrow=c(1,1), las=2 );
		  with( taxtmp_percent, barplot( as.matrix(t(taxtmp_percent)), legend=colnames(taxtmp_percent), args.legend=list(x="right", inset=c(-0.06,0), xpd=T), col=palette()) );
		  dev.off();
    } else {
      taxtmp_percent_molten <- reshape2::melt(taxtmp_percent, id.vars=c(variable), value="Relative abundance", variable.name=taxlevelsngl)
      file <- paste0( basename, ".", taxlevelpl, "_", variable, ".svg" );
      palette( RColorBrewer::brewer.pal(10, "Set3") );
    }
	}
names(outlist) <- taxlevels
return(outlist)

}

