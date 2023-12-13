#' @export plotNMDS

plotNMDS <- function( d, data, variable, color.variable, varpart.variable, basename, statistics=T, palette='default', permu=999, width=3.5, height=3.5 ){

  message("Plotting NMDS")

  	explvar <- NULL
	permanovaf <- NULL
	ord <- vegan::metaMDS( d, k=2, try=100, trymax=200 )
	title <- ""
	if( !is.null(varpart.variable) ){
	  message("Performing variance partitioning")
		f1 <- as.formula(paste0(" ~ ", variable))
		f2 <- as.formula(paste0(" ~ ", varpart.variable))
		ordvarpart <- vegan::varpart(d, f1, f2, data=data)
		explvar <- paste0(sprintf("%.2f", 100 * (ordvarpart$part$indfract[1,3])), '%');
		title <- paste0("\nVariance explained: ", explvar );
	}

	if( as.logical(statistics) ) {
	  message("Performing PERMANOVA and betadisper analyses")
		ordbetadisper <- vegan::betadisper( d, data[[variable]], bias.adjust=T)
		message("Betadisper done")
		ordbetadisperanova <- anova( ordbetadisper, permu=permu )
		message("Betadisper tested")
		betadisperf <- sprintf( "%.2f", ordbetadisperanova$"F value"[1] );
		betadisperp <- sprintf( "%.2e", ordbetadisperanova$"Pr(>F)"[1] );
		f <- as.formula(paste0("d ~ ", variable))
		ordpermanova <- vegan::adonis2(f, data=data, permu=999)
		message("PERMANOVA done")
		permanovaf <- sprintf( "%.2f", ordpermanova$F[1] );
		permanovap <- sprintf( "%.3f", ordpermanova$"Pr(>F)"[1] );
		title <- paste0("ANOVA F = ", permanovaf, ", p = ", permanovap, "\nbetadisper F = ", betadisperf, ", p = ", betadisperp, title)
	}
	file <- paste0(basename, ".nmds.svg")
	svg( file, width=width, height=height, pointsize=6 );
	par(mar=c(4,4,4,2)+0.1, mfrow=c(1,1), las=2);
	vegan::ordiplot(ord, type='none', main=title);
	points( ord, pch=21+as.numeric(as.factor(data[[color.variable]])), col='black', bg=as.numeric(as.factor(data[[variable]])) );
	vegan::ordiellipse( ord, groups=data[[variable]], label=T);
	legend('topright', legend=levels(as.factor(data[[color.variable]])), title=color.variable, pch=21+sort(unique(as.numeric(as.factor(data[[color.variable]])))));
	legend('topleft', legend=levels(as.factor(data[[variable]])), title=variable, pch=21, pt.bg=sort(unique(as.numeric(as.factor(data[[variable]])))));
	dev.off();

}

