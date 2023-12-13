#' @export plotRDA

plotRDA <- function(asvtable, data, variable, color.variable, varpart.variable, basename, statistics=T, palette='default', permu=999, width=3.5, height=3.5){
	explvar <- NULL
	permanovaf <- NULL
	f <- as.formula(paste0("d ~ ", variable))
	ord <- vegan::rda( f, data=data )
	title <- ""
	if( is.null(varpart.variable) ){
	  message("Performing variance partitioning")
		f1 <- as.formula(paste0(" ~ ", variable))
		f2 <- as.formula(paste0(" ~ ", varpart.variable))
		ordvarpart <- vegan::varpart(d, f1, f2, data=data)
		explvar <- paste0(sprintf("%.2f", 100 * (ordvarpart$part$indfract[1,3]), '%'));
		title <- paste0("\nVariance explained: ", explvar );
	}

	if( as.logical(statistics) ) {
	  message("Testing model and performing betadisper analysis")
		ordbetadisper <- vegan::betadisper( d, data[[variable]], bias.adjust=T)
		ordbetadisperanova <- anova( ordbetadisper, permu=permu )
		betadisperf <- sprintf( "%.2f", ordbetadisperanova$"F value"[1] );
		betadisperp <- sprintf( "%.2e", ordbetadisperanova$"Pr(>F)"[1] );
		ordrdaanova <- anova(ord, permu=permu)
		rdaanovaf <- sprintf( "%.2f",ordrdaanova$F[1] );
		rdaanovap <- sprintf( "%.3f", ordrdaanova$"Pr(>F)"[1] );
		title <- paste0("ANOVA F = ", rdaanovaf, ", p = ", rdaanovap, "\nbetadisper F = ", betadisperf, ", p = ", betadisperp, title)
	}

	file <- paste0(basename, ".rda.svg")
	svg( file, width=width, height=height, pointsize=6 );
	par(mar=c(4,4,4,2)+0.1, mfrow=c(1,1), las=2);
	vegan::ordiplot(ord, type='none', main=title);
	points( ord, pch=21+as.numeric(as.factor(data[[color.variable]])), col='black', bg=as.numeric(as.factor(data[[variable]])) );
	vegan::ordiellipse( ord, groups=data[[variable]], label=T);
	legend('topright', legend=levels(as.factor(data[[color.variable]])), title=color.variable, pch=21+sort(unique(as.numeric(as.factor(data[[color.variable]])))));
	legend('topleft', legend=levels(as.factor(data[[variable]])), title=variable, pch=21, pt.bg=sort(unique(as.numeric(as.factor(data[[variable]])))));
	dev.off();

}


