#' @export plotAlphaDiversity

plotAlphaDiversity <- function(
    asvtab,
    data,
    variable,
    reflevel=NULL,
    annotation='A',
    measures=c('H', 'E', 'sobs'),
    basename,
    cols=2,
    width=3.5,
    height=3.5,
    palette='default',
    comps=NULL) { # comparisons list, if NULL all possible comparisons are analyzed
	stopifnot(variable %in% colnames(data), identical(rownames(asvtab), rownames(data)), levels(data[[variable]]) > 1 )
message("Alpha-diversity analysis")
	H <- NULL
	E <- NULL
	simpson <- NULL
	invsimpson <- NULL
	sobs <- NULL
	chao <- NULL
	ace <- NULL

	bdiv <- data.frame( row.names=rownames(asvtab), variable=data[[variable]]);
	for( meas in measures ) {
		if( meas == 'H' ){
			bdiv$H <- vegan::diversity(asvtab)
		} else if( meas == 'simpson' ) {
			bdiv$simpson <- vegan::diversity(asvtab, index='simpson')
		} else if( meas == 'invsimpson' ) {
			bdiv$invsimpson <- vegan::diversity(asvtab, index='invsimpson')
		} else if( meas == 'E' ) {
			bdiv$E <- vegan::diversity(asvtab)/log(vegan::specnumber(asvtab))
		} else if( meas == 'sobs' ) {
			bdiv$sobs <- vegan::specnumber(asvtab)
		} else if( meas == 'chao' ) {
			bdiv$chao <- t(vegan::estimateR(asvtab))[,2]
		} else if( meas == 'ace' ) {
			bdiv$ace <- t(vegan::estimateR(asvtab))[,4]
		}
	}


	xlabel <- gsub("_", " ", variable)



	for( meas in measures ){
		if( length(levels(as.factor(bdiv$variable))) == 2 ) {
			testtype = 'wilcox.test'
			general_p.value = NA
		} else {
			testtype = 'kruskal.test'
			general_p.value <- kruskal.test(bdiv[[meas]], bdiv$variable)$p.value
			print( paste(meas, general_p.value ) )
		}

		if( meas == 'H' ){
			title = "Diversity (Shannon's H')";
			ylab = "Shannon's H'";
			maxy = ceiling(max(bdiv[[meas]]) + 0.5 * max(bdiv[[meas]]))
			p_shannon <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='black', add='jitter', add.params=list(size=0.5), width=0.5, size=0.1, ggtheme=theme_pubr(base_size=5), font.label=list(size=5) ) + ylim(0,maxy) + stat_compare_means(method=testtype, size=1.5, label.y=0.95*maxy) + ylab(ylab) + xlab(xlabel) + ggtitle(title)
			if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_shannon <- p_shannon + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=1, comparisons=comps)
			}
		} else if( meas == 'simpson' ) {
			title = "Diversity (Simpson)"
			ylab = "Simpson's index"
			maxy = ceiling(max(bdiv[[meas]]) + 0.5 * max(bdiv[[meas]]))
			p_simpson <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='black', add="jitter", add.params=list(size=0.5), width=0.5, size=0.1, ggtheme=theme_pubr(base_size=5), font.label=list(size=5) ) + ylim(0,maxy) + stat_compare_means(method=testtype, size=1.5, label.y=0.95*maxy) + ylab(ylab) + xlab(xlabel) + ggtitle(title)
			if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_simpson <- p_simpson + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=1, comparisons=comps)
			}
		} else if(meas == 'invsimpson') {
			title = "Diversity (inv. Simpson)"
			ylab = "Inv. Simpson's index"
			maxy = ceiling(max(bdiv[[meas]]) + 0.5 * max(bdiv[[meas]]))
			p_invsimpson <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='black', add="jitter", add.params=list(size=0.5), width=0.1, size=0.5, ggtheme=theme_pubr(base_size=5), font.label=list(size=5) ) + ylim(0,maxy) + stat_compare_means(method=testtype, size=1.5, label.y=0.95*maxy) + ylab(ylab) + xlab(xlabel) + ggtitle(title)
			if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_invsimpson <- p_invsimpson + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=1, comparisons=comps)
			}

		} else if( meas == 'E' ){
			title = "Evenness (Shannon's E)";
			ylab = "Shannon's E";
			maxy = max(bdiv[[meas]]) + 0.5 * max(bdiv[[meas]])
			p_e <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='black', add="jitter", add.params=list(size=0.5), width=0.5, size=0.1, ggtheme=theme_pubr(base_size=5), font.label=list(size=5) ) + ylim(0,maxy) + stat_compare_means(method=testtype, size=1.5, label.y=0.95*maxy) + ylab(ylab) + xlab(xlabel) + ggtitle(title)
			if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_e <- p_e + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=1, comparisons=comps)
			}

		} else if( meas == 'sobs'){
			title = "Species richness";
			ylab = "Obs. no. ASVs";
			maxy = ceiling(max(bdiv[[meas]]) + 0.5 * max(bdiv[[meas]]))
			p_sobs <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='black', add="jitter", add.params=list(size=0.5), width=0.5, size=0.1, ggtheme=theme_pubr(base_size=5), font.label=list(size=5) ) + ylim(0,maxy) + stat_compare_means(method=testtype, size=1.5, label.y=0.95*maxy) + ylab(ylab) + xlab(xlabel) + ggtitle(title)
			if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_sobs <- p_sobs + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=1, comparisons=comps)
			}

		} else if( meas == 'chao' ) {
			title = "Chao1"
			ylab = "Est. no. ASVs"
			maxy = ceiling(max(bdiv[[meas]]) + 0.5 * max(bdiv[[meas]]))
			p_chao <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='black', add="jitter", add.params=list(size=0.5), width=0.5, size=0.1, ggtheme=theme_pubr(base_size=5), font.label=list(size=5) ) + ylim(0,maxy) + stat_compare_means(method=testtype, size=1.5, label.y=0.95*maxy) + ylab(ylab) + xlab(xlabel) + ggtitle(title)
			if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_chao <- p_chao + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=1, comparisons=comps)
			}

		} else if( meas == 'ace' ) {
			title = "ACE"
			ylab = "Est. no. ASVs"
			maxy = ceiling(max(bdiv[[meas]]) + 0.5 * max(bdiv[[meas]]))
			p_ace <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='black', add="jitter", add.params=list(size=0.5), width=0.5, size=0.1, ggtheme=theme_pubr(base_size=5), font.label=list(size=5) ) + ylim(0,maxy) + stat_compare_means(method=testtype, size=1.5, label.y=0.95*maxy) + ylab(ylab) + xlab(xlabel) + ggtitle(title)
			if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_ace <- p_ace + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=1, comparisons=comps)
			}

		}
	}
	file <- paste0(basename, ".", variable, ".alpha_diversity.svg")
	final_plot <- (p_sobs + p_e)/(p_shannon + plot_spacer()) + patchwork::plot_annotation(tag_levels=annotation)

	svg(file, width=width, height=height, pointsize=3);
	print(final_plot);
	dev.off();
  return(bdiv)
}

