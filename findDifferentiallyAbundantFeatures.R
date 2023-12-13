#' @export findDifferentiallyAbundantFeatures

findDifferentiallyAbundantFeatures <- function(
    formula,
    asvtab,
    data,
    nproc = 1,
    basename = "difffeatures.",
    annot=NULL,
    mergeby=0,
    qthreshold=0.05,
    log2FCthreshold=1,
    write.single=F,
    write.final=F ){

  var1 <- all.vars(as.formula(formula))[1];
  stopifnot( length(all.vars(as.formula(formula))) == 1, all.vars(as.formula(formula)) %in% colnames(data), length(levels(as.factor(data[[var1]]))) >= 2, identical(rownames(asvtab), rownames(data)) );

  message(paste0("Finding differentially abundant features"))

  asvtabnzv <- mixOmics::nearZeroVar(asvtab);
  if( nrow(asvtabnzv$Metrics) > 0){
    asvtablv <- asvtab[ , -asvtabnzv$Position ];
  } else {
    asvtablv <- asvtab
  }

  message("Large variance features selected");
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=as.matrix(t(asvtablv)+1), colData=data, design=as.formula(formula) );
  message("DDS object created");
  contr <- combn( levels(as.factor(data[[all.vars(formula)[1]]])), 2 );
  filenameaccu <- paste0(basename, "_", var1, ".diff.csv");
  print( filenameaccu );
  #
  if( nproc > 1 ){
    dds <- DESeq2::DESeq(dds, parallel=T, BPPARAM=BiocParallel::MulticoreParam(nproc), fitType='mean');
    for( i in 1:ncol(contr) ){
      con <- contr[,i];
      constr <- paste0(con[1], "_", con[2]);
      message(paste0("Analyzing ", constr))
      filename <- paste0(basename, ".", var1, ".", constr, ".diff.csv");
      if( !is.null(annot) ){
        row <- rep( paste0(con[1], "_", con[2]), ncol(annot)+7);
      } else {
        row <- rep( paste0(con[1], "_", con[2]), 7);
      }
      con <- c(var1, con);
      ddsres <- DESeq2::results(dds, contrast=con, parallel=T, BPPARAM=BiocParallel::MulticoreParam(nproc));
      ddsres.df <- as.data.frame(ddsres)[ !is.na(ddsres$padj) & ddsres$padj < qthreshold & abs(ddsres$log2FoldChange) > log2FCthreshold, ];
      ddsres.df <- ddsres.df[ order(ddsres.df$log2FoldChange), ];
      if( !is.null(annot) ){
        ddsres.df <- merge(ddsres.df, annot, by.x=0, by.y=mergeby, all.x=T);
      }
      if( i == 1 ){
        accuddsres.df <- ddsres.df;
        cnames <- colnames(ddsres.df)
        names(row) <- cnames
        accuddsres.df <- rbind(row, accuddsres.df);
        colnames(accuddsres.df) <- cnames
      } else {
        names(row) <- colnames(accuddsres.df)
        accuddsres.df <- rbind(accuddsres.df, row);
        colnames(accuddsres.df) <- cnames
        accuddsres.df <- rbind(accuddsres.df, ddsres.df);
        colnames(accuddsres.df) <- cnames
      }
      if( write.single ){
        write.table(ddsres.df, filename, sep="\t", col.names=NA);
      }
    }
  } else {
    dds <- DESeq2::DESeq(dds);
    for( i in 1:ncol(contr) ){
      con <- contr[,i];
      constr <- paste0(con[1], "_", con[2]);
      filename <- paste0(basename, ".", var1, ".", constr, ".diff.csv");
      if( !is.null(annot) ){
        row <- rep( paste0(con[1], "_", con[2]), ncol(annot)+7);
      } else {
        row <- rep( paste0(con[1], "_", con[2]), 7);
      }
      con <- c(var1, con);
      ddsres <- DESeq2::results(dds, contrast=con);
      ddsres.df <- as.data.frame(ddsres)[ !is.na(ddsres$padj) & ddsres$padj < qthreshold & abs(ddsres$log2FoldChange) > log2FCthreshold, ];
      ddsres.df <- ddsres.df[ order(ddsres.df$log2FoldChange), ];
      if( !is.null(annot) ){
        ddsres.df <- merge(ddsres.df, annot, by.x=0, by.y=mergeby, all.x=T);
      }

      if( i == 1 ){
        accuddsres.df <- ddsres.df;
        cnames <- colnames(ddsres.df)
        names(row) <- cnames
        accuddsres.df <- rbind(row, accuddsres.df);
        colnames(accuddsres.df) <- cnames
      } else {
        names(row) <- colnames(accuddsres.df)
        accuddsres.df <- rbind(accuddsres.df, row);
        colnames(accuddsres.df) <- cnames
        accuddsres.df <- rbind(accuddsres.df, ddsres.df);
        colnames(accuddsres.df) <- cnames
      }
 #     write.table(ddsres.df, filename, sep="\t", col.names=NA);
    }
  }
  if( write.final ){
    write.table(accuddsres.df, filenameaccu, sep="\t", col.names=NA);
  }
  return(accuddsres.df);
}
