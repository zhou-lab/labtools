convert.betas.txt2rda <- function(txt) {
  betas <- read.table(txt,sep='\t',stringsAsFactor=F,check.names=F);
  basename <- strsplit(txt,'\\.')[[1]][1];
  rda <- paste0(basename,'.rda');
  cat(sprintf('Read %d probes and %d samples.\n', dim(betas)[1], dim(betas)[2]));
  save(betas,file=rda);
  cat(sprintf('Write to file %s.\n', rda));
}

omit.allna.rows <- function(df) {
  ## it is slow to use apply
  df2 <- df[rowSums(is.na(df)) == 0,];
  cat(sprintf('Removed %d rows from %d.\n', dim(df)[1]-dim(df2)[1], dim(df)[1]));
  return(df2);
}

uniformly.methylated <- function(df, thres=0.8) {
  select <- apply(df,1, function(x){all(x>thres)});
  dfm <- df[select & !is.na(select),];
  cat(sprintf('Selected %d uniformly methylated probes.\n', dim(dfm)[1]));
  return(dfm);
}

uniformly.unmethylated <- function(df, thres=0.2) {
  select <- apply(df, 1, function(x){all(x<thres)});
  dfu <- df[select & !is.na(select),];
  cat(sprintf('Selected %d uniformly unmethylated probes.\n', dim(dfm)[1]));
  return(dfu);
}

data.load.samples <- function(cancertypes, probes=NULL) {
  if(length(cancertypes)==0) return;
  data.home <- '/Users/wandingzhou/projects/hs-tcga/data/2015_03_05_TCGA_450/rda/'
  allbetas <- c();
  for(cancertype in cancertypes) {
    load(paste0(data.home,cancertype,'.rda'));
    allbetas <- cbind(allbetas, betas);
  }
  cat(sprintf('Retrieved %d probes and %d samples.\n', dim(allbetas)[1], dim(allbetas)[2]));
  return(allbetas);
}


data.load.blood <- function() {
  load('/Users/wandingzhou/projects/hs-tcga/data/2015_04_10_sorted_cell_population/blood_beta.rda');
  blood.names <- read.table('/Users/wandingzhou/projects/hs-tcga/data/2015_04_10_sorted_cell_population/Sorted_Blood/sample_sheet_IDAT.csv.unix.tsv',sep='\t',stringsAsFactors = F,comment.char='',header=T, row.names='barcode');
  colnames(betas)  <- paste(blood.names[colnames(betas),'Type'],colnames(betas),sep='_');
  return(betas);
}

split.tumor.normal <- function(df) {
  r <- list();
  # as.matrix prevent reduction of matrix to array
  r$tumor <- as.matrix(df[,substring(colnames(df),14,14)=='0']);
  r$normal <- as.matrix(df[,substring(colnames(df),14,14)=='1']);
  r$cellline <- as.matrix(df[,substring(colnames(df),14,14)=='2']);
  return(r);
}

positive_diff <- function(b1, b2, thres=0.4) {
  return(na.omit(rownames(b1)[(b1-b2)>thres]))
}







