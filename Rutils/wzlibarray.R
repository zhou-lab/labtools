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

GetTCGA <- function(
  category=NULL,
  cancer.type=NULL,
  tcga.id.map.fn='/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/merged_mapping',
  barcodes=NULL,
  barcodes20=NULL,
  base.dir='/primary/projects/laird/projects/2016_04_05_TCGA_pancan_renormalization/450k/betas',
  probes=NULL,
  ## '/Volumes/projects_primary/laird/projects/2016_01_29_NIH_3T3_run2/IDAT_merge/betas_bycancertype/'
  ## '/primary/projects/laird/projects/2016_01_29_NIH_3T3_run2/magetab/merged_mapping'
  target=c('betas','signalset','rnaseq'),
  mc=FALSE, mc.cores=8,
  n.max=NULL, nprob.max=NULL) {

  ## category can be tumor, normal or cellline
  target <- match.arg(target)

  ## read meta information
  tcga.id.map <- read.table(
    tcga.id.map.fn,
    col.names=c('cancertype','barcode','idatname'), stringsAsFactors=FALSE)
  tcga.id.map$catgry <- as.factor(sapply(
    substr(tcga.id.map$barcode,14,14), 
    function(x) switch(x, '2'='cellline','0'='tumor','1'='normal')))

  ## filtering subsets
  if (!is.null(cancer.type))
    tcga.id.map <- subset(tcga.id.map, cancertype==cancer.type)

  if (!is.null(category))
    tcga.id.map <- subset(tcga.id.map, catgry==category)

  if (!is.null(n.max))
    tcga.id.map <- tcga.id.map[1:min(nrow(tcga.id.map),n.max),,drop=FALSE]

  if (!is.null(barcodes))
    tcga.id.map <- subset(tcga.id.map, barcode %in% barcodes)

  if (!is.null(barcodes20))
    tcga.id.map <- subset(tcga.id.map, substr(barcode,1,20) %in% barcodes20)

  message('Loading ',nrow(tcga.id.map),' sample(s).')
  uniq.idats <- unique(tcga.id.map$idatname)
  message('There are ', length(uniq.idats), ' unique idats.')
  rda.fns <- list.files(base.dir, pattern='*.rda')
  rda.fns <- rda.fns[substr(rda.fns,1,nchar(rda.fns)-4) %in% uniq.idats]
  ## retrieve target
  if (target == 'signalset') {

    .processTarget <- function(rda.fn) {
      message('\r',rda.fn, '.', appendLF=FALSE)
      load(file.path(base.dir, rda.fn))
      if (is.null(probes)) {
        return(sset)
      } else {
        sset <- sset[probes]
        return(sset)
      }
    }
    if (mc) {
      library(parallel)
      ssets <- mclapply(rda.fns, .processTarget, mc.cores=mc.cores, mc.preschedule=FALSE)
    } else {
      ssets <- lapply(rda.fns, .processTarget)
    }
    names(ssets) <- substr(rda.fns,1,nchar(rda.fns)-4)
    ## attr(ssets, 'barcode') <- tcga.id.map$barcode[match(names(ssets), tcga.id.map$idatname)]
    lapply(seq_along(ssets), function(i) 
      attr(ssets[[i]],'barcode') <<- tcga.id.map$barcode[match(names(ssets)[i], tcga.id.map$idatname)])
    message('\nLoaded ', length(ssets), ' signalsets')
    return(ssets)
    
  } else if (target == 'betas') {

    all.betas <- sapply(rda.fns, function(rda.fn) {
      message('\r', rda.fn, '.', appendLF=FALSE)
      load(file.path(base.dir, rda.fn))
      if (is.null(probes)) {
        return(betas)
      } else {
        return(betas[probes])
      }
    })
    colnames(all.betas) <- substr(rda.fns,1,nchar(rda.fns)-4)
    attr(all.betas, 'barcode') <- tcga.id.map$barcode[match(colnames(all.betas),
                                                            tcga.id.map$idatname)]
    message('\nLoaded ', length(all.betas), ' betas')
    return(all.betas)

  } else {
    stop('target not implemented')
  }
}

# TCGA-07-0227-20A-01D-A418-05
