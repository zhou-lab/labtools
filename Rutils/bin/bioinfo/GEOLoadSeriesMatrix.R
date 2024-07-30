#!/usr/bin/env Rscript
library(stringr)
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
sname = args[1]
system(sprintf("zcat %s | awk '!/^!/ && length($0)>0' >tmp_%s", sname, sname))
system(sprintf("zcat %s | awk '/^!Sample/' | wzmanip transpose - >%s_samplesheets.tsv", sname, sname))

samples <- read.table(sprintf("%s_samplesheets.tsv", sname), stringsAsFactors=F, header=T, sep='\t')

select.cols <- c(1,2)
select.cols <- sort(unique(c(select.cols, grep('platform.id', colnames(samples)))))
select.cols <- sort(unique(c(select.cols, grep('characteristics', colnames(samples)))))
## chr1.cols <- sort(unique(c(select.cols, grep('ch1', colnames(samples)))))
select.cols <- sort(unique(c(select.cols, grep('source_name', colnames(samples)))))

samples0 <- samples
samples <- samples0[,select.cols]
cat('Expect:', nrow(samples), 'samples\n')

## if there is a column with all sex, then remove sex and that column title
for (i in 1:ncol(samples)) {
  if (all(grepl('^Sex: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('Sex: ','', samples[,i])
    colnames(samples)[i] <- 'Sex'
  }
  if (all(grepl('^gender: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('gender: ','', samples[,i])
    colnames(samples)[i] <- 'sex'
  }
  if (all(grepl('^age: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('age: ','', samples[,i])
    colnames(samples)[i] <- 'age'
  }
  if (all(grepl('^race: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('race: ','', samples[,i])
    colnames(samples)[i] <- 'race'
  }
  if (all(grepl('^cell line: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('cell line: ','', samples[,i])
    colnames(samples)[i] <- 'cellline'
  }
  if (all(grepl('^cell type: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('cell type: ','', samples[,i])
    colnames(samples)[i] <- 'celltype'
  }
  if (all(grepl('^cell population: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('cell population: ','', samples[,i])
    colnames(samples)[i] <- 'celltype'
  }
  if (all(grepl('^disease state: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('disease state: ','', samples[,i])
    colnames(samples)[i] <- 'diseasestate'
  }
  if (all(grepl('^smoking status: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('smoking status: ','', samples[,i])
    colnames(samples)[i] <- 'smoking'
  }
  if (all(grepl('^differentiation stage: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('differentiation stage: ','', samples[,i])
    colnames(samples)[i] <- 'differentiationstage'
  }
  if (all(grepl('^age \\(y\\): ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('age \\(y\\): ','', samples[,i])
    colnames(samples)[i] <- 'age'
  }
  if (all(grepl('^source: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('source: ','', samples[,i])
    colnames(samples)[i] <- 'source'
  }
  if (all(grepl('^plate: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('plate: ','', samples[,i])
    colnames(samples)[i] <- 'plate'
  }
  if (all(grepl('^ethnicity: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('ethnicity: ','', samples[,i])
    colnames(samples)[i] <- 'race'
  }
  if (all(grepl('^tissue: ', samples[,i], ignore.case = TRUE))) {
    samples[,i] <- sub('tissue: ','', samples[,i])
    colnames(samples)[i] <- 'tissue'
  }

  fk <- str_match(samples[,i], '^([^:]*):\\s+([^ ]*)')
  if ((!any(is.na(fk))) && ncol(fk)==3 && length(unique(fk[,2])) == 1) {
    samples[,i] <- fk[,3]
    colnames(samples)[i] <- unique(fk[,2])
  }

  samples[,i] <- gsub('\\s+','.',samples[,i])
  ## samples[,i] <- gsub('[-]+','.',samples[,i])
  cat(colnames(samples)[i],'\t')
  if (all(!is.na(as.numeric(samples[,i])))) {
    cat('numeric')
    samples[,i] <- as.numeric(samples[,i])
  }
  cat('\n')
}
colnames(samples) <- gsub('\\s+','.', colnames(samples))
colnames(samples)[colnames(samples) == 'X.Sample_geo_accession'] <- 'geo'
colnames(samples)[colnames(samples) == "X.Sample_platform_id"] <- "platform"
colnames(samples)[colnames(samples) == 'X.Sample_title'] <- 'title'
colnames(samples)[colnames(samples) == 'X.Sample_source_name_ch1'] <- 'sourceName'
samples = cbind(samples[,c("geo","platform")], samples)
extra <- cbind(
	samples[,!(colnames(samples) %in% c("geo","platform","title","sourceName"))],
	description=samples0[,grepl("Sample_description", colnames(samples0))])

extra <- do.call(paste, c(as.data.frame(extra), sep="|||"))
samples <- cbind(samples[,c("geo","platform")],
	title=samples0[,grepl("Sample_title", colnames(samples0))], 
	sourceName=samples[,"sourceName"], extra=extra)
colnames(samples) <- sub("X.Sample_description","description",colnames(samples))

if (length(args) > 1) {
  if (args[2] == 'usetitle') {
    rownames(samples) <- paste0(samples$geo, '.', samples$title)
  } else if (args[2] == 'usegeo') {
    rownames(samples) <- samples$geo
  } else if (args[2] == 'usesource') {
    rownames(samples) <- paste0(samples$geo, '.', samples$sourceName)
  }
} else {
  rownames(samples) <- samples$geo
}

suppressWarnings(suppressPackageStartupMessages(library(data.table)))
a <- fread(sprintf('tmp_%s', sname), sep='\t', header=TRUE)
betas <- as.matrix(a[,2:ncol(a)])
rownames(betas) <- a$ID_REF

if (!is.numeric(betas)) {
  aa <- as.numeric(betas)
  dim(aa) <- dim(betas)
  dimnames(aa) <- dimnames(betas)
  betas <- aa
}
write.table(samples, file=sprintf('%s_samples.tsv', sname), quote=F, sep='\t', row.names=FALSE)
saveRDS(betas, file=sprintf('%s_betas_GEOLoadSeriesMatrix.rds', sname))

system(sprintf('rm -f tmp_%s', sname))
cat('Read:', ncol(betas), 'samples.\n')
cat('Meta:', colnames(samples), '\n')
cat('Sname:', rownames(samples)[1:10], '\n')
cat('#Probes:', nrow(betas), '\n')

## system('rm -f 1')
