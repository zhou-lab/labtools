
printf <- function(...) cat(sprintf(...));

wz <- function(df) {
    df.dim <- dim(df);
    print(df[1:min(10,df.dim[1]),1:min(10,df.dim[2])]);
    printf("%d rows %d columns\n", df.dim[1], df.dim[2]);
}

pbsgen <- function(
  dest=NULL,
  commands=NULL, args=NULL,
  jobname='WandingJob',
  queue='default',
  ppn=5,
  memG=2,
  walltime=12,
  submit=FALSE) {
  
  sink(paste0(dest,'.pbs'))
  cat(sprintf("
#!/bin/bash
###
#PBS -S /bin/bash
#PBS -N %s
#PBS -e %s.pbs.stderr
#PBS -o %s.pbs.stdout
#PBS -q %s
#PBS -l nodes=1:ppn=%s,mem=%sgb,walltime=%s:00:00
Rscript -e \"load('%s.rda'); myfun(args);\"
", jobname, dest, dest, queue, ppn, memG, walltime, dest))
  sink()

  myfun <- commands
  args <- args
  save(myfun, args, file=paste0(dest,'.rda'))

  cat('pbsfile: ', dest, '.pbs\n', sep='')
  if(submit) {
    system(paste0('qsub ', dest, '.pbs'))
  }
}

rowMax <- function(x) {apply(x,1,max)}

suppressWarnings(library(ggplot2))
suppressWarnings(library(reshape2))
suppressWarnings(suppressPackageStartupMessages(library(devtools)))
suppressWarnings(suppressPackageStartupMessages(library(wheatmap)))
suppressWarnings(suppressPackageStartupMessages(library(sesame)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr, quiet=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(GenomicRanges, quiet=TRUE)))
suppressWarnings(suppressPackageStartupMessages(library(limma)))
theme_wz <- theme_set(theme_classic(15))
theme_wz <- theme_update(
  axis.line.x = element_line(colour = "black"),
  axis.line.y = element_line(colour = "black"),
  axis.text = element_text(colour='black'))

wzbind.list <- function(x, cat.names=NULL) {
  if (!is.null(cat.names))
    names(x) <- cat.names
  if (is.null(names(x)))
    stop('list has no names, please provide cat.names=')
  data.frame(x=do.call(c, x), cat=rep(names(x), sapply(x,length)))
}

wzbind <- function(..., names=NULL) {
  x <- list(...)
  if (is.null(names))
    names <- paste0('V',1:length(x))
  data.frame(x=do.call(c, x), cat=rep(names, sapply(x, length)))
}

wzvenn <- function(..., dnames=NULL) {
  data <- list(...)
  if (is.null(dnames))
    dnames <- 1:length(data)
  names(data) <- dnames
  library(VennDiagram)
  g <- venn.diagram(data, filename = NULL)
  grid.newpage()
  pushViewport(viewport(unit(0.1,'npc'),unit(0.1,'npc'),unit(0.8,'npc'),unit(0.8,'npc'), just=c('left','bottom')))
  grid.draw(g)
}

