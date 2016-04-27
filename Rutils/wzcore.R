
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
  queue='shortq',
  ppn=5,
  memG=2,
  walltime=12,
  submit=FALSE) {
  
  sink(paste0(f$dest,'.pbs'))
  cat(sprintf("
#!/bin/bash
###
#PBS -S /bin/bash
#PBS -N %s
#PBS -e %s.stderr
#PBS -o %s.stdout
#PBS -q %s
#PBS -l nodes=1:ppn=%s,mem=%sgb,walltime=%s:00:00
Rscript -e \"load('%s.rda'); myfun(args);\"
", f$jobname, f$dest, f$dest, f$queue,
              f$ppn, f$memG, f$walltime, f$dest))
  sink()

  myfun <- f$commands
  args <- f$args
  save(myfun, args, file=paste0(f$dest,'.rda'))

  cat('pbsfile: ', f$dest, '.pbs\n', sep='')
  if(submit) {
    system(paste0('qsub ', f$dest, '.pbs'))
  }
}
