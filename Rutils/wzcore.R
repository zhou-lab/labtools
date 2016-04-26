
printf <- function(...) cat(sprintf(...));

wz <- function(df) {
    df.dim <- dim(df);
    print(df[1:min(10,df.dim[1]),1:min(10,df.dim[2])]);
    printf("%d rows %d columns\n", df.dim[1], df.dim[2]);
}

PBS <- function(
  jobname=NULL,
  dest=NULL,
  queue='shortq',
  ppn=5,
  memG=2,
  walltime=12,
  commands=NULL) {
  pbsf <- list(jobname=jobname, dest=dest, queue=queue,
               ppn=ppn,memG=memG,walltime=walltime,command=command)
  class(pbsf) <- 'PBSFile'
  pbsf
}

pbsgen <- function(f, submit=FALSE) {
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
%s
", f$jobname, f$dest, f$dest, f$queue,
              f$ppn, f$memG, f$walltime, f$command))
  sink()
  cat('pbsfile: ', f$dest, '.pbs\n')
  if(submit) {
    system(paste0('qsub ', f$dest, '.pbs'))
  }
}
