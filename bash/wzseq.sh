#!/bin/bash

function auto_setup_methlevelaverage() {
  # usage: auto_setup_methlevelaverage data/methlevelaverage /data/sequencing/analysis/H5MW5BGXX/results/H5MW5BGXX/H5MW5BGXX_1_PL430BS1 mm10

  local basedir=$1;
  if [[ ! -d $basedir ]]; then mkdir -p $basedir; fi

  sampledir=$2;
  genome=$3;
  mla_fns=($(find $sampledir -name *.$genome.fa.bam.MethLevelAverages.metric.txt));

  if [[ ${#mla_fns[@]} -eq 1 ]]; then
    mla_fn=${mla_fns[0]};
  else
    echo ${#mla_fns[@]}" methlevelaverage files. Abort.";
    printf "%s\n" "${mla_fns[@]}";
    return 1;
  fi
  
  if [[ $mla_fn =~ ResultCount_([^.]*).$genome.fa.bam.MethLevelAverages.metric.txt ]]; then
    scode=${BASH_REMATCH[1]};
  else
    echo "Unmatched methlevelaverage name: "$mla_fn;
    echo "Abort";
    return 1;
  fi
  ln -s $mla_fn $basedir/$scode".txt" && echo "Linked "$mla_fn;

}

function auto_setup_links() {
  # usage: auto_setup_links data /data/sequencing/analysis/H55TVBGXX/

  local basedir=$1;
  
  rootdir=$2;
  flowcell=$(basename $rootdir); #H55TVBGXX
  samples=($(/bin/ls $rootdir/results/$flowcell)); # /data/sequencing/analysis/H5MW5BGXX/results/H5MW5BGXX/H5MW5BGXX_1_PL430BS1
  
  [[ -d $basedir ]] || mkdir -p $basedir;
  [[ -e $basedir/root ]] || ln -s $rootdir $basedir/root;

  genome="";

  for sample in ${samples[@]}; do

    echo "Found sample "$sample
    samplef=$rootdir/results/$flowcell/$sample;
    bams=($(find $samplef -name *.fa.realign.mdups.recal.bam));
    if [[ ${#bams[@]} -eq 1 ]]; then
      bam=${bams[0]};
    else
      echo ${#bams[@]}" bam files exist for the sample "$sample
      [[ ${#bams[@]} -gt 1 ]] && printf "%s\n" "${bams[@]}";
      echo "Skip."
      continue;
    fi

    if [[ $bam =~ ResultCount_([^.]*).([^.]*).fa.realign.mdups.recal.bam ]]; then
      ln -s $bam $basedir/${BASH_REMATCH[1]}.bam && echo "Linked "$bam;
      ln -s ${bam%m}i $basedir/${BASH_REMATCH[1]}.bam.bai;
      _genome=${BASH_REMATCH[2]};
    else
      echo "Unmatched bam file name: "$bam
      echo "Abort.";
      return 1;
    fi

    if [[ $genome && $genome != $_genome ]]; then
      echo "Discrepant reference genome "$genome" vs " $_genome
      echo "Abort"
      return 1;
    fi
    [[ -z $genome ]] && genome=$_genome

    auto_setup_methlevelaverage $basedir/methlevelaverages $samplef $genome
  done
}

function adaptor() {

  local base=$1;
  
  [[ -d adaptor ]] || mkdir adaptor;

  parallel wzadaptor {} '>' adaptor/{/.}.txt ::: $base/root/*R1*.fastq.gz

  for f in adaptor/*.txt; do
    tail -2 $f | cut -d":" -f2 | cut -d" " -f2 | paste -s -d"\t" ;
  done | awk 'BEGIN{print "adaptorC\tadaptorC2T"}{print $1,$2/$1}' > merged

}

function merge_methlevelaverages() {

  # usage: merge_methlevelaverages data/methlevelaverage

  local basedir=$1;
  for f in $basedir/*.txt; do
    sed -e 's/://g' -e 's/%//' $f | awk -f wanding.awk -e 'BEGIN{split("",k);split("",n);split("",v);}(length($0)>0){k[length(k)+1]=$1;n[length(n)+1]=$2;v[length(v)+1]=$3/100;}END{for(i=1;i<=length(k);++i){kn[i]=k[i]"n"};print(join(k,1,length(k),"\t")"\t"join(kn,1,length(kn),"\t"));print(join(v,1,length(v),"\t")"\t"join(n,1,length(n),"\t"))}' > ${f%.txt}.processed
  done

  wzmanip concat -f $basedir/*.processed | awk -f wanding.awk -e '{n=NF;print;}END{repeat("NA", n, rep); print joina(rep,"\t");}' > $basedir/merge;

}

function wzqsub() {
  find ~/pbs/ -type f | LC_ALL=C sort | sed -n "${1},${2}"p | xargs -I {} qsub {};
}

function methpipe_bwameth_bissnp() {
  scode=$1;
  fastq1=$2;
  fastq2=$3;
  reference=$4;

  bamdir=$(pwd)"/bam";
  [[ -d $bamdir ]] || mkdir -p $bamdir;

  base=$(pwd);
  pbsgen clean;
  pbsgen one "cd $base; bwameth --reference $reference $fastq1 $fastq2 --prefix $bamdir/$scode"
  jobid=$(wzqsub 1 1);
  echo $jobid;
}
