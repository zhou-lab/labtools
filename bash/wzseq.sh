#!/bin/bash

export WZSEQ_REF_BASE=/home/uec-00/shared/production/genomes/
export WZSEQ_REF_MM10=$WZSEQ_REF_BASE/mm10/mm10.fa
export WZSEQ_REF=$WZSEQ_REF_MM10

####### auto setup links ######

function auto_setup_links_methlevelaverage() {
  # usage: auto_setup_links_methlevelaverage . /data/sequencing/analysis/H5MW5BGXX/results/H5MW5BGXX/H5MW5BGXX_1_PL430BS1 mm10

  local base=$1
  local mla_dir=$base/data/methlevelaverage
  if [[ ! -d $mla_dir ]]; then mkdir -p $mla_dir; fi

  sampledir=$2;
  local genome=$3;
  mla_fns=($(find $sampledir -name *.$genome.fa.bam.MethLevelAverages.metric.txt));

  if [[ ${#mla_fns[@]} -eq 1 ]]; then
    mla_fn=${mla_fns[0]};
  else
    echo '[Error] '${#mla_fns[@]}" methlevelaverage files for sample "$sample;
    printf "%s\n" "${mla_fns[@]}";
    return 1;
  fi
  
  if [[ $mla_fn =~ ResultCount_([^.]*).$genome.fa.bam.MethLevelAverages.metric.txt ]]; then
    scode=${BASH_REMATCH[1]};
  else
    echo "[Error] Unmatched methlevelaverage name: "$mla_fn;
    echo "Abort";
    return 1;
  fi
  ln -s $mla_fn $mla_dir/$scode".txt" && echo "Linked "$mla_fn;

}

function auto_setup_links_fastq() {

  # auto_setup_links_fastq basedir rootdir
  # auto_setup_links_fastq . /data/sequencing/analysis/H5MW5BGXX/
  
  local base=$(readlink -f $1)
  local rootdir=$(readlink -f $2)

  echo "Linking fastq..."
  fastqdir=$base/data/fastq
  [[ -d $fastqdir ]] || mkdir -p $fastqdir
  fastqs=($(find $rootdir -maxdepth 1 -name *.fastq.gz));
  for fastq in ${fastqs[@]}; do
    ln -s $fastq $fastqdir/$(basename $fastq) && echo "Linked "$fastq;
  done
  
}

function auto_setup_links_bam() {

  local base=$(readlink -f $1)
  local sample=$(readlink -f $2)
  local genome=$3

  bamdir=$base/data/bam
  [[ -d $bamdir ]] || mkdir -p $bamdir;
  
  # bams=($(find $sample -name *.fa.realign.mdups.recal.bam));
  bams=($(find $sample -name *.$genome.fa.mdups.bam));
  if [[ ${#bams[@]} -eq 1 ]]; then
    bam=${bams[0]};
    # if [[ $bam =~ ResultCount_([^.]*).([^.]*).fa.realign.mdups.recal.bam ]]; then
    if [[ $bam =~ ResultCount_([^.]*).$genome.fa.mdups.bam ]]; then
      ln -s $bam $bamdir/${BASH_REMATCH[1]}.bam && echo "Linked "$bam
      # ln -s ${bam%m}i $bamdir/${BASH_REMATCH[1]}.bam.bai;
      ln -s $bam.bai $bamdir/${BASH_REMATCH[1]}.bam.bai
    else
      echo "[Error] unmatched bam file name: "$bam
      return 1
    fi
  else
    echo "[Error] "${#bams[@]}" bam files exist for the sample "$sample
    [[ ${#bams[@]} -gt 1 ]] && printf "%s\n" "${bams[@]}";
    return 1
  fi

}

function auto_setup_links() {
  # usage: auto_setup_links -b . -r /data/sequencing/analysis/H55TVBGXX/ -g mm10

  local OPTIND opt base genome root
  while getopts "g:r:b:" opt; do
    case $opt in
      g) genome=$OPTARG ;;
      r) rootdir=$(readlink -f $OPTARG) ;;
      b) base=$(readlink -f $OPTARG) ;;
      \?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
      :) echo "Option -$OPTARG requires an argument." >&2; exit1;;
    esac
  done

  if [[ -z ${rootdir+x} ]]; then
    echo "[Error] Please provide root dir via -r";
    return 1;
  fi
  [[ -z ${base+x} ]] && base=$(pwd);

  ## link fastq
  auto_setup_links_fastq $base $rootdir
  echo

  flowcell=$(basename $rootdir);
  samples=($(/bin/ls $rootdir/results/$flowcell));
  echo "Found "${#samples[@]}" samples:"
  printf "\t%s\n" ${samples[@]}

  local datadir=$base/data
  [[ -d $datadir ]] || mkdir -p $datadir;
  [[ -e $datadir/root ]] || ln -s $rootdir $datadir/root;

  echo "Linking bams..."
  for sample in ${samples[@]}; do
    samplef=$rootdir/results/$flowcell/$sample
    auto_setup_links_bam $base $samplef $genome
  done
  echo

  echo "Linking methlevelaverages..."
  for sample in ${samples[@]}; do
    samplef=$rootdir/results/$flowcell/$sample
    auto_setup_links_methlevelaverage $base $samplef $genome
  done
  echo
}

######## methpipe #########

function methpipe_run() {

  local base
  [[ $# -eq 1 ]] && base=$(readlink -f $1) || base=$(pwd)

  methpipe_align $base
  # capture dependencies
  methpipe_fastqc $base
  methpipe_qualimap $base
  methpipe_merge_methlevelaverages $base
  methpipe_adaptor_c2t $base
  
}

function methpipe_adaptor() {

  local base=$1;
  
  [[ -d adaptor ]] || mkdir adaptor;

  parallel wzadaptor {} '>' adaptor/{/.}.txt ::: $base/root/*R1*.fastq.gz

  for f in adaptor/*.txt; do
    tail -2 $f | cut -d":" -f2 | cut -d" " -f2 | paste -s -d"\t" ;
  done | awk 'BEGIN{print "adaptorC\tadaptorC2T"}{print $1,$2/$1}' > merged

}

function methpipe_merge_methlevelaverages() {
  # usage: methpipe_merge_methlevelaverages data/methlevelaverage

  local basedir=$1;
  for f in $basedir/*.txt; do
    sed -e 's/://g' -e 's/%//' $f | awk -f wanding.awk -e 'BEGIN{split("",k);split("",n);split("",v);}(length($0)>0){k[length(k)+1]=$1;n[length(n)+1]=$2;v[length(v)+1]=$3/100;}END{for(i=1;i<=length(k);++i){kn[i]=k[i]"n"};print(join(k,1,length(k),"\t")"\t"join(kn,1,length(kn),"\t"));print(join(v,1,length(v),"\t")"\t"join(n,1,length(n),"\t"))}' > ${f%.txt}.processed
  done

  wzmanip concat -f $basedir/*.processed | awk -f wanding.awk -e '{n=NF;print;}END{repeat("NA", n, rep); print joina(rep,"\t");}' > $basedir/merge;

}

function methpipe_align() {
  # usage: requires base/samples

  local base
  [[ $# -eq 1 ]] && base=$(readlink -f $1) || base=$(pwd)
  while read samplecode fastq1 fastq2; do
    methpipe_bwameth -b $base -s $samplecode base/fastq/$fastq1 base/fastq/$fastq2
  done < $base/samples

}  

function absolute_path() {
  local in=$1
  out=()
  for x in $@; do
    out+=($(readlink -f $x));
  done
  echo "${out[@]}"
}

function methpipe_bwameth() {
  # usage: methpipe_bwameth -s HF5FGBGXX_1_PL150528WGBS2 -b <base> R1.fastq.gz R2.fastq.gz
  
  local OPTIND opt base reference sample_code
  while getopts "s:r:b:" opt; do
    case $opt in
      s) sample_code=$OPTARG ;;
      r) reference=$OPTARG ;;
      b) base=$(readlink -f $OPTARG) ;;
      \?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
      :) echo "Option -$OPTARG requires an argument." >&2; exit1;;
    esac
  done

  if [[ -z ${sample_code+x} ]]; then
    echo "Please provide sample code via -s. Abort."
    return 1;
  fi
  [[ -z ${reference+x} ]] && reference=$WZSEQ_REF;
  [[ -z ${base+x} ]] && base=$(pwd);

  shift $(( OPTIND - 1 ))
  
  bamdir=$base"/bam";
  [[ -d $bamdir ]] || mkdir -p $bamdir;

  pbsdir=$base/pbs
  [[ -d $pbsdir ]] || mkdir -p $pbsdir
  jobname=$sample_code"_bwameth"
  local fastqs=$(absolute_path $@)
  pbsgen one "cd $base; bwameth --reference $reference $fastqs --prefix $bamdir/$sample_code" -name $jobname -dest $pbsdir/$jobname
  jobid=$(qsub $pbsdir/$jobname);
  echo $jobid;
}

function methpipe_fastqc() {
  # methpipe_qc <base_dir>

  [[ $# -eq 1 ]] && base=$(readlink -f $1) || base=$(pwd);
  fastqcdir=$base"/fastqc"
  [[ -d $fastqcdir ]] || mkdir -p $fastqcdir
  parallel "fastqc -f bam {} -o fastqcdir/" ::: $base/data/fastq/*.fastq.gz

}

function methpipe_qualimap() {
  # methpipe_qualimap <base_dir>
  
  [[ $# -eq 1 ]] && base=$(readlink -f $1) || base=$(pwd);
  qualimapdir=$base"/qualimap"
  [[ -d $qualimap ]] || mkdir -p $qualimapdir
  parallel "qualimap --java-mem-size=10G bamqc -nt 10 -bam {} -outdir $qualimapdir/{/.} -c" ::: $base/bam/*.bam  
  
}

function methpipe_mergebam() {
  # methpipe_mergebam <base> merged.bam bam1 bam2
  
  local base=$1
  local dest=$2
  shift 2

  [[ -d $base/bam ]] || mkdir -p $base/bam

  :>$base/bam/$dest.rg.txt
  for x in $@; do
    samtools view -H $x | grep '^@RG' >>$base/bam/$dest.rg.txt
  done

  samtools merge -h $base/bam/$dest.rg.txt $base/bam/$dest.bam $@
  
}
