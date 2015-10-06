#!/bin/bash

# export WZSEQ_REF_BASE=/primary/home/wandingzhou/genomes #/data/reference #/home/uec-00/shared/production/genomes/
# export WZSEQ_REF_MM10=$WZSEQ_REF_BASE/mm10/mm10.fa
# export WZSEQ_REF=$WZSEQ_REF_MM10
# export WZSEQ_REF_CGIBED=$WZSEQ_REF_BASE/mm10/cpg_island/cpgIsland.bed
# export WZSEQ_REFVERSION=mm10

# export WZSEQ_TOOLS=/primary/home/wandingzhou/tools
export PICARD=$WZSEQ_TOOLS/picard/picard-tools-1.135/picard.jar

######################
## auto setup links ##
######################

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

# auto_setup_links_fastq [basedir] [targetdir]
function auto_setup_links_fastq() {

  if [[ $# -eq 0 ]]; then
    echo " auto_setup_links_fastq basedir rootdir"
    echo " auto_setup_links_fastq . /data/sequencing/analysis/H5MW5BGXX/"
    return 1;
  fi
  
  local base=$(readlink -f $1)
  local rootdir=$(readlink -f $2)

  echo "Linking fastq..."
  fastqdir=$base/fastq
  [[ -d $fastqdir ]] || mkdir -p $fastqdir
  fastqs=($(find $rootdir -maxdepth 1 -name *.fastq.gz));
  for fastq in ${fastqs[@]}; do
    ln -s $fastq $fastqdir/$(basename $fastq) && echo "Linked "$fastq;
  done
  
}

# auto_setup_links_bam [targetdir]
function auto_setup_links_bam() {
  base=$(pwd)
  [[ -d bam ]] || mkdir -p bam
  find $1 -name *.bam |
    while read f; do
      ln -s $(readlink -f $f) bam/
      [[ -s ${f/.bam/.bam.bai} ]] && ln -s $(readlink -f ${f/.bam/.bam.bai}) bam/
      [[ -s ${f/.bam/.bai} ]] && ln -s $(readlink -f ${f/.bam/.bam.bai}) bam/
    done
}

function auto_setup_links_bam_usc() {

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

#########################################
### Whole Genome Bisulfite Sequencing ###
#########################################

##### alignment #####

# wgbs_biscuit_index_reference
# run mm10generic to setup $WZ_BISCUIT_INDEX
# use ~/pbs
function wgbs_biscuit_index_reference {
  base=$(pwd);
  [[ -d ~/pbs ]] || mkdir ~/pbs
  cmd="
biscuit index $WZSEQ_BISCUIT_INDEX
"
  jobname="biscuit_index"
  pbsfn=~/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 20 -ppn 1
  [[ $1 == "do" ]] && qsub $pbsfn
}

function wgbs_biscuit_align {
  # requires base/samples
  # format:
  # sample_code fastq1 fastq2

  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  while read sname sread1 sread2; do
    # while read samplecode fastq1 fastq2 _junk_; do
    cmd="
biscuit align $WZSEQ_BISCUIT_INDEX $base/fastq/$sread1 $base/fastq/$sread2 |samtools sort - $base/bam/${sname}
samtools index $base/bam/${sname}.bam
"
    jobname="biscuit_align_$sname"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
    [[ $1 == "do" ]] && qsub $pbsfn

    # biscuit_bwameth -b -s $samplecode -j jid_bwameth $base/fastq/$fastq1 $base/fastq/$fastq2
    # biscuit_mdup -j jid_mdup -d $jid_bwameth -i bam/$samplecode
  done < samples
  return
}

# function wgbs_biscuit_run {
#   local base
#   [[ $# -eq 1 ]] && base=$(readlink -f $1) || base=$(pwd)
#   biscuit_align $base
#   biscuit_merge_methlevelaverages $base
#   biscuit_adaptor_c2t $base
# }

# bwameth index reference
function wgbs_bwameth_index_reference {
  [[ -d ~/pbs ]] || mkdir ~/pbs
  cmd="
bwameth.py index $WZSEQ_BWAMETH_INDEX
"
  jobname="bwameth_index"
  pbsfn=~/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 20 -ppn 1
  [[ $1 == "do" ]] && qsub $pbsfn
}

function wgbs_bwameth() {

  if [[ $# -eq 0 ]]; then
    echo "biscuit_bwameth -s HF5FGBGXX_1_PL150528WGBS2 R1.fastq.gz R2.fastq.gz"
    echo "PBS: biscuit_bwameth -s HF5FGBGXX_1_PL150528WGBS2 -j jid R1.fastq.gz R2.fastq.gz"
    return 1
  fi
  
  local OPTIND opt base reference sample_code jobid _jobid depend
  while getopts "s:r:b:d:j:" opt; do
    case $opt in
      s) sample_code=$OPTARG ;;
      r) reference=$(readlink -f $OPTARG) ;;
      b) base=$(readlink -f $OPTARG) ;;
      j) jobid=$OPTARG ;;
      d) depend=$OPTARG ;;
      \?) echo "Invalid option: -$OPTARG" >&2; return 1;;
      :) echo "Option -$OPTARG requires an argument." >&2; return 1;;
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

  local fastqs=$(absolute_path $@)
  cmds="bwameth --reference $reference $fastqs -t 20 --prefix $bamdir/$sample_code"
  echo $cmds
  if [[ -z ${jobid+x} ]]; then
    $cmds;
  else
    pbsdir=$base/pbs
    [[ -d $pbsdir ]] || mkdir -p $pbsdir
    jobname=$sample_code"_bwameth"
    pbsgen one "$cmds" -name $jobname -dest $pbsdir/$jobname -ppn 20 -memG 20
    _jobid=$(qsub $pbsdir/$jobname);
    eval $jobid=$_jobid;
    echo "Submitted "$_jobid;
  fi
}

function wgbs_adaptor() {

  local base=$1;
  
  [[ -d adaptor ]] || mkdir adaptor;

  parallel wzadaptor {} '>' adaptor/{/.}.txt ::: $base/root/*R1*.fastq.gz

  for f in adaptor/*.txt; do
    tail -2 $f | cut -d":" -f2 | cut -d" " -f2 | paste -s -d"\t" ;
  done | awk 'BEGIN{print "adaptorC\tadaptorC2T"}{print $1,$2/$1}' > merged

}

function wgbs_merge_methlevelaverages() {
  # usage: biscuit_merge_methlevelaverages data/methlevelaverage

  local basedir=$1;
  for f in $basedir/*.txt; do
    sed -e 's/://g' -e 's/%//' $f | awk -f wanding.awk -e 'BEGIN{split("",k);split("",n);split("",v);}(length($0)>0){k[length(k)+1]=$1;n[length(n)+1]=$2;v[length(v)+1]=$3/100;}END{for(i=1;i<=length(k);++i){kn[i]=k[i]"n"};print(join(k,1,length(k),"\t")"\t"join(kn,1,length(kn),"\t"));print(join(v,1,length(v),"\t")"\t"join(n,1,length(n),"\t"))}' > ${f%.txt}.processed
  done

  wzmanip concat -f $basedir/*.processed | awk -f wanding.awk -e '{n=NF;print;}END{repeat("NA", n, rep); print joina(rep,"\t");}' > $basedir/merge;

}

function absolute_path() {
  local in=$1
  out=()
  for x in $@; do
    out+=($(readlink -f $x));
  done
  echo "${out[@]}"
}

function biscuit_mdup() {

  if [[ $# -eq 0 ]]; then
    echo "biscuit_mdup -i bam/H75J7BGXX_1_Undetermined.bam"
    echo "PBS: biscuit_mdup -j jid -d depend -i bam/H75J7BGXX_1_Undetermined.bam"
    return 1
  fi

  local OPTIND opt base jobid _jobid depend bam1 bam2 duplog tmp
  while getopts "b:i:d:j:" opt; do
    case $opt in
      b) base=$(readlink -f $OPTARG) ;;
      i) bam1=$(readlink -f $OPTARG) ;; # required
      j) jobid=$OPTARG ;;
      d) depend=$OPTARG ;;
      \?) echo "Invalid option: -$OPTARG" >&2; return 1;;
      :) echo "Option -$OPTARG requires an argument." >&2; return 1;;
    esac
  done

  [[ -z ${base+x} ]] && base=$(pwd);
  [[ -z ${depend+x} ]] && depend="" || depend="-depend "$depend;
  
  bam2=${bam1%.bam}.mdup.bam
  duplog=${bam1%.bam}.mdup_metrics
  tmp=$base/tmp
  [[ -d $tmp ]] || mkdir -p $tmp
  cmds="java -Xmx7g -jar $PICARD MarkDuplicates CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=$duplog READ_NAME_REGEX=null INPUT=$bam1 OUTPUT=$bam2 TMP_DIR=$tmp"
  echo $cmds
  if [[ -z ${jobid+x} ]]; then
    $cmds
  else
    pbsdir=$base/pbs
    [[ -d $pbsdir ]] || mkdir -p $pbsdir
    jobname=$(basename $bam1)"_mdup"
    pbsgen one "$cmds" -name $jobname -dest $pbsdir/$jobname $depend
    _jobid=$(qsub $pbsdir/$jobname)
    eval $jobid=$_jobid;
    echo "Submitted "$_jobid;
  fi
}

function biscuit_pileup() {

  local OPTIND opt base reference verbose

  verbose=false
  while getopts "b:r:v" opt; do
    case $opt in
      r) reference=$(readlink -f $OPTARG) ;;
      b) base=$(readlink -f $OPTARG) ;;
      v) verbose=true ;;
      \?) echo "Invalid option: -$OPTARG" >&2; return 1;;
      :) echo "Option -$OPTARG requires an argument." >&2; return 1;;
    esac
  done

  [[ -z ${base+x} ]] && base=$(pwd)
  [[ -z ${reference+x} ]] && reference=$WZSEQ_REF
  
  echo "[$(date)] Pileup reference: $reference"
  [[ -d $base/pileup ]] || mkdir -p $base/pileup
  for bam in $base/data/bam/*.bam; do
    sample_code=$(basename $bam .bam);
    echo "[$(date)] Pileing up $bam"
    cmd="pileup_cytosine -r $reference -i $bam -q 10 | bgzip > $base/pileup/$sample_code.pileup.gz;"

    # for bsmap bam file, make sure use NM filter "-n 2".
    # for bwa-meth bam file, mapping quality filter is used by default.
    [[ $verbose == true ]] && echo "[$(date)] $cmd"
    eval $cmd
    if [[ $? -ne 0 ]]; then
      echo "[$(date)] Pileup failure"
      return 1
    fi
    echo "[$(date)] Done"
  done

  # we need a relatively new gnu sort >= 8.23 that supports parallel sorting
  temp=$base/temp
  [[ -d $temp ]] || mkdir -p $temp

  # sort and index pileup file
  for pileup in $base/pileup/*.pileup.gz; do
    echo "[$(date)] Indexing" $pileup
    tabix -p bed $pileup
    echo "[$(date)] Done"
  done

  echo "[$(date)] All done"

}

function decho() {
  echo "[$(date)] "$@ >&2
}

function biscuit_diffmeth() {

  # usage: biscuit_diffmeth -t pileup1 -n pileup2 [-b base] [-c mincov]
  local OPTARG OPTIND opt base pileup1 pileup2 mincov analysis

  while getopts "b:t:n:c:" opt; do
    case $opt in
      t) pileup1=$(readlink -f $OPTARG) ;;
      n) pileup2=$(readlink -f $OPTARG) ;;
      b) base=$(readlink -f $OPTARG) ;;
      c) mincov=$OPTARG ;;
      \?) echo "Invalid option: -$OPTARG" >&2; return 1;;
      :) echo "Option -$OPTARG requires an argument." >&2; return 1;;
    esac
  done
  
  local max_hyper_len=500
  local max_hypo_len=2000
  local minhypercnt=6
  local minhypocnt=6
  
  [[ -z ${base+x} ]] && base=$(pwd)
  [[ -z ${mincov+x} ]] && mincov=3

  analysis=$base/diffmeth
  [[ -d $analysis ]] || mkdir -p $analysis
  local contrast=$analysis/$(basename $pileup1 .pileup.gz)"_vs_"$(basename $pileup2 .pileup.gz)

  echo "[$(date)] Generating "$contrast.diff
  bedtools intersect -a <(zcat $pileup1 | awk -v mincov=$mincov '$8!="."&&$9!="."&&$8+$9>=mincov{print $1,$2,$3,$4,$6,$7,$8,$9}') -b <(zcat $pileup2 | awk -v mincov=$mincov '$8!="."&&$9!="."&&$8+$9>=mincov{print $1,$2,$3,$4,$6,$7,$8,$9}') -sorted -wo | awk -f wanding.awk -e '$5~/[ATCG]CG/{bt=$7/($7+$8); bn=$15/($15+$16); print joinr(1,16),bt,bn,bt-bn}' > $contrast.diff
  echo "[$(date)] Done"

  echo "[$(date)] Plotting delta beta distribution"
  echo "[$(date)] Generating "$contrast.diff.dist.png
  cut -f19 $contrast.diff | wzplot hist -c 1 -o $contrast.diff.dist.png --xlabel "delta_beta" --ylabel "#CpG"
  echo "[$(date)] Done"
  
  echo "[$(date)] Generating "$contrast.bw
  cut -f1,2,3,19 $contrast.diff > $contrast.bedGraph
  bedGraphToBigWig $contrast.bedGraph ~/genomes_link/mm10/mm10.fa.fai $contrast.bw
  rm -f $contrast.bedGraph
  echo "[$(date)] Done"
  
  echo "[$(date)] Generating "$contrast.diff.window
  perl -alne 'BEGIN{my @window; my @poses; my @chrs;}{if (scalar @window == 10) {$p1 = scalar grep {$_>0.3} @window; $n1 = scalar grep {$_<-0.3} @window; print $chrs[0]."\t".$poses[0]."\t".$poses[9]."\t".$p1."\t".$n1; shift(@window); shift(@poses); shift(@chrs);} if (scalar @window>0 && $chrs[scalar @chrs-1] ne $F[0]) {@window=();@chrs=();@poses=();} push(@window,$F[18]); push(@chrs, $F[0]); push(@poses, $F[1]);}' $contrast.diff > $contrast.diff.window
  echo "[$(date)] Done"

  local windowsizeplot=$contrast"_diffmeth_windowsize_dist.png"
  echo "[$(date)] Plotting differential methylation window size distribution"
  echo "[$(date)] Generating "$windowsizeplot
  awk '{print $3-$2}' $contrast.diff.window | wzplot hist -c 1 -o $windowsizeplot --xlabel "window size (bp)" --ylabel "count"
  echo "[$(date)] Done"

  decho "Plotting hyper-methylation window size distribution."
  awk -v p=$minhypercnt '$4>p{print $3-$2}' $contrast.diff.window | wzplot hist -o $contrast"_diff_hyper_window_size_dist.png" --xlabel "window length" --ylabel "#windows" --xlog

  decho "Plotting hypo-methylation window size distribution."
  awk -v p=$minhypocnt '$5>p{print $3-$2}' $contrast.diff.window | wzplot hist -o $contrast"_diff_hypo_window_size_dist.png" --xlabel "window length" --ylabel "#windows" --xlog

  decho "Combining hyper-methylation segment."
  awk -v m=$max_hyper_len -v p=$minhypercnt '(($3-$2)<m) && ($4>p)' $contrast.diff.window | bedtools merge -i - -c 4 -o count > $contrast.hyper.bed
  decho "Combining hypo-methylation segment."
  awk -v m=$max_hypo_len -v p=$minhypocnt '(($3-$2)<m) && ($5>p)' $contrast.diff.window | bedtools merge -i - -c 4 -o count > $contrast.hypo.bed

  decho "Intersecting CpG Island."
  
  bedtools intersect -a $contrast.hyper.bed -b $WZSEQ_REF_CGIBED > $contrast.hyper.cgi.bed
  bedtools intersect -a $contrast.hypo.bed -b $WZSEQ_REF_CGIBED > $contrast.hypo.cgi.bed
  decho "There are "$(wc -l $contrast.hyper.cgi.bed | cut -d" " -f1)" hyper-methylated CGIs and "$(wc -l $contrast.hypo.cgi.bed | cut -d" " -f1)" hypo-methylated CGIs."

  decho "TransVar annotation."
  awk '{print $1":"$2"_"$3}' $contrast.hyper.bed | transvar ganno -l - --refversion $WZSEQ_REFVERSION --ccds >$contrast.hyper.transvar
  awk '{print $1":"$2"_"$3}' $contrast.hypo.bed | transvar ganno -l - --refversion $WZSEQ_REFVERSION --ccds >$contrast.hypo.transvar
  
  decho "All done"
}

function biscuit_cpgisland() {

  local cgipileup
  [[ $# -eq 1 ]] && base=$(readlink -f $1) || base=$(pwd)
  [[ -d $base/cpgisland ]] || mkdir -p $base/cpgisland
  for pileup in pileup/*.pileup.gz; do
    decho "Averaging "$pileup
    pileupname=$(basename $pileup .pileup.gz)
    cgipileup=$base/cpgisland/$pileupname.cgi
    # wzcpgisland.py methlevelaverage -p $pileup -c $WZSEQ_REF_CGIBED > $cgipileup.bed
    decho "Plotting "$cgipileup.png
    awk '$12!="NA"' $cgipileup.bed | wzplot hexbin -x 5 -y 12 --nolog --bins 20 -o $cgipileup.png --xlabel "CGI size" --ylabel "methlevelaverage"
  done

  decho 'Done'
}

#######################
## ChIP-seq pipeline ##
#######################

function chippipe_bcp() {

  [[ $# -eq 3 ]] || { echo "not enough arguments"; return 1; }
  targetbed=$1
  inputbed=$2
  results=$3
  BCP_HM -1 $targetbed -2 $inputbed -f 200 -w 200 -p 0.001 -3 $results

}

###########################
##### RNA-seq pipeline ####
###########################

# flowchat
# rnaseq_tophat2 => rnaseq_cufflinks => rnaseq_cuffmerge => rnaseq_cuffquant (optional) => rnaseq_cuffdiff

# rnaseq_tophat2
function rnaseq_tophat2() {
  if [[ ! -s samples ]]; then
    echo "No [samples] file. Abort."
    return 1;
  fi
  samplefn=samples
  base=$(pwd);
  [[ -d bam ]] || mkdir -p bam
  [[ -d pbs ]] || mkdir -p pbs
  while read sname sread1 sread2; do
    sfile1=$(readlink -f fastq/$sread1);
    [[ $sread2 == "." ]] && sfile2="" || sfile2=$(readlink -f fastq/$sread2); # in case of single ended.
    odir=$(readlink -f bam/$sname);
    # customize the following
    # default: no novel junction, 28 threads
    cmd="
tophat2 -p 28 -G $WZSEQ_GTF --library-type fr-unstranded -o $odir --no-novel-juncs $WZSEQ_BOWTIE2_INDEX $sfile1 $sfile2
samtools index $odir/accepted_hits.bam
samtools flagstat $odir/accepted_hits.bam > $odir/accepted_hits.bam.flagstat
"
    jobname="tophat_$sname"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
    [[ $1 == "do" ]] && qsub $pbsfn
  done <$samplefn
}

# rnaseq_tophat2_firststrand samplefn
function rnaseq_tophat2_firststrand() {
  if [[ ! -s samples ]]; then
    echo "file: samples missing. Abort."
    return 1;
  fi
  base=$(pwd);
  [[ -d bam ]] || mkdir -p bam
  [[ -d pbs ]] || mkdir -p pbs
  while read sname sread1 sread2; do
    sfile1=$(readlink -f fastq/$sread1);
    sfile2=$(readlink -f fastq/$sread2);
    odir=$(readlink -f bam/$sname);
    cmd="
tophat2 -p 28 -G /primary/vari/genomicdata/genomes/hg19/tophat/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf --library-type fr-firststrand -o $odir --no-novel-juncs ~/references/hg19/bowtie2/hg19 $sfile1 $sfile2
samtools index $odir/accepted_hits.bam
samtools flagstat $odir/accepted_hits.bam > $odir/accepted_hits.bam.flagstat
cd $base/bam
ln -s $sname/accepted_hits.bam $sname.bam
ln -s $sname/accepted_hits.bam.bai $sname.bam.bai
"
    jobname="tophat_$sname"
    pbsgen one "$cmd" -name $jobname -dest $base/pbs/$jobname.pbs -hour 24 -memG 250 -ppn 28
  done <samples
}

# rnaseq_cufflinks <do>
function rnaseq_cufflinks() {
  [[ -d cufflinks ]] || mkdir -p cufflinks;
  base=$(pwd);
  for sample in bam/*.bam; do
    sample=$(basename $sample .bam)
    # [[ -s bam/$sample.bam ]] || continue
    cmd="cufflinks -p 28 -g /primary/vari/genomicdata/genomes/hg19/tophat/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf -o $base/cufflinks/$sample $base/bam/$sample.bam -q"
    jobname="cufflinks_$sample"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 100 -ppn 28
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

# rnaseq_cuffmerge <do>
# require: cuffmerge, gtf_to_sam, cuffcompare
#
# example of cuffmerge/assemblies.txt
# cufflinks/82_SL121326/transcripts.gtf
# cufflinks/83_SL121327/transcripts.gtf
# cufflinks/84_SL121328/transcripts.gtf
# cufflinks/85_SL121329/transcripts.gtf
# cufflinks/BC/transcripts.gtf
function rnaseq_cuffmerge() {
  [[ -d cuffmerge ]] || mkdir -p cuffmerge;
  base=$(pwd)
  if [[ ! -s cuffmerge/assemblies.txt ]]; then
    echo "cuffmerge/assemblies.txt nonexistent. Create new."
    :>cuffmerge/assemblies.txt
    for f in $base/cufflinks/*; do
      [[ -s $f/transcripts.gtf ]] && echo $f/transcripts.gtf >> cuffmerge/assemblies.txt;
    done
  fi
  cmd="
cd $base/cuffmerge/
cuffmerge -g $WZSEQ_GTF -s $WZSEQ_REFERENCE -p 10 $base/cuffmerge/assemblies.txt
"
  jobname='cuffmerge'
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 1 -memG 2 -ppn 10
  [[ $1 == "do" ]] && qsub $pbsfn
}

# run cuffquant before this
# cuffquant <do>
function rnaseq_cuffquant() {

  base=$(pwd)
  gtf=$(readlink -f cuffmerge/merged_asm/merged.gtf)
  if [[ ! -s $gtf ]]; then
    echo "cuffmerge/merged_asm/merged.gtf missing. fall back to $WZSEQ_GTF"
    return 1
  fi

  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    bf=$(basename $f .bam)
    cmd="
cuffquant $gtf $fn -o $base/cuffmerge/cuffquant_$bf -p 8 -q
"
    jobname="cuffquant_$bf"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 8
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

# cuffdiff <do>
# cuffdiff/sample_sheet file example:
# cond1\tcond2\tbams1\tbams2
# e.g.
# aza     PBS     bam/aza.bam     bam/PBS.bam
# for cuffdiff >2.2 bams can be replaced by cxb files.
# bams1 and bams2 can be comma separated if containing multiple samples
function rnaseq_cuffdiff() {
  base=$(pwd);
  [[ -d cuffdiff ]] || mkdir -p cuffdiff;

  if [[ ! -s cuffdiff/sample_sheet ]]; then
    echo "file: cuffdiff/sample_sheet missing"
    return 1
  fi

  while read cond1 cond2 bams1 bams2; do
    # use merged gtf if available otherwise, use back-up gtf
    gtf=$base/cuffmerge/merged_asm/merged.gtf
    [[ -s $gtf ]] || gtf=$WZSEQ_GTF

    cmd="
cd $base
cuffdiff -o $base/cuffdiff/${cond1}_vs_${cond2} -q -b $WZSEQ_REFERENCE -p 8 -L $cond1,$cond2 -u $gtf $bams1 $bams2
"
    jobname="cuffdiff_${cond1}_vs_${cond2}"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 8
    [[ $1 == "do" ]] && qsub $pbsfn
  done < cuffdiff/sample_sheet
}

# TODO: cuffnorm
function rnaseq_cuffnorm() {
  return 1
}

# RSeQC
# need to define WZSEQ_RSEQ_GENE_BED
function rnaseq_rseqc() {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bam)

    cmd="
cd $base

# compute and plot SAM quality score
echo `date` 'running read_quality.py ...' 1>&2
read_quality.py -i $f -o rseqc/${bfn}_read_quality

# computes nucleotide composition tables and make plots
echo `date` 'running read_NVC.py ...' 1>&2
read_NVC.py -i $f -o rseqc/${bfn}_nucleotide_composition

# number of reads mapped to each genomic feature, CDS_Exons, 3'UTR_Exons, etc
echo `date` 'running read_distribution.py ...' 1>&2
read_distribution.py -i $fn -r $WZSEQ_RSEQC_GENE_BED >rseqc/${bfn}_read_distribution

# 5'->3' gene coverage metrics for each sample (coverage uniformity)
echo `date` 'running geneBody_coverage.py ...' 1>&2
geneBody_coverage.py -i $fn -r $WZSEQ_RSEQC_GENE_BED -o rseqc/${bfn}_genebody_coverage
"
    jobname="rseqc_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done

}

# subjunc, subread (if only differential gene expression is concerned, faster), here I
# only use subjunc, which is supposed to be more optimal
# need to build subread index first (for subjunc)
# bin/subread-buildindex -o ~/references/hg19/subread/hg19 ~/references/hg19/hg19.fa
# other options: -S fr -d <min_insert_size> -D <max_insert_size>
function rnaseq_subjunc() {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  while read sname sread1 sread2; do
    sfile1=$(readlink -f fastq/$sread1);
    sfile2=$(readlink -f fastq/$sread2);
    mkdir -p bam/${sname}_subjunc;
    cmd="
cd $base
/primary/home/wanding.zhou/tools/subread/subread-1.4.6-p5-Linux-x86_64/bin/subjunc -T 28 -I 16 -i $WZSEQ_SUBREAD_INDEX -r fastq/$sread1 -R fastq/$sread2 --gzFASTQinput -o bam/${sname}_subjunc/$sname.bam --BAMoutput
[[ -e bam/${sname}.bam ]] || ln -s ${sname}_subjunc/$sname.bam bam/${sname}.bam
"
    jobname="subjunc_${sname}"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done <samples
}

# TODO: featureCounts
function rnaseq_featurecounts() {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bam)
    cmd="
featureCounts -T 28 -t exon -g gene_id -a annotation.gtf -o featureCounts/$bfn_counts.txt $fn
"
    jobname="featurecounts_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done

}

# TODO: limma

# TODO: RSEM
# ~/tools/rsem/rsem-1.2.22/rsem-prepare-reference
# rsem-calculate-expression -p 20 --calc-ci --ci-memory 12294 --bowtie-chunkmbs 2000 --paired-end --bowtie-path $WZSEQ_BOWTIE1 --rsem-index RSEM

# TODO: EBSeq

# STAR
# make a genome index first
# STAR --runThreadN 28 --runMode genomeGenerate --genomeDir $WZSEQ_STAR_INDEX --genomeFastaFiles $WZSEQ_REFERENCE --sjdbGTFfile $WZSEQ_GTF --sjdbOverhang 49
function rnaseq_star() {
  if [[ ! -s samples ]]; then
    echo "file: samples missing. Abort"
    return 1
  fi
  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  while read sname sread1 sread2; do
    cmd="
mkdir $base/bam/$sname
STAR --runThreadN 28 --genomeDir $WZSEQ_STAR_INDEX --readFilesIn $base/fastq/$sread1 $base/fastq/$sread2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $base/bam/$sname/$sname
samtools index $base/bam/$sname/$snameAligned.sortedByCoord.out.bam
"
    jobname="star_$sname"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
    [[ $1 == "do" ]] && qsub $pbsfn
  done <samples
}


# mapsplice
# input fastq must be uncompressed and with no space
function rnaseq_mapsplice() {

  if [[ ! -s samples ]]; then
    echo "file: samples missing. Abort"
    return 1
  fi
  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  while read sname sread1 sread2; do
    cmd="
zcat $base/fastq/$sread1 | sed 's/ /_/g' >$base/fastq/_${sread1}_tmp
zcat $base/fastq/$sread2 | sed 's/ /_/g' >$base/fastq/_${sread2}_tmp
python /home/wanding.zhou/tools/mapsplice/MapSplice-v2.2.0/mapsplice.py -p 28 -o $base/bam/${sname}_mapsplice --bam -c $WZSEQ_BOWTIE1_INDEX -x $WZSEQ_BOWTIE1_INDEX/mm10 -1 $base/fastq/_${sread1}_tmp -2 $base/fastq/_${sread2}_tmp --gene-gtf $WZSEQ_GTF_ENSEMBL
rm -f $base/fastq/_${sread1}_tmp $base/fastq/_${sread2}_tmp
"
    jobname="mapsplice_$sname"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
    [[ $1 == "do" ]] && qsub $pbsfn
  done <samples

}

#####################
### other utility ###
#####################

function wzseq_sra_to_fastq() {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d fastq ]] || mkdir fastq
  for f in sra/*.sra; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .sra)
    # if you want to append 1/2 to read name use -I. But this
    # usually causes trouble
    cmd="
fastq-dump --split-files $fn -O $base/fastq --gzip
"
    jobname="sra2fastq_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 10 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

# wzseq_merge_fastq target1.fastq.gz source1.fq.gz,source2.fq.gz target2.fq.gz source3.fq.gz,source4.fq.gz ... <do>
function wzseq_merge_fastq {
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  until [[ -z $2 ]]; do
    target=$1
    sources=${2//,/ }
    targetdir=$(dirname $target)
    cmd="
cd $base
[[ -d $targetdir ]] || mkdir -p $targetdir
zcat $sources | gzip -c > $target
echo 'read count for '$target
zcat $target | wc -l
echo 
for source in $sources; do 
echo \$source
zcat \$source | wc -l
done
"
    jobname="merge_fastq_${target//\//_}"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
    shift 2
  done
}

function wzseq_merge_bam {
  # wzseq_mergebam merged.bam bam1 bam2...
  # if $dest.rg.txt does not exist, will create by inferring from file name

  base=$(pwd)
  local dest=$1
  shift

  [[ -d $base/bam ]] || mkdir -p $base/bam

  if [[ ! -s $base/bam/$dest.rg.txt ]]; then
    :>$base/bam/$dest.rg.txt
    for x in $@; do
      samtools view -H $x | grep '^@RG' >>$base/bam/$dest.rg.txt
    done
  fi

  samtools merge -h $base/bam/$dest.rg.txt $base/bam/$dest.bam $@
  samtools index $base/bam/$dest.bam

  rm -f $base/bam/$dest.rg.txt;
  
  echo "Done"
  
}


# wzseq_index_flagstat_bam <do>
function wzseq_index_flagstat_bam() {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  find bam -name *.bam |
    while read f; do
      fn=$(readlink -f $f)
      cmd="
samtools index $fn;
samtools flagstat $fn > $fn.flagstat;
"
      jobname="bamindex_"${f//\//_}
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 1 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function wzseq_fastqc() {
  base=$(pwd);
  [[ -d fastqc ]] || mkdir fastqc
  [[ -d pbs ]] || mkdir pbs
  for f in fastq/*.+(fastq|fq|fastq.gz|fq.gz); do
    fn=$(readlink -f $f)
    bfn=$(basename $f)
    bfn=${bfn%.fastq.gz}
    bfn=${bfn%.fastq}
    bfn=${bfn%.fq.gz}
    bfn=${bfn%.fq}
    cmd="
[[ -d $base/fastqc/$bfn ]] || mkdir -p $base/fastqc/$bfn
fastqc -f fastq $fn -o $base/fastqc/$bfn
"
    jobname="fastqc_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 5 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

# wzseq_qualimap <do>
function wzseq_qualimap() {
  
  base=$(pwd);
  [[ -d qualimap ]] || mkdir -p qualimap
  qualimapdir=$base/qualimap
  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bam)     # fn might be a symbolic link
    cmd="
qualimap --java-mem-size=10G bamqc -nt 10 -bam $fn -outdir $qualimapdir/$bfn -c
"
    jobname="qualimap_"${f//\//_}
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 1 -memG 50 -ppn 10
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

################
# trim adaptor #
################

function wzseq_trimmomatic {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  if [[ ! -s samples ]]; then
    echo "file: samples missing"
    return 1;
  fi
  
  [[ -d fastq/trimmed ]] || mkdir -p fastq/trimmed
  [[ -d fastq/untrimmed ]] || mkdir -p fastq/untrimmed
  while read sname sread1 sread2; do
    sread1base=$(basename $sread1 .fastq.gz)
    sread2base=$(basename $sread2 .fastq.gz)
    cmd="
java -jar /primary/home/wanding.zhou/tools/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 $base/fastq/$sread1 $base/fastq/$sread2 $base/fastq/trimmed/$sread1 $base/fastq/trimmed/${sread1base}_unpaired.fastq.gz $base/fastq/trimmed/$sread2 $base/fastq/trimmed/${sread2base}_unpaired.fastq.gz ILLUMINACLIP:/primary/home/wanding.zhou/tools/trimmomatic/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
mv $base/fastq/$sread1 $base/fastq/untrimmed
mv $base/fastq/$sread2 $base/fastq/untrimmed
cd $base/fastq
ln -s trimmed/$sread1 .
ln -s trimmed/$sread2 .
"
    jobname="trimmomatic_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done < samples

}
