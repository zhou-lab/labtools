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

###############
### biscuit ###
###############

function biscuit_run() {

  local base
  [[ $# -eq 1 ]] && base=$(readlink -f $1) || base=$(pwd)

  biscuit_align $base
  # capture dependencies
  biscuit_fastqc $base
  # biscuit_qualimap $base
  biscuit_merge_methlevelaverages $base
  biscuit_adaptor_c2t $base
  
}

function biscuit_adaptor() {

  local base=$1;
  
  [[ -d adaptor ]] || mkdir adaptor;

  parallel wzadaptor {} '>' adaptor/{/.}.txt ::: $base/root/*R1*.fastq.gz

  for f in adaptor/*.txt; do
    tail -2 $f | cut -d":" -f2 | cut -d" " -f2 | paste -s -d"\t" ;
  done | awk 'BEGIN{print "adaptorC\tadaptorC2T"}{print $1,$2/$1}' > merged

}

function biscuit_merge_methlevelaverages() {
  # usage: biscuit_merge_methlevelaverages data/methlevelaverage

  local basedir=$1;
  for f in $basedir/*.txt; do
    sed -e 's/://g' -e 's/%//' $f | awk -f wanding.awk -e 'BEGIN{split("",k);split("",n);split("",v);}(length($0)>0){k[length(k)+1]=$1;n[length(n)+1]=$2;v[length(v)+1]=$3/100;}END{for(i=1;i<=length(k);++i){kn[i]=k[i]"n"};print(join(k,1,length(k),"\t")"\t"join(kn,1,length(kn),"\t"));print(join(v,1,length(v),"\t")"\t"join(n,1,length(n),"\t"))}' > ${f%.txt}.processed
  done

  wzmanip concat -f $basedir/*.processed | awk -f wanding.awk -e '{n=NF;print;}END{repeat("NA", n, rep); print joina(rep,"\t");}' > $basedir/merge;

}

function biscuit_align() {
  # usage: requires base/samples
  # format:
  # sample_code fastq1 fastq2

  local base
  [[ $# -eq 1 ]] && base=$(readlink -f $1) || base=$(pwd)
  while read samplecode fastq1 fastq2 _junk_; do
    biscuit_bwameth -b $base -s $samplecode -j jid_bwameth $base/fastq/$fastq1 $base/fastq/$fastq2
    biscuit_mdup -j jid_mdup -d $jid_bwameth -i bam/$samplecode
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

function biscuit_bwameth() {

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

function biscuit_fastqc() {
  # biscuit_qc <base_dir>

  [[ $# -eq 1 ]] && base=$(readlink -f $1) || base=$(pwd);
  fastqcdir=$base"/fastqc"
  [[ -d $fastqcdir ]] || mkdir -p $fastqcdir
  parallel "fastqc -f bam {} -o fastqcdir/" ::: $base/fastq/*.fastq.gz

}

function biscuit_mergebam() {
  # biscuit_mergebam <base> merged.bam bam1 bam2...
  
  local base=$1
  local dest=$2
  shift 2

  [[ -d $base/bam ]] || mkdir -p $base/bam

  :>$base/bam/$dest.rg.txt
  for x in $@; do
    samtools view -H $x | grep '^@RG' >>$base/bam/$dest.rg.txt
  done

  samtools merge -h $base/bam/$dest.rg.txt $base/bam/$dest.bam $@
  samtools index $base/bam/$dest.bam

  rm -f $base/bam/$dest.rg.txt;
  
  echo "Done"
  
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
###############
## ChIP pipe ##
###############

function chippipe_bcp() {

  [[ $# -eq 3 ]] || { echo "not enough arguments"; return 1; }
  targetbed=$1
  inputbed=$2
  results=$3
  BCP_HM -1 $targetbed -2 $inputbed -f 200 -w 200 -p 0.001 -3 $results

}

###########################
##### RNA-seq pipeline ####
#
# rnaseq_tophat2
# rnaseq_cufflinks
# rnaseq_cuffmerge <do>
# rnaseq_cuffquant <do>
###########################

# rnaseq_tophat2
function rnaseq_tophat2() {
  if [[ -z samples ]]; then
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
  if [[ -z samples ]]; then
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
ln -s $odir/accepted_hits.bam $base/bam/$sname.bam
ln -s $odir/accepted_hits.bam.bai $base/bam/$sname.bam.bai
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
##### cuffmerge/assemblies.txt #########
# cufflinks/82_SL121326/transcripts.gtf
# cufflinks/83_SL121327/transcripts.gtf
# cufflinks/84_SL121328/transcripts.gtf
# cufflinks/85_SL121329/transcripts.gtf
# cufflinks/BC/transcripts.gtf
########################################
function rnaseq_cuffmerge() {
  [[ -d cuffmerge ]] || mkdir -p cuffmerge;
  base=$(pwd)
  if [[ -z cuffmerge/assemblies.txt ]]; then
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

  gtf=$(readlink -f cuffmerge/merged_asm/merged.gtf)
  if [[ -z $gtf ]]; then
    echo "cuffmerge/merged_asm/merged.gtf missing. Abort"
    return 1
  fi

  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    bf=$(basename $f .bam)
    cmd="
cuffquant $WZSEQ_GTF $fn -o cuffmerge/cuffquant_$bf -p 8
"
    jobname="cuffquant_$bf"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 8
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

# cuffdiff
function rnaseq_cuffdiff() {
  base=$(pwd);
  cond1=$1
  cond2=$2
  bams1=$(readlink -f $3)       # comma separated, if multiple
  bams2=$(readlink -f $4)       # comma separated, if multiple
  [[ -d cuffdiff ]] || mkdir -p cuffdiff;
  
  # use merged gtf if available otherwise, use back-up gtf
  gtf=$base/cuffmerge/merged_asm/merged.gtf
  [[ -s $gtf ]] || gtf=$WZSEQ_GTF
  
  cmd="
cuffdiff -o $base/cuffdiff/${cond1}_vs_${cond2} -q -b $WZSEQ_REFERENCE -p 28 -L $cond1,$cond2 -u $gtf $bams1 $bams2
"
  jobname="cuffdiff_${cond1}_vs_${cond2}"
  pbsgen one "$cmd" -name $jobname -dest $base/pbs/$jobname.pbs -hour 12 -memG 100 -ppn 28
}

# cuffnorm
function rnaseq_cuffnorm() {
  return 1
}

# STAR
# make a genome index first
# STAR --runThreadN 28 --runMode genomeGenerate --genomeDir $WZSEQ_STAR_INDEX --genomeFastaFiles $WZSEQ_REFERENCE --sjdbGTFfile $WZSEQ_GTF --sjdbOverhang 49
function rnaseq_star() {
  if [[ -z samples ]]; then
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

#############
### other ###
#############

function wzseq_sra_to_fastq() {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d fastq ]] || mkdir fastq
  for f in sra/*.sra; do
    fn=$(readlink -f $f)
    bfn=$(basename $f)
    cmd="
fastq-dump -I --split-files $fn -O $base/fastq
gzip $base/fastq/*.fastq
"
    jobname="sra2fastq_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 10 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
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

# wzseq_qualimap <do>
function wzseq_qualimap() {
  
  base=$(pwd);
  [[ -d qualimap ]] || mkdir -p qualimap
  qualimapdir=$base/qualimap
  find bam -name *.bam |
    while read f; do
      [[ $f == *unmapped* ]] && continue
      fn=$(readlink -f $f)
      bfn=${f#bam/}
      bfn=${bfn%.bam}
      bfn=${bfn//\//_}
      cmd="
qualimap --java-mem-size=10G bamqc -nt 10 -bam $fn -outdir $qualimapdir/$bfn -c
"
      jobname="qualimap_"${f//\//_}
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 1 -memG 50 -ppn 10
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}
