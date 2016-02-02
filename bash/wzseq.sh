#!/bin/bash

# tips:
# pbsgen one "$cmds" -name $jobname -dest $pbsdir/$jobname $depend

################################################################################
# DNA-seq, mutation calling etc.
################################################################################

# samtools mpileup call mutation #
# check http://samtools.sourceforge.net/mpileup.shtml

# mpileup one sample
# usage: wzseq_mpileup1 bam/sname.bam
# one could use "vcfutils.pl varFilter -d10" to control vcf minimum depth
function wzseq_mpileup1 {

  base=$(pwd)
  bam=$1;
  sname=$(basename $bam .bam)
  [[ -d mpileup ]] || mkdir -p mpileup
  cmd="
cd $base
samtools mpileup -q 10 -uf $WZSEQ_REFERENCE $bam | bcftools call -cv -O b -o mpileup/$sname.raw.bcf
bcftools view mpileup/$sname.raw.bcf | vcf-sort -t mpileup/ | bgzip -c > mpileup/$sname.sorted.vcf.gz
"
  [[ ${!#} == "do" ]]  && qsub $pbsfn
}

# GATK best practice
function wzseq_seqtk_trimfq1 {
  [[ -d trimfq ]] || mkdir -p trimfq;
  [[ -d pbs ]] || mkdir -p pbs;
  $fq=$1;
  $fqname=${fq#fastq};
  $fqname=${fqname%.*};
  cmd="
seqtk trimfq $1 fastq/$fqname.fastq >trimfq/${fqname}_trimmed.fastq
"
  jobname="trimfq_$fqname"
  pbsfn="pbs/$jobname.pbs"
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 2 -ppn 1
  [[ ${!#} == "do" ]] && qsub $pbsfn
}

# wzseq_bwa_mem
function wzseq_bwa_mem {

  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmd="
cd $base
bwa mem -M -R \"@RG\tLB:$WZSEQ_REFVERSION\tID:${sname}\tPL:Illumina\tPU:hiseq2500\tSM:${sname}\" -t 28 $WZSEQ_BWA_INDEX $sread1 $sread2 | samtools sort -O bam -T bam/$sname.tmp -o bam/$sname.bam -
samtools index bam/$sname.bam
samtools flagstat bam/${sname}.bam > bam/$sname.bam.flagstat
"
      jobname="bwamem_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# indel realignment
function wzseq_GATK_indelrealign {

  [[ -d pbs ]] || mkdir -p pbs
  [[ -d indelrealn ]] || mkdir -p indelrealn
  while read sname sread1 sread2; do
    nthreads=28;
    interval=logs/${sname}_gatk_indels.intervals
    cmd="
# create interval
logs/$(sample_id)_gatk_indels.intervals:bams/$(sample_id)_rg_dedup.bam.bai bams/$(sample_id)_rg_dedup.bam
java -Xmx60g -Djava.io.tmpdir=./indelrealn/ -jar ~/software/GATK/GATK-3.3.0/GenomeAnalysisTK.jar --fix_misencoded_quality_scores --num_threads $nthreads -I bam/$bam -R $WZSEQ_REFERENCE -T RealignerTargetCreator --filter_mismatching_base_and_quals -o $interval --known $WZSEQ_GATK_KNOWN_INDEL

# actual realignment
java -Xmx2g -Djava.io.tmpdir=./indelrealn/ -jar ~/software/GATK/GATK-3.3.0/GenomeAnalysisTK.jar --fix_misencoded_quality_scores -o indelrealn/${sname}.indelrealn.bam -I bam/$bam -R $WZSEQ_REFERENCE -T IndelRealigner --filter_mismatching_base_and_quals -rf BadCigar -targetIntervals $interval --maxReadsForRealignment 200000 -known ${WZSEQ_GATK_KNOWN_INDEL}
"
    jobname="indelrealn_$sname"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
  # -known $(broad_files)/Mills_and_1000G_gold_standard.indels.hg19.vcf -known $(broad_files)/1000G_phase1.indels.hg19.vcf
}

################################################################################
# Whole Genome Bisulfite Sequencing
################################################################################

function examplepipeline_wgbs {
  cat <<- EOF
=== pipeline 2015-10-01 ===
[o] wgbs_adaptor => [o] wzseq_fastqc

 => (+) wgbs_biscuit_align / (+) wgbs_bwameth / (+) wgbs_bismark_bowtie1 / (+) wgbs_bismark_bowtie2 / wgbs_bsmap

 => (+) wgbs_biscuit_align_lambdaphage => (+) wgbs_biscuit_pileup_lambdaphage (TODO: exclude human reads)

 => (+) wzseq_GATK_realign (TODO: wgbs_indel_realign) => (+) wzseq_picard_markdup (TODO: wgbs_biscuit_markdup) => (+) wzseq_clean_intermediate => (+) TODO: wgbs_basequal_recal

 => [o] wzseq_merge_bam => (+) wzseq_qualimap => (defunct) wzseq_picard_WGSmetrics => (+) wzseq_bam_coverage

 => (+) wgbs_methpipe => (+) wgbs_methpipe_methylome => (+) wgbs_methpipe_allele

 => (+) wgbs_biscuit_pileup => (+) wgbs_vcf2tracks => (+) wgbs_cpgcoverage

 => (+) wgbs_methylKit_summary

 => (+) wgbs_diffmeth_simple => (+) wgbs_methylKit_diffmeth => (+) wgbs_methpipe_diff
EOF
}

# check adaptor conversion rate
function wgbs_adaptor() {
  local base=$(pwd);
  [[ -d adaptor ]] || mkdir adaptor;
  cmd="
cd $base
parallel -j 4 ~/wzlib/pyutils/wzadapter.py {} '>' adaptor/{/.}.txt ::: $base/fastq/*R1*.fastq.gz
for f in adaptor/*.txt; do
  tail -2 \$f | cut -d\":\" -f2 | cut -d\" \" -f2 | awk -v fn=\$f '{a[NR]=\$1}END{print fn,a[1],a[2]}'
done | awk 'BEGIN{print \"filename\tadaptorC\tadaptorT\tadaptorC2T\"}{print \$1,\$2,\$3,\$3/\$2*100\"%\"}' > adaptor/adaptor.stats
"
  jobname="adaptor_analysis"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 8 -ppn 4
  [[ ${!#} == "do" ]] && qsub $pbsfn
}

########################
## section1: alignment
########################

# biscuit alignment
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
  # requires base/samples with format:
  # sample_code fastq1 fastq2

  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      # while read samplecode fastq1 fastq2 _junk_; do
      cmd="
biscuit align $WZSEQ_BISCUIT_INDEX -t 28 $base/fastq/$sread1 $base/fastq/$sread2 | samtools sort -T $base/bam/${sname} -O bam -o $base/bam/${sname}.bam
samtools index $base/bam/${sname}.bam
samtools flagstat $base/bam/${sname}.bam > $base/bam/${sname}.bam.flagstat
"
      jobname="biscuit_align_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 100 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn

      # biscuit_bwameth -b -s $samplecode -j jid_bwameth $base/fastq/$fastq1 $base/fastq/$fastq2
      # biscuit_mdup -j jid_mdup -d $jid_bwameth -i bam/$samplecode
    done
  return
}

function wgbs_biscuit_align_lambdaphage {
  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmd="
biscuit align $WZSEQ_BISCUIT_INDEX_LAMBDAPHAGE -t 28 $base/fastq/$sread1 $base/fastq/$sread2 | samtools view -h -F 0x4 - | samtools sort -T $base/bam/${sname} -O bam -o $base/bam/${sname}_lambdaphage.bam
samtools index $base/bam/${sname}_lambdaphage.bam
samtools flagstat $base/bam/${sname}_lambdaphage.bam > $base/bam/${sname}_lambdaphage.bam.flagstat
"
      jobname="biscuit_align_lambdaphage_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# BWA-meth
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

  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmds="
cd $base
bwameth --reference $WZSEQ_BWAMETH_INDEX fastq/$sread1 fastq/$sread2 -t 28 --prefix bam/${sname}_bwameth
samtools flagstat bam/${sname}_bwameth.bam > bam/${sname}_bwameth.bam.flagstat
"
      jobname="bwameth_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmds" -name $jobname -dest $pbsfn -hour 24 -ppn 28 -memG 250
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# bismark with bowtie1
function wgbs_bismark_bowtie1_prepare_reference {
  base=$(pwd)
  cmd="
export PATH=~/tools/bismark/default:~/tools/bowtie1/default:$PATH
bismark_genome_preparation --bowtie1 $WZSEQ_BISMARK_BT1_INDEX
"
  jobname="bismark_bt1_prepare_reference"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
  [[ $1 == "do" ]] && qsub $pbsfn
}

function wgbs_bismark_bowtie1 {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmd="
export PATH=~/tools/bismark/default:~/tools/bowtie1/default:$PATH
cd $base
bismark $WZSEQ_BISMARK_BT1_INDEX --chunkmbs 2000 -1 $base/fastq/$sread1 -2 $base/fastq/$sread2
"
      jobname="bismark_bt1_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# bismark with bowtie2
function wgbs_bismark_bowtie2_prepare_reference {
  base=$(pwd)
  cmd="
export PATH=~/tools/bismark/default:~/tools/bowtie2/default:$PATH
bismark_genome_preparation --bowtie2 --verbose $WZSEQ_BISMARK_BT2_INDEX
"
  jobname="bismark_bt2_prepare_reference"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
  [[ $1 == "do" ]] && qsub $pbsfn
}

function wgbs_bismark_bowtie2 {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      base1=${sread1%.fastq.gz}
      base2=${sread2%.fastq.gz}
      base1=${base1%.fq.gz}
      base2=${base2%.fq.gz}
      cmd="
export PATH=~/tools/bismark/default:~/tools/bowtie2/default:$PATH
cd $base

# trim reads
# this produce fastq/${base1}_val_1* and fastq/${base2}_val_2*
if [[ ! -d fastq/${sname}_trim_galore/ ]]; then
  mkdir -p fastq/${sname}_trim_galore/
  ~/software/trim_galore/default/trim_galore --paired fastq/$sread1 fastq/$sread2 -o fastq/${sname}_trim_galore/
fi

bismark $WZSEQ_BISMARK_BT2_INDEX --bowtie2 --chunkmbs 2000 -p 4 -o $base/bam/${sname}_bismark_bt2/ -1 $base/fastq/${sname}_trim_galore/${base1}_val_1* -2 $base/fastq/${sname}_trim_galore/${base2}_val_2 --temp_dir $base/bam/${sname}_bismark_bt2/*
"
      jobname="bismark_bt2_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      # it seems that bismark with bowtie2 can never reach full potential of parallelization
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 60 -memG 50 -ppn 8
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# bsmap
function wgbs_bsmap {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam

  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmd="
cd $base
~/tools/bsmap/default/bsmap -a fastq/$sread1 -b fastq/$sread2 -d $WZSEQ_BSMAP_INDEX -o bam/${sname}_bsmap_unsorted.bam -p 28 -s 16 -v 10 -q 2 -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
samtools sort -O bam -o bam/${sname}_bsmap.bam -T bam/${sname}_bsmap.tmp bam/${sname}_bsmap_unsorted.bam
rm -f bam/${sname}_bsmap_unsorted.bam
samtools index bam/${sname}_bsmap.bam
samtools flagstat bam/${sname}_bsmap.bam > bam/${sname}_bsmap.bam.flagstat
"
      jobname="bsmap_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

#### rmapbs (doesn't allow gaps) ###
# rmapbs -c hg18 -o Human_NHFF.mr Human_NHFF.fastq
# rmapbs-pe -c hg18 -o Human_ESC.mr Human_ESC_1.fastq Human_ESC_2.fastq


########################
## section 2: analysis
########################

# summary
function wgbs_methylKit_summary {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d methylKit ]] || mkdir methylKit
  for f in pileup/*.vcf.gz; do
    bfn=$(basename $f .vcf.gz)
    cmd="
cd $base
mkdir methylKit/$bfn
biscuit vcf2bed -t cg -u -c $f | pybiscuit.py to_methylKit >methylKit/$bfn/$bfn.methylKit
~/wzlib/Rutils/bin/bioinfo/methylKit.r summary methylKit/$bfn/$bfn.methylKit -o methylKit/$bfn/
"
    jobname="methylKit_summary_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 50 -ppn 3
    [[ $1 == "do" ]] && qsub $pbsfn
  done

}

## TODO bissnp ####

function wgbs_merge_methlevelaverages {
  # usage: biscuit_merge_methlevelaverages data/methlevelaverage

  local basedir=$1;
  for f in $basedir/*.txt; do
    sed -e 's/://g' -e 's/%//' $f | awk -f wanding.awk -e 'BEGIN{split("",k);split("",n);split("",v);}(length($0)>0){k[length(k)+1]=$1;n[length(n)+1]=$2;v[length(v)+1]=$3/100;}END{for(i=1;i<=length(k);++i){kn[i]=k[i]"n"};print(join(k,1,length(k),"\t")"\t"join(kn,1,length(kn),"\t"));print(join(v,1,length(v),"\t")"\t"join(n,1,length(n),"\t"))}' > ${f%.txt}.processed
  done

  wzmanip concat -f $basedir/*.processed | awk -f wanding.awk -e '{n=NF;print;}END{repeat("NA", n, rep); print joina(rep,"\t");}' > $basedir/merge;
}

# note that I use my own version of to-mr which is agnostic of input bam type
# /home/wanding.zhou/software/methpipe/default/to-mr
# ~/software/methpipe/default/to-mr -o methpipe/$bfn -m bsmap 
# ~/software/methpipe/default/to-mr -o methpipe bam/Undetermined_bismark_bt2/Undetermined_L000_R1_001_val_1.fq.gz_bismark_bt2_pe.bam -m bismark
# ~/software/methpipe/default/to-mr -o methpi2 bam/Undetermined_bsmap.bam -m bsmap
# duplicate-remover -S 
function wgbs_methpipe {

  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d methpipe ]] || mkdir methpipe
  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bam)
    obfn=methpipe/$bfn
    cmd="
cd $base
pybiscuit.py to_mr -i $fn -o $obfn.mr
LC_ALL=C sort -T methpipe/ -k1,1 -k2,2n -k3,3n -k6,6 -o $obfn.mr.sorted_start $obfn.mr

echo \"[\$(date)] Remove duplicates...\"
~/tools/methpipe/default/bin/duplicate-remover -S ${obfn}_dremove_stat.txt -o $obfn.mr.dremove.untrim $obfn.mr.sorted_start

echo \"[\$(date)] Remove random contigs and haplotype...\"
cut -f1 $obfn.mr.dremove.untrim | uniq | while read chrm; do [[ -s $WZSEQ_REFERENCE_SPLIT/\$chrm.fa ]] || echo \"^\"\$chrm; done >$obfn.mr.uchr
grep -v -f $obfn.mr.uchr $obfn.mr.dremove.untrim > $obfn.mr.dremove

echo \"[\$(date)] Estimate bisulfite conversion rate...\"
~/tools/methpipe/default/bin/bsrate -c $WZSEQ_REFERENCE_SPLIT -o $obfn.bsrate $obfn.mr.dremove
LC_ALL=C sort -T methpipe/ -k1,1 -k3,3n -k2,2n -k6,6 -o $obfn.mr.sorted_end_first $obfn.mr.dremove

echo \"[\$(date)] Compute single-site methylation levels...\"
~/tools/methpipe/default/bin/methcounts -c $WZSEQ_REFERENCE_SPLIT -o $obfn.meth $obfn.mr.sorted_end_first

echo \"[\$(date)] Merge C,G in CpG context...\"
~/tools/methpipe/default/bin/symmetric-cpgs -m -o $obfn.CpG.meth $obfn.meth

echo \"[\$(date)] Get statistics of methylation levels...\"
~/tools/methpipe/default/bin/levels -o $obfn.levels $obfn.meth

echo \"[\$(date)] Clean up intermediate files...\"
rm -f $obfn.mr $obfn.mr.sorted_start $obfn.mr.dremove.untrim $obfn.mr.uchr 
# not sure if these 2 should be removed as intermediate $obfn.mr.dremove $obfn.mr.sorted_end_first
echo \"[\$(date)] Done.\"
"
    jobname="methpipe_basic_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function wgbs_methpipe_methylome {
  base=$(pwd)
  [[ -d methpipe ]] || mkdir methpipe
  for f in methpipe/*.mr; do
    bfn=$(basename $f .mr)
    obfn=methpipe/$bfn
    cmd="
cd $base

echo \"[\$(date)] Hypomethylated (primarily) and hypermethylated region...\"
~/tools/methpipe/default/bin/hmr -p $obfn.hmr.params -o $obfn.hmr $obfn.CpG.meth

# use trained parameter
# echo \"[\$(date)] Train parameters...\"
# ~/tools/methpipe/default/bin/hmr -p $obfn.hmr.params -o $obfn.hmr $obfn.CpG.meth
# echo \"[\$(date)] Use trained parameter to call HMR...\"
# ~/tools/methpipe/default/bin/hmr -P $obfn.hmr.params -o $obfn.hmr $obfn.CpG.meth

echo \"[\$(date)] Partially methylated regions (PMR)...\"
# not sure what this is
~/tools/methpipe/default/bin/hmr -partial -o $obfn.hmr $obfn.CpG.meth

# echo \"[\$(date)] Hypermethylated region in plant (HyperMR)...\"
# ~/tools/methpipe/default/bin/hypermr -o $obfn.hypermr $obfn.CpG.meth

echo \"[\$(date)] Partially methylated domain (PMD, 1kb bin level segmentation).\"
~/tools/methpipe/default/bin/pmd -o $obfn.pmd $obfn.CpG.meth

echo \"[\$(date)] Done.\"
"
    jobname="methpipe_methylome_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function wgbs_methpipe_allele {
  base=$(pwd)
  [[ -d methpipe ]] || mkdir methpipe
  for f in methpipe/*.mr; do
    bfn=$(basename $f .mr)
    obfn=methpipe/$bfn
    cmd="
cd $base

echo \"[\$(date)] Make read distribution, (.epiread)\"
~/tools/methpipe/default/bin/methstates -c $WZSEQ_REFERENCE_SPLIT -o ${obfn}.epiread ${obfn}.mr.dremove

echo \"[\$(date)] Single site ASM scoring (.allelic)\"
~/tools/methpipe/default/bin/allelicmeth -c $WZSEQ_REFERENCE_SPLIT -o ${obfn}.allelic ${obfn}.epiread

echo \"[\$(date)] Allelically methylated region (AMRs)\"
~/tools/methpipe/default/bin/amrfinder -o $obfn.amr -c $WZSEQ_REFERENCE_SPLIT $obfn.epiread

# echo \"[\$(date)] Test allele-specific methylation at promoter of imprinted region\"
# amrtester -o ${obfn}.amr -c hg19 target_interval.bed ${obfn}.epiread

echo \"[\$(date)] Calculate meth-entropy\"
# -v verbose
~/tools/methpipe/default/bin/methentropy -w 5 -v -o $obfn.entropy $WZSEQ_REFERENCE_SPLIT $obfn.epiread

echo \"[\$(date)] Done.\"
"
    jobname="methpipe_allele_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function wgbs_methpipe_diff {
  # assuming the corresponding methpipe run (mr file) has finished
  # for each listed vcf file in the corresponding methpipe folder
  base=$(pwd)
  awk '/^\[/{p=0}/\[diffmeth\]/{p=1;next} p&&!/^$/' samples |
    while read sname vcfpath1 vcfpath2; do
      echo 1
    done
}

# TODO: the following gives seg fault now
function wgbs_biscuit_markdup {

  # this differentiate strands (parent/daughter) when applied to WGBS.
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam/before_mdup ]] || mkdir -p bam/before_mdup
  for f in bam/*.bam; do
    bfn=$(basename $f .bam)
    cmd="
cd $base
mv $f bam/before_mdup/$bfn.bam
[[ -e $f.bai ]] && mv $f.bai bam/before_mdup/$bfn.bam.bai
[[ -e $f.flagstat ]] && mv $f.flagstat bam/before_mdup/$bfn.bam.flagstat
biscuit markdup -q bam/before_mdup/$bfn.bam $f 2>bam/before_mdup/$bfn.mdup.stats
samtools index $f
samtools flagstat $f
"
    jobname="biscuit_markdup_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 60 -memG 50 -ppn 5
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function wgbs_biscuit_pileup() {
  # wgbs_biscuit_pileup [-nome] do

  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d pileup ]] || mkdir pileup
  [[ $1 == "-nome" ]] && nome="-N" || nome=""; # whether to pileup using the nomeseq mode
  for f in bam/*.bam; do
    [[ $f =~ "lambdaphage" ]] && continue;
    fn=$(readlink -f $f)
    bfn=$(basename $f .bam)
    cmd="
cd $base
biscuit pileup -r $WZSEQ_REFERENCE $nome -i $fn -o pileup/$bfn.vcf -q 28
bgzip pileup/$bfn.vcf
tabix -p vcf pileup/$bfn.vcf.gz
"
    jobname="biscuit_pileup_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 28
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done

}

function wgbs_biscuit_pileup_lambdaphage() {

  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d pileup ]] || mkdir pileup
  for f in bam/*_lambdaphage.bam; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bam)
    cmd="
cd $base
biscuit pileup -r $WZSEQ_REFERENCE_LAMBDAPHAGE -i $fn -o pileup/$bfn.vcf -q 28
bgzip pileup/$bfn.vcf
tabix -p vcf pileup/$bfn.vcf.gz
"
    jobname="biscuit_lambdaphage_pileup_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 28
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

# wgbs_vcf2tracks [-nome] do
function wgbs_vcf2tracks {

  # compute 1) base-pair resolution methylation track 2) mean methylation in window track
  [[ -d tracks ]] || mkdir tracks
  [[ -d qc ]] || mkdir qc
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ $1 == "-nome" ]] && items=(hcg gch) || items=(cg)
  for f in pileup/*.vcf.gz; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .vcf.gz)
    cmd="
cd $base

for pt in ${items[@]}; do

  echo processing pileup type \$pt >&2

  biscuit vcf2bed -t \${pt} pileup/${bfn}.vcf.gz | LC_ALL=C sort -k1,1 -k2,2n -T tracks > tracks/${bfn}.\${pt}.bedg
  bedGraphToBigWig tracks/${bfn}.\${pt}.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}.\${pt}.bw

  for i in 100000,100k 1000,1k 100,100; do 
    IFS=\",\"; set \$i; 
    echo processing \$i >&2
    bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w \$1 | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ | bedtools intersect -a - -b tracks/${bfn}.\${pt}.bedg -wo | bedtools groupby -i - -g 1-3 -c 7 -o mean >tracks/${bfn}.\${pt}.window\$2.bedg;

    bedGraphToBigWig tracks/${bfn}.\${pt}.window\$2.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}.\${pt}.window\$2.bw
    rm -f tracks/${bfn}.\${pt}.window\$2.bedg
    unset IFS
  done
done
"
    jobname="vcf2tracks_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done

}

# obsolete make tracks from USC format beds
function wgbs_methcallbed_to_tracks {

  [[ -d tracks ]] || mkdir tracks;
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  for f in bed/*.bed; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bed)
    cmd="
cd $base
awk 'NR>1&&\$8>3{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$7}' $fn | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ >tracks/${bfn}_cov3.bedg
bedGraphToBigWig tracks/${bfn}_cov3.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}_cov3.bw

bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w 100000 | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ | bedtools intersect -a - -b tracks/${bfn}_cov3.bedg -wo | bedtools groupby -i - -grp 1-3 -c 7 -o mean >tracks/${bfn}_cov3_window100k.bedg;
bedGraphToBigWig tracks/${bfn}_cov3_window100k.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}_cov3_window100k.bw
rm -f tracks/${bfn}_cov3_window100k.bedg

bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w 100 | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ | bedtools intersect -a - -b tracks/${bfn}_cov3.bedg -sorted -wo | bedtools groupby -i - -grp 1-3 -c 7 -o mean >tracks/${bfn}_cov3_window100.bedg;
bedGraphToBigWig tracks/${bfn}_cov3_window100.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}_cov3_window100.bw
rm -f tracks/${bfn}_cov3_window100.bedg

bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w 1000 | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ | bedtools intersect -a - -b tracks/${bfn}_cov3.bedg -sorted -wo | bedtools groupby -i - -grp 1-3 -c 7 -o mean >tracks/${bfn}_cov3_window1k.bedg;
bedGraphToBigWig tracks/${bfn}_cov3_window1k.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}_cov3_window1k.bw
rm -f tracks/${bfn}_cov3_window1k.bedg
"
    jobname="methcallbed_to_tracks_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function wgbs_cpgcoverage {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d cpg ]] || mkdir cpg
  for f in pileup/*.vcf.gz; do
    bfn=$(basename $f .vcf.gz)
    cmd="
cd $base
biscuit vcf2bed -t cg -k 0 -c $f > cpg/$bfn.cg.bedg
wzplot cumhist -t cpg/$bfn.cg.bedg -c 5 --xlabel \"coverage\" --ylabel \"cumulative distribution\" -o cpg/cpg_coverage.$bfn.cumhist.png
bedtools intersect -a cpg/$bfn.cg.bedg -b $WZSEQ_CGIBED | wzplot cumhist -t - -c 5 --xlabel \"coverage\" --ylabel \"cumulative distribution\" -o cpg/cpg_coverage_cpgisland.$bfn.cumhist.png

wzplot hist --maxline 1000000000 -t cpg/$bfn.cg.bedg -c 4 --xlabel \"beta values\" -o cpg/beta_dist.$bfn.hist.png
bedtools intersect -a cpg/$bfn.cg.bedg -b $WZSEQ_CGIBED | wzplot hist --maxline 1000000000 -t - -c 4 --xlabel \"beta values\" -o cpg/beta_dist_cpgisland.$bfn.hist.png

awk '\$5>5' cpg/$bfn.cg.bedg | bedtools intersect -a - -b $WZSEQ_CGIBED -wo | sort -k6,6 -k7,7n | bedtools groupby -i - -g 6-8 -c 4 -o mean | wzplot hist --maxline 1000000000 -t - -c 4 --xlabel \"methlevelaverage per CGI\" -o cpg/methlevelaverage_b5_dist_cpgisland.$bfn.hist.png
rm -f cpg/$bfn.cg.bedg
"
    jobname="cpg_coverage_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wgbs_methylKit_diffmeth {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d methylKit ]] || mkdir methylKit
  awk '/^\[/{p=0}/\[diffmeth\]/{p=1;next} p&&!/^$/' samples |
    while read sname vcfpath1 vcfpath2; do
      echo $sname $vcfpath1
      base1=$(basename $vcfpath1 .vcf.gz) # tumor
      base2=$(basename $vcfpath2 .vcf.gz) # normal
      [[ $base1 == $base2 ]] && base2=$base2"_1"
      dir1=$(dirname $vcfpath1)
      dir2=$(dirname $vcfpath2)
      dir1=${dir1%/pileup}
      dir2=${dir2%/pileup}
      cmd="
cd $base
mkdir methylKit/$sname
biscuit vcf2bed -t cg -u -c $vcfpath1 | pybiscuit.py to_methylKit >methylKit/$sname/$base1.methylKit
biscuit vcf2bed -t cg -u -c $vcfpath2 | pybiscuit.py to_methylKit >methylKit/$sname/$base2.methylKit
~/wzlib/Rutils/bin/bioinfo/methylKit.r diff -t 1,0 -b $base1,$base2 methylKit/$sname/$base1.methylKit,methylKit/$sname/$base2.methylKit -o methylKit/$sname -g $WZSEQ_CGIBED_METHYLKIT
rm -f methylKit/$sname/$base1.methylKit
rm -f methylKit/$sname/$base2.methylKit
"
      jobname="methylKit_diff_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 50 -ppn 3
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function wgbs_diffmeth_simple {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d diffmeth ]] || mkdir diffmeth
  awk '/^\[/{p=0}/\[diffmeth\]/{p=1;next} p&&!/^$/' samples |
    while read sname vcfpath1 vcfpath2; do
      cmd="
cd $base
[[ -d diffmeth/$sname ]] || mkdir -p diffmeth/$sname
~/wzlib/bash/wzmethdiff.sh -t $vcfpath1 -n $vcfpath2 -b $base/diffmeth/$sname -v $WZSEQ_REFVERSION -g $WZSEQ_CGIBED -r $WZSEQ_REFERENCE -s $WZSEQ_TSSBED 2>$base/diffmeth/$sname/run.log
"
      jobname="simple_methdiff_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 5 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

################################################################################
# ChIP-seq pipeline
################################################################################

function examplepipeline_chipseq {
  cat <<- EOF
=== pipeline 2015-01-28 ===
 (+) wzseq_bwa_aln_se / (+) wzseq_bwa_mem =>

 (+) chipseq_bcp => (+) chipseq_macs2 => (+) chipseq_gem => (+) chipseq_sissrs
 (+) chipseq_HOMER => (+) chipseq_HOMER_peak

EOF
}

# BWA-aln single-ended
function wzseq_bwa_aln_se {

  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam;
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread _junk_; do
      sfile=$(readlink -f fastq/$sread);
      cmd="
cd $base
bwa aln -t 10 $WZSEQ_BWA_INDEX $sfile >bam/${sname}.sai
bwa samse $WZSEQ_BWA_INDEX bam/${sname}.sai $sfile | samtools view -bS - | samtools sort -T $base/bam/$sname -O bam -o bam/${sname}.bam
samtools index bam/${sname}.bam
samtools flagstat bam/${sname}.bam > bam/${sname}.bam.flagstat
rm -f bam/${sname}.sai
"
      jobname="bwa_aln_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 10
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# http://homer.salk.edu/homer/
# got to install seqlogo from webLogo
# http://weblogo.berkeley.edu/
# the sequences can be obtained from
# perl configureHomer.pl -install mm10
# sequence and annotation gets downloaded to where you install HOMER
function chipseq_HOMER {
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam;
  [[ -d HOMER ]] || mkdir HOMER;
  for f in bam/*.bam; do
    bfn=$(basename $f .bam)
    cmd="
cd $base
echo Creating Tags for HOMER/$bfn
~/software/HOMER/default/bin/makeTagDirectory HOMER/$bfn $f
"
    jobname="HOMER_$(basename $base)"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 2 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function chipseq_HOMER_peak {
  ## make sure you have seqlogo, gs, blat and samtools available from command line.
  ## targettype: factor, histone, super, groseq, tss, dnase, mC
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[peak\]/{p=1;next} p&&!/^$/' samples |
    while read targetname controlname targettype; do

      if [[ $targettype == "factor" ]]; then
        targetsize=100;
      else
        targetsize="given"
      fi

      cmd="
cd $base
~/software/HOMER/default/bin/findPeaks HOMER/$targetname -style $targettype -o auto -i HOMER/$controlname
export PATH=$PATH:~/software/HOMER/default/bin/:~/software/webLogo/weblogo/
~/software/HOMER/default/bin/findMotifsGenome.pl HOMER/$targetname/peaks.txt $WZSEQ_REFVERSION HOMER/${targetname}/Motif -size $targetsize # -mask? 
~/software/HOMER/default/bin/annotatePeaks.pl HOMER/$targetname/peaks.txt $WZSEQ_REFVERSION > HOMER/$targetname/peaks.txt.annotation
"
      # TODO: there are some more plotting methods to add
      # see http://homer.salk.edu/homer/ngs/quantification.html
      jobname="HOMER_peak_$targetname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 4 -ppn 2
      [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function chipseq_bcp() {

  base=$(pwd);
  [[ -d bcp ]] || mkdir bcp
  [[ -d pbs ]] || mkdir pbs

  awk '/^\[/{p=0}/\[peak\]/{p=1;next} p&&!/^$/' samples |
    while read targetname controlname targettype; do

      tfile=$(readlink -f bam/$targetname);
      cfile=$(readlink -f bam/$controlname);

      mkdir bcp/$targetname;
      tbed=bcp/$targetname/$targetname.bcp.bed
      cbed=bcp/$targetname/$controlname.bcp.bed

      if [[ $targettype == "factor" ]]; then
        runcmd1="~/tools/bcp/BCP_v1.1/BCP_TF -1 $tbed -2 $cbed -3 bcp/$targetname/peaks.tfbs"
      fi

      if [[ $targettype == "histone" ]]; then
        runcmd1="~/tools/bcp/BCP_v1.1/BCP_HM -1 $tbed -2 $cbed -3 bcp/$targetname/peaks.hm";
      fi
      
      cmd="
cd $base
wzbam bed6 -bam bam/$targetname.bam -o $tbed
wzbam bed6 -bam bam/$controlname.bam -o $cbed

$runcmd1

rm -f $tbed $cbed
"
      jobname="bcp_$targetname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function chipseq_macs2 {

  base=$(pwd);
  [[ -s macs2 ]] || mkdir macs2
  [[ -d pbs ]] || mkdir pbs

  awk '/^\[/{p=0}/\[peak\]/{p=1;next} p&&!/^$/' samples |
    while read targetname controlname targettype; do

      cmd="
cd $base
macs2 callpeak -t bam/$targetname.bam -c bam/$controlname.bam -f BAM -g $WZSEQ_MACS_SHORT -n macs2/$targetname -B
"
      jobname="macs2_$targetname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function chipseq_gem {

  base=$(pwd);
  [[ -d gem ]] || mkdir -p gem
  [[ -d pbs ]] || mkdir -p pbs

  awk '/^\[/{p=0}/\[peak\]/{p=1;next} p&&!/^$/' samples |
    while read targetname controlname targettype; do

      [[ -d gem/$targetname ]] || mkdir gem/$targetname;
      tbed=gem/$targetname/$targetname.bed
      cbed=gem/$targetname/$controlname.bed
      cmd="
cd $base
wzbam bed6 -bam bam/$targetname.bam -o $tbed
wzbam bed6 -bam bam/$controlname.bam -o $cbed

java -Xmx100G -jar ~/tools/gem/gem/gem.jar --d ~/tools/gem/gem/Read_Distribution_default.txt --g $WZSEQ_REFERENCE.fai --genome $WZSEQ_REFERENCE_SPLIT --s 2000000000 --expt $tbed --ctrl $cbed --f BED --out gem/$targetname --k_min 6 --k_max 13

rm -f $tbed $cbed
"
      jobname="gem_$targetname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 50 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# http://dir.nhlbi.nih.gov/papers/lmi/epigenomes/sissrs/SISSRs-Manual.pdf
# from Keji Zhao lab, allows both with and without background
# WZSEQ_REFERENCE_SIZE can come from the last row of .fai
function chipseq_sissrs {

  [[ -d sissrs ]] || mkdir sissrs
  [[ -d pbs ]] || mkdir pbs
  base=$(pwd)
  awk '/^\[/{p=0}/\[peak\]/{p=1;next} p&&!/^$/' samples |
    while read targetname controlname targettype; do

      [[ -d sissrs/$targetname ]] || mkdir sissrs/$targetname;
      tbed=sissrs/$targetname/$targetname.bed
      cbed=sissrs/$targetname/$controlname.bed
      cmd="
cd $base
wzbam bed6 -bam bam/$targetname.bam -o $tbed
wzbam bed6 -bam bam/$controlname.bam -o $cbed
~/software/SISSERs/v1.4/sissrs.pl -i $tbed -b $cbed -o sissrs/$targetname/$targetname.bsites -s $WZSEQ_REFERENCE_SIZE
rm -f $tbed $cbed
"
      jobname="sissrs_$targetname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 20 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

################################################################################
# HiC pipeline
################################################################################
# TODO: http://homer.salk.edu/homer/interactions/HiCtagDirectory.html

################################################################################
# RNA-seq pipeline
################################################################################

function examplepipeline_rnaseq {
cat<<- EOF
=== pipeline 2015-10-01 ===
(+) rnaseq_tophat2_firststrand / rnaseq_tophat2 => (+) rnaseq_cufflinks => (+) rnaseq_cuffmerge => (+) rnaseq_cuffquant (optional) => (+) rnaseq_cuffdiff
EOF

cat<<- EOF
=== pipeline 2016-01-12 ===
 (+) wzseq_fastqc

 => (+) edit samples [alignment]; (+) editrnaseq_tophat2_firststrand / (+) rnaseq_tophat2 / (+) rnaseq_STAR / (+) rnaseq_gsnap / (+) rnaseq_subjunc

 => (+) wzseq_qualimap rnaseq_se/rnaseq_pe_stranded

 => [o] rnaseq_splitallele (for ASE) => (+) wzseq_bam_coverage

 => (+) rnaseq_cufflinks => (+) rnaseq_cuffmerge => (+) edit samples [diffexp] ; rnaseq_cuffdiff

 => (+) rnaseq_edgeR => (+) rnaseq_DESeq2

 => (+) rnaseq_allelomePro
EOF
}

#########################
## section 1: alignment
#########################
# rnaseq_tophat2
function rnaseq_tophat2() {
  # single-end reads use "." as placeholder
  if [[ ! -s samples ]]; then
    echo "No [samples] file. Abort."
    return 1;
  fi
  base=$(pwd);
  [[ -d bam ]] || mkdir -p bam
  [[ -d pbs ]] || mkdir -p pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      sfile1=$(readlink -f fastq/$sread1);
      [[ $sread2 == "." ]] && sfile2="" || sfile2=$(readlink -f fastq/$sread2); # in case of single ended.
      odir=$(readlink -f bam/$sname);
      # customize the following
      # default: no novel junction, 28 threads
      cmd="
cd $base
tophat2 -p 28 -G $WZSEQ_GTF --library-type fr-unstranded -o $odir --no-novel-juncs $WZSEQ_BOWTIE2_INDEX $sfile1 $sfile2
cd $base/bam
ln -s $sname/accepted_hits.bam $sname.bam
samtools index $sname.bam
samtools flagstat $sname.bam >$sname.bam.flagstat
"
      jobname="tophat_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# rnaseq_tophat2_firststrand samplefn
# first strand
function rnaseq_tophat2_firststrand() {
  if [[ ! -s samples ]]; then
    echo "file: samples missing. Abort."
    return 1;
  fi
  base=$(pwd);
  [[ -d bam ]] || mkdir -p bam
  [[ -d pbs ]] || mkdir -p pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      sfile1=$(readlink -f fastq/$sread1);
      sfile2=$(readlink -f fastq/$sread2);
      odir=$(readlink -f bam/$sname);
      cmd="
tophat2 -p 28 -G $WZSEQ_GTF --library-type fr-firststrand -o $odir --no-novel-juncs $WZSEQ_BOWTIE2_INDEX $sfile1 $sfile2
cd $base/bam
ln -s $sname/accepted_hits.bam $sname.bam
samtools index $sname.bam
samtools flagstat $sname.bam >$sname.bam.flagstat
"
      jobname="tophat_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# STAR
# make a genome index first
# STAR --runThreadN 28 --runMode genomeGenerate --genomeDir $WZSEQ_STAR_INDEX --genomeFastaFiles $WZSEQ_REFERENCE --sjdbGTFfile $WZSEQ_GTF
# STAR --runThreadN 28 --runMode genomeGenerate --genomeDir . --genomeFastaFiles mm10.fa --sjdbGTFfile ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming [--sjdbOverhang 49 (default 100)]
# note that the GTF must be unzipped!!
# 
# mapping quality:
# 255 = uniquely mapped reads
# 3 = read maps to 2 locations
# 2 = read maps to 3 locations
# 1 = reads maps to 4-9 locations
# 0 = reads maps to 10 or more locations
function rnaseq_STAR {

  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmd="
cd $base
mkdir $base/bam/$sname
~/software/STAR/default/source/STAR --runThreadN 28 --genomeDir $WZSEQ_STAR_INDEX --readFilesIn $base/fastq/$sread1 $base/fastq/$sread2 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $base/bam/$sname/$sname --outBAMsortingThreadN 6
cd bam
ln -s $sname/${sname}Aligned.sortedByCoord.out.bam $sname.bam
samtools index $sname.bam
samtools flagstat $sname.bam >$sname.bam.flagstat
"
      jobname="STAR_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# gsnap require some prefix setting at compilation, e.g., ./configure --prefix /home/wanding.zhou/software/gsnap/gmap-2015-12-31/
# indexing: ~/software/gsnap/default/util/gmap_build -D ~/references/mm10/gsnap -k 15 -d mm10 ~/references/mm10/bowtie1/*.fa
# note: -d gives the name of the reference, -D gives the location where one looks for the reference with the name (supplied from -d)
# 
# make splice sites and intron sites
# cat ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming | ~/software/gsnap/default/bin/gtf_splicesites > ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gsnap.splicesites
# cat ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming | ~/software/gsnap/default/bin/gtf_introns > ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gsnap.introns
# 
# make splice sites and intron sites from UCSC refGene
# gunzip -c refGene.txt.gz | psl_splicesites -s 1 > foo.splicesites # -s 1 means skip 1 column from the left
# gunzip -c refGene.txt.gz | psl_introns -s 1 > foo.introns

# make .iit files
# cat ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gsnap.introns ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gsnap.splicesites | ~/software/gsnap/default/bin/iit_store -o ~/references/mm10/gsnap/mm10/mm10.maps/Mus_musculus.GRCm38.82.gtf.gsnap.splicesites.iit
# make sure those files are stored within the directory of mm10/gsnap/mm10/mm10.maps/
# WZSEQ_GSNAP_SPLICE=Mus_musculus.GRCm38.82.gtf.gsnap.splicesites

function rnaseq_gsnap {
  # GMAP, the original program is designed for mapping cDNA and EST sequences.
  # GSNAP is the derivative for mapping shot-gun sequencing reads
  # recommend nthreads not exceed 4, otherwise performance is compromised. program design problem
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      
      cmd="
cd $base
~/software/gsnap/default/src/gsnap -A sam -t 4 --gunzip --novelsplicing=1 --use-splicing=$WZSEQ_GSNAP_SPLICE -D $WZSEQ_GSNAP_INDEX -d $WZSEQ_REFVERSION $base/fastq/$sread1 $base/fastq/$sread2 -o bam/$sname.sam
samtools sort -O bam -T bam/$sname.tmp -o bam/$sname.bam bam/$sname.sam
samtools index bam/$sname.bam
#rm -f bam/$sname.sam
"
      jobname="GSNAP_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 80 -memG 50 -ppn 4
      [[ $1 == "do" ]] && qsub $pbsfn
    done
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
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
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
    done
}

#######################################
## section 2: differential expression
#######################################

# rnaseq_cufflinks <do>
function rnaseq_cufflinks() {
  [[ -d cufflinks ]] || mkdir -p cufflinks;
  base=$(pwd);
  for sample in bam/*.bam; do
    sample=$(basename $sample .bam)
    # [[ -s bam/$sample.bam ]] || continue
    cmd="cufflinks -p 28 -g $WZSEQ_GTF -o $base/cufflinks/$sample $base/bam/$sample.bam -q"
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
  jobname='cuffmerge_'$(basename base)
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

# require [diffexp] section in "samples"
# format: cond1\tcond2\tbams1\tbams2
# bams1 and bams2 can be comma separated if containing multiple samples
# e.g.
# [diffexp]
# aza     PBS     bam/aza.bam     bam/PBS.bam
# 
# for cuffdiff >2.2 bams can be replaced by cxb files.
function rnaseq_cuffdiff() {
  base=$(pwd);
  [[ -d cuffdiff ]] || mkdir -p cuffdiff;

  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
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
    done
}

# TODO: cuffnorm
function rnaseq_cuffnorm() {
  return 1
}

# TODO: kallisto
function rnaseq_edgeR {
  base=$(pwd)
  [[ -d edgeR ]] || mkdir edgeR
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      # use merged gtf if available otherwise, use back-up gtf
      gtf=$base/cuffmerge/merged_asm/merged.gtf
      [[ -s $gtf ]] || gtf=$WZSEQ_GTF
      cmd="
cd $base
~/wzlib/Rutils/bin/bioinfo/edgeR.r -g $WZSEQ_REFVERSION -G $WZSEQ_GTF_ENSEMBL -a $cond1 -b $cond2 -A $bams1 -B $bams2 -o edgeR/${cond1}_vs_${cond2}_diffexp.tsv 2> edgeR/${cond1}_vs_${cond2}_diffexp.log
"
      jobname="edgeR_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 20 -ppn 2
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function rnaseq_DESeq2 {
  base=$(pwd)
  [[ -d DESeq2 ]] || mkdir DESeq2
  [[ -d pbs ]] || mkdir pbs
  grep '\[experiment\] stranded' samples && stranded="" || stranded="--ignoreStrand"
  grep '\[experiment\] single-end' samples && singleEnd="--singleEnd" || singleEnd=""
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      # use merged gtf if available otherwise, use back-up gtf
      gtf=$base/cuffmerge/merged_asm/merged.gtf
      [[ -s $gtf ]] || gtf=$WZSEQ_GTF
      cmd="
cd $base
~/wzlib/Rutils/bin/bioinfo/DESeq2.r $stranded $singleEnd -g $WZSEQ_REFVERSION -G $WZSEQ_GTF_ENSEMBL -a $cond1 -b $cond2 -A $bams1 -B $bams2 -o DESeq2/${cond1}_vs_${cond2}_diffexp.tsv 2> DESeq2/${cond1}_vs_${cond2}_diffexp.log
"
      jobname="DESeq2_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 36 -memG 100 -ppn 14
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# TODO: DESeq2, EBSeq, Voom

# RSeQC
# need to define WZSEQ_RSEQ_GENE_BED
function rnaseq_rseqc() {
  base=$(pwd);
  [[ -d rseqc ]] || mkdir rseqc
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
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do

      if [[ -s bam/$sname.bam ]]; then
        echo "Bam bam/$sname.bam exists. skip"
        continue
      fi
      
      cmd="
cd $base
[[ -d bam/$sname ]] || mkdir bam/$sname
/primary/home/wanding.zhou/tools/subread/subread-1.4.6-p5-Linux-x86_64/bin/subjunc -T 28 -I 16 -i $WZSEQ_SUBREAD_INDEX -r fastq/$sread1 -R fastq/$sread2 --gzFASTQinput -o bam/$sname/$sname.unsrtd.bam --BAMoutput
samtools sort -O bam -o bam/${sname}.bam -T bam/$sname/${sname}.tmp bam/$sname/${sname}.unsrtd.bam
samtools index bam/${sname}.bam
samtools flagstat bam/${sname}.bam >bam/${sname}.flagstat
#rm -f bam/$sname/$sname.unsrtd.bam
"
      jobname="subjunc_${sname}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function rnaseq_count_rmsk {
  bam=$1
  cmd=$(cat<<EOF
cd $(pwd)
bedtools coverage -a ~/references/hg19/annotation/rmsk.txt.bed -b bam/${bam}.bam -sorted -split -counts > repeatmasker/$bam.rmsk
awk '{a[\$5]+=\$8;b[\$6]+=\$8;c[\$7]+=\$8}END{for(i in a){print "1\t"i"\t"a[i]}; for(i in b){print "2\t"i"\t"b[i]}; for(i in c){print "3\t"i"\t"c[i]}}' repeatmasker/${bam}.rmsk | sort -k1,1 -k2,2nr > repeatmasker/${bam}.rmsk.cat
EOF
     )

  jobname="rmsk_${bam}"
  pbsfn=pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 2 -ppn 1
}

# TODO: featureCounts (from subread package)
function rnaseq_featurecounts {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bam)
    cmd="
featureCounts -T 28 -t exon -g gene_id -a annotation.gtf -o featureCounts/$bfn_counts.txt $fn
~/tools/subread/subread-1.4.6-p5-Linux-x86_64/bin/featureCounts -a repeatmasker/rmsk7.bed -F SAF -o repeatmasker/aza bam/aza.bam -T 28
~/tools/subread/subread-1.4.6-p5-Linux-x86_64/bin/featureCounts -t exon -g gene_id -a ~/references/hg19/gtf/Homo_sapiens.GRCh37.75.gtf -o 1 bam/aza.bam
"
    jobname="featurecounts_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

# TODO: MISO

# TODO: HTseq-count

# TODO: RSEM
# ~/tools/rsem/rsem-1.2.22/rsem-prepare-reference
# rsem-calculate-expression -p 20 --calc-ci --ci-memory 12294 --bowtie-chunkmbs 2000 --paired-end --bowtie-path $WZSEQ_BOWTIE1 --rsem-index RSEM

##########################################
# section 3: allele-specific expression
##########################################

# wzbam splitallele -bam bam/Dot1LE10.bam -snp ~/references/mm10/snp/pairwise/BL6_vs_JF1.snp.mm10.bed -name BL6,JF1 -o testsplit

function rnaseq_splitallele {
  # [split_allele]
  # PL_female_wt.bam        /primary/vari/genomicdata/genomes/mm10/snp/pairwise/BL6_vs_JF1.snp.mm10.bed     BL6,JF1
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam_allele ]] || mkdir bam_allele
  awk '/^\[/{p=0}/\[split_allele\]/{p=1;next} p&&!/^$/' samples |
    while read bamfile snpfile background; do
      bfn=$(basename $bamfile .bam)
      cmd="
cd $base
wzbam splitallele -bam bam/$bamfile -snp $snpfile -name $background -o bam_allele/$bfn 2>bam_allele/$bfn.splitallele.log
IFS=',' read -r -a samples <<< \"$background\"
sub1=bam_allele/$bfn.\${samples[0]}.bam
sub2=bam_allele/$bfn.\${samples[1]}.bam
samtools index \$sub1
samtools index \$sub2
samtools flagstat \$sub1 >\$sub1.flagstat
samtools flagstat \$sub2 >\$sub2.flagstat

# make coverage tracks
minmapq=10
for sub in \$sub1 \$sub2; do
  subbase=\$(basename \$sub .bam)
  bedtools genomecov -ibam \$sub -g ${WZSEQ_REFERENCE}.fai -bga -split | LC_ALL=C sort -k1,1 -k2,2n -T bam_allele/  >bam_allele/\${subbase}.coverage.bedg
  bedGraphToBigWig bam_allele/\${subbase}.coverage.bedg ${WZSEQ_REFERENCE}.fai bam_allele/\${subbase}.coverage.bw
  rm -f bam_allele/\${subbase}.coverage.bedg

  samtools view -q \$minmapq -b \$sub | bedtools genomecov -ibam stdin -g ${WZSEQ_REFERENCE}.fai -bga -split | LC_ALL=C sort -k1,1 -k2,2n -T bam_allele/  >bam_allele/\${subbase}.coverage.q10.bedg
  bedGraphToBigWig bam_allele/\${subbase}.coverage.q10.bedg ${WZSEQ_REFERENCE}.fai bam_allele/\${subbase}.coverage.q\$minmapq.bw
  rm -f bam_allele/\${subbase}.coverage.q10.bedg
done
"
      jobname="split_allele_$bfn"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 3 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function rnaseq_allelomePro {

  # see wzlib/other/allelomePro.config for example of config file
  # abbreviations:
  # MAT: maternal bias
  # PAT: paternal bias
  # strain1 bias, strain2 bias
  # BAE: Biallelic expression
  # NI: non-informative (low SNP coverage)
  # NS: no SNP
  # I_score: imprinted score, S_score: strain bias score
  
  # RPSM: (10^6*A)/(B*C), A: number of mappable reads at the given single nucleotide position
  # B: number of all mappable reads in the sample, C: read length
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[allelomePro\]/{p=1;next} p&&!/^$/' samples |
    while read config; do
      cmd="
cd $base
mkdir -p allelomePro/$(basename $config .config)
~/tools/allelomepro/Allelome_PRO/allelome_pro.sh -c allelomePro/$config
"
      jobname="allelomePro_${config}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

################################################################################
# other utility
################################################################################

function wzseq_liftbw {
  input=$1
  chain=$2
  output=$3
  chromsize=$4
  bigWigToBedGraph $input $input.bedg.tmp
  liftOver $input.bedg.tmp $chain $output.bedg.tmp
  bedGraphToBigWig $output.bedg.tmp $chromsize $output
}

function wzseq_picard_markdup {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam/before_mdup ]] || mkdir -p bam/before_mdup
  [[ -d tmp ]] || mkdir tmp
  for f in bam/*.bam; do
    bfn=$(basename $f .bam)
    cmd="
cd $base
mv $f bam/before_mdup/$bfn.bam
[[ -e $f.bai ]] && mv $f.bai bam/before_mdup/$bfn.bam.bai
[[ -e $f.flagstat ]] && mv $f.flagstat bam/before_mdup/$bfn.bam.flagstat
java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=bam/before_mdup/$bfn.mdup.stats READ_NAME_REGEX=null INPUT=bam/before_mdup/$bfn.bam OUTPUT=$f TMP_DIR=tmp
(cd bam; ln -s $bfn.bai $bfn.bam.bai;)
samtools flagstat $f >$f.flagstat
" # other options: REMOVE_DUPLICATES=true
    jobname="picard_markdup_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 60 -memG 50 -ppn 5
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function wzseq_clean_intermediate {

  if [[ -d bam/before_mdup ]]; then
    echo "Remove mark duplicate intermediate"
    [[ -n $(compgen -G bam/before_mdup/*.bam) ]] && rm -i bam/before_mdup/*.bam
    [[ -n $(compgen -G bam/before_mdup/*.bai) ]] && rm -i bam/before_mdup/*.bai
  fi
}

function wzseq_check_failure {
  grep 'job killed' pbs/*.stderr
}

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
  awk '/^\[/{p=0}/\[merge_fastq\]/{p=1;next} p&&!/^$/' samples |
    while read target sources; do
      sources=${sources//,/ }
      cmd="
cd $base/fastq
zcat $sources | pigz -p 4 -c > $target
"
      jobname="merge_fastq_${target//\//_}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 4
      [[ ${!#} == "do" ]] && qsub $pbsfn
    done
}

function wzseq_merge_bam {

  # "samtools merge -r" only add read group to each record, not to
  # the header, reheader helped insert read groups to header
  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[merge\]/{p=1;next} p&&!/^$/' samples | 
    while read merged sourcebams; do
      cmd="
cd $base
samtools merge -r bam/$merged.tmp.bam \$(echo $sourcebams | tr ',' ' ')
~/wzlib/pyutils/addRGtoSAMHeader.py -H -i bam/$merged.tmp.bam -o - | samtools reheader - bam/$merged.tmp.bam >bam/$merged.bam
rm -f bam/$merged.tmp.bam
samtools index bam/$merged.bam
samtools flagstat bam/$merged.bam >bam/$merged.bam.flagstat
"
      jobname="bam_merge_$merged"
      pbsfn=pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function wzseq_merge_bam_picard {

  # this properly handles read groups, though not sure how useful read groups are...
  # this add RG to original source bam and output temporary files, merge works on
  # those temporary files
  base=$(pwd)
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[merge\]/{p=1;next} p&&!/^$/' samples |
    while read merged sourcebams; do
      cmd="
cd $base
rm -rf bam/tmp
mkdir -p bam/tmp
sourceid=0
mergeinput=\"\"
for sourcebam in \$(echo $sourcebams | tr ',' ' '); do
  sourceid=\$((\$sourceid + 1))
  rgname=\$(basename \$sourcebam .bam)
  java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar AddOrReplaceReadGroups I=\$sourcebam O=bam/tmp/\${sourceid}.bam RGID=\$rgname RGLB=NA RGPL=illumina RGPU=NA RGSM=\$rgname
  mergeinput=\$mergeinput\" I=bam/tmp/\${sourceid}.bam\"
done

java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar MergeSamFiles \$mergeinput O=bam/$merged.bam ASSUME_SORTED=true CREATE_INDEX=true
samtools flagstat bam/$merged.bam >bam/$merged.bam.flagstat
rm -rf bam/tmp
"
      jobname="picard_bam_merge_$merged"
      pbsfn=pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

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
  for f in fastq/*.+\(fastq|fq|fastq.gz|fq.gz\); do
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

# create coverage track, unique and nonunique mapping
function wzseq_bam_coverage {

  base=$(pwd);
  [[ -d tracks ]] || mkdir tracks
  [[ -d pbs ]] || mkdir pbs
  [[ -d qc ]] || mkdir qc
  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bam)
    cmd="
cd $base

# coverge to bw tracks
bedtools genomecov -ibam $fn -g ${WZSEQ_REFERENCE}.fai -bga -split | LC_ALL=C sort -k1,1 -k2,2n -T tracks/  >tracks/${bfn}.coverage.bedg
bedGraphToBigWig tracks/${bfn}.coverage.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}.coverage.bw
rm -f tracks/${bfn}.coverage.bedg

# unique reads
minmapq=10
samtools view -q \$minmapq -b $fn | bedtools genomecov -ibam stdin -g ${WZSEQ_REFERENCE}.fai -bga -split | LC_ALL=C sort -k1,1 -k2,2n -T tracks/  >tracks/${bfn}.coverage.q10.bedg
bedGraphToBigWig tracks/${bfn}.coverage.q10.bedg ${WZSEQ_REFERENCE}.fai tracks/${bfn}.coverage.q\$minmapq.bw
rm -f tracks/${bfn}.coverage.q10.bedg

# coverage statistics
bedtools genomecov -ibam $fn -g ${WZSEQ_REFERENCE}.fai -max 100 >qc/$bfn.coverage_stats.tsv
grep '^genome' qc/$bfn.coverage_stats.tsv | ~/wzlib/Rutils/bin/basics/wzplot.r barplot -c 5 -n 2 -o qc/$bfn.coverage_stats.barplot.pdf --xlab BaseCoverage --ylab BaseFraction
"
    jobname="bam_coverage_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 2 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wzseq_qualimap() {
  # wzseq_qualimap
  # or wzseq_qualimap rnaseq_se
  # or wzseq_qualimap rnaseq_pe_stranded
  base=$(pwd);
  [[ -d qualimap ]] || mkdir -p qualimap
  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .bam)     # fn might be a symbolic link
    cmd="
cd $base
qualimap --java-mem-size=10G bamqc -nt 10 -bam $fn -outdir qualimap/$bfn -c
"
    if [[ $1 == "rnaseq_se" ]]; then
      cmd="$cmd
qualimap --java-mem-size=10G rnaseq -bam $fn -gtf $WZSEQ_GTF -outdir qualimap/rnaseq_$bfn
"
    fi
    if [[ $1 == "rnaseq_pe_stranded" ]]; then
      cmd="$cmd
qualimap --java-mem-size=10G rnaseq -bam $fn -gtf $WZSEQ_GTF -outdir qualimap/rnaseq_$bfn --paired -p strand-specific-forward
"
    fi
    jobname="qualimap_"${f//\//_}
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 50 -ppn 10
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wzseq_picard_index_fasta {
  dictfn=${$WZSEQ_REFERENCE%.fa}
  cmd="
java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar CreateSequenceDictionary REFERENCE=$WZSEQ_REFERENCE OUTPUT=${WZSEQ_REFERENCE}.dict
"
  jobname="picard_index_reference"
  pbsfn=~/pbs/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
  [[ $1 == "do" ]] && qsub $pbsfn
}

function wzseq_picard_WGSmetrics {
  # right now, I keep getting this error
  # Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: 14945
  base=$(pwd);
  [[ -d qc ]] || mkdir qc
  for f in bam/*.bam; do
    bfn=$(basename $f .bam)
    cmd="
cd $base
java -Xmx5g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar CollectWgsMetrics INPUT=$fn OUTPUT=qc/$bfn.wgsmetrics REFERENCE_SEQUENCE=$WZSEQ_REFERENCE
"
    jobname="picard_WGSmetrics_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

# trimming of adaptor and quality 
function wzseq_trimmomatic {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  if [[ ! -s samples ]]; then
    echo "file: samples missing"
    return 1;
  fi
  
  [[ -d fastq/trimmed ]] || mkdir -p fastq/trimmed
  [[ -d fastq/untrimmed ]] || mkdir -p fastq/untrimmed
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
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
    done
}

# basic summary of the folder
function wzseq_basic_summary() {
  # for RNAseq
  for f in bam/*.bam; do b=${f%.bam}; b=${b#bam/}; sec=$(cat ${f%.bam}/accepted_hits.bam.flagstat | grep 'secondary' | awk '{match($1, /([0-9]*) \+/, a); print a[1];}'); tot=$(cat ${f%.bam}/accepted_hits.bam.flagstat | grep 'total (' | awk '{match($1, /([0-9]*) \+/, a); print a[1];}'); echo ${b%.bam} $(($tot-$sec)); done
  echo -e "\nraw read counts (paired-end, single-end)"
  for f in fastq/*.fastq.gz; do c=$(zcat $f | lc); echo $f $(($c / 2)) $(($c / 4)); done
}

function absolute_path() {
  local in=$1
  out=()
  for x in $@; do
    out+=($(readlink -f $x));
  done
  echo "${out[@]}"
}

#############
# RepEnrich
#
# python RepEnrich_setup.py ~/references/hg19/repeatmasker/hg19_repeatmasker.txt ~/references/hg19/hg19.fa
#   ~/references/hg19/RepEnrich/
# https://github.com/nerettilab/RepEnrich
#############

################################################################################
## auto setup links 
################################################################################

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


# an example of "dependency coding" (not recommended)
# function wgbs_picard_mdup {

#   if [[ $# -eq 0 ]]; then
#     echo "biscuit_mdup -i bam/H75J7BGXX_1_Undetermined.bam"
#     echo "PBS: biscuit_mdup -j jid -d depend -i bam/H75J7BGXX_1_Undetermined.bam"
#     return 1
#   fi

#   local OPTIND opt base jobid _jobid depend bam1 bam2 duplog tmp
#   while getopts "b:i:d:j:" opt; do
#     case $opt in
#       b) base=$(readlink -f $OPTARG) ;;
#       i) bam1=$(readlink -f $OPTARG) ;; # required
#       j) jobid=$OPTARG ;;
#       d) depend=$OPTARG ;;
#       \?) echo "Invalid option: -$OPTARG" >&2; return 1;;
#       :) echo "Option -$OPTARG requires an argument." >&2; return 1;;
#     esac
#   done

#   [[ -z ${base+x} ]] && base=$(pwd);
#   [[ -z ${depend+x} ]] && depend="" || depend="-depend "$depend;

#   bam2=${bam1%.bam}.mdup.bam
#   duplog=${bam1%.bam}.mdup_metrics
#   tmp=$base/tmp
#   [[ -d $tmp ]] || mkdir -p $tmp
#   cmds="java -Xmx7g -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar MarkDuplicates CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=$duplog READ_NAME_REGEX=null INPUT=$bam1 OUTPUT=$bam2 TMP_DIR=$tmp"
#   echo $cmds
#   if [[ -z ${jobid+x} ]]; then
#     $cmds
#   else
#     pbsdir=$base/pbs
#     [[ -d $pbsdir ]] || mkdir -p $pbsdir
#     jobname=$(basename $bam1)"_mdup"
#     pbsgen one "$cmds" -name $jobname -dest $pbsdir/$jobname $depend
#     _jobid=$(qsub $pbsdir/$jobname)
#     eval $jobid=$_jobid;
#     echo "Submitted "$_jobid;
#   fi
# }
