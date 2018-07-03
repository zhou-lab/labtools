#!/bin/bash
shopt -s extglob
shopt -s expand_aliases
# Remarks:
# pbsgen one "$cmds" -name $jobname -dest $pbsdir/$jobname $depend
# note PBS scripts won't expand aliases
# to enforce aliases, use:
# shopt -s expand_aliases

#########################
## Usage
#########################
function pipeline_template {
  cat <<'EOF'
source ~/wzlib/bash/wzseq.sh
wzref_hg19
pipeline_prepare

while read sname; do
  jump_comments

done <<EOM
sample1
sample2
EOM
EOF
}

###########################################
## Helper functions
###########################################

_pipeline_prepare=$(cat <<'EOF'
callargarray=("$@"); base=$(pwd); [[ -d pbs ]] || mkdir pbs; leveljobids=();
# echo ${#callargarray[@]},0
if [[ ${#callargarray[@]} > 0 && ${callargarray[-1]} == "do" ]]; then
  pipeline_submit=true
  unset 'callargarray[${#callargarray[@]}-1]'
else
  pipeline_submit=false
fi
# echo ${#callargarray[@]},1

if [[ ${#callargarray[@]} > 0 && ${callargarray[-1]} =~ ^[0-9]+$ ]]; then
  pipeline_select=${callargarray[-1]}
elif [[ ${#callargarray[@]} > 0 && "${callargarray[-1]}" =~ ([0-9]*)-([0-9]*) ]]; then
  pipeline_select=${callargarray[-1]}
  pipeline_range=true
  pipeline_range_start=${BASH_REMATCH[1]}
  pipeline_range_end=${BASH_REMATCH[2]}
else
  pipeline_select="all"
fi
# echo ${#callargarray[@]},2
submitted_jobids=()
level_jobids=()

# echo ${#callargarray[@]},3
# echo "submit:" $pipeline_submit
# echo "select:" $pipeline_select

# pipeline_submit=false
hour=24; memG=10; ppn=1; queue=default
EOF
)

alias pipeline_prepare='eval "$_pipeline_prepare"'
alias jump_comments='sname_re="^#"; [[ "$sname" =~ $sname_re ]] && continue; depend="";'

## define pipeline_select and pipeline_submit
function pipeline_eval {
  pipeline_component=$1
  if [[ "$pipeline_select" == "all" || "$pipeline_component" == "$pipeline_select" || ( "${pipeline_select: -1}" == "+" && "$pipeline_component" -ge "${pipeline_select::-1}" ) || ($pipeline_range && $pipeline_component -ge $pipeline_range_start && $pipeline_component -le $pipeline_range_end) ]]; then
    echo "submit:" $pipeline_submit
    echo "select:" $pipeline_select
    echo "component:" $pipeline_component

    $2
    pbsfn=$base/pbs/${jobname}.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour $hour -memG $memG -ppn $ppn -queue $queue

    ## whether to submit
    if $pipeline_submit; then
      jobid=$(qsub $depend $pbsfn)
      echo "depend: " $depend
      echo "submitted jobid: " $jobid
      echo

      # by default, let next job depend on this one.
      depend="-W depend=afterok:$jobid"

      # collect job id
      submitted_jobids[$pipeline_component]=$jobid
      
      # collect all jobids on the level
      if [[ -z ${leveljobids[$pipeline_component]} ]]; then
        leveljobids[$pipeline_component]=$jobid
      else
        leveljobids[$pipeline_component]="${leveljobids[$pipeline_component]}:$jobid"
      fi
    fi
  fi
}

function pipeline_depend {
  depended_component=$1

  if [[ $depended_component == 'none' ]]; then
    depend=''
  fi

  # if dependency is not submitted, assume it's done
  if [[ -z ${submitted_jobids[$depended_component]} ]]; then
    depend=""
  else
    depend="-W depend=afterok:"${submitted_jobids[$depended_component]}
  fi
}

function pipeline_dependlevel {
  if [[ $# -lt 1 ]]; then
    depend=""
  elif [[ -z "${leveljobids[$1]}" ]]; then
    depend=""
  else
    depend="-W depend=afterok:${leveljobids[$1]}"
  fi
}

################################################################################
# reference configuration
# on a new machine, you just need to update this section
################################################################################

##########
## mouse
##########
###### mouse mm10 #####
function wzref_mm10 {
  export WZSEQ_REFVERSION=mm10
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/mm10/mm10.fa
  export WZSEQ_REFERENCE_SIZE=2785373478

  export WZSEQ_REFERENCE_LAMBDAPHAGE=/home/wanding.zhou/references/lambdaphage/biscuit/NC_001416.fa
  export WZSEQ_GTF=/primary/vari/genomicdata/genomes/mm10/tophat/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf
  export WZSEQ_GTF_ENSEMBL=/primary/vari/genomicdata/genomes/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gz

  # $ zcat ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.gz | awk '$1~/^#/{print}$1~/^[0-9XYM]/{if ($1=="MT"){$1="chrM"}else{$1="chr"$1};print $0}' | gzip -c >~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming.gz
  export WZSEQ_GTF_ENSEMBL_UCSCNAMING=/primary/vari/genomicdata/genomes/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming

  # [~/references/mm10/gtf]$ python ~/.Renv/versions/3.2.3/lib64/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py Mus_musculus.GRCm38.82.gtf.UCSCnaming Mus_musculus.GRCm38.82.gtf.UCSCnaming.DEXSeq.gff
  export WZSEQ_GTF_DEXSEQ=/home/wanding.zhou/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming.DEXSeq.gff

  export WZSEQ_BOWTIE2_INDEX=/home/wanding.zhou/references/mm10/bowtie2/mm10
  export WZSEQ_BOWTIE1_INDEX=/primary/vari/genomicdata/genomes/mm10/bowtie1/mm10
  export WZSEQ_BWA_INDEX=/home/wanding.zhou/references/mm10/bwa/mm10.fa
  export WZSEQ_STAR_INDEX=/primary/vari/genomicdata/genomes/mm10/STAR
  export WZSEQ_GSNAP_INDEX=/primary/vari/genomicdata/genomes/mm10/gsnap
  # needs $WZSEQ_GSNAP_INDEX/$WZSEQ_GSNAP_SPLICE.iit
  export WZSEQ_GSNAP_SPLICE=Mus_musculus.GRCm38.82.gtf.gsnap.splicesites
  export WZSEQ_SUBREAD_INDEX=/primary/vari/genomicdata/genomes/mm10/subread/mm10

  export WZSEQ_KALLISTO_INDEX=/primary/vari/genomicdata/genomes/mm10/kallisto/mm10.kallisto

  # WGBS indices
  export WZSEQ_BISCUIT_INDEX_LAMBDAPHAGE=/home/wanding.zhou/references/lambdaphage/biscuit/NC_001416.fa
  export WZSEQ_BISCUIT_INDEX=/home/wanding.zhou/references/mm10/biscuit/mm10.fa
  export WZSEQ_BISCUIT_QC_SETUP=/home/wanding.zhou/tools/biscuit/development/biscuit/test/QC_assets/mm10_QC_assets/setup.sh
  export WZSEQ_BWAMETH_INDEX=/home/wanding.zhou/references/mm10/bwameth/mm10.fa

  export WZSEQ_BISMARK_BT1_INDEX=/home/wanding.zhou/references/mm10/bismark_bt1
  export WZSEQ_BISMARK_BT2_INDEX=/home/wanding.zhou/references/mm10/bismark_bt2
  export WZSEQ_BSMAP_INDEX=/home/wanding.zhou/references/mm10/bsmap/mm10.fa

  export WZSEQ_REFERENCE_SPLIT=/primary/vari/genomicdata/genomes/mm10/bowtie1/

  export WZSEQ_MACS_SHORT=mm
  export WZSEQ_CPGBED=/primary/vari/genomicdata/genomes/mm10/annotation/cpg.bed
  export WZSEQ_CGIBED=/home/wanding.zhou/references/mm10/annotation/cpgisland/cpgIslandExt.bed
  export WZSEQ_CGIBED_METHYLKIT=/home/wanding.zhou/references/mm10/annotation/cpgisland/cpgIslandExt.methylKit.bed
  export WZSEQ_TSSBED=/home/wanding.zhou/references/mm10/annotation/TSS/mm10.refseq.tss.bed

  # UCSC table
  export WZSEQ_UCSC_REFSEQ=/home/wanding.zhou/references/mm10/annotation/mm10_RefSeq_FromUCSC.bed
  export WZSEQ_UCSC_CGIBED=/home/wanding.zhou/references/mm10/annotation/mm10_cpgIsland_FromUCSC.bed

  # rmsk
  export WZSEQ_RMSK=/primary/vari/genomicdata/genomes/mm10/annotation/rmsk/rmsk.bed
  # build the following using UCSC table builder
  export WZSEQ_RMSK_GTF=/home/wanding.zhou/references/mm10/annotation/rmsk/rmsk.mm10.gtf
}

#############
## Human
#############

###### human hg19 ####
function wzref_hg19 {
  export WZSEQ_REFVERSION=hg19
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/hg19/hg19.fa
  export WZSEQ_REFERENCE_SIZE=3199901561

  export WZSEQ_GTF=/primary/vari/genomicdata/genomes/hg19/tophat/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
  export WZSEQ_EXON=/primary/vari/genomicdata/genomes/hg19/annotation/gtf_exon_merged.bed
  export WZSEQ_BWA_INDEX=/home/wanding.zhou/references/hg19/bwa/hg19.fa
  export WZSEQ_BOWTIE2_INDEX=/home/wanding.zhou/references/hg19/bowtie2/hg19
  export WZSEQ_BOWTIE1_INDEX=/home/wanding.zhou/references/hg19/bowtie1/hg19
  export WZSEQ_STAR_INDEX=/primary/vari/genomicdata/genomes/hg19/STAR
  export WZSEQ_RSEQC_GENE_BED=/primary/vari/genomicdata/genomes/hg19/rseqc/hg19_GENCODE_GENE_V19_comprehensive.bed

  # awk '!/^#/{if($1~/^[0-9XY]*$/) $1="chr"$1; if($1=="MT") $1="chrM"; print $0}/^#/' gtf/Homo_sapiens.GRCh37.75.gtf >gtf/Homo_sapiens.GRCh37.75.gtf.UCSCnaming
  export WZSEQ_GTF_ENSEMBL_UCSCNAMING=/primary/vari/genomicdata/genomes/hg19/gtf/Homo_sapiens.GRCh37.75.gtf.UCSCnaming
  export WZSEQ_BISMARK_BT2_INDEX=/home/wanding.zhou/references/hg19/bismark_bt2
  export WZSEQ_BWAMETH_INDEX=/home/wanding.zhou/references/hg19/bwameth/hg19.fa
  export WZSEQ_SUBREAD_INDEX=/primary/vari/genomicdata/genomes/hg19/subread/hg19
  export WZSEQ_REFERENCE_SPLIT=/primary/vari/genomicdata/genomes/hg19/tophat/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes
  # export WZSEQ_HISAT2_INDEX=/primary/vari/genomicdata/genomes/hg19/hisat/genome
  export WZSEQ_HISAT2_INDEX=/primary/vari/genomicdata/genomes/hg19_noContig/hisat2/genome
  export WZSEQ_EXOME_CAPTURE=/primary/vari/genomicdata/genomes/hg19/annotation/hg19.exomes.bed
  export WZSEQ_KALLISTO_INDEX=/primary/vari/genomicdata/genomes/hg19/kallisto/hg19.kallisto

  # bedtools makewindows -w 100 -g hg19.fa.fai >windows/hg19.windows100bp.bed
  # seqtk comp -r windows/hg19.windows100bp.bed ~/references/hg19/hg19.fa | awk '$4+$5+$6+$7>0{print $1,$2,$3,($5+$6)/($4+$5+$6+$7)}' >windows/hg19.windows100bp.gc_content.bed
  # sort -k4,4n windows/hg19.windows100bp.gc_content.bed >windows/hg19.windows100bp.gc_content.srted.bed
  # head -3000000 windows/hg19.windows100bp.gc_content.srted.bed | sortbed >windows/hg19.windows100bp.gc_content.top10p.bed
  # tail -3000000 windows/hg19.windows100bp.gc_content.srted.bed | sortbed >windows/hg19.windows100bp.gc_content.bot10p.bed
  export WZSEQ_TOPGC_BED=/primary/vari/genomicdata/genomes/hg19/windows/hg19.windows100bp.gc_content.top10p.bed
  export WZSEQ_BOTGC_BED=/primary/vari/genomicdata/genomes/hg19/windows/hg19.windows100bp.gc_content.bot10p.bed

  # WGBS
  export WZSEQ_BISCUIT_INDEX=/home/wanding.zhou/references/hg19/biscuit/hg19.fa
  export WZSEQ_BSMAP_INDEX=/home/wanding.zhou/references/hg19/bsmap/hg19.fa
  export WZSEQ_CPGBED=/primary/vari/genomicdata/genomes/hg19/annotation/cpg.bed
  export WZSEQ_CGIBED=/primary/vari/genomicdata/genomes/hg19/annotation/cpgisland/cpgIslandExt.bed
  export WZSEQ_TSSBED=/primary/vari/genomicdata/genomes/hg19/annotation/hg19.refseq.tss.bed
  export WZSEQ_MACS_SHORT=hs
  export WZSEQ_BISCUIT_QC_SETUP=/home/wanding.zhou/tools/biscuit/development/biscuit/test/QC_assets/hg19_QC_assets/setup.sh

  # rmsk
  export WZSEQ_RMSK=/primary/vari/genomicdata/genomes/hg19/annotation/rmsk/rmsk.txt.bed
  # build the following using the UCSC table builder
  export WZSEQ_RMSK_GTF=/home/wanding.zhou/references/hg19/annotation/rmsk/rmsk.hg19.gtf

  export WZSEQ_HICPRO_CONFIG=/home/wanding.zhou/wzprojects/2017/2017_05_13_CTCF_EZH2_hicpro_config_template
}

###### human hg19_noContig ####
function wzref_hg19_noContig {
  export WZSEQ_REFVERSION=hg19
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/hg19_noContig/hg19_noContig.fa

  export WZSEQ_BISMARK_BT2_INDEX=/home/wanding.zhou/references/hg19_noContig/bismark_bt2
  export WZSEQ_BWAMETH_INDEX=/home/wanding.zhou/references/hg19_noContig/bwameth/hg19_noContig.fa

  # WGBS
  export WZSEQ_BISCUIT_INDEX=/home/wanding.zhou/references/hg19_noContig/biscuit/hg19_noContig.fa
  export WZSEQ_BSMAP_INDEX=/home/wanding.zhou/references/hg19_noContig/bsmap/hg19_noContig.fa
}


function wzref_hg38 {
  export WZSEQ_REFVERSION=hg38
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/hg38/hg38.fa
  # export WZSEQ_REFERENCE_SIZE=3199901561

  # export WZSEQ_GTF=/primary/vari/genomicdata/genomes/hg19/tophat/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
  # export WZSEQ_BWA_INDEX=/home/wanding.zhou/references/hg19/bwa/hg19.fa
  # export WZSEQ_BOWTIE2_INDEX=/home/wanding.zhou/references/hg19/bowtie2/hg19
  # export WZSEQ_BOWTIE1_INDEX=/home/wanding.zhou/references/hg19/bowtie1/hg19
  # export WZSEQ_STAR_INDEX=/primary/vari/genomicdata/genomes/hg19/STAR
  # export WZSEQ_RSEQC_GENE_BED=/primary/vari/genomicdata/genomes/hg19/rseqc/hg19_GENCODE_GENE_V19_comprehensive.bed
  export WZSEQ_BISMARK_BT2_INDEX=/home/wanding.zhou/references/hg38/bismark_bt2
  export WZSEQ_BISCUIT_INDEX=/primary/vari/genomicdata/genomes/hg38/biscuit/hg38.fa
  export WZSEQ_BISCUIT_QC_SETUP=/home/wanding.zhou/tools/biscuit/development/biscuit/test/QC_assets/hg38_QC_assets/setup.sh
  # export WZSEQ_SUBREAD_INDEX=/primary/vari/genomicdata/genomes/hg19/subread/hg19
  # export WZSEQ_REFERENCE_SPLIT=/primary/vari/genomicdata/genomes/hg19/tophat/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes

  # WGBS
  # export WZSEQ_BISCUIT_INDEX=/home/wanding.zhou/references/hg19/biscuit/hg19.fa
  # export WZSEQ_CGIBED=/primary/vari/genomicdata/genomes/hg19/annotation/cpgisland/cpgIslandExt.bed
  # export WZSEQ_MACS_SHORT=hs

  # rmsk
  # export WZSEQ_RMSK=/primary/vari/genomicdata/genomes/hg19/annotation/rmsk/rmsk.txt.bed
}

###### human hg19 rCRS ####
function wzref_hg19rCRS {
  export WZSEQ_REFVERSION=hg19
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/hg19-rCRS/hg19_rCRS.fa
  export WZSEQ_DBSNP=/primary/vari/genomicdata/genomes/hg19-rCRS/dbsnp_137.hg19.vcf
  export WZSEQ_EXOME_CAPTURE=/primary/vari/genomicdata/genomes/hg19/annotation/hg19.exomes.bed
}


##################
## lancelet
##################

function wzref_braFlo1 {
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/braFlo1/braFlo1.fa
  export WZSEQ_BISCUIT_INDEX=/primary/vari/genomicdata/genomes/braFlo1/biscuit/braFlo1.fa
}

function wzref_braFlo2 {
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/braFlo2/braFlo2.fa
  export WZSEQ_BISCUIT_INDEX=/primary/vari/genomicdata/genomes/braFlo2/biscuit/braFlo2.fa
}

##################
## Elephant Shark
##################

function wzref_calMil1 {
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/calMil1/calMil1.fa
  export WZSEQ_BISCUIT_INDEX=/primary/vari/genomicdata/genomes/calMil1/biscuit/calMil1.fa
}

############
## Lamprey
############

function wzref_petMar2 {
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/petMar2/petMar2.fa
  export WZSEQ_BISCUIT_INDEX=/primary/vari/genomicdata/genomes/petMar2/biscuit/petMar2.fa
}

##############
## Sea Squirt
##############

function wzref_ci2 {
  export WZSEQ_REFERENCE=/primary/vari/genomicdata/genomes/ci2/ci2.fa
}

################################################################################
# DNA-seq, mutation calling etc.
################################################################################

function examplepipeline_dna {
  cat <<- EOF
=== pipeline 2016-02-05 ===

=> (+) wzseq_pindel => (+) wzseq_pindel_somatic

EOF
}

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

function __wzseq_bwa_mem_PE {

  cmd='
cd '$base'
mkdir -p bam;
bwa mem -M -R \"@RG\tLB:$WZSEQ_REFVERSION\tID:${sname}\tPL:Illumina\tPU:nextseq500\tSM:${sname}\" -t 28 $WZSEQ_BWA_INDEX fastq/$sread1 fastq/$sread2 | samtools sort -O bam -T bam/$sname.tmp -o bam/$sname.bam -
samtools index bam/$sname.bam
samtools flagstat bam/${sname}.bam > bam/$sname.bam.flagstat
'
  jobname="bwamem_$sname"
}

function __wzseq_bwa_mem_SE {
  cmd='
cd '$base'
mkdir -p bam;
bwa mem -M -R "@RG\tLB:'$WZSEQ_REFVERSION'\tID:'$sname'\tPL:Illumina\tPU:nextseq500\tSM:'$sname'" -t '$ppn' '$WZSEQ_BWA_INDEX' '$fastq' | samtools sort -O bam -T bam/'$sname'.tmp -o bam/'$sname'.bam -
samtools index bam/'$sname'.bam
samtools flagstat bam/'$sname'.bam > bam/'$sname'.bam.flagstat
'
  jobname="bwamem_"$sname
}

# wzseq_bwa_mem
function wzseq_bwa_mem {

  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      __wzseq_bwa_mem
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

function wzseq_pindel {

  # call pindel_somatic directly if you want somatic indels
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d pindel ]] || mkdir pindel

  awk '/^\[/{p=0}/\[pindel\]/{p=1;next} p&&!/^$/' samples |
    while read sname bam; do
      cmd="
cd $base
[[ -d pindel/$sname ]] || mkdir pindel/$sname
echo -e $(readlink -f $bam)\"\\t250\\t\"$sname >pindel/$sname/pindel.config
~/software/pindel/default/pindel -i pindel/$sname/pindel.config -f $WZSEQ_REFERENCE -o pindel/$sname/$sname -c ALL -T 28 2>&1 > pindel/$sname/$sname.pindel_log
"
      jobname="pindel_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done

}

function wzseq_pindel_somatic {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d pindel ]] || mkdir pindel
  # sname normal tumor
  awk '/^\[/{p=0}/\[pindel_somatic\]/{p=1;next} p&&!/^$/' samples |
    while read sname sname1 sname2; do
      cmd="
cd $base
[[ -d pindel/$sname ]] || mkdir pindel/$sname
echo -e bam/$sname1.bam\"\\t250\\t\"$sname1 >pindel/$sname/pindel.config
echo -e bam/$sname2.bam\"\\t250\\t\"$sname2 >>pindel/$sname/pindel.config
~/software/pindel/default/pindel -i pindel/$sname/pindel.config -f $WZSEQ_REFERENCE -o pindel/$sname/$sname -c ALL -T 28 2>&1 > pindel/$sname/$sname.pindel_log

# filter somatic events
grep \"ChrID\" pindel/$sname/${sname}_D >pindel/$sname/head_DSI
grep \"ChrID\" pindel/$sname/${sname}_SI >>pindel/$sname/head_DSI
cat <<EOT >pindel/$sname/somatic.config
indel.filter.input = pindel/$sname/head_DSI
indel.filter.vaf = 0.1
indel.filter.cov = 5
indel.filter.hom = 2
indel.filter.pindel2vcf = /primary/vari/software/pindel/default/pindel2vcf
indel.filter.reference = $WZSEQ_REFERENCE
indel.filter.referencename = NA
indel.filter.referencedate = NA
indel.filter.output = pindel/$sname/$sname.indel.filter
EOT
perl ~/software/pindel/default/somatic_filter/somatic_indelfilter.pl pindel/$sname/somatic.config
"
      jobname="pindel_somatic_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# compute base coverage for both tumor and normal
function wzseq_basecov_somatic {
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d basecov ]] || mkdir basecov
  awk '/^\[/{p=0}/\[pindel_somatic\]/{p=1;next} p&&!/^$/' samples |
    while read sname sname1 sname2; do
      cmd="
cd $base
samtools mpileup -s bam/$sname1.bam bam/$sname2.bam | perl -alne 'BEGIN{\$b=0}{\$c=0; foreach(split //, \$F[6]){if (ord(\$_)>53) {\$c++;}} \$d=0; foreach(split //, \$F[10]){if(ord(\$_)>53){\$d++;}} if(\$c>4 && \$d>4){\$b++;}}END{print \$b}' > basecov/${sname1}_vs_${sname2}.cov;
"
      jobname="basecov_somatic_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function wzseq_lofreq {

  base=$(pwd)
  [[ -d lofreq ]] || mkdir lofreq
  [[ -d pbs ]] || mkdir pbs

  for bam in bam/*.bam; do
    sname=$(basename $bam .bam)
    cmd="
cd $base
~/software/lofreq/default/bin/lofreq call -f $WZSEQ_REFERENCE -s -S $WZSEQ_DBSNP -o lofreq/$sname.vcf $bam -b 1
transvar ganno --vcf lofreq/$sname.vcf --ccds >lofreq/$sname.vcf.transvar
"
    jobname="lofreq_$sname"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 60 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wzseq_vardict {
  base=$(pwd)
  [[ -d vardict ]] || mkdir vardict
  [[ -d pbs ]] || mkdir pbs
  bedregion=""
  for bam in bam/*.bam; do
    sname=$(basename $bam .bam)
    cmd="
cd $base
export PATH=/primary/vari/software/vardict/VarDict:$PATH
## AF_THR=0.001 # minimum allele frequency
# vardict -G $WZSEQ_REFERENCE -f 0.001 -N ${sname} -b $bam -c 1 -S 2 -E 3 -g 4 $WZSEQ_EXOME_CAPTURE > vardict/${sname}_out.raw
# cat vardict/${sname}_out.raw | teststrandbias.R >vardict/${sname}_out.strandbias
# cat vardict/${sname}_out.strandbias | var2vcf_valid.pl -N sample_name -E -f 0.001 >vardict/${sname}_out.vcf
transvar ganno --vcf vardict/${sname}_out.vcf --ccds >vardict/${sname}_out.vcf.transvar
"
    jobname="vardict_$sname"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
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

 => (+) wgbs_biscuit_pileup [-nome] => (+) wgbs_vcf2tracks [-nome] => (+) wgbs_cpgcoverage => (+) wgbs_repeat => (+) wgbs_repeat_diff

 => (+) wgbs_methylKit_summary

 => (+) wgbs_diffmeth_simple => (+) wgbs_methylKit_diffmeth => (+) wgbs_methpipe_diff => (+) wgbs_metilene
EOF
}

# check adaptor conversion rate
function wgbs_adaptor() {
  local base=$(pwd);
  [[ -d adaptor ]] || mkdir adaptor;
  [[ -d pbs ]] || mkdir pbs
  cmd="
cd $base
parallel -j 4 ~/wzlib/pyutils/wzadapter.py {} '>' adaptor/{/.}.txt ::: $base/fastq/*R1*.fastq.gz
for f in adaptor/*.txt; do
  tail -2 \$f | cut -d\":\" -f2 | cut -d\" \" -f2 | awk -v fn=\$f '{a[NR]=\$1}END{print fn,a[1],a[2]}'
done | awk 'BEGIN{print \"filename\tadaptorC\tadaptorT\tadaptorC2T\"}{print \$1,\$2,\$3,\$3/\$2*100\"%\"}' > adaptor/adaptor.stats
"
  jobname="adaptor_analysis_"$(basename $base)
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 8 -ppn 4
  [[ ${!#} == "do" ]] && qsub $pbsfn
}

function __wzseq_clean_bam {

  # input: bamfn
  cmd="
set -xe
cd $base
rm -f bam/$bamfn
"
  jobname="clean_bam_$bamfn"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 1 -memG 1 -ppn 1
}

function __wzseq_clean_fastq {

  # input: sname, sread1, sread2
  cmd="
set -xe
cd $base
[[ -f fastq/$sread1 ]] && rm -f fastq/$sread1
[[ -f fastq/$sread2 ]] && rm -f fastq/$sread2
"
  jobname="clean_fastq_$sname"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 1 -memG 1 -ppn 1
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

function __wgbs_biscuit_align_SE {
  cmd='
cd '$base'
mkdir -p bam
~/tools/biscuit/master/biscuit/biscuit align '$WZSEQ_BISCUIT_INDEX' -t '$ppn' '$fastq' | samtools sort -T '$output_bam'_tmp -O bam -o '$output_bam'
'
  jobname="biscuit_align_"$sname"_SE"
}

function __wgbs_biscuit_align_PE {
  cmd='
cd '$base'
mkdir -p bam
~/bin/biscuit align '$WZSEQ_BISCUIT_INDEX' -t '$ppn' '$fastq1' '$fastq2' | samtools sort -T '$output_bam'_tmp -O bam -o '$output_bam'
'
  jobname="biscuit_align_"$sname"_PE"
}

function __wgbs_biscuit_align_PE_both {
  cmd='
cd '$base'
mkdir -p bam
# ~/tools/biscuit/master/biscuit/biscuit
~/bin/biscuit align '$WZSEQ_BISCUIT_INDEX' -b 1 -t '$ppn' '$fastq1' '$fastq2' >'$output_bam'.sam
samtools sort -T '$output_bam'_tmp -O bam -o '$output_bam' '$output_bam'.sam
#samtools index '$output_bam';
#samtools flagstat '$output_bam' > '$output_bam'.flagstat
'
  jobname="biscuit_align_"${sname}"_PE_both"
}

function __wgbs_biscuit_align_PE_Walid_lib {
  cmd='
cd '$base'
mkdir -p bam
# ~/tools/biscuit/master/biscuit/biscuit
~/bin/biscuit align '$WZSEQ_BISCUIT_INDEX' -b 1 -t '$ppn' -J AGATCGGAAGAGCGGTTCAGCA -K TCTAGCCTTCTCGCAGCACA '$fastq1' '$fastq2' >'$output_bam'.sam
samtools sort -T '$output_bam'_tmp -O bam -o '$output_bam' '$output_bam'.sam
#samtools index '$output_bam';
#samtools flagstat '$output_bam' > '$output_bam'.flagstat
'
  jobname="biscuit_align_"${sname}"_PE_both"
}

function __wgbs_biscuit_align_PE_POETIC {
  cmd='
cd '$base'
mkdir -p bam
~/bin/biscuit align '$WZSEQ_BISCUIT_INDEX' -t '$ppn' -J AGATCGGAAGAGC -K AAATCAAAAAAAC '$fastq1' '$fastq2' >'$output_bam'.sam
samtools sort -T '$output_bam'_tmp -O bam -o '$output_bam' '$output_bam'.sam
# samtools index '$output_bam';
# samtools flagstat '$output_bam' > '$output_bam'.flagstat
'
  jobname="biscuit_align_"${sname}"_PE_POETIC"
}

function __wgbs_biscuit_markdup {
  cmd='
cd '$base'
~/tools/biscuit/development/biscuit/biscuit markdup '$input_bam' '$output_bam'_unsrt.bam 2>'$output_bam'_markdup_report.txt
samtools sort -T '$output_bam'_unsrt.tmp -o '$output_bam' '$output_bam'_unsrt.bam
rm -f '$output_bam'_unsrt.bam
mkdir -p multiqc/raw/biscuit/
ln -sf `readlink -f '$output_bam'_markdup_report.txt` multiqc/raw/biscuit/
'
  jobname='biscuit_markdup_'$sname
}

function __wgbs_biscuit_QC {
  cmd='
cd '$base'
~/tools/biscuit/development/biscuit/scripts/QC.sh -v '$input_vcf' '$WZSEQ_BISCUIT_QC_SETUP' '$sname' '$input_bam'
mkdir -p multiqc/raw/BISCUITqc/'$sname'
ln -sf `readlink -f BISCUITqc` multiqc/raw/BISCUITqc
'
  jobname='biscuit_QC_'$sname
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
if [[ \"$sread2\" == \".\" ]]; then
  bwameth --reference $WZSEQ_BWAMETH_INDEX fastq/$sread1 -t 28 --prefix bam/${sname}_bwameth
else
  bwameth --reference $WZSEQ_BWAMETH_INDEX fastq/$sread1 fastq/$sread2 -t 28 --prefix bam/${sname}_bwameth
fi
samtools flagstat bam/${sname}_bwameth.bam > bam/${sname}_bwameth.bam.flagstat
"
      jobname="bwameth_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmds" -name $jobname -dest $pbsfn -hour 8 -ppn 28 -memG 250
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
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 100 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

## TRIM GALORE!
##################
function __wzseq_trim_galore_SE {
  cmd='
set -xe
cd '$base'
[[ -d '$trim_galore_dir' ]] && rm -rf '$trim_galore_dir'
mkdir -p '$trim_galore_dir'
~/software/trim_galore/default/trim_galore --fastqc '$fastq' --gzip -o '$trim_galore_dir'
mkdir -p multiqc/raw/trim_galore/
ln -fs `readlink -f '$trim_galore_dir'` multiqc/raw/trim_galore/
'
  jobname="trim_galore_SE_"$sname
}

function __wzseq_trim_galore_PE {
  cmd='
set -xe
cd '$base'
[[ -d '$trim_galore_dir' ]] && rm -rf '$trim_galore_dir'
mkdir -p '$trim_galore_dir'
~/software/trim_galore/default/trim_galore --fastqc --paired '$fastq1' '$fastq2' -o '$trim_galore_dir'
mkdir -p multiqc/raw/trim_galore/
ln -s `readlink -f '$trim_galore_dir'` multiqc/raw/trim_galore/
'
  jobname="trim_galore_PE_"$sname
}

function __wzseq_customize {
  cmd="
set -xe
cd $base
$customized_command
"
  jobname="wzseq_customize_$sname"
}

## bismark with bowtie2
###########################
function wgbs_bismark_bowtie2_prepare_reference {
  base=$(pwd)
  cmd="
export PATH=~/tools/bismark/default:~/tools/bowtie2/default:$PATH
bismark_genome_preparation --bowtie2 --verbose $WZSEQ_BISMARK_BT2_INDEX
"
  jobname="bismark_bt2_prepare_reference"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 10 -ppn 1
  [[ $1 == "do" ]] && qsub $pbsfn
}

function __wgbs_bismark_bowtie2_SE {
  : '
fastq=fastq/${sname}_trim_galore/${sname}_merged.fq.gz
direction="--non_directional"
bismark_bt2_dir=bam/${sname}_bismark_bt2
bismark_bt2_bam_unsorted=bam/${sname}_bismark_bt2/${sname}_merged.fq.gz_bismark_bt2.bam
bismark_bt2_bam_final=bam/${sname}_bismark_bt2.bam
hour=200; memG=180; ppn=28
'
  cmd='
export PATH=~/tools/bismark/default:~/tools/bowtie2/default:$PATH
set -xe
cd '$base'
mkdir -p bam
rm -rf '$bismark_bt2_dir'
bismark '$direction' '$WZSEQ_BISMARK_BT2_INDEX' --bowtie2 --chunkmbs 2000 -p '$ppn' -o '$bismark_bt2_dir' '$fastq' --temp_dir '$bismark_bt2_dir'/tmp
samtools sort -O bam -o '$bismark_bt2_bam_final' -T '${bismark_bt2_bam_unsorted}'_tmp '${bismark_bt2_bam_unsorted}'
samtools index '${bismark_bt2_bam_final}'
samtools flagstat '${bismark_bt2_bam_final}' >'${bismark_bt2_bam_final}'.flagstat
mkdir -p multiqc/raw/bismark/
ln -sf `readlink -f '$bismark_bt2_dir'/[PS]E_report.txt` multiqc/raw/bismark/; done
'
  jobname="bismark_bt2_SE_"$sname
}

function __wgbs_bismark_bowtie2_PE {
  # fastq1=fastq/${sname}_trim_galore/${sname}_R1_val_1.fq.gz
  # fastq2=fastq/${sname}_trim_galore/${sname}_R2_val_2.fq.gz
  # bismark_bt2_dir=bam/${sname}_bismark_bt2
  # bismark_bt2_bam_unsorted=bam/${sname}_bismark_bt2/${sname}_pe.bam
  # bismark_bt2_bam_final=bam/${sname}_bismark_bt2.bam
  # hour=200; memG=180; ppn=28
  cmd='
export PATH=~/tools/bismark/default:~/tools/bowtie2/default:$PATH
set -xe
cd '$base'
mkdir -p bam
rm -rf '$bismark_bt2_dir'
bismark '$WZSEQ_BISMARK_BT2_INDEX' --bowtie2 --chunkmbs 2000 -p '$ppn' -o '$bismark_bt2_dir' -B '$sname' -1 '$fastq1' -2 '$fastq2' -temp_dir '$bismark_bt2_dir'/tmp
samtools sort -O bam -o '$bismark_bt2_bam_final' -T '$bismark_bt2_bam_unsorted'_tmp '$bismark_bt2_bam_unsorted'
samtools index '$bismark_bt2_bam_final'
samtools flagstat '$bismark_bt2_bam_final' >'$bismark_bt2_bam_final'.flagstat
mkdir -p multiqc/raw/bismark
ln -s `readlink -f '$bismark_bt2_dir'` multiqc/raw/bismark/
'
  jobname="bismark_bt2_PE_$sname"
}

function __wgbs_bismark_deduplicate {
  # Please note that this function is very memory-consuming
  # input_bam=bam/${sname}_bismark_bt2.bam
  # hour=48; memG=50; ppn=2
  # library="-p"
  ## paired-end
  # library="-p"
  ## single-end
  # library="-s"
  cmd='
export PATH=~/tools/bismark/default:~/tools/bowtie2/default:$PATH
cd '$base'
mkdir -p multiqc/raw/bismark/
deduplicate_bismark '$library' --bam '$input_bam'
ln -sf `readlink -f '${input_bam%.bam}'.deduplication_report.txt` multiqc/raw/bismark/
'
  jobname="bismark_deduplicate_"$sname
}

function __wgbs_bismark_methylextraction {
  : '
# for paired-end
library="--no_overlap"
# for single-end
library=""

input_bam=bam/${sname}_bismark_bt2.deduplicated.bam
hour=48; memG=10; ppn=1
'
  cmd='
export PATH=~/tools/bismark/default:~/tools/bowtie2/default:$PATH
cd '$base'
mkdir -p bismark_methylextraction
bismark_methylation_extractor '$library' --gzip --bedGraph '$input_bam' -o bismark_methylextraction
zcat bismark_methylextraction/'$(basename $input_bam .bam)'.bismark.cov.gz | awk '\''{print $1,$2-1,$3,$4/100,$5+$6}'\'' | sortbed | biscuit mergecg '$WZSEQ_REFERENCE' - >bismark_methylextraction/'$sname'.cpg_methylation.bed
ln -sf `readlink -f bismark_methylextraction/'$(basename $input_bam)'_splitting_report.txt` multiqc/raw/bismark/
'
  jobname="bismark_methylation_extraction_"$sname
}

## bsmap
###########
function wgbs_bsmap {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam

  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      cmd="
cd $base
if [[ \"$sread2\" == \".\" ]]; then
  ~/tools/bsmap/default/bsmap -a fastq/$sread1 -d $WZSEQ_BSMAP_INDEX -o bam/${sname}_bsmap_unsorted.bam -p 28 -s 16 -v 10 -q 2 -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
else
  ~/tools/bsmap/default/bsmap -a fastq/$sread1 -b fastq/$sread2 -d $WZSEQ_BSMAP_INDEX -o bam/${sname}_bsmap_unsorted.bam -p 28 -s 16 -v 10 -q 2 -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
fi

samtools sort -O bam -o bam/${sname}_bsmap.bam -T bam/${sname}_bsmap.tmp bam/${sname}_bsmap_unsorted.bam
rm -f bam/${sname}_bsmap_unsorted.bam
samtools index bam/${sname}_bsmap.bam
samtools flagstat bam/${sname}_bsmap.bam > bam/${sname}_bsmap.bam.flagstat
"
      jobname="bsmap_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 8 -memG 250 -ppn 28
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

#### rmapbs (doesn't allow gaps) ###
# rmapbs -c hg18 -o Human_NHFF.mr Human_NHFF.fastq
# rmapbs-pe -c hg18 -o Human_ESC.mr Human_ESC_1.fastq Human_ESC_2.fastq

# methylKit summary
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

#############
## MethPIPE
#############

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
# not sure if these 2 should be removed as intermediate: $obfn.mr.dremove $obfn.mr.sorted_end_first
echo \"[\$(date)] Done.\"
"
    jobname="methpipe_basic_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

## methpipe
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
  for f in methpipe/*.mr.dremove; do
    obfn=${f%.mr.dremove}
    cmd="
cd $base

echo \"[\$(date)] Make read distribution, (.epiread)\"
~/tools/methpipe/default/bin/methstates -c $WZSEQ_REFERENCE_SPLIT -o ${obfn}.epiread ${obfn}.mr.dremove

echo \"[\$(date)] Single site ASM scoring (.allelic)\"
~/tools/methpipe/default/bin/allelicmeth -c $WZSEQ_REFERENCE_SPLIT -o ${obfn}.allelic ${obfn}.epiread

echo \"[\$(date)] Allelically methylated region (AMRs)\"
~/tools/methpipe/default/bin/amrfinder -o $obfn.amr -c $WZSEQ_REFERENCE_SPLIT $obfn.epiread.tmp
sort -k1,1 -k2,2n -T methpipe/ $obfn.epiread.tmp > $obfn.epiread

# echo \"[\$(date)] Test allele-specific methylation at promoter of imprinted region\"
# amrtester -o ${obfn}.amr -c hg19 target_interval.bed ${obfn}.epiread

echo \"[\$(date)] Calculate meth-entropy\"
# -v verbose
~/tools/methpipe/default/bin/methentropy -w 5 -v -o $obfn.entropy $WZSEQ_REFERENCE_SPLIT $obfn.epiread

echo \"[\$(date)] Done.\"
"
    jobname="methpipe_allele_"$(basename $obfn)
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

####################
## BISCUIT
####################
function __wgbs_biscuit_pileup {
  : '
input_bam=bam/${sname}_markdup.bam
output_vcf=pileup/${sname}.vcf.gz
hour=24; memG=100; ppn=10
pipeline_depend bamjob
'
  cmd='
cd '$base'
mkdir -p pileup
~/tools/biscuit/development/biscuit/biscuit pileup -q '$ppn' '$WZSEQ_REFERENCE' '$input_bam' >'${output_vcf%.gz}'
bgzip '${output_vcf%.gz}'
tabix -p vcf '$output_vcf'

~/tools/biscuit/master/biscuit/biscuit vcf2bed -k 0 -c -t cg '$output_vcf' | cut -f1-5 | LC_ALL=C sort -k1,1 -k2,2n -T pileup | gzip -c >'$output_vcf'_cg_bed.gz

# make a bw file, arbitrarily chose cov5
~/tools/biscuit/master/biscuit/biscuit vcf2bed -k 5 -t cg '$output_vcf' | LC_ALL=C sort -k1,1 -k2,2n -T pileup >'$output_vcf'_cg_cov5.bed
bedGraphToBigWig '$output_vcf'_cg_cov5.bed '$WZSEQ_REFERENCE'.fai '$output_vcf'_cg_cov5.bw
'
  jobname="biscuit_pileup_"$sname
}

function __wgbs_biscuit_pileup_nome {
  # input bamfn
  cmd="
set -xe
cd $base
mkdir -p pileup
biscuit pileup -r $WZSEQ_REFERENCE -N -i bam/$bamfn.bam -o pileup/$bamfn.vcf -q 28
bgzip pileup/$bamfn.vcf
tabix -p vcf pileup/$bamfn.vcf.gz
biscuit vcf2bed -k 0 -c -t hcg pileup/${bamfn}.vcf.gz | cut -f1-5 | LC_ALL=C sort -k1,1 -k2,2n -T pileup | gzip -c >pileup/${bamfn}.hcg.bed.gz
biscuit vcf2bed -k 10 -t hcg pileup/${bamfn}.vcf.gz | LC_ALL=C sort -k1,1 -k2,2n -T pileup >pileup/${bamfn}.hcg.cov10.bedg
biscuit vcf2bed -k 10 -t gch pileup/${bamfn}.vcf.gz | LC_ALL=C sort -k1,1 -k2,2n -T pileup >pileup/${bamfn}.gch.cov10.bedg
bedGraphToBigWig pileup/${bamfn}.hcg.cov10.bedg ${WZSEQ_REFERENCE}.fai pileup/${bamfn}.hcg.cov10.bw
bedGraphToBigWig pileup/${bamfn}.gch.cov10.bedg ${WZSEQ_REFERENCE}.fai pileup/${bamfn}.gch.cov10.bw
"
  jobname="biscuit_pileup_$bamfn"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 8 -memG 200 -ppn 28
}

function wgbs_biscuit_pileup {
  # wgbs_biscuit_pileup [-nome] do

  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d pileup ]] || mkdir pileup
  for bamfn0 in bam/*.bam; do
    bamfn=$(basename $bamfn0 .bam);
    # awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    # while read bamfn sread1 sread2; do
    __wgbs_biscuit_pileup
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

function __wgbs_vcf2tracks {

  # input: bamfn, items

  cmd="
set -xe
cd $base
mkdir -p tracks

for pt in ${items[@]}; do

  echo processing pileup type \$pt >&2

  biscuit vcf2bed -k 0 -c -t \${pt} pileup/${bamfn}.vcf.gz | cut -f1-5 | LC_ALL=C sort -k1,1 -k2,2n -T tracks | gzip -c >tracks/${bamfn}.\${pt}.bed.gz

  biscuit vcf2bed -k 10 -t \${pt} pileup/${bamfn}.vcf.gz | LC_ALL=C sort -k1,1 -k2,2n -T tracks >tracks/${bamfn}.\${pt}.bedg
  bedGraphToBigWig tracks/${bamfn}.\${pt}.bedg ${WZSEQ_REFERENCE}.fai tracks/${bamfn}.\${pt}.bw

  for i in 100000,100k 1000,1k 100,100; do 
    IFS=\",\"; set \$i; 
    echo processing \$i >&2
    bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w \$1 | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ | bedtools intersect -a - -b tracks/${bamfn}.\${pt}.bedg -wo | bedtools groupby -i - -g 1-3 -c 7 -o mean >tracks/${bamfn}.\${pt}.window\$2.bedg;

    bedGraphToBigWig tracks/${bamfn}.\${pt}.window\$2.bedg ${WZSEQ_REFERENCE}.fai tracks/${bamfn}.\${pt}.window\$2.bw
#rm -f tracks/${bamfn}.\${pt}.window\$2.bedg
    unset IFS
  done
done
"
  jobname="vcf2tracks_$bamfn"
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 200 -ppn 1
}

# wgbs_vcf2tracks [-nome] do
function wgbs_vcf2tracks {

  # compute 1) base-pair resolution methylation track 2) mean methylation in window track
  [[ -d tracks ]] || mkdir tracks
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ $1 == "-nome" ]] && items=(hcg gch) || items=(cg)
  for f in pileup/*.vcf.gz; do
    fn=$(readlink -f $f)
    bamfn=$(basename $f .vcf.gz)
    __wgbs_vcf2tracks
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wgbs_vcf2chtracks {

  # compute 1) base-pair resolution methylation track 2) mean methylation in window track
  [[ -d tracks ]] || mkdir tracks
  [[ -d pbs ]] || mkdir pbs
  items=(cg ch)
  base=$(pwd)
  for f in pileup/*.vcf.gz; do
    fn=$(readlink -f $f)
    bfn=$(basename $f .vcf.gz)
    cmd="
cd $base
for pt in ${items[@]}; do
  echo processing pileup type \$pt >&2
  biscuit vcf2bed -k 5 -t \${pt} pileup/${bfn}.vcf.gz | LC_ALL=C sort -k1,1 -k2,2n -T tracks > tracks/${bfn}.\${pt}.bedg
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
    jobname="vcf2chtracks_$bfn"
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
~/wzlib/Rutils/bin/bioinfo/methylKit.r diff -t 1,0 -b $base1,$base2 methylKit/$sname/$base1.methylKit,methylKit/$sname/$base2.methylKit -o methylKit/$sname/${base1}_vs_${base2}.tsv -g $WZSEQ_CGIBED_METHYLKIT
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

function wgbs_metilene {

  # output format:
  # chr start stop
  # q-value mean-methylation-difference #CpGs
  # p (Mann-Whitney U) p (2D KS) mean-g1 mean-g2
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d metilene ]] || mkdir metilene
  awk '/^\[/{p=0}/\[diffmeth\]/{p=1;next} p&&!/^$/' samples |
    while read sname vcf1 vcf2; do
      sname1=$(basename $vcf1 .vcf.gz)
      sname2=$(basename $vcf2 .vcf.gz)
      cmd="
cd $base

echo \"\$(date) Converting to bedgraph\"
biscuit vcf2bed -t cg $vcf1 >metilene/$sname1.bedg
biscuit vcf2bed -t cg $vcf2 >metilene/$sname2.bedg

echo \"\$(date) Merge\"
~/software/metilene/default/metilene_input.pl -in1 metilene/$sname1.bedg -in2 metilene/$sname2.bedg --h1 $sname1 --h2 $sname2 --out metilene/$sname.metilene
~/software/metilene/default/metilene_linux64 -t 4 -a $sname1 -b $sname2 metilene/$sname.metilene >metilene/$sname.metilene.result
sort -k4,4n metilene/$sname.metilene.result | head -n 1000 > metilene/$sname.metilene.result.top1k.tsv

echo \"\$(date) Clean\"
rm -f metilene/$sname1.bedg
rm -f metilene/$sname2.bedg
rm -f metilene/$sname.metilene
"
      jobname="metilene_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 4
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

## wgbs_repeat [-nome]
function wgbs_repeat {
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d rmsk ]] || mkdir rmsk
  [[ $1 == "-nome" ]] && items=(hcg gch) || items=(cg)
  for bfn0 in bam/*.bam; do
    bfn=$(basename $bfn0 .bam);
    # awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    # while read bfn; do
    cmd="
cd $base

for pt in ${items[@]}; do
  biscuit vcf2bed -t \${pt} -k 5 -c pileup/$bfn.vcf.gz > rmsk/$bfn.\${pt}.bedg
  bedtools map -a $WZSEQ_RMSK -b rmsk/$bfn.\${pt}.bedg -c 4 -o count,mean,collapse >rmsk/$bfn.\${pt}.bedg.rmsk
done
"
    jobname="rmsk_meth_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 4
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wgbs_repeat_diff {
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d rmsk ]] || mkdir rmsk
  awk '/^\[/{p=0}/\[diffmeth\]/{p=1;next} p&&!/^$/' samples |
    while read sname vcf1 vcf2; do
      base=$(pwd)
      [[ -d pbs ]] || mkdir pbs
      sname1=$(basename $vcf1 .vcf.gz)
      sname2=$(basename $vcf2 .vcf.gz)
      cmd="
cd $base
mkdir -p rmsk/$sname
biscuit vcf2bed -t cg -k 5 $vcf1 > rmsk/$sname/$sname1.cg.bedg
bedtools map -a $WZSEQ_RMSK -b rmsk/$sname/$sname1.cg.bedg -c 4 -o count,mean,collapse >rmsk/$sname/$sname1.cg.bedg.rmsk
biscuit vcf2bed -t cg -k 5 $vcf2 > rmsk/$sname/$sname2.cg.bedg
bedtools map -a $WZSEQ_RMSK -b rmsk/$sname/$sname2.cg.bedg -c 4 -o count,mean,collapse >rmsk/$sname/$sname2.cg.bedg.rmsk
paste rmsk/$sname/$sname1.cg.bedg.rmsk rmsk/$sname/$sname2.cg.bedg.rmsk >rmsk/$sname/merged.rmsk
awk -f wanding.awk -e '\$9!=\".\"&&\$19!=\".\"{print joinr(1,9)\"\\t\"\$10\"\\t\"\$19\"\\t\"\$20}' rmsk/$sname/merged.rmsk | wzstats Utest -c1 10 -c2 12 - -p 1-9,11 --outfc >rmsk/$sname/$sname.tsv

# awk '\$11<0.01&&\$12>1' rmsk/$sname/$sname.tsv | sort -k12,12nr > rmsk/$sname/$sname.top_hyper.tsv
# awk '\$11<0.01&&\$12<-1' rmsk/$sname/$sname.tsv | sort -k12,12n > rmsk/$sname/$sname.top_hypo.tsv

rm -f rmsk/$sname/$sname1.cg.bedg
rm -f rmsk/$sname/$sname1.cg.bedg.rmsk
rm -f rmsk/$sname/$sname2.cg.bedg
rm -f rmsk/$sname/$sname2.cg.bedg.rmsk
rm -f rmsk/$sname/merged.rmsk
"

      # output format
      # 1-7 basic repeat info
      # 8: mean beta of sample 1
      # 9: mean beta of sample 2
      # 10: p-value
      # 11: log2(fold change)
      # chr1 3000000 3002128 - L1_Mus3 LINE L1  5   0.9028  0.8968  0.917  -0.010
      # chr1 3003152 3003994 - L1Md_F  LINE L1  8   0.90775 0.94025 0.345   0.051
      # chr1 3004270 3005001 + L1_Rod  LINE L1  1   0.895   0.859   0.317   -0.059
      # chr1 3005570 3006764 + Lx9     LINE L1  3   0.95    0.957   0.663   0.011
      jobname="rmsk_diffmeth_$sname"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
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
 => (+) wzseq_bam_coverage

 (+) chipseq_bcp => (+) chipseq_macs2 => (+) chipseq_gem => (+) chipseq_sissrs
 (+) chipseq_HOMER => (+) chipseq_HOMER_peak

  => (+) chipseq_diffbind_window_with_internal_control
  => (+) chipseq_bam2track_with_internal_control
EOF
}

function pipeline_20170917_chipseqPE {

  pipeline_start=${!#}
  base=$(pwd);
  mkdir -p pbs;
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do

      pipeline_init

      pipeline_eval 1 __wzseq_bwa_mem

      fn=bam/${sname}.bam
      pipeline_eval 2 __wzseq_bam_coverage

#[[ ${!#} == "clean" ]] && pipeline_eval 3 __wzseq_clean_bam

    done
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
~/software/HOMER/default/bin/findPeaks HOMER/$targetname -style $targettype -o HOMER/$targetname/peaks.txt -i HOMER/$controlname
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

function __chipseq_macs2_noCtrl {
  : '
target_bam=bam/
hour=12; memG=10; ppn=1
pipeline_eval x __chipseq_macs2
'
  cmd='
cd '$base'
mkdir -p macs2
macs2 callpeak -t '$target_bam' -f BAM -g '$WZSEQ_MACS_SHORT' -n macs2/'$sname' -B --broad
'
  jobname="macs2_"$sname
}

function __chipseq_macs2 {
  : '
target_bam=bam/
control_bam=bam/
hour=12; memG=10; ppn=1
pipeline_eval x __chipseq_macs2
'
  cmd='
cd '$base'
mkdir -p macs2
macs2 callpeak -t '$target_bam' -c '$control_bam' -f BAM -g '$WZSEQ_MACS_SHORT' -n macs2/'$sname' -B --broad
'
  jobname="macs2_"$sname
}

function chipseq_macs2 {

  base=$(pwd);
  [[ -s macs2 ]] || mkdir macs2
  [[ -d pbs ]] || mkdir pbs

  awk '/^\[/{p=0}/\[peak\]/{p=1;next} p&&!/^$/' samples |
    while read targetname controlname targettype; do

      cmd="
cd $base
macs2 callpeak -t bam/$targetname.bam -c bam/$controlname.bam -f BAM -g $WZSEQ_MACS_SHORT -n macs2/$targetname -B --broad
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

function chipseq_diffbind_window_with_internal_control {
  [[ -d pbs ]] || mkdir pbs
  [[ -d diffbind ]] || mkdir diffbind

  # read in internal control
  declare -A sname2control
  while read sname controlval; do
    sname2control[$sname]=$controlval
  done <<EOF
$(awk '/^\[/{p=0}/\[internalcontrol\]/{p=1;next} p&&!/^$/' samples)
EOF

  base=$(pwd)
  awk '/^\[/{p=0}/\[diffbind\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 sname1 sname2; do
      odir=diffbind/${cond1}_vs_${cond2}
      # bam is usually unordered, would be very memory-consuming
      # here bam is first converted to bam and sorted
      [[ ${sname2control[$sname1]+x} ]] && control1=${sname2control[$sname1]} || control1=1
      [[ ${sname2control[$sname2]+x} ]] && control2=${sname2control[$sname2]} || control2=1
      cmd="
cd $base
[[ -d $odir ]] || mkdir $odir

bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w 100 -s 50 | sort -k1,1 -k2,2n -T $odir | bedtools intersect -a - -b <(bedtools bamtobed -i bam/$sname1.bam | awk '\$5>=20' | sort -k1,1 -k2,2n -T $odir) -wao -sorted | awk '{if (\$10==\".\") \$10=0; print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$10}' | bedtools groupby -i - -g 1-3 -c 4 -o sum >$odir/$sname1.bed
wc -l $odir/$sname1.bed
bedtools makewindows -g ${WZSEQ_REFERENCE}.fai -w 100 -s 50 | sort -k1,1 -k2,2n -T $odir | bedtools intersect -a - -b <(bedtools bamtobed -i bam/$sname2.bam | awk '\$5>=20' | sort -k1,1 -k2,2n -T $odir) -wao -sorted | awk '{if (\$10==\".\") \$10=0; print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$10}' | bedtools groupby -i - -g 1-3 -c 4 -o sum >$odir/$sname2.bed
wc -l $odir/$sname2.bed

# control fold change and maximum of the two
paste $odir/$sname1.bed $odir/$sname2.bed | awk -f wanding.awk -e '{rawma=max(\$4,\$8); n1=(\$4+1); n2=(\$8+1)*$control1/$control2; ma=max(n1,n2); mi=min(n1,n2); strand=n1<n2?\"+\":\"-\"; if (rawma>1000 && (ma/mi)>5) print \$1\"\t\"\$2\"\t\"\$3\"\t.\t0\t\"strand\"\t\"log(n2/n1)/log(2)\"\t\"\$4\"\t\"\$8}' | sort -k1,1 -k2,2n -T $odir | bedtools merge -s -i - -delim \";\" -c 7,7,8,9 -o mean,count,mean,mean >$odir/${cond1}_vs_${cond2}_diffbind.tsv
"
      jobname="diffbind_${cond1}_vs_${cond2}_with_internal_control"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 80 -ppn 7
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

##################################################
## normalize chipseq signal with internal control
##################################################
function chipseq_bam2track_with_internal_control {

  [[ -d pbs ]] || mkdir pbs
  [[ -d tracks ]] || mkdir tracks

  # read in internal control
  declare -A sname2control
  while read sname controlval; do
    sname2control[$sname]=$controlval
  done <<EOF
$(awk '/^\[/{p=0}/\[internalcontrol\]/{p=1;next} p&&!/^$/' samples)
EOF

  minmapq=10

  base=$(pwd)
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname fastq1 fastq2; do
      # bam is usually unordered, would be very memory-consuming
      # here bam is first converted to bam and sorted
      [[ ${sname2control[$sname]+x} ]] && control=${sname2control[$sname]} || control=1
      cmd="
cd $base
[[ -d tracks ]] || mkdir tracks

bedtools genomecov -ibam bam/$sname.bam -g ${WZSEQ_REFERENCE}.fai -bga -split | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4/$control}' | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ >tracks/${sname}.coverage.bedg
bedGraphToBigWig tracks/${sname}.coverage.bedg ${WZSEQ_REFERENCE}.fai tracks/${sname}.coverage.bw

# rm -f tracks/${sname}.coverage.bedg

samtools view -q $minmapq -b bam/$sname.bam | bedtools genomecov -ibam stdin -g ${WZSEQ_REFERENCE}.fai -bga -split | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4/$control}' | LC_ALL=C sort -k1,1 -k2,2n -T tracks/ >tracks/${sname}.coverage.q10.bedg
bedGraphToBigWig tracks/${sname}.coverage.q10.bedg ${WZSEQ_REFERENCE}.fai tracks/${sname}.coverage.q${minmapq}.bw

# rm -f tracks/${sname}.coverage.q10.bedg
"
      jobname="chipseq_bam2track_with_internal_control_${sname}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 5 -ppn 2
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}


################################################################################
# HiC pipeline
################################################################################
# TODO: http://homer.salk.edu/homer/interactions/HiCtagDirectory.html

function examplepipeline_hic {
cat<<- EOF
=== pipeline ===
(+) hic_hicpro
EOF
}

## generate fragments Hind-III A^AGCTT, Mbol: ^GATC
## utils/digest_genome.py -r A^AGCTT -o mm9_hindiii.bed mm9.fasta
## [~/software/HiC_Pro/default/annotation]$ ../bin/utils/digest_genome.py -r ^GATC -o Mbol_resfrag_hg19.bed ~/references/hg19_noContig/hg19_noContig.fa
## 
## java -Xmx8g -jar ~/software/juicer/juicer_tools_linux_0.8.jar pre hicpro/O3_57_CTCF_HiCHP-41521568_output/hic_results/data/O3_57_CTCF_HiCHP-41521568/O3_57_CTCF_HiCHP-41521568_allValidPairs.forjuicer hicpro/O3_57_CTCF_HiCHP-41521568_output/hic_results/data/O3_57_CTCF_HiCHP-41521568/O3_57_CTCF_HiCHP-41521568_allValidPairs.forjuicer.hic hg19
function hic_hicpro {

  base=$(pwd)
  [[ -d pbs ]] || mkdir -p pbs;
  [[ -d hicpro ]] || mkdir -p hicpro;
  awk '/^\[/{p=0}/\[hicpro\]/{p=1;next} p&&!/^$/' samples |
    while read sname fastqs; do
      ## need config file with name $sname.hicpro
      mkdir -p hicpro/${sname}_input/$sname;

      ## check whether to overwrite existing output
      if [[ -d hicpro/${sname}_output ]]; then
        read -p "Do you wish to replace existing output: hicpro/${sname}_output [yn]? " yn </dev/tty
        case $yn in
          [Yy]* ) rm -rf hicpro/${sname}_output; ;;
          [Nn]* ) echo "Skip $sname."; continue;;
              * ) echo "Please answer yes or no."; exit;;
        esac
      fi

      ## hicpro_config file, sample-specific config file has higher precedence
      for fastq in ${fastqs//,/ }; do
        ln -sf `rf fastq/$fastq` hicpro/${sname}_input/$sname/$fastq;
      done
      configfile=hicpro/${sname}_hicpro_config
      if [[ ! -e $configfile ]]; then configfile=hicpro/hicpro_config; fi;

      ## write pbs file
    cmd="
cd $base
/home/wanding.zhou/software/HiC_Pro/default/bin/HiC-Pro -c $configfile -i hicpro/${sname}_input -o hicpro/${sname}_output
mkdir -p hicpro/hic
awk '{gsub(/chr/,\"\",\$2); gsub(/chr/,\"\",\$5); print 0,\$2,\$3,0,0,\$5,\$6,1,60}' hicpro/${sname}_output/hic_results/data/${sname}/${sname}_allValidPairs | sort -k1,1 -k6,6 >hicpro/hic/${sname}.juicer
java -Xmx8g -jar ~/software/juicer/juicer_tools_linux_0.8.jar pre hicpro/hic/${sname}.juicer hicpro/hic/${sname}.hic $WZSEQ_REFVERSION
rm -f hicpro/hic/${sname}.juicer
"
      jobname="hicpro_${sname}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 100 -ppn 28
      [[ ${!#} == "do" ]] && qsub $pbsfn
    done
}

function __hic_hicpro {

  ## sname, fastq
  ## need config file with name $sname.hicpro
  ## note that in the config file we specify N_CPU=20, this is hard-coded
  ## hour=48; memG=100; ppn=20; queue=shortq
  mkdir -p hicpro/${sname}_input/$sname;

  sed "s/WZREPLACE_N_CPU/"$ppn"/" $WZSEQ_HICPRO_CONFIG >hicpro/config
  
  ## check whether to overwrite existing output
  rm -rf hicpro/${sname}_output
  ## hicpro_config file, sample-specific config file has higher precedence
  for fastq in ${fastqs//,/ }; do
    ln -sf `readlink -f fastq/$fastq` hicpro/${sname}_input/$sname/
  done
  
  ## write pbs file
  cmd='
cd '$base'
/home/wanding.zhou/software/HiC_Pro/default/bin/HiC-Pro -c hicpro/config -i hicpro/'${sname}'_input -o hicpro/'${sname}'_output
mkdir -p hicpro/hic
awk '\''{gsub(/chr/,"",$2); gsub(/chr/,"",$5); print 0,$2,$3,0,0,$5,$6,1,60}'\'' hicpro/'${sname}'_output/hic_results/data/'${sname}'/'${sname}'_allValidPairs | sort -k1,1 -k6,6 >hicpro/hic/'${sname}'.juicer
java -Xmx8g -jar ~/software/juicer/juicer_tools_linux_0.8.jar pre hicpro/hic/'${sname}'.juicer hicpro/hic/'${sname}'.hic '$WZSEQ_REFVERSION'
rm -f hicpro/hic/'${sname}'.juicer
'
  jobname="hicpro_"${sname}
  pbsfn=$base/pbs/$jobname.pbs
}

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

 => (+) edit samples [alignment]; (+) editrnaseq_tophat2_firststrand / (+) rnaseq_tophat2 / (+) rnaseq_STAR / (+) rnaseq_gsnap / (+) rnaseq_subjunc / (+) rnaseq_mapsplice / (+) rnaseq_hisat2

 => (+) rnaseq_kallisto => rnaseq_kallisto_diff

 => (+) wzseq_qualimap rnaseq_se/rnaseq_pe_stranded

 => [o] rnaseq_splitallele (for ASE) => (+) wzseq_bam_coverage
 => (+) rnaseq_allelomePro

 => (+) rnaseq_cufflinks => (+) rnaseq_cuffmerge => (+) edit samples [diffexp] ; rnaseq_cuffdiff

 => (+) rnaseq_featureCounts => (+) rnaseq_edgeR => (+) rnaseq_DESeq2 => (+) rnaseq_edgeR_rmsk

 => (+) rnaseq_splitstrand_pe / rnaseq_splitstrand_se => (+) rnaseq_count_rmsk_stranded => (+) rnaseq_count_rmsk_stranded_edgeR
 => (+) rnaseq_count_rmsk_unstranded => (+) rnaseq_count_rmsk_unstranded_edgeR

 => (+) rnaseq_dexseq

EOF
}

function pipeline_20170818_rnaseq_SE_stranded {

  pipeline_start=${!#}
  base=$(pwd)
  mkdir -p pbs
  all_aln_jobs=""
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
  while read sname sread1 sread2; do

    pipeline_init

    pipeline_eval 1 __rnaseq_hisat2
    all_aln_jobs="$all_aln_jobs:$jobid"

    pipeline_eval 2 __rnaseq_splitstrand_se

    pipeline_eval 3 __rnaseq_count_rmsk_stranded

  done 

  allbams="bam/*.bam"
  stranded="-s 1"
  pairEnd=""
  depend="-W depend=afterok:$all_aln_jobs"
  echo $depend
  pipeline_eval 4 __rnaseq_featureCounts
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
# fr-firststrand:dUTP, NSR, NNSR Same as above except we enforce the rule that
# the right-most end of the fragment (in transcript coordinates) is the first
# sequenced (or only sequenced for single-end reads).
# F2R1 for forward transcript
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

function __rnaseq_hisat2_SE {

  # hour=24; memG=30; ppn=14; queue=shortq
  cmd='
set -xe
cd '$base'
mkdir -p bam

export PATH='$PATH':/primary/vari/software/hisat/default
hisat2 -x '$WZSEQ_HISAT2_INDEX' -U '$fastq' -p '$ppn' | samtools view -bS - | samtools sort -T bam/'$sname'.tmp -O bam -o bam/'$sname'.bam
samtools index bam/'$sname'.bam
samtools flagstat bam/'$sname'.bam >bam/'$sname'.bam.flagstat
'
  jobname="hisat_${sname}_SE"
}

function __rnaseq_hisat2_PE {

  cmd='
set -xe
cd '$base'
mkdir -p bam

export PATH=$PATH:/primary/vari/software/hisat/default
hisat2 -x '$WZSEQ_HISAT2_INDEX' -1 '$fastq1' -2 '$fastq2' -p '$ppn' | samtools view -bS - | samtools sort -T bam/'$sname'.tmp -O bam -o bam/'$sname'.bam
samtools index bam/'$sname'.bam
samtools flagstat bam/'$sname'.bam >bam/'$sname'.bam.flagstat
'
  jobname="hisat_${sname}_PE"
}

function rnaseq_hisat2 {
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read sname sread1 sread2; do
      fastq=fastq/$sread1
      ppn=28
      queue=shortq
      __rnaseq_hisat2_SE
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 240 -ppn $ppn
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
samtools flagstat bam/$sname.bam >bam/$sname.bam.flagstat
rm -f bam/$sname.sam
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

# subjunc, subread (if only differential gene expression is concerned, faster), here I
# only use subjunc, which is supposed to be more optimal
# need to build subread index first (for subjunc)
# bin/subread-buildindex -o ~/references/hg19/subread/hg19 ~/references/hg19/hg19.fa
# other options: -S fr -d <min_insert_size> -D <max_insert_size>
function rnaseq_subjunc() {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam ]] || mkdir bam

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

#################################################
## section 2: differential expression
#################################################

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

###### kallisto
###### index for mm10
# make transcript fasta using the tophat's gtf_to_fasta and some postprocessing
# ~/software/tophat2/default/gtf_to_fasta ~/references/mm10/gtf/Mus_musculus.GRCm38.82.gtf.UCSCnaming ~/references/mm10/mm10.fa ~/references/mm10/kallisto/mm10.transcripts.fa
# awk 'match($0,/>(\S+)\s+(\w+)/,a){print ">"a[2];next;}1' ~/references/mm10/kallisto/mm10.transcripts.fa | gzip -c >~/references/mm10/kallisto/mm10.transcripts.fa.gz
# rm ~/references/mm10/kallisto/mm10.transcripts.fa
# ~/software/kallisto/kallisto_linux-v0.42.4/kallisto index ~/references/mm10/kallisto/mm10.transcripts.fa.gz -i ~/references/mm10/kallisto/mm10.kallisto
# 
###### index for human
# cd ~/references/hg19
# ~/software/tophat2/default/gtf_to_fasta gtf/Homo_sapiens.GRCh37.75.gtf.UCSCnaming hg19.fa kallisto/hg19.transcripts.fa
# awk 'match($0,/>(\S+)\s+(\w+)/,a){print ">"a[2];next;}1' kallisto/hg19.transcripts.fa | gzip -c >kallisto/hg19.transcripts.fa.gz
# ~/software/kallisto/kallisto_linux-v0.42.4/kallisto index kallisto/hg19.transcripts.fa.gz -i kallisto/hg19.kallisto
function rnaseq_kallisto {
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d kallisto ]] || mkdir kallisto

  while read sname fastq1 fastq2; do
    if [[ "$fastq2" == "." ]]; then
      input2=""
      ## the following is requried, read length is 75 sd 10. this needs to be tuned for each library
      additional_option="--single -l 75 -s 10"
    else
      input2="fastq/$fastq2"
      additional_option=""
    fi
    cmd="
cd $base
[[ -d kallisto/$sname ]] || mkdir kallisto/$sname
~/software/kallisto/kallisto_linux-v0.42.4/kallisto quant -i $WZSEQ_KALLISTO_INDEX $additional_option -o kallisto/$sname fastq/$fastq1 $input2 -t 4
~/wzlib/pyutils/wzseqtk.py ensembl2name --transcript -g $WZSEQ_GTF_ENSEMBL_UCSCNAMING -i kallisto/$sname/abundance.tsv -o kallisto/$sname/abundance.anno.tsv
"
    jobname="kallisto_$sname"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 4
    [[ $1 == "do" ]] && qsub $pbsfn
  done < <(awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples)
}

function rnaseq_kallisto_diff {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      snames1=""
      n1=0
      for f in $(echo $bams1 | tr ',' ' '); do
        [[ -d kallisto/$(basename $f .bam) ]] || continue;
        [[ ${#snames1} == 0 ]] || snames1=$snames1",";
        snames1=$snames1$(basename $f .bam);
        ((n1++))
      done
      
      snames2=""
      n2=0
      for f in $(echo $bams2 | tr ',' ' '); do
        [[ -d kallisto/$(basename $f .bam) ]] || continue;
        [[ ${#snames2} == 0 ]] || snames2=$snames2",";
        snames2=$snames2$(basename $f .bam);
        ((n2++))
      done

      if [[ $n1 -gt 1 ]] && [[ $n2 -gt 1 ]]; then # voom requires more than 1 replica
        cmd="
cd $base
~/wzlib/Rutils/bin/bioinfo/kallisto_voom.r -a $cond1 -b $cond2 -A $snames1 -B $snames2 -o kallisto/diff_${cond1}_vs_${cond2}.tsv
~/wzlib/pyutils/wzseqtk.py ensembl2name --transcript -g $WZSEQ_GTF_ENSEMBL_UCSCNAMING -i kallisto/diff_${cond1}_vs_${cond2}.tsv -o kallisto/diff_${cond1}_vs_${cond2}.anno.tsv
"
        jobname="kallisto_diff_${cond1}_vs_${cond2}"
        pbsfn=$base/pbs/$jobname.pbs
        pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 3 -ppn 1
        [[ $1 == "do" ]] && qsub $pbsfn
      fi
    done
}

# TODO: sailfish
# TODO: salmon

function __rnaseq_featureCounts {

  # input: base, pairEnd, stranded, allbams
  prog=~/tools/subread/default/bin/featureCounts

  # -p paired-end
  # -g goup by gene_id in GTF
  # -t exon, together with -g, it means count exon (as feature) and gene_id (as meta_feature)
  # -T: number of threads
  # -Q: min mapping quality 20

  # optional:
  # -P -d 50 -D 600 set insert size range to [50,600]
  # -O allowMultiOverlap
  # -M allow multi-mapping
  # -F SAF or GTF, format of annotation file

  cmd="
set -xe
cd $base
mkdir -p featureCounts

# count genes
$prog $pairEnd $stranded -T 5 -t exon -g gene_id -a $WZSEQ_GTF_ENSEMBL_UCSCNAMING -o featureCounts/genes.tsv $allbams --primary -Q 20 --ignoreDup

~/wzlib/pyutils/wzseqtk.py cnt2rpkm -i featureCounts/genes.tsv | ~/wzlib/pyutils/wzseqtk.py ensembl2name --gene -g $WZSEQ_GTF_ENSEMBL_UCSCNAMING -H >featureCounts/genes.rpkm.tsv

# count repeat loci, -f suppresses meta-feature counts
$prog $pairEnd $stranded -T 5 -t exon -g gene_id -f -a $WZSEQ_RMSK_GTF -o featureCounts/rmsk_loci.tsv $allbams --primary -Q 20 --ignoreDup

~/wzlib/pyutils/wzseqtk.py cnt2rpkm -i featureCounts/rmsk_loci.tsv >featureCounts/rmsk_loci.rpkm.tsv

awk -f wanding.awk -e 'NR==FNR{\$2+=1; a[\$1\":\"\$2\"-\"\$3]=joinr(4,NF)}NR!=FNR{if(\$1==\"ID\") print \"strand\tfam1\tfam2\tfam3\t\"\$0; else {ak=\$2\":\"\$3\"-\"\$4; print a[ak],\$0;}}' $WZSEQ_RMSK featureCounts/rmsk_loci.rpkm.tsv >featureCounts/rmsk_loci.rpkm.fam.tsv
"
  jobname="featurecounts_"$(basename $base)
  pbsfn=$base/pbs/$jobname.pbs
  pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 5
}

# R-based methods
# this doesn't work very well for repeats, since they lost the family information
# the most generic version
function rnaseq_featureCounts {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  grep '\[experiment\] single-end' samples && pairEnd="" || pairEnd="-P"
  # grep '\[experiment\] unstranded' samples && stranded="-s 0" || stranded="-s 1"
  # stranded="-s 0" # unstranded
  # stranded="-s 1" # for forward-stranded
  # stranded="-s 2" # for reverse-stranded (most of Minmin's library)
  [[ -z "$pairEnd" ]] && stranded=""

  allbams="bam/*.bam"
  [[ -d bam_allele ]] && allbams=$allbams" bam_allele/*.bam"
  __rnaseq_featureCounts
  [[ ${!#} == "do" ]] && qsub $pbsfn
}

function rnaseq_featureCounts_SE_stranded {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  pairEnd=""
  stranded="-s 1"
  allbams="bam/*.bam"
  [[ -d bam_allele ]] && allbams=$allbams" bam_allele/*.bam"
  __rnaseq_featureCounts
  [[ ${!#} == "do" ]] && qsub $pbsfn
}


function rnaseq_featureCounts_SE_revstranded {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  pairEnd=""
  stranded="-s 2"
  allbams="bam/*.bam"
  [[ -d bam_allele ]] && allbams=$allbams" bam_allele/*.bam"
  __rnaseq_featureCounts
  [[ ${!#} == "do" ]] && qsub $pbsfn
}

function rnaseq_edgeR {
  # 2016-02-16 update:
  # switch from GenomicAlignment to featureCounts
  ## original command: ~/wzlib/Rutils/bin/bioinfo/edgeR.r -g $WZSEQ_REFVERSION -G $WZSEQ_GTF_ENSEMBL -a $cond1 -b $cond2 -A $bams1 -B $bams2 -o edgeR/${cond1}_vs_${cond2}_diffexp.tsv 2> edgeR/${cond1}_vs_${cond2}_diffexp.log
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
~/wzlib/Rutils/bin/bioinfo/edgeR_featureCounts.r -G $WZSEQ_GTF_ENSEMBL_UCSCNAMING -a $cond1 -b $cond2 -A $bams1 -B $bams2 -o edgeR/${cond1}_vs_${cond2}_diffexp.tsv 2> edgeR/${cond1}_vs_${cond2}_diffexp.log
~/wzlib/pyutils/wzseqtk.py ensembl2name --gene -g $WZSEQ_GTF_ENSEMBL_UCSCNAMING -i edgeR/${cond1}_vs_${cond2}_diffexp.tsv -o edgeR/${cond1}_vs_${cond2}_diffexp.anno.tsv
"
      jobname="edgeR_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 20 -ppn 2
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

function rnaseq_edgeR_rmsk {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d featureCounts ]] || mkdir featureCounts
  [[ -d edgeR/rmsk ]] || mkdir -p edgeR/rmsk
  [[ -d edgeR/rmsk_categories ]] || mkdir -p edgeR/rmsk_categories
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      cmd="
cd $base
~/wzlib/Rutils/bin/bioinfo/edgeR_featureCounts.r -s 7 -F featureCounts/rmsk_loci.tsv -a $cond1 -b $cond2 -A $bams1 -B $bams2 -o edgeR/rmsk/${cond1}_vs_${cond2}.tsv 2> edgeR/rmsk/${cond1}_vs_${cond2}.log
~/wzlib/Rutils/bin/bioinfo/edgeR_featureCounts.r -s 3 -F featureCounts/rmsk_categories.tsv -a $cond1 -b $cond2 -A $bams1 -B $bams2 -o edgeR/rmsk_categories/${cond1}_vs_${cond2}.tsv 2> edgeR/rmsk_categories/${cond1}_vs_${cond2}.log
"
      jobname="edgeR_rmsk_${cond1}_vs_${cond2}"
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

# TODO: EBSeq, Voom
##############################
## section 3: repeats
##############################
function rnaseq_splitstrand_pe {
  # split stranded RNAseq into 2 strands
  # that only works for paired-end
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d stranded ]] || mkdir stranded
  for bam in bam/*.bam; do
    sname=$(basename $bam .bam)
    cmd="
cd $base
minmapq=10

# first read, positive strand
samtools view -H bam/${sname}.bam > stranded/${sname}_p.sam
samtools view -q \$minmapq -f 0x50 -F 0x120 bam/${sname}.bam >> stranded/${sname}_p.sam
# second read, reverse strand
samtools view -q \$minmapq -f 0xa0 -F 0x110 bam/${sname}.bam >> stranded/${sname}_p.sam
samtools view -b stranded/${sname}_p.sam | samtools sort -o stranded/${sname}_p.bam -O bam -T stranded/${sname}_tmp

# first read, reverse strand
samtools view -H bam/${sname}.bam > stranded/${sname}_r.sam
samtools view -q \$minmapq -f 0x60 -F 0x110 bam/${sname}.bam >> stranded/${sname}_r.sam
# second read, positive strand
samtools view -q \$minmapq -f 0x90 -F 0x120 bam/${sname}.bam >> stranded/${sname}_r.sam
samtools view -b stranded/${sname}_r.sam | samtools sort -o stranded/${sname}_r.bam -O bam -T stranded/${sname}_tmp

bedtools genomecov -ibam stranded/${sname}_p.bam -g ${WZSEQ_REFERENCE}.fai -bga -split | LC_COLLATE=C sort -k1,1 -k2,2n -T stranded/ >stranded/${sname}_p.bedg
bedGraphToBigWig stranded/${sname}_p.bedg ${WZSEQ_REFERENCE}.fai stranded/${sname}_p.bw

bedtools genomecov -ibam stranded/${sname}_r.bam -g ${WZSEQ_REFERENCE}.fai -bga -split | LC_COLLATE=C sort -k1,1 -k2,2n -T stranded/ | awk -F\"\\t\" -v OFS=\"\t\" '{print \$1,\$2,\$3,-\$4}' >stranded/${sname}_r.bedg
bedGraphToBigWig stranded/${sname}_r.bedg ${WZSEQ_REFERENCE}.fai stranded/${sname}_r.bw
rm -f stranded/${sname}_p.sam stranded/${sname}_r.sam
rm -f stranded/${sname}_p.bam stranded/${sname}_r.bam
"
    jobname="splitstrand_${sname}"
    pbsfn=pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function __rnaseq_splitstrand_se {
  # hour=12; memG=10; ppn=1; queue=shortq
  # input: base, sname
  cmd='
set -xe
cd '$base'
minmapq=10
mkdir -p stranded

# first read, positive strand
samtools view -H bam/'${sname}'.bam > stranded/'${sname}'_p.sam
samtools view -q $minmapq -F 0x110 bam/'${sname}'.bam >> stranded/'${sname}'_p.sam
samtools view -b stranded/'${sname}'_p.sam | samtools sort -o stranded/'${sname}'_p.bam -O bam -T stranded/'${sname}'_tmp

# first read, reverse strand
samtools view -H bam/'${sname}'.bam > stranded/'${sname}'_r.sam
samtools view -q $minmapq -f 0x10 -F 0x100 bam/'${sname}'.bam >> stranded/'${sname}'_r.sam
samtools view -b stranded/'${sname}'_r.sam | samtools sort -o stranded/'${sname}'_r.bam -O bam -T stranded/'${sname}'_tmp

totalreads=$(awk -F" " '\''NR==1{print $1}'\'' bam/'${sname}'.bam.flagstat)

bedtools genomecov -ibam stranded/'${sname}'_p.bam -g '${WZSEQ_REFERENCE}'.fai -bga -split | LC_COLLATE=C sort -k1,1 -k2,2n -T stranded/ >stranded/'${sname}'_p.bedg
awk -v totalreads=$totalreads '\''{print $1,$2,$3,$4/totalreads*10000000}'\'' stranded/'${sname}'_p.bedg >stranded/'${sname}'_p_norm.bedg
bedGraphToBigWig stranded/'${sname}'_p_norm.bedg '${WZSEQ_REFERENCE}'.fai stranded/'${sname}'_p_norm.bw
bedGraphToBigWig stranded/'${sname}'_p.bedg '${WZSEQ_REFERENCE}'.fai stranded/'${sname}'_p.bw

## negative strand always have negative counts
bedtools genomecov -ibam stranded/'${sname}'_r.bam -g '${WZSEQ_REFERENCE}'.fai -bga -split | LC_COLLATE=C sort -k1,1 -k2,2n -T stranded/ | awk -F"\t" -v OFS="\t" '\''{print $1,$2,$3,-$4}'\'' >stranded/'${sname}'_r.bedg
awk -v totalreads=$totalreads '\''{print $1,$2,$3,$4/totalreads*10000000}'\'' stranded/'${sname}'_r.bedg >stranded/'${sname}'_r_norm.bedg
bedGraphToBigWig stranded/'${sname}'_r_norm.bedg '${WZSEQ_REFERENCE}'.fai stranded/'${sname}'_r_norm.bw
bedGraphToBigWig stranded/'${sname}'_r.bedg '${WZSEQ_REFERENCE}'.fai stranded/'${sname}'_r.bw

rm -f stranded/'${sname}'_p.sam stranded/'${sname}'_r.sam
rm -f stranded/'${sname}'_p.bam stranded/'${sname}'_r.bam
'
  jobname="splitstrand_${sname}_se"
}

function rnaseq_splitstrand_se {
  # split stranded RNAseq into 2 strands
  # that only works for single-end
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  for bam in bam/*.bam; do
    sname=$(basename $bam .bam)
    __rnaseq_splitstrand_se
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

## the following restrict to intergenic ERVs
## for f in rmsk/*.tsv; do echo $f; bedtools intersect -a <(sed -n '2,$p' $f) -b ERVIntergenic.bed -sorted -wo | cut -f1-13 | cat <(echo -e "chrm\tbeg\tend\tstr\tclass1\tclass2\tclass3\ttlen\tBaseNeg\tBasePos\tRPKMneg\tRPKMpos\tRPKM") - >rmsk/$(basename $f).LTRintergenic.bed; done
## for f in rmsk/*.tsv; do echo $f; bedtools intersect -a <(sed -n '2,$p' $f) -b ERVIntergenic.bed -sorted -wo | cut -f1-13 | awk -f wanding.awk -e '$9<0 && $10>0 && min(-$9,$10)/max(-$9,$10)>=0.5' | cat <(echo -e "chrm\tbeg\tend\tstr\tclass1\tclass2\tclass3\ttlen\tBaseNeg\tBasePos\tRPKMneg\tRPKMpos\tRPKM") - >rmsk/$(basename $f).LTRintergenic.doubleStranded.bed; done

function __rnaseq_count_rmsk_stranded {

  # input: sname, base
  # hour=24; memG=20; ppn=2;
  cmd='
set -xe
cd '$base'
mkdir -p rmsk

bedtools intersect -a '$WZSEQ_RMSK' -b stranded/'${sname}'_r.bedg -wao -sorted | awk -f wanding.awk -e '\''{print joinr(1,7)"\t"$11*$12;}'\'' | bedtools groupby -g 1-7 -c 8 -o sum > rmsk/'${sname}'_r.tsv

## find all cummulative base counts
## all bases mapped to positive strand
all_p=$(awk '\''{a+=($3-$2)*$4}END{print a}'\'' stranded/'${sname}'_p.bedg)
## all bases mapped to negative strand
all_r=$(awk '\''{a+=($3-$2)*$4}END{print a}'\'' stranded/'${sname}'_r.bedg)

nmap=$(($all_p-$all_r))

echo "cumulative base cnt: $nmap"

bedtools intersect -a '$WZSEQ_RMSK' -b stranded/'${sname}'_p.bedg -wao -sorted | awk -f wanding.awk -e '\''{print joinr(1,7)"\t"$11*$12;}'\'' | bedtools groupby -g 1-7 -c 8 -o sum > rmsk/'${sname}'_p.tsv

paste rmsk/'${sname}'_p.tsv rmsk/'${sname}'_r.tsv | awk '\''$2==$10 && $5==$13'\'' | cut -f1-7,8,16 | awk -v alln=$nmap '\''BEGIN{print "chrm\tbeg\tend\tstrand\tcat1\tcat2\tcat3\ttlen\tposBaseCnt\tnegBaseCnt\tposRPKM\tnegRPKM\tRPKM"}{tlen=$3-$2; pp=$8/tlen*1000/alln*1000000; rr=$9/tlen*1000/alln*1000000; print $0"\t"tlen"\t"pp"\t"rr"\t"pp-rr}'\'' >rmsk/'$sname'.tsv

## count category
awk -v alln=$nmap '\''{n=$8-$9; a[$5]+=n; b[$6]+=n; c[$7]+=n;} END{print "Genome\t0\t"alln"\t1.0"; for (i in a) {print i"\t1\t"a[i]"\t"a[i]/alln} for(i in b){print i"\t2\t"b[i]"\t"b[i]/alln} for(i in c){print i"\t3\t"c[i]"\t"c[i]/alln}}'\'' rmsk/'$sname'.tsv | sort -k2,2n -k1,1 >rmsk/'$sname'.tsv.categories

rm -f rmsk/'${sname}'_p.tsv rmsk/'${sname}'_r.tsv
'
  # for locating double-strand transcription
  # awk -v OFS="\t" -f wanding.awk -e 'max($14,-$15)>0 && min($14,-$15)/max($14,-$15)>0.5 && min($14,-$15)>100{print $_"\t"($14-$15)/($8-$9+10)}' merged.rmsk.bed | sort -k16,16nr >merged.rmsk.bed.double.up
  jobname="rmsk_${sname}"
}

function rnaseq_count_rmsk_stranded {
  # require stranded/, this actually gives the cumulative
  # base counts, not the read counts
  # equivalent to RPKM because base counts are normalized by total base counts
  # both numerator and denominator are scaled by the read length
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  for bedg_r in stranded/*_r.bedg; do
    sname=$(basename $bedg_r _r.bedg)
    __rnaseq_count_rmsk_stranded
    pbsfn=$base/pbs/$jobname.pbs
    hour=24; memG=20; ppn=2
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour $hour -memG $memG -ppn $ppn
    [[ $1 == "do" ]] && qsub $pbsfn
  done
}

function rnaseq_count_rmsk_stranded_edgeR {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d rmsk ]] || mkdir rmsk
  [[ -d rmsk/diff ]] || mkdir rmsk/diff
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      # skip non-replicated designs
      # Remark: At non-replicated design, we would sometimes have
      # "f() values at end points not of opposite sign" error
      # for "deviance" estimate of common dispersion
      # othertimes, it may work (usually for genes but not rmsk)
      # [[ -n $(grep -o "," <<< "$bams1") ]] || continue
      # [[ -n $(grep -o "," <<< "$bams2") ]] || continue

      # make sure the tsv all exists
      allexist=1
      tsv1=""
      for i in ${bams1//,/ }; do
        _tsv=rmsk/$(basename $i .bam).tsv 
        if [[ ! -e $_tsv ]]; then
          allexist=0;
          break;
        fi
        [[ -n $tsv1 ]] && tsv1=$tsv1","
        tsv1=$tsv1""$_tsv;
      done

      tsv2=""
      for i in ${bams2//,/ }; do
        _tsv=rmsk/$(basename $i .bam).tsv
	      if [[ ! -e $_tsv ]]; then
          allexist=0;
          break;
        fi
        [[ -n $tsv2 ]] && tsv2=$tsv2","
        tsv2=$tsv2""$_tsv;
      done

      [[ $allexist == 0 ]] && continue;

      # compare rmsk
      cmd="
cd $base
~/wzlib/Rutils/bin/bioinfo/edgeR_rmsk.r -a $cond1 -b $cond2 -A $tsv1 -B $tsv2 -o rmsk/diff/${cond1}_vs_${cond2}.diff.tsv
"
      jobname="rmsk_diff_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn

    done
}

function rnaseq_count_rmsk_unstranded {

  # this count the cumulative base coverage
  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d rmsk ]] || mkdir rmsk
  for bam in bam/*.bam; do
    sname=$(basename $bam .bam)
    cmd="
cd $base
minmapq=10
samtools view -q \$minmapq -b $bam | bedtools genomecov -ibam - -g ${WZSEQ_REFERENCE}.fai -bga -split | LC_COLLATE=C sort -k1,1 -k2,2n -T rmsk/ >rmsk/${sname}.bedg

bedGraphToBigWig rmsk/${sname}.bedg ${WZSEQ_REFERENCE}.fai rmsk/${sname}.bw

all=\$(awk '{a+=(\$3-\$2)*\$4}END{print a}' rmsk/${sname}.bedg)

bedtools intersect -a $WZSEQ_RMSK -b rmsk/${sname}.bedg -wao -sorted | awk -f wanding.awk -e '{print joinr(1,7)\"\\t\"\$11*\$12;}' | bedtools groupby -g 1-7 -c 8 -o sum | awk -v alln=\$all 'BEGIN{print \"chrm\tbeg\tend\tstrand\tcat1\tcat2\tcat3\ttlen\tbaseCnt\tRPKM\"}{tlen=\$3-\$2; rpkm=\$8/tlen*1000/alln*1000000; print \$0\"\t\"tlen\"\t\"rpkm}' > rmsk/${sname}.tsv

## count category 
awk -v all=\$all '{a[\$5]+=\$8; b[\$6]+=\$8; c[\$7]+=\$8;} END{print \"Genome\t0\t\"all\"\t1.0\"; for (i in a) {print i\"\t1\t\"a[i]\"\t\"a[i]/all} for(i in b){print i\"\t2\t\"b[i]\"\t\"b[i]/all} for(i in c){print i\"\t3\t\"c[i]\"\t\"c[i]/all}}' rmsk/$sname.tsv | sort -k2,2n -k1,1 >rmsk/$sname.tsv.categories

rm -f rmsk/${sname}.bedg
"
    # The following only count reads
    # bedtools coverage -a $WZSEQ_RMSK -b bam/${bam}.bam -sorted -split -counts > repeatmasker/$bam.rmsk
    # awk '{a[\$5]+=\$8;b[\$6]+=\$8;c[\$7]+=\$8}END{for(i in a){print "1\t"i"\t"a[i]}; for(i in b){print "2\t"i"\t"b[i]}; for(i in c){print "3\t"i"\t"c[i]}}' repeatmasker/${bam}.rmsk | sort -k1,1 -k2,2nr > repeatmasker/${bam}.rmsk.cat
    jobname="rmsk_${sname}"
    pbsfn=pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 2 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function rnaseq_count_rmsk_unstranded_edgeR {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d rmsk ]] || mkdir rmsk
  [[ -d rmsk/diff ]] || mkdir rmsk/diff
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      # skip non-replicated designs
      # Remark: At non-replicated design, we would sometimes have
      # "f() values at end points not of opposite sign" error
      # for "deviance" estimate of common dispersion
      # othertimes, it may work (usually for genes but not rmsk)
      # [[ -n $(grep -o "," <<< "$bams1") ]] || continue
      # [[ -n $(grep -o "," <<< "$bams2") ]] || continue

      # make sure the tsv all exists
      allexist=1
      tsv1=""
      for i in ${bams1//,/ }; do
        _tsv=rmsk/$(basename $i .bam).tsv 
        if [[ ! -e $_tsv ]]; then
          allexist=0;
          break;
        fi
        [[ -n $tsv1 ]] && tsv1=$tsv1","
        tsv1=$tsv1""$_tsv;
      done

      tsv2=""
      for i in ${bams2//,/ }; do
        _tsv=rmsk/$(basename $i .bam).tsv
	      if [[ ! -e $_tsv ]]; then
          allexist=0;
          break;
        fi
        [[ -n $tsv2 ]] && tsv2=$tsv2","
        tsv2=$tsv2""$_tsv;
      done

      [[ $allexist == 0 ]] && continue;

      # compare rmsk
      cmd="
cd $base
~/wzlib/Rutils/bin/bioinfo/edgeR_rmsk.r -U -a $cond1 -b $cond2 -A $tsv1 -B $tsv2 -o rmsk/diff/${cond1}_vs_${cond2}.diff.tsv
"
      jobname="rmsk_diff_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn

    done
}

###################################
# section 4: alternative splicing
###################################

function rnaseq_dexseq {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d DEXSeq ]] || mkdir DEXSeq
  awk '/^\[/{p=0}/\[diffexp\]/{p=1;next} p&&!/^$/' samples |
    while read cond1 cond2 bams1 bams2; do
      # skip non-replicated designs
      # Remark: I haven't figured out how to make DEX-seq
      # work without replica
      [[ -n $(grep -o "," <<< "$bams1") ]] || continue
      [[ -n $(grep -o "," <<< "$bams2") ]] || continue
      outdir=DEXSeq/${cond1}_vs_${cond2}
      cmd="
cd $base
[[ -d $outdir ]] || mkdir $outdir

quant1=\"\"
for b in $(echo $bams1 | sed 's/,/ /g'); do
  bb=\$(basename \$b .bam)
  python ~/.Renv/versions/3.2.3/lib64/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam $WZSEQ_GTF_DEXSEQ \$b $outdir/\$bb.quant.txt
  [[ -z \"\$quant1\" ]] || quant1=\$quant1\",\"
  quant1=\$quant1\"$outdir/\$bb.quant.txt\"
done

quant2=\"\"
for b in $(echo $bams2 | sed 's/,/ /g'); do
  bb=\$(basename \$b .bam)
  python ~/.Renv/versions/3.2.3/lib64/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam $WZSEQ_GTF_DEXSEQ \$b $outdir/\$bb.quant.txt
  [[ -z \"\$quant2\" ]] || quant2=\$quant2\",\"
  quant2=\$quant2\"$outdir/\$bb.quant.txt\"
done

echo \"input quants:\" \$quant1
echo \"input quants:\" \$quant2

~/wzlib/Rutils/bin/bioinfo/DEXSeq.r -G $WZSEQ_GTF_DEXSEQ -a mut -b wt -A \$quant1 -B \$quant2 -o $outdir/${cond1}_vs_${cond2}.tsv

~/wzlib/pyutils/wzseqtk.py ensembl2name --gene -g $WZSEQ_GTF_ENSEMBL_UCSCNAMING -c 2 -i $outdir/${cond1}_vs_${cond2}.tsv -o $outdir/${cond1}_vs_${cond2}.anno.tsv
"
      jobname="dexseq_${cond1}_vs_${cond2}"
      pbsfn=$base/pbs/$jobname.pbs
      pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 24 -memG 10 -ppn 1
      [[ $1 == "do" ]] && qsub $pbsfn
    done
}

# TODO: MATS 

# TODO: MISO (pretty old)

# TODO: HTseq-count

# TODO: RSEM
# ~/tools/rsem/rsem-1.2.22/rsem-prepare-reference
# rsem-calculate-expression -p 20 --calc-ci --ci-memory 12294 --bowtie-chunkmbs 2000 --paired-end --bowtie-path $WZSEQ_BOWTIE1 --rsem-index RSEM

# TODO: ALEXA-Seq http://www.alexaplatform.org/alexa_seq/

# TODO: Trinity for RNA-seq
# de novo assembly

##########################################
# section 5: allele-specific expression
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


###########################
# RSeQC
###########################
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

##########################################
##### thresholding bed files by coverage
##########################################


#################
## SRA-toolkits
#################

function __wzseq_fastq_dump_SE {
  # input: base, sname, srr_ids
  cmd='
set -xe
mkdir -p '$base'/fastq
cd '$base'/fastq
rm -f '${sname}'.fastq.gz
for f in '$srr_ids'; do
  ~/software/sra-toolkit/default/bin/fastq-dump --gzip $f;
  cat ${f}.fastq.gz >>'${sname}'.fastq.gz
  rm -f ${f}.fastq.gz
done
'
  jobname='fastqdump_'$sname
}

function __wzseq_fastq_dump_PE {
  cmd='
set -xe
mkdir -p '$base'/fastq
cd $base/fastq
rm -f '${sname}'_R1.fastq.gz '${sname}'_R2.fastq.gz
for f in '$srr_ids'; do
  ~/software/sra-toolkit/default/bin/fastq-dump --split-files --gzip $f;
  cat ${f}_1.fastq.gz >>'${sname}'_R1.fastq.gz
  [[ -e ${f}_2.fastq.gz ]] && cat ${f}_2.fastq.gz >>'${sname}'_R2.fastq.gz
  rm -f ${f}_1.fastq.gz ${f}_2.fastq.gz
done
'
  jobname="fastqdump_"$sname
}

# download using fastq-dump with accession numbers
# section
# [sra]
# HUES64_derived_CD56_Mesoderm  PE	SRR1067566,SRR1067568,SRR1067569,SRR1067570
function wzseq_fastq_dump() {
  base=$(pwd);
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[sra\]/{p=1;next} p&&!/^$/' samples |
    while read sname pe_or_se srr_ids; do
      srr_ids=${srr_ids//,/ };
      pipeline_init
      hour=48; memG=10; ppn=1
      if [[ $pe_or_se == "PE" ]]; then
        pipeline_eval 1 __wzseq_fastq_dump_PE
      else
        pipeline_eval 1 __wzseq_fastq_dump_SE
      fi
    done
}

# convert downloaded sra to fastq
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

##########
## fastqc
##########

function __wzseq_fastqc() {
  # hour=12; memG=5; ppn=1
  # note that fastq_sname include _R1/R2
  cmd='
set -xe
cd '$base'
mkdir -p fastqc/'$fastq_sname'
fastqc -f fastq '$fastq' -o fastqc/'$fastq_sname'
mkdir -p multiqc/raw/fastqc/
ln -sf `readlink -f fastqc` multiqc/raw/fastqc/
'
  jobname='fastqc_'$fastq_sname
}

####################################
## base/cpg coverage and uniformity
####################################

function __wzseq_uniformity_1M() {
  cmd='
  bedtools makewindows -w 1000000 -g '${WZSEQ_REFERENCE}'.fai | grep -v random | grep -v chrUn | grep -v hap | sortbed | bedtools coverage -a - -b <(samtools view -O BAM -q 40 '$input_bam') -sorted >uniformity/'${sname}'_1Mb.bed
'
  jobname='uniformity_1m_'$sname
}

# create coverage track, unique and nonunique mapping
function wzseq_bam_coverage {

  base=$(pwd);
  [[ -d tracks ]] || mkdir tracks
  [[ -d pbs ]] || mkdir pbs
  [[ -d qc ]] || mkdir qc
  for f in bam/*.bam; do
    fn=$(readlink -f $f)
    __wzseq_bam_coverage
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

function wgbs_cpgcoverage_OBSOLETE {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d cpg ]] || mkdir cpg
  for f in pileup/*.vcf.gz; do
    bfn=$(basename $f .vcf.gz)
    cmd="
cd $base

## coverage
~/tools/biscuit/master/biscuit/biscuit vcf2bed -t cg -k 0 -c $f > cpg/${bfn}_cg.bedg
bedtools intersect -a $WZSEQ_CPGBED -b cpg/${bfn}_cg.bedg -sorted -loj | awk '{if(\$11==\".\")\$11=0;print \$1,\$2,\$3,\$11}' | bedtools groupby -g 1-3 -c 4 -o sum >cpg/cpgCoverage_${bfn}.bedg
wzplot cumhist -t cpg/cpgCoverage_${bfn}.bedg -c 4 --xlabel \"Coverage\" --ylabel \"Cumulative Distribution\" -o cpg/cpgCoverage_wholeGenome_${bfn}.png
bedtools intersect -a cpg/cpgCoverage_${bfn}.bedg -b $WZSEQ_CGIBED | wzplot cumhist -t - -c 4 --xlabel \"Coverage\" --ylabel \"Cumulative Distribution\" -o cpg/cpgCoverage_CGI_${bfn}.png

## beta value distribution
wzplot hist --maxline 1000000000 -t cpg/${bfn}_cg.bedg -c 4 --xlabel \"Beta Values\" -o cpg/betaDist_wholeGenome_${bfn}.png
bedtools intersect -a cpg/${bfn}_cg.bedg -b $WZSEQ_CGIBED | wzplot hist --maxline 1000000000 -t - -c 4 --xlabel \"Beta Values\" -o cpg/betaDist_CGI_${bfn}.png

## Methylation average in CGI
awk '\$5>5' cpg/${bfn}_cg.bedg | bedtools intersect -a $WZSEQ_CGIBED -b - -wo -sorted | bedtools groupby -i - -g 1-5 -c 14,14 -o mean,count > cpg/CGImethAverage_cov5_${bfn}.bed
wzplot hist --maxline 1000000000 -t cpg/CGImethAverage_cov5_${bfn}.bed -c 6 --xlabel \"CGImethAverage\" -o cpg/CGImethAverage_cov5_${bfn}.png

# bedtools closest -a $WZSEQ_TSSBED -b cpg/$bfn.cgi.bed -d >cpg/CGImethAverage_cov5_${bfn}_tss.bed
# rm -f cpg/$bfn.cg.bedg
"
    jobname="cpgCoverage_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 12 -memG 10 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

##########################
# reverse bam to fastq
##########################
#
# bam2fastq
#
# wzseq_bam2fastq
# unmark secondary (some bam have only secondary mapping but primary mapping missing)
# samtools view -h in.bam | awk '!/^@/{$2=and($2,compl(0x100)); print $0}/^@/' | samtools view -bo tmp0.bam

function __wzseq_bam2fastq {

  cmd='
set -xe
# group reads by read names, collate is faster than sort
cd '$base';
mkdir -p bam/collate;
i=1;

rm -f fastq/'$sname'.pe1.fq.gz
rm -f fastq/'$sname'.pe2.fq.gz
rm -f fastq/'$sname'.se.fq.gz
rm -f fastq/'$sname'.paired_nolabel.fq.gz
for sourcebam in '${sourcebams//,/ }'; do

  ## NOTE: I extracted only primary mapping
  samtools view -h $sourcebam | awk -F"\t" -v OFS="\t" '\''!/^@/{$2=and($2,compl(0x100)); print $0}/^@/'\'' > bam/collate/'$sname'_${i}_tmp1.sam
  samtools collate -u bam/collate/'$sname'_${i}_tmp1.sam bam/collate/'$sname'_${i}_tmp2
  ## NOTE: use the following to skip primary mapping extraction
  ## samtools collate -u $sourcebam collate/${sample}_\$i.bam;
  
  ## NOTE: sometimes I used -O too. But -O can cause seg-fault on some malformed bams
  samtools fastq -n -0 fastq/'$sname'_$i.paired_nolabel.fq -1 fastq/'$sname'_$i.pe1.fq -2 fastq/'$sname'_$i.pe2.fq -s fastq/'$sname'_$i.se.fq bam/collate/'$sname'_${i}_tmp2.bam;
  gzip -c fastq/'$sname'_$i.paired_nolabel.fq >>fastq/'$sname'.paired_nolabel.fq.gz
  gzip -c fastq/'$sname'_$i.pe1.fq >>fastq/'$sname'.pe1.fq.gz
  gzip -c fastq/'$sname'_$i.pe2.fq >>fastq/'$sname'.pe2.fq.gz
  gzip -c fastq/'$sname'_$i.se.fq >>fastq/'$sname'.se.fq.gz 
  rm -f fastq/'$sname'_$i.paired_nolabel.fq fastq/'$sname'_$i.pe1.fq fastq/'$sname'_$i.pe2.fq fastq/'$sname'_$i.se.fq
  rm -f bam/collate/'$sname'_$i_tmp{1,2}*;
  i=$((i+1))
done
'
  jobname='bam2fastq_'$sname
}

# remove suffix _1, _2, mark read info to flag and
#  samtools view -h bam/SRR1029055.bam chr19 | awk '!/^@/{inpair=substr($1,length($1),1);$1=substr($1,1,length($1)-2);if(inpair==1) {$2=or($2,0x40);} else {$2=or($2,0x80);} $2=or($2,0x1); if (!(and($2, 0x100))) print $0}/^@/' | samtools collate -uO - SRR1029055tmp >SRR1029055.collate.bam
#  samtools view -h SRR1029055.collate.bam | awk 'BEGIN{key=""; line=""}!/^@/{if (key==$1) {print line; print $0;} key=$1; line=$0;}/^@/' | samtools fastq - -1 fastq_chr19/read1.fastq -2 fastq_chr19/read2.fastq -0 fastq_chr19/unpaired.fastq
function wzseq_bam2fastq {
  base=$(pwd);
  [[ -d bam ]] || mkdir bam
  [[ -d pbs ]] || mkdir pbs
  awk '/^\[/{p=0}/\[bam2fastq\]/{p=1;next} p&&!/^$/' samples |
    while read sample sourcebams; do
      cmd="
set -xe
# group reads by read names, collate is faster than sort
cd $base;
mkdir -p bam/collate;
i=1;
for sourcebam in ${sourcebams//,/ }; do
  samtools view -h \$sourcebam | awk -F\"\\t\" -v OFS=\"\t\" '!/^@/{\$2=and(\$2,compl(0x100)); print \$0}/^@/' > bam/collate/${sample}_\${i}_tmp1.sam
  samtools collate -u bam/collate/${sample}_\${i}_tmp1.sam bam/collate/${sample}_\${i}_tmp2
  # samtools collate -u \$sourcebam collate/${sample}_\$i.bam;
  samtools fastq -On -0 fastq/${sample}_\$i.paired_nolabel.fq -1 fastq/${sample}_\$i.pe1.fq -2 fastq/${sample}_\$i.pe2.fq -s fastq/${sample}_\$i.se.fq bam/collate/${sample}_\${i}_tmp2.bam;
  gzip -c fastq/${sample}_\$i.paired_nolabel.fq >>fastq/${sample}.paired_nolabel.fq.gz
  gzip -c fastq/${sample}_\$i.pe1.fq >>fastq/${sample}.pe1.fq.gz
  gzip -c fastq/${sample}_\$i.pe2.fq >>fastq/${sample}.pe2.fq.gz
  gzip -c fastq/${sample}_\$i.se.fq >>fastq/${sample}.se.fq.gz 
  rm -f fastq/${sample}_\$i.paired_nolabel.fq fastq/${sample}_\$i.pe1.fq fastq/${sample}_\$i.pe2.fq fastq/${sample}_\$i.se.fq
  rm -f bam/collate/${sample}_\${i}_tmp{1,2}*;
  i=\$((i+1))
done
"
## the following cleaning creates piping error, better check the size of output fastq manually
# [[ -z \$(gunzip -c fastq/${sample}.paired_nolabel.fq.gz | head -c1) ]] && rm -f fastq/${sample}.paired_nolabel.fq.gz
# [[ -z \$(gunzip -c fastq/${sample}.pe1.fq.gz | head -c1) ]] && rm -f fastq/${sample}.pe1.fq.gz
# [[ -z \$(gunzip -c fastq/${sample}.pe2.fq.gz | head -c1) ]] && rm -f fastq/${sample}.pe2.fq.gz
# [[ -z \$(gunzip -c fastq/${sample}.se.fq.gz | head -c1) ]] && rm -f fastq/${sample}.se.fq.gz
    jobname="bam2fastq_$sample"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 48 -memG 5 -ppn 1
    [[ ${!#} == "do" ]] && qsub $pbsfn
  done
}

#####################
## Qualimap
#####################

function __wzseq_qualimap_bamqc() {
  # NOTE: be sure to set memG correctly, otherwise, it would be very very slow and weird
  # hour=24; memG=200; ppn=28; queue=longq
  # input_bam=
  # output_sname=
  cmd='
cd '$base'
qualimap --java-mem-size='$memG'G bamqc -nt '$ppn' -bam '$input_bam' -outdir qualimap/'$output_sname' -c
mkdir -p multiqc/raw/qualimap
ln -sf `readlink -f qualimap/'$sname'` multiqc/raw/qualimap/
'
  jobname='qualimap_bamqc_'$sname
}

function __wzseq_qualimap_rnaseqSE() {
  # hour=24; memG=50; ppn=10
  cmd='
cd '$base'
qualimap --java-mem-size=10G rnaseq -bam '$input_bam' -gtf '$WZSEQ_GTF' -outdir qualimap/rnaseq_'$sname'
'
  jobname='qualimap_rnaseqSE_'$sname
}

function __wzseq_qualimap_rnaseqPEstranded() {
  # hour=24; memG=50; ppn=10
  cmd='
cd '$base'
qualimap --java-mem-size=10G rnaseq -bam '$input_bam' -gtf '$WZSEQ_GTF' -outdir qualimap/rnaseq_'$sname' --paired -p strand-specific-forward
'
  jobname='qualimap_rnaseqPEst_'$sname
}

############
## Picard
############

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

function wzseq_picard_markdup {

  base=$(pwd)
  [[ -d pbs ]] || mkdir pbs
  [[ -d bam/before_mdup ]] || mkdir -p bam/before_mdup
  [[ -d tmp ]] || mkdir tmp
  awk '/^\[/{p=0}/\[alignment\]/{p=1;next} p&&!/^$/' samples |
    while read bfn sread1 sread2; do
      f=bam/$bfn.bam
      cmd="
cd $base
if [[ ! -e bam/before_mdup/$bfn.bam ]]; then
  mv $f bam/before_mdup/$bfn.bam
  [[ -e $f.bai ]] && mv $f.bai bam/before_mdup/$bfn.bam.bai
  [[ -e $f.flagstat ]] && mv $f.flagstat bam/before_mdup/$bfn.bam.flagstat
fi
[[ -d bam/picard_mdup ]] || mkdir bam/picard_mdup
java -Xmx10g -Djava.io.tmpdir=./tmp/ -jar /home/wanding.zhou/software/picard/picard-tools-1.138/picard.jar MarkDuplicates CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=bam/picard_mdup/$bfn.mdup.stats READ_NAME_REGEX=null INPUT=bam/before_mdup/$bfn.bam OUTPUT=bam/picard_mdup/$bfn.bam TMP_DIR=tmp
samtools flagstat bam/picard_mdup/$bfn.bam >bam/picard_mdup/$bfn.bam.flagstat

cd bam; 
[[ -h $bfn.bam ]] || [[ ! -e $bfn.bam ]] && ln -sf picard_mdup/$bfn.bam .
[[ -h $bfn.bam.bai ]] || [[ ! -e $bfn.bam.bai ]] && ln -sf picard_mdup/$bfn.bai $bfn.bam.bai
[[ -h $bfn.bam.flagstat ]] || [[ ! -e $bfn.bam.flagstat ]] && ln -sf picard_mdup/$bfn.bam.flagstat .

" # other options: REMOVE_DUPLICATES=true
    jobname="picard_markdup_$bfn"
    pbsfn=$base/pbs/$jobname.pbs
    pbsgen one "$cmd" -name $jobname -dest $pbsfn -hour 60 -memG 50 -ppn 5
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

################
## trimmomatic
################
# trimming of adaptor and quality
# see http://www.usadellab.org/cms/?page=trimmomatic for explanation
function __wzseq_trimmomatic_SE {
  # hour=12; memG=10; ppn=10
  cmd='
mkdir trimmomatic
java -jar /primary/home/wanding.zhou/tools/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads '$ppn' -phred33 '$input_fastq' '$output_fastq' ILLUMINACLIP:/primary/home/wanding.zhou/tools/trimmomatic/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>'$input_fastq'_trimmomatic_report
mkdir -p multiqc/raw/trimmomatic/
ln -sf `readlink -f '$input_fastq'_trimmomatic_report` multiqc/raw/trimmomatic/
'
  jobname='trimmomatic_SE_'$sname
}

function __wzseq_trimmomatic_PE {
  # hour=12; memG=10; ppn=10
  cmd='
mkdir -p multiqc/raw/trimmomatic/
java -jar /primary/home/wanding.zhou/tools/trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads '$ppn' -phred33 '$input_fastq1' '$input_fastq2' '$output_fastq1' '$output_fastq1'_unpaired.fastq.gz '$output_fastq2' '$output_fastq2'_unpaired.fastq.gz ILLUMINACLIP:/primary/home/wanding.zhou/tools/trimmomatic/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>multiqc/raw/trimmomatic/'$sname'_trimming_report
'
  jobname='trimmomatic_PE_'$sname
}

#########
## others
#########

## assume chrm, beg, end, beta, coverage
function wzseq_cov5 {
  f=$1
  zcat $f | awk '$5>=5' | gzip -c >${f%.bed.gz}.cov5.bed.gz
  echo `zcat ${f%.bed.gz}.cov5.bed.gz | wc -l` "CpGs covered 5X"
}

## assume chrm, beg, end, beta, coverage
function wzseq_cov10 {
  f=$1
  zcat $f | awk '$5>=10' | gzip -c >${f%.bed.gz}.cov10.bed.gz
  echo `zcat ${f%.bed.gz}.cov10.bed.gz | wc -l` "CpGs covered 10X"
}

function wzseq_liftbw {
  # liftOver bigwig file
  # Usage: wzseq_liftbw input.bigWig ~/tools/liftover/mm9ToMm10.over.chain.gz output.bigWig ~/references/mm10/mm10.fa.fai
  input=$1
  chain=$2
  output=$3
  chromsize=$4
  echo "[$(date)] Converting bigwig to bedgraph.."
  bigWigToBedGraph $input $input.bedg.tmp
  echo "[$(date)] Lifting over.."
  liftOver $input.bedg.tmp $chain $output.bedg.tmp $output.tmp.unmapped
  echo "  Mapped:   $(wc -l $output.bedg.tmp)"
  echo "  Unmapped: $(wc -l $output.tmp.unmapped)"
  echo "[$(date)] Sorting bedGraph and skip overlapping ..."
  sortbed $output.bedg.tmp | wzbedtools deoverlap -i - -o $output.bedg.tmp.sorted
  echo "  Before skipping: "$(wc -l $output.bedg.tmp)
  echo "  After skipping:  "$(wc -l $output.bedg.tmp.sorted)
  echo "[$(date)] Converting bedGraph to bigWig .."
  bedGraphToBigWig $output.bedg.tmp.sorted $chromsize $output
  echo "[$(date)] Cleaning"
  rm -f $input.bedg.tmp $output.tmp.unmapped $output.bedg.tmp $output.bedg.tmp.sorted
  echo "[$(date)] Done."
}

# fetch all SRR using SRX 
# wzseq_srx2srr SRX306253
function wzseq_srx2srr {
  curl "https://www.ncbi.nlm.nih.gov/sra?term=$1" | grep -o "SRR[[:digit:]]*" | sort | uniq | paste -d, -s
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

function __wzseq_index_bam() {
  : '
bam=bam/${sname}.bam
hour=5; memG=5; ppn=1
pipeline_eval x __wzseq_index_bam
'
  cmd='
samtools index '$bam';
samtools flagstat '$bam' > '$bam'.flagstat
mkdir -p multiqc/raw/flagstats/
ln -sf `readlink -f '$bam'.flagstat` multiqc/raw/flagstats/
'
  jobname="bamindex_"${bam//\//_}
}

# basic summary of the folder
function wzseq_basic_summary() {
  # for RNAseq
  for f in bam/*.bam; do b=${f%.bam}; b=${b#bam/}; sec=$(cat ${f%.bam}/accepted_hits.bam.flagstat | grep 'secondary' | awk '{match($1, /([0-9]*) \+/, a); print a[1];}'); tot=$(cat ${f%.bam}/accepted_hits.bam.flagstat | grep 'total (' | awk '{match($1, /([0-9]*) \+/, a); print a[1];}'); echo ${b%.bam} $(($tot-$sec)); done
  echo -e "\nraw read counts (paired-end, single-end)"
  for f in fastq/*.fastq.gz; do c=$(zcat $f | lc); echo $f $(($c / 2)) $(($c / 4)); done
}

function __wzseq_bam_mapq() {
  cmd='
cd '$base'
mkdir -p multiqc/raw/mapq
samtools view -F 0x100 -f 0x4 '$input_bam' | wc -l | cat <(echo -ne "unmapped\t") - >'$input_bam'_mapq_table
samtools view -F 0x104 '$input_bam' | awk '\''{cnt[$5]+=1}END{for(mapq in cnt) {print mapq"\t"cnt[mapq];}}'\'' | sort -k1,1n >>'$input_bam'_mapq_table
mkdir -p multiqc/raw/mapq/
ln -fs `readlink -f '$input_bam'_mapq_table` multiqc/raw/mapq/
'
  jobname='mapqdist_'$sname
}

