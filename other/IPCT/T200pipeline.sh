function T200pipeline {

  if [[ $# -ne 4 ]]; then
      echo "T200pipeline <samplename> <fastq_dir> <bam_dir> <reference>"
      echo "need PICARD_MARKDUP set up"
      return;
  fi
  sample=$1
  fastq_dir=$2
  bam_dir=$3
  reference=$4

  fastq1_fns=($(/bin/ls $fastq_dir/*R1*));
  fastq2_fns=($(/bin/ls $fastq_dir/*R2*));
  if (( ${#fastq1_fns} != ${#fastq1_fns} )); then
    echo "Paired fastq names unmatched. Abort."
    exit;
  fi
  smallbams=()
  for ((i=0; i<${#fastq1_fns[@]}; ++i)); do
    fq1=${fastq1_fns[$i]}
    fq2=${fastq2_fns[$i]}
    fq1sai=$bam_dir/$(basename $fq1).sai
    fq2sai=$bam_dir/$(basename $fq2).sai
    bwa aln -t 4 $reference $fq1 > $fq1sai
    bwa aln -t 4 $reference $fq2 > $fq2sai
    rgprefix=${sample}_$(($i+1))
    bamprefix=$bam_dir/rgprefix
    bwa sampe $reference $fq1sai $fq2sai $fq1 $fq2 | samtools view -b -S - | samtools sort - $bamprefix
    if [[ i -eq 0 ]]; then samtools view -H $bamprefix.bam > $bam_dir/header.txt; fi
    echo -e "@RG\tID:$rgprefix\tLB:$sample\tSM:$sample" >> $bam_dir/header.txt;
    smallbams+=($bamprefix.bam)
  done

  sortedbam=$bam_dir/$sample.sorted.bam
  samtools merge -r -h $bam_dir/header.txt -f $sortedbam ${smallbams[@]}
  samtools index $sortedbam
  samtools flagstat $sortedbam $sortedbam.flagstat

  rmdupbam=$bam_dir/$sample.rmdup.bam
  java -jar PICARD_MARKDUP I=$sortedbam O=$rmdupbam M=$bam_dir/$sample.markdup.metrics AS=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
  samtools index $rmdupbam
  
}
