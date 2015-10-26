function wzcore_bam_count_total {
  awk 'BEGIN{sum=0}{sum+=$3-$2}END{print sum}' $1
}
