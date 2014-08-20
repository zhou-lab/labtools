function wz_ipct_sampleloc03 {
  for project in /projects*; do
    for flowcell in $project/flowcell*; do
      if [[ -d $flowcell/data ]]; then
	samples=($(/bin/ls $flowcell/data));
	for sample in ${samples[@]}; do
	  if [[ -d $flowcell/data/$sample/gz ]]; then
	    echo -e $(basename $flowcell)"\t"$sample"\tdqs03\t"$flowcell/data/$sample/gz;
	  fi
	done
      fi
    done
  done
}

function wz_ipct_sampleloc06 {
  for flowcell in /projects/flowcell*; do
    if [[ -d $flowcell/data ]]; then
      samples=($(/bin/ls $flowcell/data));
      for sample in ${samples[@]}; do
	if [[ -d $flowcell/data/$sample/gz ]]; then
	  echo -e $(basename $flowcell)"\t"$sample"\tdqs06\t"$flowcell/data/$sample/gz;
	fi
      done
    fi
  done
}

# Ex: wz_ipct_sample_match sampleloc
# sampleloc is generated by wz_ipct_sampleloc0{3,6}
function wz_ipct_sample_match {
  while read flowcell sample server loc; do
    if [[ $sample =~ Tumor-([0-9]*) ]]; then
      sample_id=${BASH_REMATCH[1]};
      grep Normal-$sample_id[^0-9] $1 | while read normal; do
	if [[ -n $normal ]]; then
	  echo -e $sample_id"\t"$sample"\t"$loc"\t$normal" | tr -s [:blank:] \\t;
	fi
      done
    fi
  done < $1;
}
