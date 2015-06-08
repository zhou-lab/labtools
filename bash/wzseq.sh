#!/bin/bash

function auto_setup_links() {
  # usage: auto_setup_links /data/sequencing/analysis/H55TVBGXX/

  rootdir=$1;
  flowcell=$(basename $rootdir); #H55TVBGXX
  samples=($(/bin/ls $rootdir/results/$flowcell)); # /data/sequencing/analysis/H5MW5BGXX/results/H5MW5BGXX/H5MW5BGXX_1_PL430BS1
  
  mkdir data;
  ln -s $rootdir data/root;

  bams=($(find $rootdir -name *.fa.realign.mdups.recal.bam));
  if [[ $bams ]]; then
    for bam in ${bams[@]}; do
      if [[ $bam =~ ResultCount_([^.]).([^.]).fa.realign.mdups.recal.bam ]]; then
        ln -s $bam data/${BASH_REMATCH[2]}.bam
      fi
  fi
}

