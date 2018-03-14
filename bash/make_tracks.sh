#!/usr/bin/env sh

[ -z "$1" ] && echo "prog samplename parentname rootpath suffix";

samplename=$1
parentname=$2
root=$3
suffix=$4
template=$(cat <<EOF
         track $samplename
         parent $parentname on
         type bigWig
         color 0,102,255
         maxHeightPixels 128:25:10
         viewLimits 0.0:1.0
         shortLabel $samplename
         longLabel  $samplename
         bigDataUrl $root/$samplename.$suffix
         priority 1
EOF
);

echo "$template"
