#!/usr/bin/env sh

[ -z "$1" ] && echo "prog longname parentname rootpath suffix shortname gsm";

longname=$1
parentname=$2
root=$3
suffix=$4
shortname=$5
gsm=$6
template=$(cat <<EOF
         track $longname
         parent $parentname on
         type bigWig
         color 0,102,255
         maxHeightPixels 128:25:10
         viewLimits 0.0:1.0
         shortLabel $shortname
         longLabel  $longname
         bigDataUrl $root/$gsm.$suffix
         priority 1
EOF
);

echo "$template"
