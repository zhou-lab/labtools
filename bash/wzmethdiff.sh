#!/bin/bash
# usage: biscuit_diffmeth -t pileup1 -n pileup2 [-b base] [-c mincov] -v $WZSEQ_REFVERSION -g $WZSEQ_CGIBED -r $WZSEQ_REFERENCE

function decho() {
  echo "[$(date)] "$@ >&2
}

while getopts "b:t:n:c:v:g:r:" opt; do
  case $opt in
    t) pileup1=$(readlink -f $OPTARG) ;;
    n) pileup2=$(readlink -f $OPTARG) ;;
    b) base=$(readlink -f $OPTARG) ;;
    c) mincov=$OPTARG ;;
    v) wzseq_refversion=$OPTARG ;;
    g) wzseq_cgibed=$OPTARG ;;
    r) wzseq_reference=$OPTARG ;;
    \?) echo "Invalid option: -$OPTARG" >&2; return 1;;
    :) echo "Option -$OPTARG requires an argument." >&2; return 1;;
  esac
done

max_hyper_len=500
max_hypo_len=2000
hyperwindow_minhypercnt=7
hyperwindow_maxhypocnt=2
hypowindow_minhypocnt=7
hypowindow_maxhypercnt=2

[[ -z ${base+x} ]] && base=$(pwd)/diffmeth
[[ -z ${mincov+x} ]] && mincov=3

[[ -d $base ]] || mkdir -p $base
contrast=$base/$(basename $pileup1 .pileup.gz)"_vs_"$(basename $pileup2 .pileup.gz)

######### single base analysis ############

decho "Generating "$contrast.diff
# bedtools intersect -a <(zcat $pileup1 | awk -v mincov=$mincov '$8!="."&&$9!="."&&$8+$9>=mincov{print $1,$2,$3,$4,$6,$7,$8,$9}') -b <(zcat $pileup2 | awk -v mincov=$mincov '$8!="."&&$9!="."&&$8+$9>=mincov{print $1,$2,$3,$4,$6,$7,$8,$9}') -sorted -wo | awk -f wanding.awk -e '$5~/[ATCG]CG/{bt=$7/($7+$8); bn=$15/($15+$16); print joinr(1,16),bt,bn,bt-bn}' > $contrast.diff
bedtools intersect -a <(biscuit vcf2bed -ct cg -k 5 $pileup1) -b <(biscuit vcf2bed -ct cg -k 5 $pileup2) -sorted -wo | awk '{print $1,$2,$3,$4-$10,$5,$4,$11,$10}' >$contrast.diff
decho "Done"

decho "Plotting delta beta distribution"
decho "Generating "$contrast.diff.dist.png
cut -f4 $contrast.diff | wzplot hist -c 1 -o $contrast.diff.dist.png --xlabel "delta_beta" --ylabel "#CpG"
decho "Done"

decho "Generating "$contrast.bw
cut -f1,2,3,4 $contrast.diff > $contrast.bedGraph
bedGraphToBigWig $contrast.bedGraph ${wzseq_reference}.fai $contrast.bw
rm -f $contrast.bedGraph
decho "Done"

decho "Whole CpG island analysis (chrm,beg,end,#cpg,median_beta1,median_beta2,delta_betas1-2,t-stat,p-value"
bedtools intersect -a $wzseq_cgibed -b $contrast.diff -sorted -wo | wzstats.py ttest -c1 16 -c2 18 --groupby 1-3 - --rmNA | awk -f wanding.awk -e '{if($9<0.05&&abs($7)>0.3){if($7>0){a="hyper"}else{a="hypo"}}else{if($4>7){a="nc"}else{a="na"}} print $0,a}' >$contrast.cpgi.diff

########### tiling window analysis #########

decho "Generating "$contrast.diff.window
perl -alne 'BEGIN{my @window; my @poses; my @chrs;}{if (scalar @window == 10) {$p1 = scalar grep {$_>0.3} @window; $n1 = scalar grep {$_<-0.3} @window; print $chrs[0]."\t".$poses[0]."\t".$poses[9]."\t".$p1."\t".$n1; shift(@window); shift(@poses); shift(@chrs);} if (scalar @window>0 && $chrs[scalar @chrs-1] ne $F[0]) {@window=();@chrs=();@poses=();} push(@window,$F[3]); push(@chrs, $F[0]); push(@poses, $F[1]);}' $contrast.diff > $contrast.diff.window
decho "Done"

windowsizeplot=$contrast"_diffmeth_windowsize_dist.png"
decho "Plotting differential methylation window size distribution"
decho "Generating "$windowsizeplot
awk '{print $3-$2}' $contrast.diff.window | wzplot hist -c 1 -o $windowsizeplot --xlabel "window size (bp)" --ylabel "count"
decho "Done"

decho "Plotting hyper-methylation window size distribution."
awk -v p=$minhypercnt '$4>p{print $3-$2}' $contrast.diff.window | wzplot hist -o $contrast"_diff_hyper_window_size_dist.png" --xlabel "window length" --ylabel "#windows" --xlog

decho "Plotting hypo-methylation window size distribution."
awk -v p=$minhypocnt '$5>p{print $3-$2}' $contrast.diff.window | wzplot hist -o $contrast"_diff_hypo_window_size_dist.png" --xlabel "window length" --ylabel "#windows" --xlog

decho "Combining hyper-methylation segment."
awk -v m=$max_hyper_len -v p=$hyperwindow_minhypercnt -v q=$hyperwindow_maxhypocnt '(($3-$2)<m) && ($4>=p) && $5<=q' $contrast.diff.window | bedtools merge -i - -c 4 -o count > $contrast.hyper.bed
decho "Combining hypo-methylation segment."
awk -v m=$max_hypo_len -v p=$hypowindow_minhypocnt -v q=$hypowindow_maxhypercnt '(($3-$2)<m) && ($5>=p) && $4<=q' $contrast.diff.window | bedtools merge -i - -c 4 -o count > $contrast.hypo.bed

decho "Intersecting CpG Island."
bedtools intersect -a $contrast.hyper.bed -b $wzseq_cgibed > $contrast.hyper.cgi.bed
bedtools intersect -a $contrast.hypo.bed -b $wzseq_cgibed > $contrast.hypo.cgi.bed
decho "There are "$(wc -l $contrast.hyper.cgi.bed | cut -d" " -f1)" hyper-methylated CGIs and "$(wc -l $contrast.hypo.cgi.bed | cut -d" " -f1)" hypo-methylated CGIs."

decho "TransVar annotation."
awk '{print $1":"$2"_"$3}' $contrast.hyper.bed | transvar ganno -l - --refversion $wzseq_refversion --ccds >$contrast.hyper.transvar
awk '{print $1":"$2"_"$3}' $contrast.hypo.bed | transvar ganno -l - --refversion $wzseq_refversion --ccds >$contrast.hypo.transvar

decho "All done"
