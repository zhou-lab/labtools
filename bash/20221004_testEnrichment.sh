#!/usr/bin/bash
# source ~/repo/labtools/bash/20221004_testEnrichment.sh
# testEnrichment ~/references/mm10/annotation/cpg/cpg_nocontig.bed.gz /mnt/isilon/zhou_lab/projects/20191221_references/mm10/features/EnsRegBuild.20220710.bed.gz /scr1/users/zhouw3/projects/20220609_ExpressionMethylationCorrelation/20220815_Clark/Output/EBcells/AllTestedCPG.bed /scr1/users/zhouw3/projects/20220609_ExpressionMethylationCorrelation/20220815_Clark/Output/EBcells/NegSig.5.CPGonly.bed /scr1/users/zhouw3/tmp/FeatureAggregationtest2

function set_environment {
    Reference=$1
    Feature=$2
    Universe=$3
    Query=$4
    Output=$5
}

function testEnrichment() (     # this spawn a subshell
    set_environment $1 $2 $3 $4 $5
    echo
    echo "================="
    echo "Query:     $Query"
    echo "Reference: $Reference"
    echo "Universe:  $Universe"
    echo "Feature:   $FeatureFILE"
    echo "Output:    $Output"
    echo "================="

    base=$(basename ${Query} .bed)
    TMPFDR=${Output}.tmp
    mkdir -p $TMPFDR
    rm -rf $TMPFDR/*

    # assume universe is sorted
    bedtools intersect -a ${Reference} -b ${Universe} -sorted -wo | cut -f1-3 | uniq >$TMPFDR/in_u
    # make query cpg-based if not yet
    zcat -f ${Query} | cut -f1-3 | sortbed >$TMPFDR/in_q
    # make feature cpg-based if not yet
    zcat -f ${Feature} | sortbed >$TMPFDR/in_f
    bedtools intersect -a $TMPFDR/in_u -b $TMPFDR/in_f -loj -sorted | bedtools intersect -a - -b $TMPFDR/in_q -loj -sorted |
	awk '{if($9==".") {$9=0;} else {$9=1;}print $1,$2,$3,$7,$9;}' | sort -k4,4 -k1,1 -k2,2n |
	bedtools groupby -g 1-4 -c 5 -o sum | awk '{if($5>0) $5=1; print;}' >$TMPFDR/overlaps
    n_q=$(bedtools intersect -a $TMPFDR/in_u -b $TMPFDR/in_q -sorted | wc -l)
    n_u=$(cat $TMPFDR/in_u | wc -l)
    awk -v n_q=$n_q -v n_u=$n_u '$4!="."{k=$1":"$2"_"$3; if(k!=k0 || g!=$4) {if($5==0) {cnt0[$4]+=1;} if($5==1) {cnt1[$4]+=1;} features[$4]=1;} k0=k; g=$4;}END{print "Feature\tnfmq\tnfq\tnq\tnu"; for(i in features){print i,cnt0[i],cnt1[i],n_q,n_u;}}' $TMPFDR/overlaps >$TMPFDR/stat_cnts
    Rscript ~/repo/labtools/Rutils/testFisher.R $TMPFDR/stat_cnts ${Output}

    echo "Created temporary folder $TMPFDR. Feel free to delete."
    echo "All completed."
    echo
)
