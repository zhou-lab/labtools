#!/bin/bash

export KYCG_species=mm10
export KYCG_basedir=/scr1/users/zhouw3/projects/2021124_scMeth/Annotation/

# should list the excel spreadsheet
function KYCG_listDBGroups {
  ls $KYCG_basedir/$KYCG_species
}

function KYCG_dbStats {

  local options OPTARG OPTIND
  while getopts ":m:d:" options; do
    case "${options}" in
      m) methcall_file=${OPTARG} ;;
      d) db_name=${OPTARG} ;;
    esac
  done
  shift $((OPTIND-1))

  mkdir -p tmp
  sname=$(basename $methcall_file .bed.gz)
  db_name=${db_name%.bed.gz}
  db_file=$KYCG_basedir/$KYCG_species/${db_name}.bed.gz

  sname=R106W_BSseq
  db_name=ChromHMM.20220318
  db_file=/scr1/users/zhouw3/projects/2021124_scMeth/Annotation/mm10/ChromHMM.20220318.bed.gz
  methcall_file=/mnt/isilon/zhoulab/labprojects/20220312_JoeZhou/methcall/R106W_BSseq-R1.bed.gz
  
  zcat $db_file | bedtools intersect -a - -b $methcall_file -sorted -wo | sort -k4,4 | bedtools groupby -g 4 -c 9,9 -o mean,count | awk -v SNAME=$sname '{print $0"\t"SNAME}' >>tmp/${sname}_${db_name}.txt
}

function KYCGutil_computeMean {
  # first column is the key, 2nd column is the value to compute mean on
  awk 'BEGIN{OFS="\t"}{cnt[$1]+=1; sum[$1]+=$2;}END{for(k in cnt){print k,sum[k]/cnt[k],cnt[k];}}' $1
}

# function KYCG_dbStats {
#   rm -f tmp/tmp1_*
#   parallel -j 24 'KYCG_dbStats1 -m {} CpGisland.20220306.bed.gz ChromHMM.20220318.bed.gz' ::: $@
#   cat tmp/*
# }

# ## KYCG_dbStats meth_call/*.bed.gz
function KYCG_mapToMM285 {
  methcall_file=$1
  zcat /mnt/isilon/zhoulab/labprojects/20200228_Mouse_Array_Project/20201202_multispecies_alignment/MM285/mapping/mm10.tsv.gz | awk 'NR>1&&$1!="NA"{print $1,$2,$3,$9;}' | sortbed | bedtools intersect -a - -b $methcall_file -sorted -wo | cut -f4,8
}
