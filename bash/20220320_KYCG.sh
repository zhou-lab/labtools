
export KYCG_SPECIES=mm10

# should list the excel spreadsheet
function KYCG_listDBGroups {
  base_dir=/scr1/users/zhouw3/projects/2021124_scMeth/Annotation/
  ls $base_dir/$KYCG_SPECIES
}

function KYCG_dbStats {
  db_name=$1; shift;
  db_file=$db_basedir/$KYCG_SPECIES
  methcall_files="$@"
  mkdir -p tmp
  for methcall_file1 in $methcall_files; do
    echo $methcall_file1 $db_file
    sname=$(basename $methcall_file1 .bed.gz)
    bedtools intersect -a $db_file -b $methcall_file1 -sorted -wo | sort -k4,4 | bedtools group -g 4 -c 9 -o mean >tmp/$sname
  done
}
