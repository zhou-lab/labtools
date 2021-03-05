#!/bin/bash
# ------------------------------------------------------------------
# [Ethan Jacob Moyer] tbk2track.sh
# 	Convert file /tbk file directory processed through tbmate into 
# 	bed file and bigwig file directories.   
# ------------------------------------------------------------------

VERSION=0.0.1
SUBJECT=some-unique-id
usage="$(basename "$0") [-h] [-d tbk] [-g hg19/hg38] [-t bw/bb] -- This function converts tbk files stored in given /tbk directory into either bigwig or bigbed formats for trackhub development. Currently, it can convert the tbk files into either the hg19 or hg38 reference. The /bed, /bw, /bb directories need not be created before hand. 

where:
    -h  show this help text
    -d  sets the tbk directory to be converted to bed and bigwig; use tbk and not tbk/
    -g  sets the genomic reference, either hg19 or hg38; default value is hg19
    -t  sets the type of trackhub format, either bigwig [bw] or bigbed [bw]; default value is bw
    -v  sets whether verbose is true or false; default is false"


# --- Options processing -------------------------------------------
if [[ $1 == '' ]] ; then
    echo "$usage"
    exit 1;
fi

verbose="false"
genomic_ref="hg19"
type="bw"
while getopts ':h:d:g:t:v:' option; do
  case "$option" in
    "h") echo "$usage"
       exit
       ;;
    "d") tbk_path=$OPTARG
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    "g") genomic_ref=$OPTARG
       ;;
    "t") type=$OPTARG
       ;;
    "v") verbose=$OPTARG
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

if [[ $genomic_ref == 'hg19' ]]; then
	cpg=/mnt/isilon/zhou_lab/projects/20191221_references/hg19/annotation/cpg/cpg.bed.gz

	sample_file=$(ls $tbk_path | grep ".tbk" | head -n+1)


	if tbmate header $tbk_path"/"$sample_file | grep -q "HM450.idx.gz"; then
		echo "HM450 to hg19 is currently not supported."
	elif tbmate header $tbk_path"/"$sample_file | grep -q "EPIC.idx.gz"; then
		to_cpg=/mnt/isilon/zhou_lab/projects/20191221_references/InfiniumArray/EPIC/EPIC_to_hg19.idx.gz
	elif tbmate header $tbk_path"/"$sample_file | grep -q "/mnt/isilon/zhou_lab/projects/20191221_references/hg19/annotation/cpg/idx.gz"; then
		to_cpg=/mnt/isilon/zhou_lab/projects/20191221_references/hg19/annotation/cpg/idx.gz
	elif tbmate header $tbk_path"/"$sample_file | grep -q "/mnt/isilon/zhou_lab/projects/20191221_references/hg38/annotation/cpg/idx.gz"; then
		to_cpg=/mnt/isilon/zhou_lab/projects/20191221_references/hg38/annotation/cpg/hg38_to_hg19_idx.gz
	else
		printf "There is not enough information to determine the index file via tbmate header %s\n" "$genomic_ref"
		exit 1
	fi
	
	chr=/mnt/isilon/zhou_lab/projects/20191221_references/hg19/hg19.fa.fai
elif [[ $genomic_ref == 'hg38' ]]; then
	cpg=/mnt/isilon/zhou_lab/projects/20191221_references/hg38/annotation/cpg/cpg.bed.gz
	sample_file=$(ls | grep ".tbk" | head -n+1)
	if tbmate header $tbk_path$sample_file | grep "HM450.idx.gz"; then
		echo "HM450 to hg38 is currently not supported."
	elif tbmate header $tbk_path$sample_file | grep "EPIC.idx.gz"; then
		to_cpg=/mnt/isilon/zhou_lab/projects/20191221_references/InfiniumArray/EPIC/EPIC_to_hg38.idx.gz
	elif tbmate header $tbk_path"/"$sample_file | grep -q "/mnt/isilon/zhou_lab/projects/20191221_references/hg38/annotation/cpg/idx.gz"; then
		to_cpg=/mnt/isilon/zhou_lab/projects/20191221_references/hg38/annotation/cpg/idx.gz
	elif tbmate header $tbk_path"/"$sample_file | grep -q "/mnt/isilon/zhou_lab/projects/20191221_references/hg19/annotation/cpg/idx.gz"; then
		to_cpg=/mnt/isilon/zhou_lab/projects/20191221_references/hg38/annotation/cpg/hg19_to_hg38_idx.gz
	else
		printf "There is not enough information to determine the index file via tbmate header %s" "$genomic_ref"
		exit 1
	fi
	
	chr=/mnt/isilon/zhou_lab/projects/20191221_references/hg38/hg38.fa.fai
else
	printf "illegal argument for -g: %s\n" "$genomic_ref"
	exit 0
fi


# if ! grep -q "tbk" <<< "$tbk_path"; then
# echo "Please provide a tbk directory (such as /mnt/isilon/zhou_lab/projects/20191212_GEO_datasets/GSE59157/tbk)"
# exit 0
# fi

if [[ $type == 'bw' ]]; then
	mkdir -p ${tbk_path//tbk/bw};
elif [[ $type == 'bb' ]] ; then
	mkdir -p ${tbk_path//tbk/bb};
else
	printf "illegal argument for -g: %s\n" "$type"
	exit 0
fi 
mkdir -p ${tbk_path//\/tbk/\/bed};


echo "Converting tbk files in \tbk to bed files in \bed..."
for f in $tbk_path/*; do 

	if grep -q ".gz" <<< "$f" || ! grep -q "tbk" <<< "$f"; then
		continue
	fi
	if [[ $verbose == 'true' ]]; then
		echo $f
	fi

	new_name="${f//tbk\//bed/}"
	new_name="${new_name//.tbk/.bed}"
	if ! test -f "$new_name"; then
		tbmate view -i $to_cpg $f -d | awk '{if($4 != "NA") {print($0)}}' > $new_name
	fi 
done

if [[ $type == 'bw' ]]; then
	if [[ $verbose == 'true' ]]; then
		echo "Converting bed files in \bed to bigwig files in \bw..."
	fi
	for f in ${tbk_path//tbk/bed}/*; do 
		if [[ $verbose == 'true' ]]; then
			echo $f
		fi

		new_name="${f//bed\//bw/}"
		new_name="${new_name//.bed/.bw}"

		if ! test -f "$new_name"; then
			bedGraphToBigWig $f $chr $new_name 
		fi
	done
	
elif [[ $type == 'bb' ]] ; then
	if [[ $verbose == 'true' ]]; then
		echo "Converting bed files in \bed to bigwig files in \bb..."
	fi
	for f in ${tbk_path//tbk/bed}/*; do 
		if [[ $verbose == 'true' ]]; then
			echo $f
		fi

		new_name="${f//bed\//bb/}"
		new_name="${new_name//.bed/.bb}"

		if ! test -f "$new_name"; then 
			bedToBigBed $f $chr $new_name
		fi
	done
fi

exit 0
