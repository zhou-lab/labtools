#!/bin/bash

awk -v a=$1 '$3==a||$1=="#CHROM"' $2 | cut -f10- | sed 's/:\(PASS\|LOWQUAL\):[[:graph:]]*//g' | awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}'