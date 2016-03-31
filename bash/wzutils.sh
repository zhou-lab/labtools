function sumbed {
  awk '{a+=$3-$2}END{print a}' $1
}
