#!/bin/bash

meme_matrix=$1
out_dir=$(dirname $1)

start_line=$(grep -n "letter-probability matrix" $meme_matrix | awk -v FS=":" '{print $1}')
start_line_array=( $start_line )
width_matrix=$(grep "letter-probability matrix" $meme_matrix | sed 's/.*w= \([0-9][0-9]*\) .*/\1/')
width_matrix_array=( $width_matrix )


for ((i=0; i<${#start_line_array[@]}; i++))
do
    matrix=$(tail -n +${start_line_array[i]} $meme_matrix | head -n ${width_matrix_array[i]})
    line1=$(head -1 <(echo -e "$matrix"))
    line1=( $line1 )
    nsites=${line1[7]}
    line_name_matrix=$((start_line_array[i] -2))
    name_matrix=$(sed -n "${line_name_matrix}p" $meme_matrix  | awk '{print $2}')
    pfm=$(tail -n +2 <(echo -e "$matrix")  | sed -e 's/^\s\s*//' -e 's/\s\s*/\t/g' -e 's/\s\s*$//g'  | awk -v n=$nsites -v OFS="\t" '{print int($1*n)+1,int($2*n)+1,int($3*n)+1,int($4*n)+1}' | sed "1iA\tC\tG\tT")
    header="MATRIX COUNT ASYMMETRIC $name_matrix SIMPLE"
    echo -e "$header\n$pfm" > "$out_dir/$name_matrix.pfm"
done

exit 0

