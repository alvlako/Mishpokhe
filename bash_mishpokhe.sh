#!/bin/bash
echo $1
echo $2
echo $3

# FIX: rewrite to numpy, seems to be slow

# print number of unique proteins in clusters
cl_prot=$(awk '{print $1}' "$1" | sort -u | sort -n | wc -l)
echo "$cl_prot"
# print number of unique proteins in normal search results
norm_prot=$(awk '{print $2}' "$2" | sort -u | sort -n | wc -l) 
echo "$norm_prot"

# check if files for scores (hits numbers for protein/all prot)
# in the cluster and norm search exist
if [ -f "${1}_scores" ];then
    echo "scores file exists and will be overwritten"
fi
rm -f "${1}_scores"
touch "${1}_scores"

if [ -f "${2}_scores" ];then
    echo "scores file exists and will be overwritten"
fi
rm -f "${2}_scores"
touch "${2}_scores"

if [ -f "${3}_all_score" ];then
    echo "scores file exists and will be overwritten"
fi
rm -f "${3}_all_score"
touch "${3}_all_score"

# go through cluster search file, $1 is for query protein accession
# the loop counts number of target prot matches for each query prot
for ind in $(awk '{print $1}' "$1" | sort -u | sort -n)
do
  paste <(echo "$ind") <(awk -v var="$ind" '$1==var { print $0 }' "$1" | wc -l) >> ${1}_scores
done
# Produces lines with empty 1st column FIX

# go through normal search file, $2 is for TARGET protein accession, CHANGE for later,
# now it because the $1 for the query is filled with the genomes ID because of the format
# for clustersearch. It doesnt matter now as query = target
# the loop counts number of target prot matches for each query prot
for ind in $(awk '{print $2}' "$2" | sort -u | sort -n)
do
  paste <(echo "$ind") <(awk -v var="$ind" '$2==var { print $0 }' "$2" | wc -l) >> ${2}_scores
done


# add the actual score column, dividing number of matches for prot by all matches
# separately for clustersearch and normal search

awk -v var="$cl_prot" '{print $0, $2/var}' OFS='\t' ${1}_scores > tmp1
cat tmp1 > ${1}_scores

awk -v var="$norm_prot" '{print $0, $2/var}' OFS='\t' ${2}_scores > tmp1
cat tmp1 > ${2}_scores

# merge 2 files by column to have both clustersearch and normal search scores
# in one file (for each protein)
# I dont know what is the 2nd column in file, but 1st must be the prot id,
# 3rd and 4th are clustersearch and normal search scores correspondingly

join  -t $'\t' -j 1 -o 1.1,1.2,1.3,2.3 <(sort -k1 ${1}_scores) <(sort -k1 ${2}_scores) > ${3}_all_score