#!/bin/bash

#REFERENCE=/n/dat/hg38/hg38.no_alt.fa

while getopts "r:b:o:c:R:" option; do
  case $option in
    r) REGIONFILE=$OPTARG ;;
    b) BEDFILE=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    c) CHROMOSOME=$OPTARG ;;
    R) REFERENCE=$OPTARG ;;
  esac
done

if [ -z ${REFERENCE+x} ]
then
  echo "supply a reference Fasta file with flag -R"
  exit 1
fi

if [ -z ${OUTPUT+x} ]
then
  echo "specify output with flag -o"
  exit 1
fi

if [ -z ${REGION+x} ]
then
  if [ -z ${BEDFILE+x} ]
  then
    echo "supply either a region file with flag -r or a bedfile with flag -b"
    exit 1
  else
  awk '{print $1":"$2"-"$3}' $BEDFILE > temp.regionfile
  REGIONFILE=temp.regionfile
  fi
fi

samtools faidx -n0 $REFERENCE -r $REGIONFILE > temp.fasta

if [ -z ${CHROMOSOME+x} ]
then
  python3 /n/scripts/mini_meow/filterCpG.py NONE temp.fasta $OUTPUT
else
  python3 /n/scripts/mini_meow/filterCpG.py $CHROMOSOME temp.fasta $OUTPUT
fi


rm -rf temp.*