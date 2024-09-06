#!/bin/bash

# by default calculates per position averages. Can return per read averages instead.

# to run this, activate conda environment rust_plus
# usage: bash /n/scripts/mini_meow/miniMeow.sh -b bamfile -p positions  -o outputfilename -m 5mCG_5hmCG

MODE="5mCG"
READAVG=0
REMOVENA=0

while getopts "b:p:r:o:m:Rn" option; do
  case $option in
    b) BAM=$OPTARG ;;
    p) POSITIONS=$OPTARG ;;
    r) REGION=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    m) MODE=$OPTARG ;;
    R) READAVG=1 ;;
    n) REMOVENA=1 ;;
  esac
done

if [ -z ${BAM+x} ]
then
    echo "no Bam file provided"
    exit 1
fi

if [ -z ${POSITIONS+x} ]
then
    echo "no positions file provided"
    exit 1
fi

if [ -z ${OUTPUT+x} ]
then
    echo "no output name provided"
    exit 1
fi

if [ $MODE != "5mCG" ] && [ $MODE != "5mCG_5hmCG" ]
then
    if [ $MODE != "5mCG_5hmCG_merge" ]
    then
        echo "$MODE not recognized as a mode. Use 5mCG or 5mCG_5hmCG"
        exit 1
    fi
fi

if [ $READAVG -eq 0 ]
then
    echo "running mpileup version"
    if [ -z ${REGION+x} ]
    then
        samtools mpileup --positions $POSITIONS --ff \
        SUPPLEMENTARY,UNMAP,SECONDARY,QCFAIL,DUP -Q 1 -q 1 --no-output-ins --no-output-ins \
        --no-output-del --no-output-del --no-output-ends --output-mods \
        -o temp.pileup.tsv $BAM
    else
        samtools mpileup -r $REGION --positions $POSITIONS --ff \
        SUPPLEMENTARY,UNMAP,SECONDARY,QCFAIL,DUP -Q 1 -q 1 --no-output-ins --no-output-ins \
        --no-output-del --no-output-del --no-output-ends --output-mods \
        -o temp.pileup.tsv $BAM
    fi
    

    if [ $MODE == "5mCG" ]
    then
        rustprogram="/n/scripts/mini_meow/./methyl_pileup"  
    elif [ $MODE == "5mCG_5hmCG" ]
    then
        #r10 mode splits 5hmCG frequency into unmethylated and methylated
        rustprogram="/n/scripts/mini_meow/./methyl_pileup_5hmc_split"
    elif [ $MODE == "5mCG_5hmCG_merge" ]
    then
        #merge mode reports the cumulative probability of 5hmCG and 5mCG
        rustprogram="/n/scripts/mini_meow/./methyl_pileup_merge" 
    fi

    paste <(cut -f 1,2 --output-delimiter='.' temp.pileup.tsv) <(cut -f 5 temp.pileup.tsv | $rustprogram T) \
    | sort -k 1,1 -k 2,2n | tr '.' '\t' | sort -k1,1 -k 2,2n > temp.output.tsv

    if [ $REMOVENA -eq 1 ]
    then
        awk '{if($3!="NA"){print $0}}' temp.output.tsv | tr ' ' '\t' > $OUTPUT
    else
        cat temp.output.tsv > $OUTPUT
    fi

    rm -rf temp.output.tsv
    rm -rf temp.pileup.tsv
else
    if [ -z ${REGION+x} ]
    then
        echo "Per read averages only supported within chromosomes currently. Please specify a chromosome with flag -r"
        exit 1
    else
        # get readname, flags, chromosome, left most position
        samtools view --no-header $BAM | awk -v region=$REGION '{if($3==region){print $1,$2,$5,$3,$4}}' | tr ' ' '\t' > temp.reads.tsv
        # get some of likelihood and then averages
        #sums=( $(samtools view --no-header $BAM | awk -v region=$REGION '{if($3==region){print $0}}' | tr '\t' '\n' | grep ML: | sed 's/ML:B:C,/ /g' |\
        #sed 's/ML:B:C/ /g' | tr ',' '+') )
        
        sums=( $(samtools view --no-header $BAM | awk -v region=$REGION '{if($3==region){print $0}}' | awk '{if(NF == 25 || NF == 24){print $24}else{if(NF==23){print $23}else{print "500,500"}}}' | \
        tr '\t' '\n' | sed 's/ML:B:C,/ /g' | sed 's/ML:B:C/500+500/g' | tr ',' '+') )

        avgs=( $( \
        for sum in ${sums[@]};\
        do \
        subsum=$(echo $sum | grep -o "+" | wc -l); \
        tot=$(echo "$subsum + 1" | bc); \
        numerator=$(echo $sum | bc); \
        echo "scale=4;$numerator / $tot / 255.0 " | bc; \
        done ) )

        # get cigars and calculate reference length
        cigars=( $( samtools view --no-header $BAM | awk -v region=$REGION '{if($3==region){print $6}}' ) )
        reflengths=( $( for cigar in ${cigars[@]};\
        do \
            operations=( $(echo $cigar | tr '[:digit:]' ' ' | tr -s ' ') ); \
            lengths=( $( echo $cigar | tr '[:alpha:]' ' ') ); \
            total=0; \
            for i in $(seq 1 ${#operations[@]}); \
            do \
                op=${operations[$i]}; \
                if [[ $op == "M" ]] || [[ $op == "D" ]] || [[ $op == "N" ]]; \
                then \
                    let total+=${lengths[$i]}; \
                fi; \
            done;\
            echo $total;\
        done ))
        #echo "Sample,type,readName,bitwiseflags,MAPQ,chromosome,refStart,refEnd,meanMethylationFrequency" > $OUTPUT
        paste <( paste <( cat temp.reads.tsv ) <( echo ${reflengths[@]} | tr ' ' '\n' ) ) <( echo ${avgs[@]} | tr ' ' '\n' ) | sed 's/1.9607/NA/g' > $OUTPUT
    fi
fi
