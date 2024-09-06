# Mini-Meow

A small utility set for ONT data that:
     - generates a positions file of CpG locations give a reference fasta file
     - pulls average probability of methylation at each site or for each read

## Installation

This utility relies on bash shell scripts and was designed to run in a Linux environment. It was tested on Ubuntu and RedHat operating systems.

Clone this repository. 
`makePositions` requires Python >= 3.9.0 and the numpy library (>=1.25.0). 
`miniMeow` is used to pull methylation frequencies and requires the following:
    - Rust >= 1.72.1
    - Samtools >= 1.17

## Usage

Run `makePositions.sh` to generate a CpG positions file given a bed file of regions of interest and a fasta file (with index).

Required flags:
    - `-R` : path to a reference fasta file
    - `-o` : output path and filename
    - `-r` or `-b`
        - `-r` : a file containing regions of interest in the format `chrX:12345-2345561`
        - `-b` : a bed file of regions of interest
Options:
    - `-c` : restrict to one chromosome specified here -- significantly speeds things up.


Run `miniMeow.sh` using your generated positions file and an aligned and indexed BAM file with methylation calls encoded in MM/ML tags to make either read-average methylation probabilities or position-average methylation probabilties.

Required flags:
    - `-b` : path to input BAM file. Must be aligned, position sorted, and indexed.
    - `-p` : path to a positions file
    - `-o` : output name
    - `-m` : CpG methylation model mode. There are two modes for interpreting `5mCG_5hmCG` methylation and one for `5mCG`.
        - `5mCG` : use this for data generated using models that predict **only** 5mCG methylation (e.g. dna_r10.4.1_e8.2_400bps_sup@v3.5.2_5mCG@v2)
        - `5mCG_5hmCG` : for later models that predict both 5mCG and 5hmCG methylation probabilies and coerce them into the same 0---255 interval (e.g. dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v1). 5hmCG probability is split between the probability of a canonical base and the probability of a 5mCG methylated base. This is consistent with how ModKit reports 5mCG methylation frequencies when ignoring 5hmCG.
        - `5mCG_5hmCG_merge` : for later models that predict both 5mCG and 5hmCG methylation probabilies and coerce them into the same 0---255 interval. Mini_meow reports the combined probability that a base is methylated with either 5mCG or 5hmCG. This is designed to be more similar to the methylation frequencies detected using bisulfite sequencing, which detects 5hmCG but does not distinguish it from 5mCG methylation.
Options:
    - `-r` : restricts analysis to only one specified chromosome, speeding things up.
    - `-n` : remove positions where no information was reported from results. By default if read depth of primary alignments excluding insertions and deletions is too low to report an average, mini_meow reports 'NA' and keeps the position in results.
    - `-R` : use read average mode. This is currently only implemented within one chromosome at a time, so must be run with flag `-r`. Reports the mean methylation probability for all reads rather than all positions. This was designed with `chrM` and numT detection in mind. 

Generally speaking, all of these tools will run faster if you sort your bed file, region file, or position input.

## Examples

### Default positions based mode

Generate positions first:

```
bash makePositions.sh -b myBedFile.bed -R hg38.fa -o myPositions.positions
```

Use the positions file with miniMeow
```
bash miniMeow.sh -b bamfile.bam -p myPositions.positions  -o positionAverages.tsv -m 5mCG_5hmCG
```

Bed file contents:
```
chr19   5236005 5236054
chr20   63216298        63216347
chr1    6781065 6781114
```
Generated positions contents:
```
chr19   5236005
chr19   5236006
chr20   63216298
chr20   63216299
chr20   63216304
chr20   63216305
chr20   63216312
chr20   63216313
chr20   63216328
chr20   63216329
chr20   63216339
chr20   63216340
chr20   63216346
chr20   63216347
chr1    6781065
chr1    6781066
chr1    6781079
chr1    6781080
chr1    6781083
chr1    6781084
```

MiniMeow output:

```
chr1    6781065 89
chr1    6781066 93
chr1    6781079 97
chr1    6781080 95
chr1    6781083 97
chr1    6781084 97
chr1    148962533       89
chr1    148962534       98
chr19   5236005 93
chr19   5236006 93
chr20   63216298        31
chr20   63216299        40
chr20   63216304        40
chr20   63216305        39
chr20   63216312        22
chr20   63216313        21
chr20   63216328        18
chr20   63216329        16
chr20   63216339        24
chr20   63216340        27
chr20   63216346        34
chr20   63216347        21
```

Frequencies are scaled to a 0-100 integer interval. 

Use these frequencies with downstream tools like the full MeOW suite or visualize results directly using your favorite R tools.

![alt text](https://github.com/mgaleyuw/mini_meow/blob/main/positionMethylation.png?raw=true)


### Read Average Mode

Make a subset of your bam file with only reads of interest before using this tool.

Supply that bam file to miniMeow and omit the positions file argument, using -r and -R instead.
```
bash miniMeow.sh -b bamfile.bam -R -r chrM  -o readAverages.tsv -m 5mCG_5hmCG
```

Example output:

```
readName,bitwiseflags,MAPQ,chromosome,refStart,refEnd,meanMethylationFrequency
bba49e20-713c-4855-89ac-7d3b8b42adee,0,60,chrM,1,16569,.2083
aac5eb1e-d995-467b-a61f-291b82a3b91f,0,60,chrM,1,16569,.1990
fddd3a68-5bdd-41ee-a979-7e68735bb036,0,60,chrM,1,16554,.1997
4efddd28-9d97-4ed5-adab-6f3ce2c8c774,0,60,chrM,1,16487,.1575
dc770049-08e9-4023-b8ef-33b14dd7ef38,0,60,chrM,1,16569,.2417
e51bd76d-2062-4260-a525-86f480df3e3c,0,60,chrM,1,16548,.2043
95dd6148-14ad-4662-89fe-a3872b412c58,0,60,chrM,1,16536,.3034
14f4ee98-64ac-406f-976f-e35a93f7ca4f,0,60,chrM,1,16569,.2612
7e9050ae-4382-4abb-ac6a-269f0ea321a1,0,60,chrM,1,16569,.1243
```

After merging outputs for different bam files and/or samples, you can visualize results in R.

![alt text](https://github.com/mgaleyuw/mini_meow/blob/main/readMethylation.png?raw=true)