# Long read Assembly process and notes for the hourglass dolphin Harua-tai-nui

<img src="images/dolphin.jpg" title="Dolphin jpg" width="200"/>  

- [Long read Assembly process and notes for the hourglass dolphin Harua-tai-nui](#long-read-assembly-process-and-notes-for-the-hourglass-dolphin-harua-tai-nui)
  - [Runs](#runs)
  - [File Conversion and Basecalling](#file-conversion-and-basecalling)
    - [Basecalling](#basecalling)
  - [Filtering](#filtering)
    - [Contaminants](#contaminants)
    - [Length and Quality](#length-and-quality)
  - [Assemblies](#assemblies)
    - [Assembly contig numbers and length](#assembly-contig-numbers-and-length)
  - [Mitochondria](#mitochondria)
  - [Busco](#busco)
  - [Purge Diploid regions and polish](#purge-diploid-regions-and-polish)
    - [Diploid regions](#diploid-regions)
    - [Medaka polish](#medaka-polish)
  - [Final purged polished](#final-purged-polished)
    - [Busco re-comparison](#busco-re-comparison)
  - [Inspector to identfy misassembly](#inspector-to-identfy-misassembly)
  - [Correcting and scaffolding and visualisation](#correcting-and-scaffolding-and-visualisation)
    - [Visualisation](#visualisation)
  - [Variant calls](#variant-calls)
    - [Clair3](#clair3)
    - [DeepVariant](#deepvariant)
    - [Variant call visualisaations](#variant-call-visualisaations)
    - [Phasing](#phasing)
  - [Annotation](#annotation)
    - [Mask Repeats](#mask-repeats)
    - [Braker3](#braker3)
    - [Annotation output](#annotation-output)
  - [Mod Calling](#mod-calling)
  - [Notes et alia](#notes-et-alia)
    - [Blobtools](#blobtools)
      - [`blobtools` is unappealing](#blobtools-is-unappealing)
    - [References (cetacean assemblies)](#references-cetacean-assemblies)
    - [Other genomes](#other-genomes)
    - [Goldrush behaviour](#goldrush-behaviour)

## Runs

Runs were from three different dates, Dec-15 and Dec-16 2022, and Apr-27 2023

## File Conversion and Basecalling

`fast5` to `pod5`

```bash
# example
pod5 convert fast5 ./fast5/*fast5  --output 2022-12-15-pod5s/ --one-to-one ./fast5/
```

### Basecalling

Guppy was done previously (6.X version). As Dorado is the focus of ONT in the future, analyses will focus on Dorado although providing `guppy` here for comparison. Same process for all three dates.

```bash
# example
# DORADO 4.1.0
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.1.0 2022-12-16-pod5s/ > dolphin-12-16-2022-dorado-sup-4-1-0.bam
```

2022-12-15: `Basecalled @ Samples/s: 4.499077e+06`  
2022-12-16: `Basecalled @ Samples/s: 4.611918e+06`

Conversion  
Transform `.bam` output from `dorado` to `fastq`

```bash
samtools fastq dolphin-12-16-2022-dorado-sup-4-1-0.bam > dolphin-12-16-2022-dorado-sup-4-1-0.fastq
```

## Filtering

### Contaminants

- below can all be simplified with the `--contam` parameter in `chopper`.

Mapped out reads against lambdaphage genome, `NC_001416.1` with `minimap2`

```bash
minimap2 -ax map-ont lambda.fasta dolphin-12-16-2022-dorado-sup-4-1-0.fastq > dolphin-12-16-2022.filt.sam
```

Removed mapped reads from `.sam` file, collected the read names, chopped the header, and selected unmapped reads from the `.fastq`. This could be simpler with `bam` to `fastq`.

```bash
samtools view -f 4 -h dolphin.lambda.sam | cut -f1 > readnames-nolambda.txt
seqkit grep -f readnames-nolambda.txt dolphin.fastq -o dolphin.nolambda.fastq
```

```bash
file                                       format  type    num_seqs          sum_len  min_len  avg_len    max_len     Q1     Q2     Q3  sum_gap    N50  Q20(%)  Q30(%)  GC(%)
dolphin-2022-12-15-dorado-sup-4-1-0.fastq  FASTQ   DNA    8,506,969   39,626,636,053        5  4,658.1  1,081,751  1,562  3,245  5,768        0  7,190   78.37   64.42   42.5
dolphin-2022-12-15.nolambda.fastq          FASTQ   DNA    8,321,380   38,894,669,060        5  4,674.1  1,081,751  1,532  3,188  5,844        0  7,289   78.37   65.08  42.48
dolphin-2022-12-16-dorado-sup-4-1-0.fastq  FASTQ   DNA    4,957,169   31,312,286,009        5  6,316.6    680,433  2,123  4,662  8,146        0  9,486   80.49   67.83  42.33
dolphin-2022-12-16.nolambda.fastq          FASTQ   DNA    4,886,144   31,032,667,713        5  6,351.2    680,433  2,099  4,720  8,210        0  9,548   80.49   68.17  42.32
dolphin-2023-04-27-dorado-sup-4-1-0.fastq  FASTQ   DNA   20,423,436  110,159,598,040        5  5,393.8    757,784  1,903  3,592  6,593        0  8,378   69.75   53.35  42.27
dolphin-2023-04-27.nolambda.fastq          FASTQ   DNA   19,721,170  107,664,477,442        5  5,459.3    757,784  1,845  3,633  6,753        0  8,585   69.71   53.61  42.24
```

### Length and Quality

2Kbp length and q10
From this point the `fastq` from all three dates have been concatenated.

```bash
# head and tailcrop at this step
cat dolphin.fastq | chopper -l 2000 -q 10 --threads 40 --headcrop 50 --tailcrop 50 > dolphin.2k.q10.fastq
```

Statistics on length and quality

```bash
# loop over all fastq
for F in *lambda.fastq; do echo $F; seqkit sample -p 0.01 $F | seqkit head -n 10000 | seqkit fx2tab -qln > ${F/fastq/sample\.qln\.txt}; done
```

<img src="images/basecall-stats.png" title="hists" width="800"/><br>
**Three dorado plots and a filtered guppy plot, 10K reads in each**


```bash
file                  format  type    num_seqs          sum_len  min_len  avg_len  max_len     Q1     Q2     Q3  sum_gap    N50  Q20(%)  Q30(%)  GC(%)
guppy.2k.q10.fastq  FASTQ   DNA   21,107,908  145,901,357,913    1,900  6,912.2  198,884  3,257  4,856  8,204        0  8,878   84.37   71.05  42.22

file                  format  type    num_seqs          sum_len  min_len  avg_len  max_len     Q1     Q2     Q3  sum_gap    N50  Q20(%)  Q30(%)  GC(%)
dorado-q10-2k.fastq  FASTQ   DNA   20,274,358  142,566,554,959    1,900  7,031.9  660,796  3,227  5,010  8,387        0  9,037   84.85   68.63  42.21
```

## Assemblies

| file                    | num_seqs | sum_len       | min_len | avg_len     | max_len    | Q1        | Q2       | Q3          | N50       |
| ----------------------- | -------- | ------------- | ------- | ----------- | ---------- | --------- | -------- | ----------- | --------- |
| hifiasm-guppy-q10       | 25,440   | 3,850,819,104 | 1,950   | 151,368.7   | 3,455,131  | 23,512.5  | 61,638.5 | 191,877     | 366,624   |
| goldrush-guppy-q10      | 40,014   | 2,572,444,687 | 1       | 64,288.6    | 28,679,029 | 743       | 2,808.5  | 9,288       | 5,168,950 |
| goldrush-guppy-2kfilter | 25,504   | 2,563,972,225 | 2,000   | 100,532.2   | 28,679,029 | 3,004     | 6,978.5  | 16,188      | 5,230,893 |
| goldrush-dorado-q10     | 29,863   | 2,443,538,685 | 1       | 81,825      | 24,473,496 | 646       | 2,000    | 5,750.5     | 6,142,064 |
| nextdenovo-guppy-q10    | 1,090    | 2,320,259,928 | 15,220  | 2,128,678.8 | 28,442,925 | 168,222   | 406,787  | 2,577,623   | 6,883,491 |
| nextdenovo-dorado-q10   | 1,098    | 2,321,300,339 | 16,390  | 2,114,116.9 | 27,479,023 | 166,555   | 420,018  | 2,472,563   | 6,549,365 |
| raven-guppy-q10         | 1,289    | 2,453,988,083 | 1,380   | 1,903,792.2 | 24,910,266 | 144,240   | 269,805  | 1,913,506   | 7,224,846 |
| raven-dorado-q10        | 1,243    | 2,443,232,195 | 1,379   | 1,965,593.1 | 39,055,960 | 141,447.5 | 258,182  | 1,955,587   | 7,919,957 |
| raven-dorado-q10-S1     | 1,138    | 2,404,292,365 | 1,384   | 2,112,734.9 | 41,189,193 | 158,031   | 421,191  | 2,416,573   | 6,424,551 |
| raven-dorado-q16        | 1,187    | 2,399,149,723 | 10,313  | 2,021,187.6 | 30,559,859 | 165,888.5 | 464,783  | 2,335,714.5 | 6,305,064 |

filtered goldrush is with contigs < 2 Kbp removed as many contigs are 1bp or so (??)

Get all length N (N10, N20, etc.)

```bash
for F in *.fasta; do echo $F; seqkit fx2tab -qln $F > ${F/fasta/q
ln\.txt}; done
```

### Assembly contig numbers and length

<img src="images/guppy-stats.png" title="NX plots" width="600"/>
**NX plots for different assemblies**

<img src="images/NX-stats-zoom.png" title="NX plots zoomed" width="600"/>
**NX plots zoomed**

<img src="images/raven-hist.png" title="raven contigs" width="600"/>
**As an example, contig lengths for Raven**

<img src="images/nextdenovo-hist.png" title="raven contigs" width="600"/>
**As an example, contig lengths for Nextdenovo**

## Mitochondria

```bash
cat mtdna-mapped.fastq | chopper --minlength 14000 --maxlength 17000 -t 8 > mtdna-filter.fastq
seqkit sample -p 0.05 mtdna-filter.fastq > mtdna-sample.fastq
```

## Busco

This is actually `compleasm.py` from Heng Li
`Hifiasm`
Version: `0.19.5-r587`

`nextDenovo`
Version: `[24589 INFO] 2023-05-15 13:46:21 version:2.5.2`

`goldrush`
Version: `goldrush v1.0.1`

`Raven`
Version: `1.8.1``

Total Laurasthenia is N:12234

| Assembler            | Single | Duplicate | Fragment | Intragenic | Missing |
| -------------------- | ------ | --------- | -------- | ---------- | ------- |
| Guppy hifiasm q10    | 8635   | 717       | 604      | 5          | 2273    |
| Guppy nextDenovoq10  | 11937  | 118       | 82       | 1          | 96      |
| Dorado nextDenovoq10 | 11936  | 116       | 78       | 0          | 104     |
| Guppy Goldrush q10   | 11795  | 210       | 98       | 0          | 131     |
| Guppy Raven q10      | 12017  | 109       | 58       | 0          | 50      |
| Dorado Raven q10     | 12024  | 107       | 55       | 0          | 48      |
| Dorado Raven q16     | 12016  | 104       | 58       | 0          | 56      |
| Dorado Goldrush q10  | 11826  | 191       | 92       | 1          | 124     |

## Purge Diploid regions and polish

### Diploid regions

```bash
## estimated commands (not tracked)
purge_haplotigs hist -b raven-dorado-no-second-sorted.bam -G ../raven-dorado-q10.fasta -t 40
purge_haplotigs cov -i raven-dorado-no-second-sorted.bam.gencov -l 9  -m 31  -h 95
purge_haplotigs purge -g ../raven-dorado-q10.fasta -c coverage_stats.csv
```

![purged-hist png](images/raven-dorado-no-second-sorted.bam.histogram.png)  

### Medaka polish

Test effects of `medaka` polishing. Make sure to use through the `mamba` `medaka` environment. `medaka` version is 1.8:

```bash
# have to run with decreased batch size using -b
medaka_consensus -i dolphin-dorado-q10-2k.fastq -d ../haplotigs/raven-dorado-curated.fasta -o medaka-polished -t 40 -m r1041_e82_400bps_sup_v4.1.a -b 40

Checking program versions
This is medaka 1.8.0
Program Version Required Pass
bcftools 1.17 1.11 True
bgzip 1.17 1.11 True
minimap2 2.26 2.11 True
samtools 1.17 1.11 True
tabix 1.17 1.11 True 
```

## Final purged polished

`contigs:895  length:2,383,719,527    min:1,396  mean:2,663,373.8  longest:39,026,300`

small contig is a contmainant from Bos taurus, remove.

check depth

```bash
minimap2 -ax map-ont --secondary=no -t 40 raven-purged-polished-clean.fasta dolphin-dorado-q10-2k.fastq > raven-purged-polished-clean.sam
nice -n 10 samtools sort -@ 20 raven-purged-polished-clean.sam > raven-purged-polished-clean-sort.bam
nice -n 10 samtools depth -@ 20 raven-purged-polished-clean-sort.bam > raven-purged-polished-clean-depth.txt
## all contigs contained in contig-names.txt
```

To get mean and median depths we just subsample the file so that the smallest contig (24kbp) should be sampled 100 times (10,000,000 lines). this is using `shuf`

```bash
shuf -n 10000000 raven-purged-polished-clean-depth.txt > raven-purged-polished-clean-depth-sample.txt
```

Turns out this is not the problem.

Try to look at the exact sequence in these contigs.

```bash
# use bed file bad-contigs.bed that contians contig names and locations from 11Kbp to 12 Kbp for each contig
seqkit subseq --bed bad-contigs.bed raven-bad-contigs.fasta > raven-bad-contigs-subseq.fasta

# get stats
seqkit stats raven-bad-contigs-subseq.fasta

# stats
# raven-bad-contigs-subseq.fasta  FASTA   DNA         48   48,000    1,000    1,000    1,000
```

These contigs are clearly repetitive. Use `jellyfish` to characterise kmer counts.

```bash
# example
for K in 13 17 25 29; do jellyfish count -m $K -s 100M -t 20 raven-bad-contigs.fasta -o raven-bad.${K}.jf; don

# histogram
for K in 13 17 25 29; do jellyfish histo raven-bad.${K}.jf > raven-bad.${K}.histo.txt; done

```

To get all counts per contig, split contigs using `seqkit` and then count for each contig, outputting a `.jf` for each. This is (currently) in the `raven-purged-polished-clean.fasta.split` dir.

### Busco re-comparison

Looks like repeating with cetartiodactyla_odb10.2019-11-20 is a good idea, 13335 markers. Also eutheria_odb10.2019-11-20 11366; mammalia 9226

```bash
compleasm.py run -a ../raven-purged-polished-sorted.fasta -o euk -t 20 -l eukaryota
compleasm.py run -a ../raven-purged-polished-sorted.fasta -o euk -t 20 -l mammalia
compleasm.py run -a ../raven-purged-polished-sorted.fasta -o euk -t 20 -l eutheria
compleasm.py run -a ../raven-purged-polished-sorted.fasta -o euk -t 20 -l cetartiodactyla
```

On polished and purged assemblies. Single copy increases, duplicated decreases by 50%, fragment decreases slightly, missing remains the same.

| Assembler            | Single | Duplicate | Fragment | Intragenic | Missing |
| -------------------- | ------ | --------- | -------- | ---------- | ------- |
| _Laurasiatheria_ raven   | 12024  | 107       | 55       | 0          | 48      |
| _Certartiodactyla_-purge-polish     | 13118  | 118       | 67       | 0          | 32      |
| _Laurasiatheria_-purge-polish   | 12034  | 99       | 53       | 0          | 48      |
| _Eutheria_-purge-polish     | 11171  | 75       | 49       | 1          | 70      |
| _Mammalia_-purge-polish     | 9088  | 63       | 34       | 0          | 41      |
| _Eukarya_-purge-polish     | 250  | 5       | 0       | 0          | 0      |

Sousa chinensis
S:92.71%, 11342
D:1.13%, 138
F:3.38%, 414
I:0.02%, 3
M:2.75%, 337

T. truncatus
S:98.14%, 12006
D:0.99%, 121
F:0.38%, 46
I:0.00%, 0
M:0.50%, 61
N:12234

## Inspector to identfy misassembly

https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02527-4#Sec12
[Inspector](https://github.com/ChongLab/Inspector)
[CRAQ](https://github.com/JiaoLaboratory/CRAQ)

## Correcting and scaffolding and visualisation

Previously (June 2023 or so) used [RagTag](https://github.com/malonge/RagTag/wiki/scaffold) to correct and scaffold against the [N. asiaeorientalis](https://www.ebi.ac.uk/ena/browser/view/GCA_026225855?show=assembly-stats) genome.

This is not good. Scaffolding should be **based off closely related species, see phylogeny below.**

<img src="images/annotated-phylogeny.png" title="phylogeny png" width="600"/><br>

Genomes to align against (see phylogeny below for relative placements):

- _Tursiops truncatus_ (bottlenose dolphin)
- _Tursiops aduncus_ (Indo-pacific)
- _Lagenorhynchus obliquidens_ (Pacific white-sided)
- _Sousa chinensis_ (Indo-pacific humpback) CNGBdb: CNP0000397
- _Lagnorhynchus acutus_ (Atlantic white-sided)
  
Summary of genomes
|file |num_seqs |sum_len |min_len |avg_len |max_len |Q1 |Q2 |Q3 |sum_gap |N50 |Q20(%) |Q30(%) |GC(%) |
|-------|-----------|-----------|-----------|-----------|-----------|---|---|---|-----------|-------|-------|-------|-------|
|_Lagenorhynchus acutus_ |105,133 |2,328,596,073 |200 |22,149 |174,168,249 |239 |296 |482 |0 |103,333,746 |0 |0 |41.2 |
|_Lagenorhynchus obliquidens_ |5,162 |2,333,909,671 |202 |452,132.8 |180,730,733 |1,436 |2,623.5 |5,825 |0 |107,447,310 |0 |0 |41.11 |
|_Sousa chinensis_ |23 |2,439,440,112 |35,025,918 |106,062,613.6 |210,152,121 |78,342,071.5 |100,036,560 |124,265,802.5 |0 |111,358,879 |0 |0 |39.57 |
|_Tursiops aduncus_ |12,471 |2,505,817,531 |904 |200,931.6 |186,145,090 |1,766.5 |2,880 |10,063.5 |0 |111,961,311 |0 |0 |39.9 |
|_Tursiops truncatus_ |362 |2,378,522,213 |1,017 |6,570,503.4 |183,742,880 |32,203 |64,819.5 |112,540|0  |108,430,135 |0|0   |41.33|

Remove _Lagnorhynchus acutus_.

- _S. chinensis_ (indo-pacific)
  - anchored ~90.7% of genome scaffolds into 22 chromosomes<br>
<img src="images/S-chinensis-hi-c.png" title="hi-c png" width="200"/><br>

First correction onto _S. chinensis_

```bash
ragtag.correct.fasta  FASTA   DNA      9,860  2,383,719,527      102  241,756.5  26,550,021  4,235.5  8,849  25,794.5        0  5,631,925       0       0  41.39
```

Onto _T. truncatus_

```bash
ragtag.correct.fasta  FASTA   DNA      1,382  2,383,719,527      154  1,724,833.2  28,376,917  48,558  200,608.5  1,695,993        0  6,845,537       0       0  41.39
```

Output and comparisons
|file|num_seqs|sum_len|min_len|avg_len|max_len|Q1|Q2|Q3|N50|GC(%)|
|----|--------|-------|-------|-------|-------|--|--|--|---|-----|
|L-obliquidens|2,109|2,383,719,527|64|1,130,260.6|28,343,814|12,713|61,795|468,697|6,845,537|41.39|
|S-chinensis|9,860|2,383,719,527|102|241,756.5|26,550,021|4,235.5|8,849|25,794.5|5,631,925|41.39|
|T-aduncus|2,875|2,383,719,527|65|829,119.8|22,069,092|7,955|34,411|151,590|6,274,499|41.39|
|T-truncatus|1,382|2,383,719,527|154|1,724,833.2|28,376,917|48,558|200,608.5|1,695,993|6,845,537|41.39|
|L-obliquidens|1,324|2,383,798,027|64|1,800,451.7|182,132,532|7,501|21,306|73,381|108,258,426|41.38|
|S-chinensis|2,067|2,384,498,827|102|1,153,603.7|180,111,841|4,172|9,767|38,069.5|105,079,198|41.37|
|T-aduncus|1,983|2,383,808,727|65|1,202,122.4|181,556,247|5,396|15,790|49,761.5|108,424,675|41.38|
|T-truncatus|591|2,383,798,627|154|4,033,500.2|182,741,771|13,517|61,399|169,946|108,574,935|41.38|

### Visualisation


## Variant calls

Used the purged, medaka-polished dorado assembly, sorted, and indexed for use in `Clair3`. This is **not** on the scaffolded contigs.

```bash
minimap2 -ax map-ont raven-purged-polished.fasta ../dolphin-dorado-q10-2k.fastq | samtools sort -o raven-purged-polished.bam
```

### Clair3

```bash
run_clair3.sh --bam_fn=raven-purged-polished.bam \
--ref_fn=raven-purged-polished.fasta --threads=40 \
--platform="ont" \
--model_path=rerio/clair3_models/r1041_e82_400bps_sup_v410 \
--output=clair3-calls --include_all_ctgs

## Clair3 final output file: ${OUTPUT_DIR}/merge_output.vcf.gz
```

Nicely, `DeepVariant` (below) comes with some handy plotting tools.

```bash
DATA_DIR="/cratch/olin/dolphin/hourglass-dolphin/dorado-assemblies/variants/clair3-calls/"
OUTPUT_DIR="/scratch/olin/dolphin/hourglass-dolphin/dorado-assemblies/variants/clair3-calls/"
sudo docker run -v "${DATA_DIR}:/input" -v "${OUTPUT_DIR}:/output" google/deepvariant:"${BIN_VERSION}" /opt/deepvariant/bin/vcf_stats_report --input_vcf /input/merge_output.vcf.gz --outfile_base /output/clair3
```

There is a clear trough of genotype quality scores at 14. We implement a cutoff there:

```bash
rtg vcffilter -Z -q 14 -i merge_output.vcf.gz -o filtered-q14.vcf.gz
```

### DeepVariant

Have to install `docker` for this.

```bash
# docker commit: c2de0811708b6d9015ed1a2c80f02c9b70c8ce7b
curl -fsSL https://get.docker.com/ | sh

## docker install
BIN_VERSION="1.5.0"
sudo docker pull google/deepvariant:"${BIN_VERSION}"

## test data
sudo docker run -v "${INPUT_DIR}":"/input" -v "${OUTPUT_DIR}":"/output"   google/deepvariant:"${BIN_VERSION}"  /opt/deepvariant/bin/run_deepvariant --model_type=ONT_R104 --ref=/input/ucsc.hg19.chr20.unittest.fasta --reads=/input/NA12878_S1.chr20.10_10p1mb.bam --output_vcf=/output/output.vcf.gz --output_gvcf=/output/output.g.vcf.gz --intermediate_results_dir /output/intermediate_results_dir --num_shards=1

## repeated the above but substitute the INPUT_DIR and OUTPUT_DIR variables with the real dirs, and substituted the names. Allowed 40 cores. Looks like it needs the full path
INPUT_DIR="/scratch/olin/dolphin/hourglass-dolphin/dorado-assemblies/variants/deepvariant-calls/input"

OUTPUT_DIR="/scratch/olin/dolphin/hourglass-dolphin/dorado-assemblies/variants/deepvariant-calls/output"

sudo docker run -v "${INPUT_DIR}":"/input" -v "${OUTPUT_DIR}":"/output" google/deepvariant:"${BIN_VERSION}" /opt/deepvariant/bin/run_deepvariant --model_type=ONT_R104 --ref=/input/raven-purged-polished.fasta --reads=/input/raven-purged-polished.bam --output_vcf=/output/output.vcf.gz --output_gvcf=/output/output.g.vcf.gz --intermediate_results_dir /output/intermediate_results_dir --num_shards=40

```

### Variant call visualisaations

Some quick plots of the results (Clair3 on top, DeepVariant on bottom)
<img src="images/cl3-variant-types.png" title="Clair3" width="600"/><br>
<img src="images/dv-variant-types.png" title="DeepVariant" width="600"/><br>

<img src="images/cl3-depth.png" title="Clair3" width="300"/> <img src="images/dv-depth.png" title="DeepVariant" width="300"/><br>

<img src="images/cl3-genotype-quality.png" title="Clair3" width="300"/> <img src="images/dv-genotype-quality.png" title="DeepVariant" width="300"/><br>

<img src="images/cl3-base-changes.png" title="Clair3" width="800"/><br>
<img src="images/dv-base-changes.png" title="DeepVariant" width="800"/><br>

<img src="images/cl3-indel-size.png" title="Clair3" width="600"/><br>
<img src="images/dv-indel-size.png" title="DeepVariant" width="600"/><br>

For `DeepVariant` there is a clear trough of genotype quality scores at 20. We implement a cutoff there:

```bash
rtg vcffilter -Z -q 20 -i output.vcf.gz -o dv-filtered-q20.vcf.gz
```

We can intersect the `Clair3` and `deepvariant` calls to get a higher confidence set of heterozygous calls.

```bash
# make a bed file from the VCF
bcftools query -f '%CHROM\t%POS0\t%POS0\n' cl3.filtered-q14.vcf.gz > cl3.filtered-q14.bed
bcftools query -f '%CHROM\t%POS0\t%POS0\n' ../deepvariant-calls/output/dv-filtered-q20.vcf.gz > dv-filtered-q20.bed

# intersect the beds
bedtools intersect -a cl3.filtered-q14.bed -b dv-filtered-q20.bed > cl3.filtered-q14.both.bed
bedtools intersect -a dv-filtered-q20.bed -b cl3.filtered-q14.bed > dv-filtered-q20.both.bed

# intersect the intersected bed with the VCFs
bedtools intersect -header -a cl3.filtered-q14.vcf.gz -b cl3.filtered-q14.both.bed > cl3.filtered-q14.both.vcf

bedtools intersect -header -a ../deepvariant-calls/output/dv-filtered-q20.vcf.gz -b dv.filtered-q20.both.bed > dv.filtered-q20.both.vcf

OUTPUT_DIR="/scratch/olin/dolphin/hourglass-dolphin/dorado-assemblies/variants"

INPUT_DIR="/scratch/olin/dolphin/hourglass-dolphin/dorado-assemblies/variants"

sudo docker run -v "${INPUT_DIR}:/input" -v "${OUTPUT_DIR}:/output" google/deepvariant:"${BIN_VERSION}" /opt/deepvariant/bin/vcf_stats_report --input_vcf /input/cl3-dv-merged.vcf.gz --outfile_base /output/cl3-dv-merge

```

The above is incorrect and inefficient. The proper way to do this is with `bcftools norm` and `bcftools isec`. Thus:

```bash
bcftools norm -D dv-filtered-q20.vcf.gz > dv-filtered-q20-norm.vcf.gz
bcftools norm -D cl3.filtered-q14.vcf.gz > cl3.filtered-q14-norm.vcf.gz
bcftools isec cl3.filtered-q14-norm.vcf.gz dv-filtered-q20-norm.vcf.gz
```
to generally count variant types use (or similar):

```bash
bcftools view -i 'ALT!="."' cl3.filtered-q14.vcf.gz | awk 'length($5)>length($4)' | wc
```

### Phasing

```bash
whatshap --version 2.0

# standard output

whatshap phase -o phased.vcf --ignore-read-groups --reference=deepvariant-calls/input/raven-purged-polished.fasta cl3-dv-merged.vcf.g
z deepvariant-calls/input/raven-purged-polished.bam                                                                                                                                                                        
This is WhatsHap 2.0 running under Python 3.10.12                                                                                                                                                                          
Working on 1 sample from 1 family                                                                                                                                                                                          
WARNING: Skipping duplicated position 3023 on chromosome 'Utg637368' Hiding further warnings of this type, use --debug to show

# Working on contig Utg634904 in individual SAMPLE
Found 5493 usable heterozygous variants (0 skipped due to missing genotypes)
Found 18964 reads covering 5491 variants
Kept 18649 reads that cover at least two variants each
Selected 2462 most phase-informative reads covering 5491 variants
Best-case phasing would result in 1 non-singleton phased block (0 singletons). 
Phasing 1 sample by solving the MEC problem ...
Largest block contains 5491 variants (100.0% of accessible variants) between position 1482 and 2363485


Maximum memory usage: 5.837 GB
Time spent reading BAM/CRAM:                 3624.7 s
Time spent parsing VCF:                        67.0 s
Time spent selecting reads:                  1698.2 s
Time spent phasing:                          18933.8 s
Time spent writing VCF:                       122.7 s
Time spent finding components:                 93.5 s
Time spent on rest:                           678.8 s
Total elapsed time:                          25218.6 s

```

Output summary is at least partly in `variants/phased-stats.tsv`.

Phasing for visualisation

```bash
bgzip phased.vcf
tabix -p vcf phased.vcf.gz

whatshap haplotag --ignore-read-groups --output-threads=40 -o haplotagged.bam --reference deepvariant-calls/input/raven-purged-polished.fasta phased.vcf.gz deepvariant-calls/input/raven-purged-polished.bam

== SUMMARY ==
Total alignments processed:                  31869302
Alignments that could be tagged:             17539735
Alignments spanning multiple phase sets:            0
Finished in 5905.2 s

```

## Annotation

### Mask Repeats

Best practice before annotation so that the gene strucuture is not abnormal. This _will not_ install from `mamba`, it hangs on `Executing transaction`. Installing from [source](https://www.repeatmasker.org/RepeatMasker/) instead.

Configured to use `HMMER3.1 & DFAM` as defaults. However, this did not seem to solve the issue, so specificed nhmmer, see below command. In addition, would not recognize artiodactyl so used mammalia.

```bash

2023-09-18
RepeatMasker -pa 10 -species mammalia -html -default_search_engine nhmmer raven-purged-polished.fasta

Search Engine: HMMER [ 3.3.2 (Nov 2020) ]

Using Master RepeatMasker Database: /home/olin/software/RepeatMasker/Libraries/RepeatMaskerLib.h5

Title    : Dfam
Version  : 3.7
Date     : 2023-01-11
Families : 19,768

Species/Taxa Search:

Mammalia [NCBI Taxonomy ID: 40674]

Lineage: root;cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata <chordates>;Vertebrata <vertebrates>;Gnathostomata <vertebrates>;Teleostomi;Euteleostomi

214 families in ancestor taxa; 10524 lineage-specific families

analyzing file raven-purged-polished.fasta
```

file name: raven-purged-polished.fasta
sequences: 895

Total length: 2,383,719,527

GC level: 41.39%

bases masked: 1,282,522,487  ( 53.80%)

|Element|Number of elements|  Length occupied|Percentage of sequence|
|----|--------|-------|-------|
|**SINEs**|935,618|111,448,831|4.68%|
|Alu/B1|101,719|9,178,606|0.39%|
|MIRs|671,880|93,430,603|3.92%|
|**LINEs**|1,438,734|704,812,271|29.57%|
|LINE1|874,495|548,769,192|23.02%|
|LINE2|481,542|133,492,629|5.60%|
|L3/CR1|614,13|15,066,173|0.63%|
|RTE|18,443|6,993,976|0.29%|
|**LTR elements**|1,421,857|236,050,105|9.90%|
|ERVL|244,158|5,6762,273|2.38%|
|ERVL-MaLRs|274,451|74,314,533|3.12%|
|ERV_classI|310,578|46,037,307|1.93%|
|ERV_classII|451,988|33,125,630|1.39%|
|**DNA elements**|1,231,135|151,923,165|6.37%|
|hAT-Charlie|586,825|71,853,271|3.01%|
|TcMar-Tigger|182,633|32,019,456|1.34%|
|Unclassified|184,336|10,917,146|0.46%|
|**Total interspersed repeats**||1,215,151,518|50.98%|
|Small RNA|668681|39,636,861|1.66%|
|Satellites|108,716|6,137,273|0.26%|
|Simple repeats|431,317|17,807,625|0.75%|
|Low complexity|70,698|3,399,449|0.14%|

==================================================<br>
The query species was assumed to be mammalia
RepeatMasker version 4.1.5 , default mode
run with nhmmscan version 3.3.2 (Nov 2020)
FamDB: HMM-Dfam_3.7

Make `.bed` file from repeat annotations for intersection with variant calls. RM2bed from [here](https://github.com/rmhubley/RepeatMasker/blob/master/util/RM2Bed.py).

```bash
python RM2Bed.py raven-purged-polished.fasta.out
### takes ~ 1m to run
```

### Braker3

Use Braker2 for fully de novo annotation followed by Augustus. Braker2 installable via `mamba`.

Get the vertyebrate OrthoDb from [here](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/)
Should cite: *Kuznetsov, D., Tegenfeldt, F., Manni, M., Seppey, M., Berkeley, M., Kriventseva, E. V., & Zdobnov, E. M. (2023). OrthoDB v11: annotation of orthologs in the widest sampling of organismal diversity. Nucleic Acids Research, 51(D1), D445-D451.*

```bash
Could not solve for environment specs
The following package could not be installed
└─ genemark   does not exist (perhaps a typo or a missing channel).
(braker) olin@agnes:/scratch/olin/dolphin/hourglass-dolphin/dorado-assemblies/variants$ braker.pl --workingdir=braker3 --genome=raven-purged-polished.fasta.masked --prot_seq=/data/databases/Vertebrata.fa --threads 40 
Unknown option: threads
#**********************************************************************************
#                               BRAKER CONFIGURATION                               
#**********************************************************************************
# BRAKER CALL: /home/olin/mambaforge/envs/braker/bin/braker.pl --workingdir=braker3 --genome=raven-purged-polished.fasta.masked --prot_seq=/data/databases/Vertebrata.fa --threads 40
# Wed Sep 27 12:51:49 2023: braker.pl version 2.1.6
# Wed Sep 27 12:51:49 2023: REMARK: Protein input detected, BRAKER will be executed in the EP mode (--epmode).
# Wed Sep 27 12:51:49 2023: Configuring of BRAKER for using external tools...
# Wed Sep 27 12:51:49 2023: Found environment variable $AUGUSTUS_CONFIG_PATH. Setting $AUGUSTUS_CONFIG_PATH to /home/olin/mambaforge/envs/braker/config/
# Wed Sep 27 12:51:49 2023: Found environment variable $AUGUSTUS_BIN_PATH. Setting $AUGUSTUS_BIN_PATH to /home/olin/mambaforge/envs/braker/bin/
# Wed Sep 27 12:51:49 2023: Found environment variable $AUGUSTUS_SCRIPTS_PATH. Setting $AUGUSTUS_SCRIPTS_PATH to /home/olin/mambaforge/envs/braker/bin/
# Wed Sep 27 12:51:49 2023: Did not find environment variable $PYTHON3_PATH
# Wed Sep 27 12:51:49 2023: Trying to guess $PYTHON3_PATH from location of python3 executable that is available in your $PATH
# Wed Sep 27 12:51:49 2023: Setting $PYTHON3_PATH to /home/olin/mambaforge/envs/braker/bin
# Wed Sep 27 12:51:49 2023: Did not find environment variable $GENEMARK_PATH  (either variable does not exist, or the path given in variable does not exist). Will try to set this variable in a different way, later.
# Wed Sep 27 12:51:49 2023: Trying to guess $GENEMARK_PATH from location of gmes_petap.pl executable that is available in your $PATH
#*********
# WARNING: Guessing the location of $GENEMARK_PATH failed.
#*********
# Wed Sep 27 12:51:49 2023: ERROR: in file /home/olin/mambaforge/envs/braker/bin/braker.pl at line 2172
$GENEMARK_PATH not set!
```

Export:
`export GENEMARK_PATH=~/software/gmes_linux_64`
`export PROTHINT_PATH=~/software/ProtHint/bin`

```bash
# Wed Sep 27 21:16:04 2023: Did not find environment variable $PROTHINT_PATH
```

**Restart** with `Braker3` not `Braker2`.
`braker.pl version 3.0.3`

`cp` gmes dir into the PrtoHint one.

Test with (make sure in `mamba braker3 env`)
`ProtHint/dependencies/GeneMarkES/check_install.bash`

errors: `Error, installation key is missing or expired`

get key here: <http://topaz.gatech.edu/GeneMark/tmp/GMtool_JrWuM/gm_key_64.gz>

```bash

braker.pl --threads 20 --genome=raven-purged-polished.fasta.masked --prot_seq=/data/databases/Vertebrata.fa 

# Thu Sep 28 09:41:09 2023: Log information is stored in file /scratch/olin/dolphin/hourglass-dolphin/dorado-assemblies/variants/braker/braker.log
#*********
# WARNING: Detected whitespace in fasta header of file /data/databases/Vertebrata.fa. This may later on cause problems! The pipeline will create a new file without spaces or "|" characters and a genome_header.map file to look up the old and new headers. This message will be suppressed from now on!
#*********
```

Errors out:

```bash
ERROR in file /home/olin/mambaforge/envs/braker3/bin/braker.pl at line 5147
Failed to execute: /home/olin/mambaforge/envs/braker3/bin/perl /home/olin//software/gmes_linux_64/gmes_petap.pl --verbose --cores=20 --ES --gc_donor 0.001 --sequence=/scratch/olin/dolphin/hourglass-dolphin/dorado-assemblies/variants/braker/genome.fa  --soft_mask auto 1>/scratch/olin/dolphin/hourglass-dolphin/dorado-assemblies/variants/braker/GeneMark-ES.stdout 2>/scratch/olin/dolphin/hourglass-dolphin/dorado-assemblies/variants/braker/errors/GeneMark-ES.stderr
The most common problem is an expired or not present file ~/.gm_key!
```

changed perl headers with script in `gmes` directory. Now running smoothly??
This is noted [here](https://github.com/Gaius-Augustus/BRAKER/issues/393)

Also started in singularity container on `gru` (error when using `singulairty` on `agnes`)
`singularity` notes are [here](https://hub.docker.com/r/teambraker/braker3)

```bash
# BRAKER CALL: /home/olin/mambaforge/envs/braker3/bin/braker.pl --threads 20 --genome=raven-purged-polished.fasta.masked --prot_seq=/data/databases/Vertebrata.fa
# Thu Sep 28 10:46:35 2023: braker.pl version 3.0.3

# 20 threads
# Thu Sep 28 13:53:49 2023] ProtHint finished.
# Thu Sep 28 16:08:58 2023: optimizing AUGUSTUS parameters

# Fri Sep 29 22:22:12 2023] Output processed
# Fri Sep 29 22:22:13 2023] ProtHint finished.

# singularity
# BRAKER CALL: /opt/BRAKER/scripts/braker.pl --threads 8 --genome=raven-purged-polished.fasta.masked --prot_seq=Vertebrata.fa
# Thu Sep 28 10:58:57 2023: braker.pl version 3.0.3

# 8 threads
# Thu Sep 28 18:47:27 2023] ProtHint finished.
# Fri Sep 29 00:25:47 2023: optimizing AUGUSTUS parameters

# STOPPED run as much slower than out of singularity

```

### Annotation output

|number|min length| max length| avg length|
|------|----------|-----------|-----------|
|25,329|4         | 4,065       |268 |

Number of exons: 129,464
Number of introns: 104,187
Mean exon length: 157
Mean intron length: 1509

intersect with vcf to find variants in each position:

```bash
## general format
bedtools intersect -a cl3-dv-merged.vcf.gz -b braker.introns1.bed -header > cl3-dv-merged.introns.vcf
```

stats:

INTRONS
SN      0       number of samples:      1
SN      0       number of records:      334697
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 303852
SN      0       number of MNPs: 0
SN      0       number of indels:       30852
SN      0       number of others:       0
SN      0       number of multiallelic sites:   1047
SN      0       number of multiallelic SNP sites:       1

EXONS
SN      0       number of samples:      1
SN      0       number of records:      25570
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 24522
SN      0       number of MNPs: 0
SN      0       number of indels:       1048
SN      0       number of others:       0
SN      0       number of multiallelic sites:   26
SN      0       number of multiallelic SNP sites:       0

TOTAL
SN      0       number of samples:      1
SN      0       number of records:      4448218
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 4020783
SN      0       number of MNPs: 0
SN      0       number of indels:       427532
SN      0       number of others:       0
SN      0       number of multiallelic sites:   15871
SN      0       number of multiallelic SNP sites:       17

total exon bp: 20,325,848
total intron bp: 157,218,183
total bp: 2,383,719,527
exon fraction: 0.00852694613
intron fraction: 0.06595

## Mod Calling

Realign basecallers to newly scaffolded genome:

```bash
dorado aligner dorado-assemblies/raven-purged-polished-sorted.fasta mod-dorado-sup-4-1-0.bam > purged-polished-mod-dorado-sup-4-1-0.bam

# onbly primary mappers
samtools view -F 0x0100 -b purged-polished-mod-dorado-sup-4-1-0.bam | samtools sort > purged-polished-mod-dorado-sup-4-1-0-sort.bam

# index
samtools index -@ 40 purged-polished-mod-dorado-sup-4-1-0-sort.bam

```

```bash

# modkit
modkit pileup purged-polished-mod-dorado-sup-4-1-0-sort.bam purged-polished-mod-dorado-sup-4-1-0.bed --log-filepath purged-polished-pileup.log

# pileup log
[src/logging.rs::54][2023-09-26 21:44:29][DEBUG] command line: modkit pileup purged-polished-mod-dorado-sup-4-1-0-sort.bam purged-polished-mod-dorado-sup-4-1-0.bed --log-filepath purged-polished-pileup.log
...
[src/pileup/subcommand.rs::786][2023-09-27 11:27:18][INFO] Done, processed 335536708 rows. Processed ~32878745 reads and skipped ~51878 reads.



```

Try to run `methylartist`

```bash
methylartist segmeth -b align-mod-dorado-sup-4-1-0.sort.bam -i ragtag.scaffold.n-asiaeorientalis.10k.bed -p 40 --ref ragtag.scaffold.n-asiaeorientalis.fasta --motif C
```

fails with error

Get rid of multimappers and reconstruct files.

```bash
samtools view -F 0x0100 -b -o primary-align-mod-dorado-sup-4-1-0.sort.bam align-mod-dorado-sup-4-1-0.sort.bam
samtools index primary-align-mod-dorado-sup-4-1-0.sort.bam 
modbam2bed -m 5mC -t 40 ragtag.scaffold.n-asiaeorientalis.fasta primary-align-mod-dorado-sup-4-1-0.sort.bam > mod-calls-scaffold-n-as.bed
```

## Notes et alia

- PSMC Li and Durbin 2011
- introgression - gene trees vs species trees
- Tajima's D or other (windowed) - is this possible with a single taxon?
- Runs of homozygosity (heterozygosity?)
- evolution Molecular versus polymorphism, methylation, synteny, diversity in sister taxa, polymorphism in other taxa, linkage in other taxa, positive selection, coding vs. non-coding
- GALBA
- Sniffles
- phase in Clair3

[jcvi](https://github.com/tanghaibao/jcvi)
[methylartist](https://github.com/adamewing/methylartist)
[clapper rail assembly](https://academic.oup.com/g3journal/article/13/8/jkad097/7148140)
[mushroom queen](https://academic.oup.com/g3journal/article/13/8/jkad102/7161644)
[hydractinia](https://academic.oup.com/g3journal/article/13/8/jkad107/7170730)
[marigold](https://www.biorxiv.org/content/10.1101/2023.07.25.550479v2.full)

### Blobtools

#### `blobtools` is unappealing

It is not clear that it is worth pursuing. The snail plots are not that informative and the contigd are all classified clearly as mammalian by `kraken2` save a small number:

`blobtoolkit v4.1.5`

Quick check with `kraken2`: `kraken2 --db /data/databases/k2_pluspf_20210517/ --output kraken2-raven-purged-polished.out --report krkane2-raven-purged-polished raven-purged-polished.fasta` - four contigs classified as fungi/bacteria. Blasting with each indicates they are mammalian.

BLAST: `export BLASTDB=/data/databases/nt/`
after download
`blastn -db nt -query raven-purged-polished.fasta -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 4 -out raven-purged-polished.blastn.out`

Stopped mid-run. Not sorth the effort after `kraken2`.`

```bash
blobtools create --fasta variants/raven-purged-polished.fasta raven-purged-polished
```

**used export LC_ALL=en_US.UTF-8 to change LANG was en_NZ.UTF-8**
see [here](https://forum.qiime2.org/t/unicodedecodeerror-utf-8-codec-cant-decode-byte-0x8b-in-position-1-invalid-start-byte/3385/6)

This **does not** resolve the error after invoking `blobtools add`:

```bash
UnicodeDecodeError: 'utf-8' codec can't decode byte 0x8b in position 1: invalid start byte
```

### References (cetacean assemblies)

- [Population genomics of finless porpoises reveal an incipient cetacean species adapted to freshwater](https://www.nature.com/articles/s41467-018-03722-x)  
Scaffold N50 values of 6.3 Mb  
- [Chromosome-Level Genome Assembly of the Rough-Toothed Dolphin (Steno bredanensis)](https://www.mdpi.com/2077-1312/11/2/418)  
N50 length of 105.53 Mb  
- [Draft genome sequencing and assembly of Risso's dolphin, Grampus griseus](https://www.jgenomics.com/v11p0009.htm)  
2.042 Mb of scaffold N50  
- [Chromosome-level Genome Provides Insights into Environmental Adaptability and Innate Immunity in the Common Dolphin (Delphinus delphis)](https://www.authorea.com/doi/full/10.22541/au.167146942.28834420)  
scaffold N50 of 108.93 Mb  
- High-quality chromosome-level genome assembly of the melon-headed whale (Peponocephala electra)  
contig N50 of 82.36 Mb
- [Gapless genome assembly of East Asian finless porpoise](https://www.nature.com/articles/s41597-022-01868-4)  
contig N50 of 84.69 Mb 
- [An Indo-Pacific Humpback Dolphin Genome Reveals Insights into Chromosome Evolution and the Demography of a Vulnerable Species](https://www.sciencedirect.com/science/article/pii/S2589004220308324)  
scaffold N50 of 27.7 Mb  

### Other genomes

_S. bridanensis_ (rough toothed)
  - scaffold N50 length of 105.53 Mb
  - BUSCO score of 90.6% and 97.3%
  - 98.17% of the assembled genome was anchored onto 22 chromosomes<br>
  - Not a great assembly
<img src="images/S-bridanensi-hi-c.webp" title="hi-c png" width="200"/><br>

- _D. delphis_ (common)
  - contig N50 of 63.85 Mb and a scaffold N50 of 108.93 Mb. Approximately 93.81% of contigs were anchored onto 22 chromosomes.<br>
  - this is still a terrible assembly
 <img src="images/D-delphinus-hi-c.png" title="hi-c png" width="200"/><br>

### Goldrush behaviour

Goldrush seems to have some odd behaviour besides the large number of small contigs. The contig lengths are very quantized. In the plots below, orange lines are draw at 1,000bp intervals starting at 1Kbp.
n.b. this is a result of rounding, which is part of the `goldrush` algorithm (see issues)
.
<img src="images/goldrush-hist.png" title="Odd distributions" width="500"/>
**NX plots zoomed**

Note: here is the Nextdenovo config file

```bash
***[General]
job_type = local # here we use SGE to manage jobs
job_prefix = nextDenovo
task = all
rewrite = yes
deltmp = yes
parallel_jobs = 2
input_type = raw
read_type = ont # clr, ont, hifi
input_fofn = input.fofn
workdir = dolphin-2k-assemble

[correct_option]
read_cutoff = 1k
genome_size = 3g # estimated genome size
sort_options = -m 50g -t 12
minimap2_options_raw = -t 12
pa_correction = 2
correction_options = -p 30

[assemble_option]
minimap2_options_cns = -t 12
nextgraph_options = -a 1
```

previous HiFiasm stats:
Stats on the 5Kbp filtered assembly:

```bash
file                               format  type  num_seqs        sum_len  min_len   avg_len  max_len      Q1      Q2      Q3  sum_gap     N50  Q20(%)  Q30(%)  GC(%)
dolphin.5k.q10.asm.bp.p_utg.fasta  FASTA   DNA    220,736  6,422,784,721    4,979  29,097.1  521,555  13,498  20,138  33,300        0  37,739       0       0  41.38
```

Comparative Genomics
Orthofinder
Orthology
