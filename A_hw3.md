# Homework 3: Pipelines

## Calculate Summaries of the Genome

1. The total number of nucleotides are 143726002.
2. The total number of N's are 1152978
3. The total number of of sequences is 1870

#### See complete output below
baseCount       143726002
nBaseCount      1152978
realBaseCount   142573024
upperBaseCount  142573024
lowerBaseCount  0
seqCount        1870
fileCount       1
meanSize        76858.8
SdSize  1382100.2
minSize 544
minSeqSize      211000022279089
maxSize 32079331
maxSeqSize      3R
medianSize      1577
nCountMean      616.6
nCountSd        6960.7
upperCountMean  76242.3
upperCountSd    1379508.4
lowerCountMean  0.0
lowerCountSd    0.0
fracMasked      0.00
fracRealMasked  0.00

#### The code below was used to verify the accuracy of the retrieved file 

``` md5sum -c <(grep "dmel-all-chromosome-r6.48.fasta.gz" md5sum.txt) ```

## Report Summarizing the Annotation

1. Total number features of each type, sorted from the most common to the least common
The code below was used to list, sort, and order features in the file from greatest to least.

``` bioawk -c gff '{print $feature}' dmel-all-r6.48.gtf.gz | sort | uniq -c | sort -nr -k1,1 ```

 190050 exon,
 163242 CDS,
  46802 5UTR,
  33738 3UTR,
  30885 start_codon,
  30825 stop_codon,
  30799 mRNA,
  17896 gene,
   3053 ncRNA,
    485 miRNA,
    365 pseudogene,
    312 tRNA,
    300 snoRNA,
    262 pre_miRNA,
    115 rRNA,
     32 snRNA,

2. Total number of genes per chromosome arm (X,Y, 2L, 2R, 3L, 3R, 4)

3515 2L      gene,
3653 2R      gene,
3489 3L      gene,
4227 3R      gene,
114 4        gene,
    
2708 X       gene,
113 Y        gene,

The code below was used to list a describe gene output only.

``` bioawk -c gff '{if($feature ~ /^gene$/)print $seqname"\t"$feature}' dmel-all-r6.48.gtf.gz | sort | uniq -c ```

      1 211000022278279 gene
      1 211000022278436 gene
      1 211000022278449 gene
      1 211000022278760 gene
      1 211000022279165 gene
      1 211000022279188 gene
      1 211000022279264 gene
      1 211000022279392 gene
      1 211000022279681 gene
      1 211000022280328 gene
      1 211000022280341 gene
      1 211000022280347 gene
      1 211000022280481 gene
      2 211000022280494 gene
      1 211000022280703 gene
   3515 2L      gene
   3653 2R      gene
   3489 3L      gene
   4227 3R      gene
    114 4       gene
     38 mitochondrion_genome    gene
     21 rDNA    gene
      2 Unmapped_Scaffold_8_D1580_D1567 gene
   2708 X       gene
    113 Y       gene
