# Homework 4: Pipelines, plotting, genome assembly

## Calculate the following for all sequences less than or equal to 100kb:


1. The total number of nucleotides are 6178042.
2. The total number of N's are 662593.
3. The total number of of sequences is 1863.

## Calculate the following for all sequences greater than 100kb:

1. The total number of nucleotides are 137547960.
2. The total number of N's are 490385.
3. The total number of of sequences is 7.

## Plots of the following for all sequences less than or equal to 100kb

![less than 100kb GC plot](C:\Users\desichase79\Desktop\Speeches and Interviews SP1\dmel6_LT_GC_plotr.png)

![less than 100kb length plot](C:\Users\desichase79\Desktop\Speeches and Interviews SP1\dmel6_LT_length_plotr.png)

![Cumulative plot](C:\Users\desichase79\Desktop\Speeches and Interviews SP1\dmelr6.lte.png)

## Plots of the following for all sequences greater than 100kb

![greater than 100kb GC plot](C:\Users\desichase79\Desktop\Speeches and Interviews SP1\dmel6_GT_GC_plotr.png)

![greater than 100kb length plot](C:\Users\desichase79\Desktop\Speeches and Interviews SP1\dmel6_GT_length_plotr.png)

![Cumulative plot](C:\Users\desichase79\Desktop\Speeches and Interviews SP1\dmelr6.gt.png)

# Genome Assembly

## Assemble a genome from MinION reads

cp /pub/jje/ee282/iso1_onp_a2_1kb.fastq /pub/desic/classrepos/EE282

basedir=~/
projname=nanopore_assembly
basedir=$basedir/$projname
raw=$basedir/$projname/data/raw
processed=$basedir/$projname/data/processed
figures=$basedir/$projname/output/figures
reports=$basedir/$projname/output/reports

createProject $projname $basedir

minimap -t 32 -Sw5 -L100 -m0 $raw/reads.fq{,} \
| gzip -1 \
> $processed/onp.paf.gz

miniasm -f $raw/reads.fq $processed/onp.paf.gz \
> $processed/reads.gfa

awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' $processed/reads.gfa \
| tee >(n50 /dev/stdin > $reports/n50.txt) \
| fold -w 60 \
> $processed/unitigs.fa


## Assembly Assessment

awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' $processed/reads.gfa \
| tee >(n50 /dev/stdin > $reports/n50.txt) \
| fold -w 60 \
> $processed/unitigs.fa

1.This gives the final answer for the N50
``` gawk ' { tot=tot+$1; print $1 "\t" tot } END { print tot } ' \
 /dfs6/pub/jje/ee282/r6.gt100kb.len_gc.txt \
| sort -k1,1rn \
| gawk ' NR == 1 { tot = $1 } NR > 1 && $2/tot >= 0.5 { print $1 } ' | head -1 ```

The N50 in the code above is 25,286,936 which is exactly the same as indicated in the link provided in the homework.

2.Assemblies compared in the code below

perl ~/bin/scaffold_to_contigs.pl -inputfile dmel-all-chromosome-r6.48.fasta -outputfile dmel_contigs.fa -minLength 0
1871 scaffolds in total
2442 contigs in total
2442 contigs over 0

``` bioawk -c fastx '{ print $name "\t" length($seq) }' dmel_contigs.fa \
| sort -k2,2rn \
| bioawk ' BEGIN { print "Name\tLength\tAssembly" }{print $0 "\tFB_Contig"}' \
> conitgs.sizes.txt ```

``` bioawk -c fastx '{ print $name "\t" length($seq) }' unitigs.fa \
| sort -k2,2rn \
| bioawk ' BEGIN { print "Name\tLength\tAssembly" }{print $0 "\tFB_Contig"}' \
> unitigs.sizes.txt ```

``` bioawk -c fastx '{ print $name "\t" length($seq) }' <dmel-all-chromosome-r6.48.fasta> \
| sort -k2,2rn | bioawk ' BEGIN { print "Name\tLength\tAssembly" }{print $0 "\tFB_Scaff"}' \
> All.chromosomes.sizes.txt ```

```  plotCDF2 *.sizes.txt plots/assemblies.png ```
``` plotCDF2 <(conitgs.sizes.txt) contigs.png ```
``` plotCDF2 tmp/{dmel-all-chromosome-r6.48.fasta,dmel_contigs.fa,truseq}_fifo /dev/stdout \
> tee r6_v_truseq.png \ ```

3. Busco scores can be found with the code below

``` busco -c 16 -i dmel-all-chromosome-r6.48.fasta -l diptera_odb10 -o nanopore -m genome -c 16 ```