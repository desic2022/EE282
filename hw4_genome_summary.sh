#!/usr/bin/env bash

#how to connect work to github
#make sure your fisrt in the right github directory

sed 's/|//g' hw4.md
git add .
git commit -m "Uploading hw4"
git push

sed 's/|//g' hw4_genome_summary.sh
git add .
git commit -m "Uploading hw4_genome_summary.sh"
git push

srun -A class-ee282 --pty --x11 bash -i

cd /pub/desic/classrepos/EE282

conda activate ee282

wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.48.fasta.gz

less dmel-all-chromosome-r6.48.fasta.gz

bioawk -c fastx ' length($seq) <= 100000 { print $name "\t" length($seq) } ' dmel-all-chromosome-r6.48.fasta.gz | less

bioawk -c fastx ' length($seq) <= 100000 { print $name "\t" length($seq) } ' dmel-all-chromosome-r6.48.fasta.gz | head |
column -t

bioawk -c fastx \
 ' length($seq) <= 100000 { print ">" $name "\n" $seq } ' \
 dmel-all-chromosome-r6.48.fasta.gz \
| head \
| less

bioawk -c fastx \
 ' length($seq) <= 100000 { print ">" $name "\n" $seq } ' \
 dmel-all-chromosome-r6.48.fasta.gz \
| less -S

bioawk -c fastx \
 ' length($seq) <= 100000 { print ">" $name "\n" $seq } ' \
 dmel-all-chromosome-r6.48.fasta.gz \
>dmelr6.lte.fa

bioawk -c fastx \
 ' length($seq) > 100000 { print ">" $name "\n" $seq } ' \
 dmel-all-chromosome-r6.48.fasta.gz \
>dmelr6.gt.fa

faSize -tab  dmel-all-chromosome-r6.48.fasta.gz

faSize dmelr6.lte.fa

faSize dmelr6.gt.fa

bioawk -c fastx \
  ' { print length($seq) "\t" gc($seq) } ' \
  dmelr6.lte.fa \
| sort -k1,1rn \
| head

bioawk -c fastx \
  ' { print length($seq) "\t" gc($seq) } ' \
  dmelr6.lte.fa \
| sort -k1,1rn \
> dmelr6.lte.lengc.txt

conda install ImageMagick

#1st plot LESS THAN: this is the cumulative distribution plot
plotCDF <(cut -f 1 dmelr6.lte.lengc.txt) dmelr6.lte.png

#this will display all png files
display *.png

#another way to do the same thing
plotCDF <(cut -f 1 dmelr6.lte.lengc.txt) /dev/stdout | display

#2nd plot GREATER THAN
bioawk -c fastx \
  ' { print length($seq) "\t" gc($seq) } ' \
  dmelr6.gt.fa \
| sort -k1,1rn \
> dmelr6.gt.lengc.txt

plotCDF <(cut -f 1 dmelr6.gt.lengc.txt) dmelr6.gt.png

display *.png

plotCDF <(cut -f 1 dmelr6.gt.lengc.txt) /dev/stdout | display

#3rd plot ALL CHROMOSOMES
bioawk -c fastx \
  ' { print length($seq) "\t" gc($seq) } ' \
  dmel-all-chromosome-r6.48.fasta.gz \
| sort -k1,1rn \
> dmelr6.all.lengc.txt

plotCDF <(cut -f 1 dmelr6.all.lengc.txt) /dev/stdout | display

#to make histograms load R in command line just type in the letter R and run
#below makes less than histograms
LT_data <- read.table("dmelr6.lte.lengc.txt",sep='\t')

head(LT_data)

colnames(LT_data) <- c("Length", "GC")

head(LT_data)

#first need to log transform column 1
LT_data[,1] <- log10(LT_data[,1])

head(LT_data)

LT_data[,2] <- log10(LT_data[,2])

head(LT_data)

library(ggplot2)

hist(LT_data$Length)

#save this for less than histogram
ggplot(LT_data, aes(x=GC)) + geom_histogram(binwidth=0.01) + xlab("GC log transformed")
ggsave("dmel6_LT_GC_plotr.png")

ggplot(LT_data, aes(x=Length)) + geom_histogram(binwidth=0.01) + xlab("Sequence Length log transformed")
ggsave("dmel6_LT_length_plotr.png")

GT_data <- read.table("dmelr6.gt.lengc.txt")
                     
GT_data <- read.table("dmelr6.gt.lengc.txt",sep='\t')

head(GT_data)
                   
GT_data[,2] <- log10(GT_data[,2])

GT_data[,1] <- log10(GT_data[,1])

> head(GT_data)

> colnames(GT_data) <- c("Length", "GC")
> head(GT_data)                    

#save this for greater than histogram
ggplot(GT_data, aes(x=GC)) + geom_histogram(binwidth=0.01) + xlab("GC log transformed")
ggsave("dmel6_GT_GC_plotr.png")

> ggplot(GT_data, aes(x=Length)) + geom_histogram(binwidth=0.01) + xlab("Sequence Length log transformed")
> ggsave("dmel6_GT_length_plotr.png")


gawk \
 ' { t = t + $1; print $1 "\t" t } ' \
  dmelr6.all.lengc.txt \
| less

faSize -tab  dmel-all-chromosome-r6.48.fasta.gz

#you can use this to find the N50, use the graph to find the halfway point
#the N50 is 25286936
less dmelr6.all.lengc.txt

#length, cumulative length, proportion: 143726002 is the total number of bases
gawk \
 ' { t = t + $1; print $1 "\t" t "\t" t / 143726002 } ' \
  dmelr6.all.lengc.txt | column -t | less

#this gives cumulative total which is 137547960
gawk ' { tot=tot+$1; print $1 "\t" tot } END { print tot }  ' /dfs6/pub/jje/ee282/r6.gt100kb.len_gc.txt

#assembly assignment
#this is now sorted in descending order with the cumulative total at the top
#rn means reverse numeric so largest to smallest

gawk ' { tot=tot+$1; print $1 "\t" tot } END { print tot } ' \
 /dfs6/pub/jje/ee282/r6.gt100kb.len_gc.txt \
| sort -k1,1rn \
| gawk ' NR == 1 { tot = $1 } NR > 1 { print $0 "\t" $2 / tot } ' | column -t

#this give everything above 50%
gawk ' { tot=tot+$1; print $1 "\t" tot } END { print tot } ' \
 /dfs6/pub/jje/ee282/r6.gt100kb.len_gc.txt \
| sort -k1,1rn \
| gawk ' NR == 1 { tot = $1 } NR > 1 && $2/tot >= 0.5 { print $1 }{ print $0 "\t" $2 / tot } ' | column -t

#this gives the final answer for the N50
gawk ' { tot=tot+$1; print $1 "\t" tot } END { print tot } ' \
 /dfs6/pub/jje/ee282/r6.gt100kb.len_gc.txt \
| sort -k1,1rn \
| gawk ' NR == 1 { tot = $1 } NR > 1 && $2/tot >= 0.5 { print $1 } ' | head -1

#the N50 is 25286936

faSize -detailed data/raw/dmel-all-chromosome-r6.48.fasta.gz \
| sort -k2, 2rn \
gawk ' { tot=tot+$1; print $1 "\t" tot } END { print tot } ' \
 /dfs6/pub/jje/ee282/r6.lte100kb.len_gc.txt \
| sort -k1,1rn \
| gawk ' NR == 1 { tot = $1 } NR > 1 && $2/tot >= 0.5 { print $1 } ' | head -1

#from class
faSize -detailed "/dfs6/pub/desic/dmel-all-chromosome-r6.48.fasta.gz" \
| sort -k2,2rn \
| gawk ' { t=t+2; print $2 "\t" t} END { print t } ' \
| sort -k1,1rn \
| gawk '
  NR == 1 { t = $1 }
  NR > 1 { print NR-1 "\t" $0 "\t" $2 / t }
  ' \
| less
 
conda install -y miniasm minimap

#Genome assembly
# Assemble a genome from MinION reads

cp /pub/jje/ee282/iso1_onp_a2_1kb.fastq /pub/desic/classrepos/EE282

basedir=~/
projname=nanopore_assembly
basedir=$basedir/$projname
raw=$basedir/$projname/data/raw
processed=$basedir/$projname/data/processed
figures=$basedir/$projname/output/figures
reports=$basedir/$projname/output/reports

createProject $projname $basedir

minimap2 -x ava-ont -t16 iso1_onp_a2_1kb.fastq iso1_onp_a2_1kb.fastq | gzip -1 > onp.paf.gz

miniasm -f iso1_onp_a2_1kb.fastq onp.paf.gz > assembly.gfa

n50 () {
  bioawk -c fastx ' { print length($seq); n=n+length($seq); } END { print n; } ' $1 \
  | sort -rn \
  | gawk ' NR == 1 { n = $1 }; NR > 1 { ni = $1 + ni; } ni/n > 0.5 { print $1; exit; } '
}

awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' assembly.gfa \
| tee >(n50 /dev/stdin > n50.txt) \
> unitigs.fa

#BUSCO Assembly Assessment
conda install -c bioconda busco=5.4.4

#looks at the quality of the assembly
busco -m genome -i unitigs.fa -o ONP_Dmel -l diptera_odb10 --cpu 4

#2. Compare your assembly to both the contig assembly and the scaffold assembly
#from the Drosophila melanogaster on FlyBase using a contiguity plot.
# run in ee282 environment
faSplitbyN dmel-all-chromosome-r6.48.fasta dmel_contigs.fa 10

perl ~/bin/scaffold_to_contigs.pl -inputfile dmel-all-chromosome-r6.48.fasta -outputfile dmel_contigs.fa -minLength 0
1871 scaffolds in total
2442 contigs in total
2442 contigs over 0

bioawk -c fastx '{ print $name "\t" length($seq) }' <dmel-all-chromosome-r6.48.fasta> \
| sort -k2,2rn | bioawk ' BEGIN { print "Name\tLength\tAssembly" }{print $0 "\tFB_Scaff"}' \
> All.chromosomes.sizes.txt

bioawk -c fastx '{ print $name "\t" length($seq) }' dmel_contigs.fa \
| sort -k2,2rn \
| bioawk ' BEGIN { print "Name\tLength\tAssembly" }{print $0 "\tFB_Contig"}' \
> conitgs.sizes.txt

bioawk -c fastx '{ print $name "\t" length($seq) }' unitigs.fa \
| sort -k2,2rn \
| bioawk ' BEGIN { print "Name\tLength\tAssembly" }{print $0 "\tFB_Contig"}' \
> unitigs.sizes.txt

plotCDF2 *.sizes.txt plots/assemblies.png
plotCDF2 <(conitgs.sizes.txt) contigs.png
plotCDF2 tmp/{dmel-all-chromosome-r6.48.fasta,dmel_contigs.fa,truseq}_fifo /dev/stdout \
> tee r6_v_truseq.png \

#Busco scores
busco -c 16 -i dmel-all-chromosome-r6.48.fasta -l diptera_odb10 -o nanopore -m genome -c 16