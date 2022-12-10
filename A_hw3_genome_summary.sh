#!/usr/bin/env bash

faSize -tab  dmel-all-chromosome-r6.48.fasta.gz

md5sum -c <(grep "dmel-all-chromosome-r6.48.fasta.gz" md5sum.txt)

srun -A class-ee282 --pty --x11 bash -i

conda activate ee282

faSize dmel-all-chromosome-r6.48.fasta.gz

wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/md5sum.txt

./hw3_genome_summary.sh

faSize -tab  dmel-all-chromosome-r6.48.fasta.gz

wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/md5sum.txt

less ./hw3_genome_summary.sh

faSize ./hw3_genome_summary.sh

faSize -tab  dmel-all-chromosome-r6.48.fasta.gz

md5sum -c md5sum.txt

md5sum.txt -c <(grep "dmel-all-chromosome-r6.48.fasta.gz" md5sum.txt)

md5sum -c <(grep "dmel-all-chromosome-r6.48.fasta.gz" md5sum.txt)

ls

less dmel-all-r6.48.gtf.gz

less -S dmel-all-r6.48.gtf.gz

bioawk -c help

bioawk -c gff '{print $feature}' | less

bioawk -c gff '{print $feature}' dmel-all-r6.48.gtf.gz | less

bioawk -c gff '{print $feature}' dmel-all-r6.48.gtf.gz | sort | uniq -c

bioawk -c gff '{print $feature}' dmel-all-r6.48.gtf.gz | sort | uniq -c | sort -nr -k1,1

bioawk -c gff '{print $seqname"\t"$feature}' dmel-all-r6.48.gtf.gz | less

bioawk -c gff '{if($feature ~ /^gene$/)print $seqname"\t"$feature}' dmel-all-r6.48.gtf.gz | less

bioawk -c gff '{if($feature ~ /^gene$/)print $seqname"\t"$feature}' dmel-all-r6.48.gtf.gz | sort | uniq -c

history > 10thNoov.txt
