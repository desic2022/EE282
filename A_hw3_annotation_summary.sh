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
