#!/bin/sh
python /home/tdfyoder/GIT/AmpUMI/AmpUMI.py Process --fastq lima_output.bc1001--bc1001.bam.Q20.fastq --fastq_out bc1001.for.fastq --umi_regex ".*TCGTCGGCAGCGTCAGATGTGTCIIIIIIIITTGGGAAAGGGACCCATGTATTC" # --write_UMI_counts

python /home/tdfyoder/GIT/AmpUMI/AmpUMI.py Process --fastq lima_output.bc1001--bc1001.bam.Q20.fastq --fastq_out bc1001.rev.fastq --umi_regex ".*CTCGTGGGCTCGGAGATGTGAGIIIIIIIIGGCTCTGTTTGGGGCTATGT" # --write_UMI_counts

grep "^\@m" bc1001.for.fastq > bc1001.for.fastq.list

grep "^\@m" bc1001.rev.fastq > bc1001.rev.fastq.list

sort -k 1 bc1001.for.fastq.list > sorted.bc1001.for.fastq.list

sort -k 1 bc1001.rev.fastq.list > sorted.bc1001.rev.fastq.list

comm -12 sorted.bc1001.for.fastq.list sorted.bc1001.rev.fastq.list > bc1001.red.list #find common read names between forwrad and reverse fastq files

grep -vf bc1001.red.list sorted.bc1001.rev.fastq.list > bc1001.non.list #subtract common read names from reverse fastq to retain only unique reads

sed "s/^.//g" bc1001.non.list > bc1001.non.list2 # reformat using sed

bash /home/tdfyoder/GIT/bbmap/filterbyname.sh in=bc1001.rev.fastq out=bc1001.dedup.rev.fastq names=bc1001.non.list2 include=t qin=33 #filter fastq file for unique reverese reads

cat bc1001.for.fastq bc1001.dedup.rev.fastq > uni.bc1001.dedup.fastq #combine forward fastq with reverse fastq

cutadapt -b TCGTCGGCAGCGTCAGATGTGTC -b CTCGTGGGCTCGGAGATGTGAG -b TTGGGAAAGGGACCCATGTATTC -b GGCTCTGTTTGGGGCTATGT  -o uni.bc1001.dedup..trim.fastq  uni.bc1001.dedup.fastq
##------------------------------------###

#ngmlr -x pacbio -i 0.30 -r /mnt/c/Export/Fang/pacBIO/ref/TRAC.fasta -q uni.bc1001.dedup.trim.fastq  -o /mnt/c/Export/Fang/PacBIO/bc1001/map/bc1001.sam

#sniffles -m sort.bc1008.bam -v bc1008.vcf -s 2 -d 6000
