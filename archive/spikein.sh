ctrl=dm6
cd trimmed
mkdir ../temp_spikein
for i in *.fq.gz
do    
    bowtie2 -p 82 -x ~/reference/$ctrl/bowtie2/genome -U ${i}  -S ../temp_spikein/${i%%.*}.sam 2>> ../qc/align_report_spikein
    samtools sort -@ 82 -o ../temp_spikein/${i%%.*}.bam  ../temp_spikein/${i%%.*}.sam
    samtools index  -@ 82 ../temp_spikein/${i%%.*}.bam
done