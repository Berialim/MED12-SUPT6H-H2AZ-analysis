##chipseq_se
source ~/.bashrc
mkdir aligned
mkdir trimmed
mkdir fastqc1
#for i in *.fq.gz
#do 
#    (fastqc ${i} -o fastqc)& 
#done
#multiqc fastqc -o preqc
##trim adapter
# move all fastq in fq dictionary
echo "usage:"
echo "bash chipseq_se.sh quality(20/30) fragment_length bed_file(for plot)"
if [ ! -n "$3" ];then
        echo "not enough parameter"
        exit 8
fi
species="hg38"
bed_file=$3
cd fq
source activate qc
for i in *.fq.gz
do
    mv $i ${i##*L0*_}
done
mkdir useless
mv ambiguous.fq.gz useless
mv undecoded.fq.gz useless
x=`ls *.fq.gz | sed 's/.fq.gz//g'`
cd ..
##trimadapter
mkdir trimmed
for line in $x
do
    trim_galore -q $1 fq/${line}.fq.gz   -o trimmed &
done
wait
cp -r trimmed keep
cd trimmed
for i in *.fq.gz
do
    mv ${i} ${i%_*}.fq.gz
done
cd ..

for i in trimmed/*.fq.gz
do 
    fastqc ${i} -o fastqc1& 
done
wait
multiqc fastqc1 -o qc

source ~/.bashrc
cd trimmed
for i in *.fq.gz
do
    echo "sample: ${i%%.*}" >> ../qc/align_report
    bowtie2 -p 16 -x ~/reference/$species/bowtie2/genome -U ${i}  -S ../aligned/${i%%.*}.sam 2>> ../qc/align_report
    samtools sort -@ 16 -o ../aligned/${i%%.*}.bam  ../aligned/${i%%.*}.sam
    picard MarkDuplicates I=../aligned/${i%%.*}.bam  O=../aligned/${i%%.*}.markdup.bam  M=../aligned/${i%%.*}.markdup.txt
    rm ../aligned/${i%%.*}.sam
    rm ../aligned/${i%%.*}.bam
    mv ../aligned/${i%%.*}.markdup.bam ../aligned/${i%%.*}.bam
    samtools index -@ 16   ../aligned/${i%%.*}.bam
done
cd ../qc
python ~/python/align_plot/alignment_se.py

# generatre bigwig file for read coverage in genome(visualize in igv)
mkdir ../bw
cd ../aligned
for i in *.bam
do
    bamCoverage --bam ${i} -o ../bw/${i%.*}.bw --binSize 10 --normalizeUsing RPGC \
     --effectiveGenomeSize 2913022398 -e $2 -p 16
done


cd ..
mv *.csv bw
cd bw
mv database.csv ../
csv=`ls *.csv|sed "s/.csv//g"`
for i in $csv
do
    mkdir $i 
    mv -rf $i.csv $i 
    cat $i.csv | while read line
    do
        name=`echo $line | sed "s/,/_/g"`
        cp ${line%,*}.bw $i/$name.bw
    done
    cd $i 
    bash ~/shell/chip/bw_plot.sh $bed_file plot 5000 *.bw 
    cd ..
done

# change sample name to treatment
cp ../database.csv ./
