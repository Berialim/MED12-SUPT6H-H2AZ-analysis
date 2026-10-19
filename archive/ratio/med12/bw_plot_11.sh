# chipseq plot with selected genes
# source ~/.bashrc
echo "parameter1 = selected.bed out_dir regtion_length filter_value bwfiles(a.bw b.bw c.bw)"
if [ ! -n "$5" ];then
    echo "need parameter"
    exit 8
fi
mkdir $2
a=${@:5}
let body=300
computeMatrix scale-regions -b 50 -a $3 --missingDataAsZero --regionBodyLength $body -R $1 -S $a -o $2/scal.gz -p 24 --outFileSortedRegions $2/scale.bed &
computeMatrix reference-point --referencePoint TES --missingDataAsZero -b 50  $3 -R $1 -S $a  -o $2/tes.gz -p 24 --outFileSortedRegions $2/tes.bed &
computeMatrix reference-point --referencePoint TSS --missingDataAsZero  -b 50 -a $3 -R $1 -S $a  -o $2/tss.gz -p 24 --outFileSortedRegions $2/tss.bed &
computeMatrix reference-point --referencePoint center --missingDataAsZero  -b 50 -a $3 -R $1 -S $a  -o $2/center.gz -p 24 --outFileSortedRegions $2/center.bed &
wait
cd $2
for i in *.gz
do
    computeMatrixOperations filterValues -m $i -o ${i%.*}_filterd.gz --max $4 --min -${4} &
done
wait
for i in *.gz
do
    plotHeatmap -m $i -out ${i%.*}.pdf --colorList "white, #960E0E" "white, #960E0E" --zMin 0 --sortUsingSamples 1 --whatToShow "heatmap and colorbar"  --plotFileFormat pdf  --heatmapHeight 22 --heatmapWidth 4 &
    plotProfile -m $i -out ${i%.*}_profile.pdf  --perGroup --plotFileFormat pdf & # --colors '#9DD0C7' '#9180AC' '#E58579' '#8AB1D2' &
    plotHeatmap --colorList "#131895, white, #960E0E" "#131895, white, #960E0E"  --whatToShow "heatmap and colorbar" --plotFileFormat pdf --zMax 1 --zMin -1 --sortUsingSamples 1 --heatmapHeight 22 --heatmapWidth 4 -m $i -out br_${i%.*}.pdf &
    plotHeatmap -m $i  --sortUsing region_length  -out ${i%.*}_len.pdf --colorList "white, #960E0E" "white, #960E0E" --zMin 0 --sortUsingSamples 1 --whatToShow "heatmap and colorbar"  --plotFileFormat pdf  --heatmapHeight 22 --heatmapWidth 4 &
    plotHeatmap --sortUsing region_length  --colorList "#131895, white, #960E0E" "#131895, white, #960E0E"  --whatToShow "heatmap and colorbar" --plotFileFormat pdf --zMax 1 --zMin -1 --sortUsingSamples 1 --heatmapHeight 22 --heatmapWidth 4 -m $i -out br_${i%.*}_len.pdf &
    plotHeatmap --sortRegions no  --colorList "#131895, white, #960E0E" "#131895, white, #960E0E"  --whatToShow "heatmap and colorbar" --plotFileFormat pdf --zMax 1 --zMin -1 --heatmapHeight 22 --heatmapWidth 4 -m $i -out br_${i%.*}_nosort.pdf &
    plotHeatmap -m $i  --sortRegions no  -out ${i%.*}_nosort.pdf --colorList "white, #960E0E" "white, #960E0E" --zMin 0  --whatToShow "heatmap and colorbar" --plotFileFormat pdf  --heatmapHeight 22 --heatmapWidth 4 &
done
wait
