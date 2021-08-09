# reset current exclusive_signal_replicate results
if [ -d "data/exclusive_signal_replicate" ]; then
    rm -r "data/exclusive_signal_replicate"
    mkdir "data/exclusive_signal_replicate"
else
    mkdir "data/exclusive_signal_replicate"
fi

# exclusive gopeaks signal
echo "gopeaks"
for bed in $(find data/intervene/ -name "1000_gopeaks.bed"); do

    method="gopeaks"
    sample=$(echo $bed | cut -d/ -f3)
    bam=$(find data/ban -name "$sample*.bam")

    cmd="cut -f1-3 $bed | bedtools intersect -C -a stdin -b $bam > data/exclusive_signal_replicate/${method}_${sample}.bed &"
    echo $cmd
    eval $cmd

done
wait

# exclusive macs2 signal
echo "macs2"
for bed in $(find data/intervene/ -name "0100_macs2.bed"); do

    method="macs2"
    sample=$(echo $bed | cut -d/ -f3)
    bam=$(find data/ban -name "$sample*.bam")

    cmd="cut -f1-3 $bed | bedtools intersect -C -a stdin -b $bam > data/exclusive_signal_replicate/${method}_${sample}.bed &"
    echo $cmd
    eval $cmd

done
wait

# exclusive seacr-relaxed signal
echo "seacr-relaxed"
for bed in $(find data/intervene/ -name "0010_seacr-relaxed.bed"); do

    method="seacr-relaxed"
    sample=$(echo $bed | cut -d/ -f3)
    bam=$(find data/ban -name "$sample*.bam")

    cmd="cut -f1-3 $bed | bedtools intersect -C -a stdin -b $bam > data/exclusive_signal_replicate/${method}_${sample}.bed &"
    echo $cmd
    eval $cmd

done
wait

# exclusive seacr-stringent signal
echo "seacr-stringent"
for bed in $(find data/intervene/ -name "0001_seacr-stringent.bed"); do

    method="seacr-stringent"
    sample=$(echo $bed | cut -d/ -f3)
    bam=$(find data/ban -name "$sample*.bam")

    cmd="cut -f1-3 $bed | bedtools intersect -C -a stdin -b $bam > data/exclusive_signal_replicate/${method}_${sample}.bed &"
    echo $cmd
    eval $cmd

done
wait
