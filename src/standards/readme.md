# GoPeaks-Compare Standards

These standards were used to assess peak calls from GoPeaks, SEACR, and MACS2, via Receiver Operating Characteristic (ROC) curves and Precision-Recall (PR) curves. 

Most files were downloaded through ENCODE with `wget` and preprocessed by filtering for peaks with logQVAL > 10 and peaks within 1kbp were merged. The following command was used:

```bash
$ zcat <standard.bed.gz> | awk '$9 > 10' | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 1000 > standard.1000.filtered.bed
```

Below is the key to each ENCODE sample.

| Standard Name  | Cell Line | Mark     | Assembly | Notes             |
| -------------- | --------- | -------- | -------- | ----------------- |
| ENCFF590NGQ    | K562      | H3K4me1  | GRCh38   | replicated peaks  |
| ENCFF256AQN    | K562      | H3K4me2  | GRCh38   | replicated peaks  |
| ENCFF246IEW    | K562      | H3K4me3  | GRCh38   | replicated peaks  |
| ENCFF403DTU    | K562      | H3K27Ac  | GRCh38   | replicated peaks  |
| ENCFF795ZOS    | K562      | H3K27me3 | GRCh38   | replicated peaks* |
| ENCFF660GHM    | K562      | CTCF     | GRCh38   | replicated peaks  |
| Kasumi_H3K27Ac | Kasumi    | H3K27Ac  | GRCh38   | More info below*  |

All the K562 marks come from the epigenome ENCSR612NLL. 

## Note

\* Kasumi_H3K27Ac comes from SRX4143063/SRX4143067, and the bed files are download through [Chip Atlas](https://chip-atlas.org/). Peaks that appear in both replicates were kept, and features were merged if they were within 1kbp. Original publication was titled [TAF1 plays a critical role in AML1-ETO driven leukemogenesis](https://www.nature.com/articles/s41467-019-12735-z) by Xu *et al.* 

\* K562 H3K27me3 broad peak used a merge distance of 3000 bp.


