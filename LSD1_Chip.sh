#Tool used: deeptools
#LSD1 bigWig track where downloaded from GEO.

bigwigCompare -b1 GSM1619041_2_H526-GSK552_LSD1_i12.bigWig -b2 GSM1619045_8_H526-GSK552_Input_i8.bigWig -p max -o H526-GSK552_LSD1_ratio.bw --operation	ratio
bigwigCompare -b1 GSM1619040_1_H526-DMSO_LSD1_i8.bigWig -b2 GSM1619044_7_H526-DMSO_Input_i7.bigWig -p max -o H526-Ctrl_LSD1_ratio.bw --operation	ratio
BLACKLIST=hg19-blacklist.v2.bed
SAMPLELIST_LSD=$(ls *ratio.bw)
Sample_coll=TE_Down_hg19.bed #hg19 coordinates of downregulated ERVs
computeMatrix scale-regions -S  $SAMPLELIST_LSD -R $Sample_coll  -o LSD1_SCLC1.mtx -bs 5 -b 5000 -a 5000 -m 1000 --smartLabels -p max --missingDataAsZero --skipZeros
computeMatrix scale-regions -S  $SAMPLELIST_LSD -R $Sample_coll  -o LSD1_SCLC2.mtx -bs 5 -b 50 -a 50 -m 500 --smartLabels -p max --missingDataAsZero --skipZeros
plotHeatmap -m LSD1_SCLC1.mtx -out LSD1_SCLC1.pdf --perGroup --plotType se
plotHeatmap -m LSD1_SCLC2.mtx -out LSD1_SCLC2.pdf --perGroup --plotType se


SAMPLELIST_LSD=$(ls *ratio.bw)
Sample_coll=LTR_LSD1.bed #hg19 coordinates of downregulated ERVs
computeMatrix scale-regions -S  $SAMPLELIST_LSD -R $Sample_coll  -o LSD1_SCLC1.mtx -bs 5 -b 5000 -a 5000 -m 1000 --smartLabels -p max --missingDataAsZero --skipZeros
computeMatrix scale-regions -S  $SAMPLELIST_LSD -R $Sample_coll  -o LSD1_SCLC2.mtx -bs 5 -b 50 -a 50 -m 500 --smartLabels -p max --missingDataAsZero --skipZeros
plotHeatmap -m LSD1_SCLC1.mtx -out LSD1_SCLC1.pdf --perGroup --plotType se
plotHeatmap -m LSD1_SCLC2.mtx -out LSD1_SCLC2.pdf --perGroup --plotType se

H526-Ctrl_H3K4
SAMPLELIST_LSD=$(ls *Ctrl_H3K4*.bw)
Sample_coll=LTR_LSD1.bed #hg19 coordinates of downregulated ERVs
computeMatrix scale-regions -S  $SAMPLELIST_LSD -R $Sample_coll  -o LSD1_SCLC1.mtx -bs 5 -b 5000 -a 5000 -m 1000 --smartLabels -p max --missingDataAsZero --skipZeros
plotHeatmap -m LSD1_SCLC1.mtx -out H3K4_SCLC1.pdf --perGroup --plotType se
