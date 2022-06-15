cd REdiscoverTE_Lung/  #work folder
SALMON_INDEX_DIR=/mnt/z/REdiscoverTE_Lung/Step_1_salmon_idx   # SALMON INDEX
METADATA_FOR_ROLLUP=/mnt/z/REdiscoverTE_Lung/Step_3_metadata_for_rollup.tsv.txt #METADATA
PERCORSOFASTQ=/mnt/z/EGA_DATA/EGAD00001000223/ #fastq files
SALMON="Salmon-0.8.2_linux_x86_64/bin/salmon"
echo "sample\tquant_sf_path" > $METADATA_FOR_ROLLUP  #FILE METADATA
 ###################
for i in `cat filenames.txt` #filenames with all _1.rnaseq.fastq.gz filenames
do
j=`basename $i _1.rnaseq.fastq.gz`
SAMPLE=$j
PAIR1=$SAMPLE"_1.rnaseq.fastq.gz"
PAIR2=$SAMPLE"_2.rnaseq.fastq.gz"
OUTPUT_SALMON=Step_2_salmon_counts/$j
SALMON_QUANT_SF=$OUTPUT_SALMON/quant.sf
echo $j
$SALMON quant --seqBias --gcBias  -i $SALMON_INDEX_DIR  -l 'A' -1 $PERCORSOFASTQ$PAIR1 -2 $PERCORSOFASTQ$PAIR2  -o $OUTPUT_SALMON -p 8
echo "$j"\t"$SALMON_QUANT_SF" >> $METADATA_FOR_ROLLUP
done
Rscript  ./rollup.R -m Step_3_metadata_for_rollup.tsv.txt --datadir=rollup_annotation --threads=8 --outdir=Step_4_rollup --assembly="hg38"
