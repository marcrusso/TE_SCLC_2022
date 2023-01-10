###ANALYSIS WITH TElocal Tool
cd /mnt/z/TE_revision
TElocal_exe="TElocal/TElocal-master/TElocal"
PATH_origin=""
TElocal_GTF="TElocal/gencode.v26.basic.annotation.gtf.ind"
TElocal_TE="TElocal/hg38_rmsk_TE.gtf.locInd"
cd $PATH_origin
for i in $(ls *_1.rnaseq.fastq.gz)
do
cd $PATH_origin
SAMPLE=$(basename $i "_1.rnaseq.fastq.gz")
EXT=".rnaseq.fastq.gz"
FOR="_1"
REV="_2"
STAR_INDEX="Genome_data/GRCh38_STARindex"
PATH_WORK="TMP/"
cp $SAMPLE$FOR$EXT $PATH_WORK$SAMPLE$FOR$EXT
cp $SAMPLE$REV$EXT $PATH_WORK$SAMPLE$REV$EXT
cd $PATH_WORK
echo $SAMPLE.trimming
#trimmomatic
trimmomatic PE -threads 8  $SAMPLE$FOR$EXT $SAMPLE$REV$EXT $SAMPLE$FOR.P.fq $SAMPLE$FOR.U.fq $SAMPLE$REV.P.fq $SAMPLE$REV.U.fq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20
###STAR
echo $SAMPLE.STAR
STAR  --readFilesIn $SAMPLE$FOR.P.fq $SAMPLE$REV.P.fq --alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8  --alignSoftClipAtReferenceEnds Yes --genomeDir $STAR_INDEX --genomeLoad NoSharedMemory --limitSjdbInsertNsj 1200000 --outFileNamePrefix $SAMPLE --outFilterIntronMotifs None  --outFilterMatchNminOverLread 0.33  --outFilterMismatchNmax 999  --outFilterMismatchNoverLmax 0.1 --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100  --outFilterScoreMinOverLread 0.33 --outFilterType BySJout  --outSAMtype BAM Unsorted --runThreadN 8 --twopassMode Basic
####SAMTOOLS
echo $SAMPLE.samtools
samtools index -@ 8 $SAMPLE"Aligned.out.bam"
####HTSEQ
echo $SAMPLE.TElocal
$TElocal_exe -b $SAMPLE"Aligned.out.bam" --GTF $TElocal_GTF --TE $TElocal_TE --project $SAMPLE.TElocal
mv $SAMPLE.TElocal* /mnt/z/TE_revision/TElocal_results/
rm -r *
done


####ANALYsis with homer
findMotifsGenome.pl LTR30.peaks hg38 LTR30 -size given -p 4 -cache 5000
findMotifsGenome.pl LTR22C.peaks hg38 LTR22C -size given -p 4 -cache 5000
findMotifsGenome.pl MER61F.peaks hg38 MER61F -size given -p 4 -cache 5000
