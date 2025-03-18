thread=$1
mem=$2

fastp -i raw/070424-LHZ-PBMC_R1.fq.gz -I raw/070424-LHZ-PBMC_R2.fq.gz \
      -o raw/PBMC_1.fastq -O raw/PBMC_2.fastq \
      -w $thread

mysamtools="/public/home/liunangroup/liangyan/software/miniconda3/bin/samtools"
myjava="/public/home/liunangroup/liangyan/software/miniconda3/envs/wgs/bin/java"
mybwa="/public/home//liunangroup/liangyan/software/bwa-0.7.17/bwa"
mypicard="/public/home//liunangroup/liangyan/software/miniconda3/envs/wgs/share/picard-2.18.10-0/picard.jar"
mygatk="/public/home/liunangroup/liangyan/software/miniconda3/envs/wgs/opt/gatk-3.8/GenomeAnalysisTK.jar"
mytabix="/public/home/liunangroup/liangyan/software/miniconda3/bin/tabix"
mybgzip="/public/home/liunangroup/liangyan/software/miniconda3/bin/bgzip"

$mysamtools --version
$mytabix --version
$mybgzip --version
$myjava -version
md5sum $mypicard
md5sum $mygatk

myrefdir="/public/home/liunangroup/liangyan/project/zzZA/20240801_WGS_liunan/ref"
myfasta=$myrefdir/Homo_sapiens_assembly38.fasta
mydict=$myrefdir/Homo_sapiens_assembly38.dict
mysnp=$myrefdir/Homo_sapiens_assembly38.dbsnp138.vcf.gz
myindel=$myrefdir/Homo_sapiens_assembly38.known_indels.vcf.gz

md5sum $myfasta
md5sum $mydict
md5sum $mysnp
md5sum $myindel

mytmpjava="tmp/GATK/"
mkdir -p $mytmpjava
[ -w $mytmpjava ] && echo "writeable" || echo "write permission denied"

fastqroot="raw"
mkdir $fastqroot
gatkroot="GATK"
mkdir $gatkroot
samplelabel="PBMC"
mkdir $gatkroot/$samplelabel

myR1=$fastqroot/$samplelabel"_1.fastq"
myR2=$fastqroot/$samplelabel"_2.fastq"
mysamplebase=$gatkroot/$samplelabel/$samplelabel

myrgID="PBMC"
myrgSM="PBMC"
myrgLB="PBMC"
myrgPL="illumina"
myrgPU="PBMC"

# time $bwa index $myfasta
time $mybwa mem \
     -aM \
     -t $thread \
     $myfasta \
     $myR1 $myR2 \
     > $mysamplebase.sam 2> $mysamplebase-bwasam.err

time $myjava -Xmx$mem -Djava.io.tmpdir=$mytmpjava -jar \
     $mypicard AddOrReplaceReadGroups \
     I=$mysamplebase.sam \
     O=$mysamplebase-sort-rg.bam \
     SORT_ORDER=coordinate \
     CREATE_INDEX=true \
     RGID=$myrgID RGSM=$myrgSM RGLB=$myrgLB RGPL=$myrgPL RGPU=$myrgPU \
     > $mysamplebase-rgsort.out 2> $mysamplebase-rgsort.err

time $myjava -Xmx$mem -Djava.io.tmpdir=$mytmpjava -jar \
     $mypicard MarkDuplicates \
     I=$mysamplebase-sort-rg.bam \
     O=$mysamplebase-sort-rg-dedup.bam \
     M=$mysamplebase-dedup-metric.txt \
     ASSUME_SORT_ORDER=coordinate \
     CREATE_INDEX=true \
     > $mysamplebase-dedup.out 2> $mysamplebase-dedup.err

time $myjava -Xmx$mem -Djava.io.tmpdir=$mytmpjava -jar \
     $mygatk -T RealignerTargetCreator \
     -R $myfasta \
     -I $mysamplebase-sort-rg-dedup.bam \
     -known $myindel \
     -o $mysamplebase-realign-targets.intervals \
     > $mysamplebase-realigntarget.out 2> $mysamplebase-realigntarget.err

time $myjava -Xmx$mem -Djava.io.tmpdir=$mytmpjava -jar \
     $mygatk -T IndelRealigner \
     -R $myfasta \
     -known $myindel \
     -targetIntervals $mysamplebase-realign-targets.intervals \
     -I $mysamplebase-sort-rg-dedup.bam \
     -o $mysamplebase-sort-rg-dedup-realigned.bam \
     --filter_bases_not_stored \
     > $mysamplebase-realignbam.out 2> $mysamplebase-realignbam.err

time $myjava -Xmx$mem -Djava.io.tmpdir=$mytmpjava -jar \
     $mypicard BuildBamIndex \
     I=$mysamplebase-sort-rg-dedup-realigned.bam \
     > $mysamplebase-realignbamindex.out 2> $mysamplebase-realignbamindex.err

time $myjava -Xmx$mem -Djava.io.tmpdir=$mytmpjava -jar \
     $mygatk -T BaseRecalibrator \
     -R $myfasta \
     -knownSites $mysnp \
     -knownSites $myindel \
     -I $mysamplebase-sort-rg-dedup-realigned.bam \
     -o $mysamplebase-recal-table.txt \
     > $mysamplebase-recalibrate.out 2> $mysamplebase-recalibrate.err

time $myjava -Xmx$mem -Djava.io.tmpdir=$mytmpjava -jar \
     $mygatk -T PrintReads \
     -R $myfasta \
     -I $mysamplebase-sort-rg-dedup-realigned.bam \
     -BQSR $mysamplebase-recal-table.txt \
     -o $mysamplebase-GATK.bam \
     > $mysamplebase-recalbam.out 2> $mysamplebase-recalbam.err

time $myjava -Xmx$mem -Djava.io.tmpdir=$mytmpjava -jar \
     $mypicard BuildBamIndex \
     I=$mysamplebase-GATK.bam \
     > $mysamplebase-GATKbamindex.out 2> $mysamplebase-GATKbamindex.err

time $myjava -Xmx$mem -Djava.io.tmpdir=$mytmpjava -jar \
     $mygatk -T HaplotypeCaller \
     -R $myfasta \
     -I $mysamplebase-GATK.bam \
     --dbsnp $mysnp \
     -o $mysamplebase-GATK-HC.vcf \
     -nct $thread