# -----------------------methylation analysis-------------------------------
cd /public/home/liunangroup/liangyan/project/zzZA/20240731_liunan-DNAmethylation/Release/Mapping

sbatch -J job -p cu -N 1 -n 1 -c 1 -o %j.out -e %j.err --wrap="python DSS-input-generate.py control1.Original_CX_report.txt.gz control1.txt"
sbatch -J job -p cu -N 1 -n 1 -c 1 -o %j.out -e %j.err --wrap="python DSS-input-generate.py patient1.Original_CX_report.txt.gz patient1.txt"

# bismark extract methylation
thread="30"
mem=$(echo "$thread * 2"|bc)G

sbatch -J job -p cu -N 1 -n 1 -o %j.out -e %j.err --cpus-per-task=${thread} --mem=${mem} --wrap="
# sambamba sort -m $mem -t $thread -n control1_P_R1_bismark_bt2_pe.bam.sort.bam -o control1.sorted.bam
~/software/Bismark-0.24.2/bismark_methylation_extractor -p --bedGraph --parallel 10 control1.sorted.bam --output_dir methylation_extract
"

sbatch -J job -p cu -N 1 -n 1 -o %j.out -e %j.err --cpus-per-task=${thread} --mem=${mem} --wrap="
# sambamba sort -m $mem -t $thread -n patient1_P_R1_bismark_bt2_pe.bam.sort.bam -o patient1.sorted.bam
~/software/Bismark-0.24.2/bismark_methylation_extractor -p --bedGraph --parallel 10 patient1.sorted.bam --output_dir methylation_extract
"

# differential methylation
cd /public/home/liunangroup/liangyan/project/zzZA/20240731_liunan-DNAmethylation/Release/Mapping/methylation_extract
# DSS package
conda activate deseq2
sbatch -J job -p cu -N 1 -n 1 -o %j.out -e %j.err --cpus-per-task=1 --wrap="
Rscript DSS-DE-analysis.r
"

sed -i 's/meanMethy1/Methylation_patient/g' DSS-DE-analysis.tab
sed -i 's/meanMethy2/Methylation_control/g' DSS-DE-analysis.tab

grep "ITGB2" gencode.v46.MANE_Select.gtf |\
awk 'BEGIN{FS="\t";OFS="\t"}{if($3=="transcript"){print $1, $4 - 1 - 200000, $5 + 200000, "ITGB2"}}' > gencode.ITGB2_extend100k.bed
tail -n+2 DSS-DE-analysis.tab | bedtools intersect -a stdin -b gencode.ITGB2_extend100k.bed -wa > DSS-DE-analysis.ITGB2.tab

# ------------------------------whole genome sequencing analysis------------------------------
# reference pipeline: https://people.duke.edu/~ccc14/duke-hts-2017/Statistics/08032017/GATK-pipeline-sample.html
cd /public/home/liunangroup/liangyan/project/zzZA/20240801_WGS_liunan

thread="32"
mem=$(echo "$thread * 5"|bc)G
sbatch -J job -p cu -N 1 -n 1 -o %j.out -e %j.err --cpus-per-task=${thread} --mem=$mem --wrap="
sh wgs_pipeline.sh $thread $mem
"

