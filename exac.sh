#!/bin/bash
#SBATCH --time=144:00:00
#SBATCH --cpus-per-task=8
#SBATCH  --mem=64G
# Specify a job name:
#SBATCH -J exac
# Specify an output file
#SBATCH -o exac.out
#SBATCH -e exac.err

module load samtools
module load bwa
module load R

myfasta="/users/aleith/scratch/broad/ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta"
mydict="/users/aleith/scratch/broad/ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.dict"
mysnp="/users/aleith/scratch/broad/ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf"
myindel="/users/aleith/scratch/broad/ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
myhapmap="/users/aleith/scratch/broad/ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf"
myomni="/users/aleith/scratch/broad/ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf"
mymills="/users/aleith/scratch/broad/ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
myhighconf="/users/aleith/scratch/broad/ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf"


for f in /users/aleith/scratch/broad/ftp.broadinstitute.org/bundle/hg19/exome_data/bottle/SRX265482/*_1.fastq
do
    echo $f
    file=$(echo "$f" | rev | cut -d"/" -f1 | rev | cut -d"_" -f1)
    samplelabel=$file
    myR1=$samplelabel"_1.fastq"
    myR2=$samplelabel"_2.fastq"
    mysamplebase=$samplelabel


time bwa mem \
         -t 8 \
         -R '@RG\tID:'$mysamplebase'\tLB:'$mysamplebase'\tSM:'$mysamplebase'\tPL:ILLUMINA' \
         $myfasta \
         $myR1 $myR2 \
         > $mysamplebase".sam"

samtools view -S -b $mysamplebase".sam" > $mysamplebase".bam"
samtools sort $mysamplebase".bam" -T temp -o $mysamplebase"_sorted.bam"
samtools index $mysamplebase"_sorted.bam"


picard MeanQualityByCycle \
    INPUT=$mysamplebase"_sorted.bam" \
    CHART_OUTPUT=$mysamplebase"_mean_quality_by_cycle.pdf" \
    OUTPUT=$mysamplebase"_read_quality_by_cycle.txt" \
    VALIDATION_STRINGENCY=LENIENT \
    REFERENCE_SEQUENCE=$myfasta

picard QualityScoreDistribution \
    INPUT=$mysamplebase"_sorted.bam" \
    CHART_OUTPUT=$mysamplebase"_mean_quality_overall.pdf" \
    OUTPUT=$mysamplebase"_read_quality_overall.txt" \
    VALIDATION_STRINGENCY=LENIENT \
    REFERENCE_SEQUENCE=$myfasta

picard CollectWgsMetrics \
    INPUT=$mysamplebase"_sorted.bam" OUTPUT=$mysamplebase"_stats_picard.txt" \
    REFERENCE_SEQUENCE=$myfasta \
    MINIMUM_MAPPING_QUALITY=20 \
    MINIMUM_BASE_QUALITY=20 \
    VALIDATION_STRINGENCY=LENIENT





picard MarkDuplicates I=$mysamplebase"_sorted.bam" \
    O=$mysamplebase"_sorted_dedup.bam" M=$mysamplebase"_dedup_metric.txt" \
    VALIDATION_STRINGENCY=LENIENT \
    REMOVE_DUPLICATES=true CREATE_INDEX=true


gatk -T RealignerTargetCreator \
    -R $myfasta \
    -I $mysamplebase"_sorted_dedup.bam" \
    -known $myindel \
    -o $mysamplebase"_realign_targets.intervals"

gatk -T IndelRealigner \
    -R $myfasta \
    -known $myindel \
    -targetIntervals $mysamplebase"_realign_targets.intervals"  \
    -I $mysamplebase"_sorted_dedup.bam" \
    -o $mysamplebase"_sorted_dedup_realigned.bam" \
    --filter_bases_not_stored

picard BuildBamIndex \
    I=$mysamplebase"_sorted_dedup_realigned.bam"


gatk -T BaseRecalibrator \
    -R $myfasta \
    -knownSites $mysnp \
    -knownSites $myindel \
    -I $mysamplebase"_sorted_dedup_realigned.bam" \
    -o $mysamplebase"_recal_table.txt" \
    -nct 8

gatk -T BaseRecalibrator \
    -R $myfasta \
    -I $mysamplebase"_sorted_dedup_realigned.bam" \
    -knownSites $mysnp \
    -knownSites $myindel \
    -BQSR $mysamplebase"_recal_table.txt" \
    -o $mysamplebase"_post_recal_table.txt"




picard CollectOxoGMetrics \
    I=$mysamplebase"_sorted_dedup_realigned.bam" \
    O=$mysamplebase"_sorted_dedup_realigned_oxoG_metrics.txt" \
    R=$myfasta \
    VALIDATION_STRINGENCY=LENIENT

picard CollectSequencingArtifactMetrics \
    I=$mysamplebase"_sorted_dedup_realigned.bam" \
    O=$mysamplebase"_sorted_dedup_realigned_artifact_metrics.txt" \
    R=$myfasta \
    VALIDATION_STRINGENCY=LENIENT

#CalculateHsMetrics is deprecated (fully removed from current picard)

picard CollectAlignmentSummaryMetrics \
    R=$myfasta \
    I=$mysamplebase"_sorted_dedup_realigned.bam" \
    O=$mysamplebase"_sorted_dedup_realigned_alignment_metrics.txt"

picard CollectInsertSizeMetrics \
    I=$mysamplebase"_sorted_dedup_realigned.bam" \
    O=$mysamplebase"_sorted_dedup_realigned_insert_size_metrics.txt" \
    H=$mysamplebase"_sorted_dedup_realigned_insert_size_histogram.pdf" \
    M=0.5





gatk -T AnalyzeCovariates \
    -R $myfasta \
    -before $mysamplebase"_recal_table.txt" \
    -after $mysamplebase"_post_recal_table.txt" \
    -plots $mysamplebase"_recalibration_plots.pdf"

gatk -T PrintReads \
    -R $myfasta \
    -I $mysamplebase"_sorted_dedup_realigned.bam" \
    -BQSR $mysamplebase"_recal_table.txt" \
    -o  $mysamplebase"_GATK.bam"

picard BuildBamIndex \
    I= $mysamplebase"_GATK.bam"

gatk -T HaplotypeCaller --disable_auto_index_creation_and_locking_when_reading_rods \
    -R $myfasta \
    -I  $mysamplebase"_GATK.bam" \
    --dbsnp  $mysnp \
    -stand_call_conf 30 \
    --minPruning 3 --maxNumHaplotypesInPopulation 200 \
    --max_alternate_alleles 3 -variant_index_parameter 128000 \
    -variant_index_type LINEAR -contamination 0.0 \
    --genotyping_mode DISCOVERY \
    --emitRefConfidence GVCF \
    -o $mysamplebase"_GATK-HC.g.vcf"
#"when running on exome data, we use -L to specify a file containing the list of exome targets corresponding to the capture kit that was used to generate the exome libraries."

done

ls -1a *GATK-HC_exac.g.vcf > vcfs_exac.list

gatk -T GenotypeGVCFs \
    -R $myfasta \
    --variant "vcfs_exac.list" \
    -o "all_GATK-HC_joined_exac.g.vcf"


picard MakeSitesOnlyVcf \
    INPUT="all_GATK-HC_joined_exac.g.vcf" \
    OUTPUT="all_GATK-HC_joined_sites_only_exac.g.vcf" \
    VALIDATION_STRINGENCY=LENIENT



gatk -T VariantRecalibrator \
    -R $myfasta \
    -input "all_GATK-HC_joined_sites_only_exac.g.vcf" \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $myhapmap \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 $myomni \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $myhighconf \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $mysnp \
    --maxGaussians 6 \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an InbreedingCoeff \
    -mode SNP \
    -allPoly -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 \
    -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 \
    -tranche 90.0 \
    -recalFile "all_recalibrate_SNP_exac.recal" \
    -tranchesFile "all_recalibrate_SNP_exac.tranches" \
    -rscriptFile "all_recalibrate_SNP_plots_exac.R"

gatk -T ApplyRecalibration \
    -R $myfasta \
    -input "all_GATK-HC_joined_sites_only_exac.g.vcf" \
    -mode SNP \
    --ts_filter_level 99.6 \
    -recalFile "all_recalibrate_SNP_exac.recal" \
    -tranchesFile "all_recalibrate_SNP_exac.tranches" \
    -o "all_recalibrated_snps_raw_indels_exac.vcf"

gatk -T VariantRecalibrator \
    -R $myfasta \
    -input "all_recalibrated_snps_raw_indels_exac.vcf" \
    -resource:mills,known=false,training=true,truth=true,prior=12.0 $mymills \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $mysnp \
    -an FS -an ReadPosRankSum -an InbreedingCoeff -an MQRankSum -an QD \
    -mode INDEL \
    -allPoly -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 \
    -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    --maxGaussians 6 \
    -recalFile "all_recalibrate_INDEL_exac.recal" \
    -tranchesFile "all_recalibrate_INDEL_exac.tranches" \
    -rscriptFile "all_recalibrate_INDEL_plots_exac.R"

gatk -T ApplyRecalibration \
    -R $myfasta \
    -input "all_recalibrated_snps_raw_indels_exac.vcf" \
    -mode INDEL \
    --ts_filter_level 95.0 \
    -recalFile "all_recalibrate_INDEL_exac.recal" \
    -tranchesFile "all_recalibrate_INDEL_exac.tranches" \
    -o "all_recalibrated_variants_exac.vcf"
