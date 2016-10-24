#!/bin/bash
#PBS -l walltime=20:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Somatic Amplicon Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_SAMPLE
version="dev"

# Directory structure required for pipeline
#
# /data
# └── results
#     └── seqId
#         ├── panel1
#         │   ├── sample1
#         │   ├── sample2
#         │   └── sample3
#         └── panel2
#             ├── sample1
#             ├── sample2
#             └── sample3
#
# Script 1 runs in sample folder, requires fastq files split by lane

phoneTrello() {
    #Call trello API
    /share/apps/node-distros/node-v4.4.7-linux-x64/bin/node \
    /data/diagnostics/scripts/TrelloAPI.js \
    "$1" "$2" #seqId & message
}

countQCFlagFails() {
    #count how many core FASTQC tests failed
    grep -E "Basic Statistics|Per base sequence quality|Per tile sequence quality|Per sequence quality scores|Per base N content" "$1" | \
    grep -v ^PASS | \
    grep -v ^WARN | \
    wc -l | \
    sed 's/^[[:space:]]*//g'
}

#load sample & pipeline variables
. *.variables
. /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel".variables

phoneTrello "$seqId" "Starting SomaticAmplicon analysis for $sampleId ..."

### Preprocessing ###

#record FASTQC pass/fail
rawSequenceQuality=PASS

#convert FASTQ to uBAM & add RGIDs
for fastqPair in $(ls "$sampleId"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do

    #parse fastq filenames
    laneId=$(echo "$fastqPair" | cut -d_ -f3)
    read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
    read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

    #trim adapters
    /share/apps/cutadapt-distros/cutadapt-1.9.1/bin/cutadapt \
    -a "$read1Adapter" \
    -A "$read2Adapter" \
    -m 50 \
    -o "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    -p "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    "$read1Fastq" \
    "$read2Fastq"

    #merge overlapping reads
    /share/apps/pear-distros/pear-0.9.10-bin-64/pear-0.9.10-bin-64 \
    -f "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    -r "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    -o "$seqId"_"$sampleId"_"$laneId"_merged.fastq \
    -j 12

    #convert fastq to ubam
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar FastqToSam \
    F1="$seqId"_"$sampleId"_"$laneId"_merged.fastq \
    O="$seqId"_"$sampleId"_"$laneId"_unaligned.bam \
    QUALITY_FORMAT=Standard \
    READ_GROUP_NAME="$seqId"_"$laneId"_"$sampleId" \
    SAMPLE_NAME="$sampleId" \
    LIBRARY_NAME="$worklistId"_"$sampleId"_"$panel" \
    PLATFORM_UNIT="$seqId"_"$laneId" \
    PLATFORM="ILLUMINA" \
    SEQUENCING_CENTER="IMG" \
    PREDICTED_INSERT_SIZE="$expectedInsertSize" \
    SORT_ORDER=queryname \
    MAX_RECORDS_IN_RAM=2000000 \
    TMP_DIR=/state/partition1/tmpdir

    #fastqc
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R1.fastq
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R2.fastq

    #check FASTQC output
    if [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R1_fastqc/summary.txt) -gt 0 ] || [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R2_fastqc/summary.txt) -gt 0 ]; then
        phoneTrello "$seqId" "$sampleId has failed raw sequence quality check for $laneId"
        rawSequenceQuality=FAIL
    fi

    #clean up
    rm "$seqId"_"$sampleId"_"$laneId"_R1.fastq "$seqId"_"$sampleId"_"$laneId"_R2.fastq "$seqId"_"$sampleId"_"$laneId"_merged.fastq

done

#merge lane bams
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MergeSamFiles \
$(ls "$seqId"_"$sampleId"_*_unaligned.bam | sed 's/^/I=/' | tr '\n' ' ') \
SORT_ORDER=queryname \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT \
USE_THREADING=true \
MAX_RECORDS_IN_RAM=2000000 \
TMP_DIR=/state/partition1/tmpdir \
O="$seqId"_"$sampleId"_unaligned.bam

#uBam2fq, map & MergeBamAlignment
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar SamToFastq \
I="$seqId"_"$sampleId"_unaligned.bam \
FASTQ=/dev/stdout \
NON_PF=true \
MAX_RECORDS_IN_RAM=2000000 \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=/state/partition1/tmpdir | \
/share/apps/bwa-distros/bwa-0.7.15/bwa mem \
-M \
-t 10 \
-p \
/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
/dev/stdin | \
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MergeBamAlignment \
ATTRIBUTES_TO_RETAIN=X0 \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM="$seqId"_"$sampleId"_unaligned.bam \
OUTPUT="$seqId"_"$sampleId"_aligned.bam \
R=/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
PAIRED_RUN=false \
SORT_ORDER="coordinate" \
IS_BISULFITE_SEQUENCE=false \
ALIGNED_READS_ONLY=false \
CLIP_ADAPTERS=false \
MAX_RECORDS_IN_RAM=2000000 \
MAX_INSERTIONS_OR_DELETIONS=-1 \
UNMAP_CONTAMINANT_READS=false \
CLIP_OVERLAPPING_READS=false \
ALIGNER_PROPER_PAIR_FLAGS=false \
ATTRIBUTES_TO_RETAIN=XS \
INCLUDE_SECONDARY_ALIGNMENTS=true \
CREATE_INDEX=true \
TMP_DIR=/state/partition1/tmpdir

#Realign soft clipped bases
/share/apps/jre-distros/jre1.8.0_101/bin/java -Xmx2g -jar /data/diagnostics/apps/AmpliconRealigner/AmpliconRealigner-1.0.0.jar \
-I "$seqId"_"$sampleId"_aligned.bam \
-O "$seqId"_"$sampleId"_amplicon_realigned.bam \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-T /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx24g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-known /data/db/human/cosmic/b37/Cosmic69Indels.vcf \
-I "$seqId"_"$sampleId"_amplicon_realigned.bam \
-o "$seqId"_"$sampleId"_indel_realigned.intervals \
-L /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip "$padding" \
-nt 12 \
-dt NONE

#Realign around indels
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-known /data/db/human/cosmic/b37/Cosmic69Indels.vcf \
-targetIntervals "$seqId"_"$sampleId"_indel_realigned.intervals \
-I "$seqId"_"$sampleId"_amplicon_realigned.bam \
-o "$seqId"_"$sampleId"_indel_realigned.bam \
-dt NONE

#soft clip PCR primers
/share/apps/jre-distros/jre1.8.0_101/bin/java -Xmx2g -jar /data/diagnostics/apps/SoftClipPCRPrimer/SoftClipPCRPrimer-1.0.0.jar \
-I "$seqId"_"$sampleId"_indel_realigned.bam \
-O "$seqId"_"$sampleId".bam \
-T /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed

### Variant calling ###

#SNPs and Indels with Illumina Pisces

### QC ###

#Convert BED to interval_list for later
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar BedToIntervalList \
I=/data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed \
O="$panel"_ROI.interval_list \
SD=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.dict 

#Calculate insert size
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar CollectInsertSizeMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_insert_metrics.txt \
H="$seqId"_"$sampleId"_insert_metrics.pdf

#HsMetrics: capture & pooling performance
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar CollectHsMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_hs_metrics.txt \
R=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
BAIT_INTERVALS="$panel"_ROI.interval_list \
TARGET_INTERVALS="$panel"_ROI.interval_list

#Generate per-base coverage: variant detection sensitivity
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx12g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_"$sampleId"_DepthOfCoverage \
-I "$seqId"_"$sampleId".bam \
-L /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed \
--countType COUNT_FRAGMENTS \
--minMappingQuality 20 \
--minBaseQuality 10 \
--omitIntervalStatistics \
-ct "$minimumCoverage" \
-nt 12 \
-dt NONE

#Calculate gene percentage coverage
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /data/diagnostics/apps/CoverageCalculator-2.0.0/CoverageCalculator-2.0.0.jar \
"$seqId"_"$sampleId"_DepthOfCoverage \
/data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_genes.txt \
/state/partition1/db/human/refseq/ref_GRCh37.p13_top_level.gff3 \
-p5 \
-d"$minimumCoverage" \
> "$seqId"_"$sampleId"_PercentageCoverage.txt

#Gather QC metrics
meanInsertSize=$(head -n8 "$seqId"_"$sampleId"_insert_metrics.txt | tail -n1 | cut -s -f5) #mean insert size
sdInsertSize=$(head -n8 "$seqId"_"$sampleId"_insert_metrics.txt | tail -n1 | cut -s -f6) #insert size standard deviation
totalReads=$(head -n8 "$seqId"_"$sampleId"_hs_metrics.txt | tail -n1 | cut -s -f6) #The total number of reads in the SAM or BAM file examine.
pctSelectedBases=$(head -n8 "$seqId"_"$sampleId"_hs_metrics.txt | tail -n1 | cut -s -f19) #On+Near Bait Bases / PF Bases Aligned.
totalTargetedUsableBases=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f2) #total number of usable bases. NB BQSR requires >= 100M, ideally >= 1B
meanOnTargetCoverage=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f3) #avg usable coverage
pctTargetBasesCt=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f7) #percentage panel covered with good enough data for variant detection

#Print QC metricss
echo -e "TotalReads\tRawSequenceQuality\tTotalTargetUsableBases\tPctSelectedBases\tpctTargetBasesCt\tMeanOnTargetCoverage\tMeanInsertSize\tSDInsertSize" > "$seqId"_"$sampleId"_qc.txt
echo -e "$totalReads\t$rawSequenceQuality\t$totalTargetedUsableBases\t$pctSelectedBases\t$pctTargetBasesCt\t$meanOnTargetCoverage\t$meanInsertSize\t$sdInsertSize" >> "$seqId"_"$sampleId"_qc.txt

#print metaline for final VCF
echo \#\#SAMPLE\=\<ID\="$sampleId",WorklistId\="$worklistId",SeqId\="$seqId",Panel\="$panel",PipelineName\=SomaticAmplicon,PipelineVersion\="$version",RawSequenceQuality\="$rawSequenceQuality",MeanInsertSize\="$meanInsertSize",SDInsertSize\="$sdInsertSize",TotalReads\="$totalReads",PctSelectedBases\="$pctSelectedBases",MeanOnTargetCoverage\="$meanOnTargetCoverage",pctTargetBasesCt\="$pctTargetBasesCt",TotalTargetedUsableBases\="$totalTargetedUsableBases",RemoteBamFilePath\=$(find $PWD -type f -name "$seqId"_"$sampleId".bam)\> > "$seqId"_"$sampleId"_meta.txt

### Clean up ###

#create final file lists
find $PWD -name "$seqId"_"$sampleId".bam >> ../FinalBams.list

#delete unused files

#check if all VCFs are written
if [ $(find .. -maxdepth 1 -mindepth 1 -type d | wc -l | sed 's/^[[:space:]]*//g') -eq $(sort ../GVCFs.list | uniq | wc -l | sed 's/^[[:space:]]*//g') ]; then
    echo -e "seqId=$seqId\npanel=$panel" > ../variables
    cp 2_SomaticAmplicon.sh .. && cd .. && qsub 2_SomaticAmplicon.sh
fi