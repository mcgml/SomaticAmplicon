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
# Script 2 runs in panel folder, requires final bams

phoneTrello() {
    #Call trello API
    /share/apps/node-distros/node-v4.4.7-linux-x64/bin/node \
    /data/diagnostics/scripts/TrelloAPI.js \
    "$1" "$2" #seqId & message
}

addMetaDataToVCF(){
    output=$(echo "$1" | sed 's/\.vcf/_meta\.vcf/g')
    grep '^##' "$1" > "$output"
    for sample in $(/share/apps/bcftools-distros/bcftools-1.3.1/bcftools query -l "$1"); do
        cat "$sample"/"$seqId"_"$sample"_meta.txt >> "$output"
    done
    grep -v '^##' "$1" >> "$output"
}

#load run & pipeline variables
. variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

### Variant calling ###

#load mono
. /opt/mono/env.sh

#SNPs and Indels with Illumina Pisces
mono /share/apps/pisces-distros/5.1.3.60/Pisces.exe \
-B $(tr '\n' ',' FinalBams.list) \
-g /data/db/human/gatk/2.8/b37 \
-i /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-f 0.01 \
-fo false \
-b 20 \
-q 100 \
-c 20 \
-a 20 \
-F 20 \
-m 20 \
-gVCF false

#rename Pisces VCF file
#TODO

#Annotate with low complexity region length using mdust
/share/apps/bcftools-distros/bcftools-1.3.1/bcftools annotate \
-a /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.mdust.v34.lpad1.bed.gz \
-c CHROM,FROM,TO,LCRLen \
-h <(echo '##INFO=<ID=LCRLen,Number=1,Type=Integer,Description="Overlapping mdust low complexity region length (mask cutoff: 34)">') \
-o "$seqId"_variants.lcr.vcf \
"$seqId"_variants.vcf

#Filter variants near to homopolymers
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.lcr.vcf \
--filterExpression "LCRLen > 8" \
--filterName "LowComplexity" 
--genotypeFilterExpression "DP < 20" \
--genotypeFilterName "LowDP" \
-L /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-o "$seqId"_filtered.vcf \
-dt NONE

#Add VCF meta data to final VCF
addMetaDataToVCF "$seqId"_filtered.vcf

### QC ###

#Variant Evaluation
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T VariantEval \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_variant_evaluation.txt \
--eval:"$seqId" "$seqId"_filtered_meta.vcf \
--comp:omni2.5 /state/partition1/db/human/gatk/2.8/b37/1000G_omni2.5.b37.vcf \
--comp:hapmap3.3 /state/partition1/db/human/gatk/2.8/b37/hapmap_3.3.b37.vcf \
--comp:cosmic78 /state/partition1/db/human/cosmic/b37/cosmic_78.b37.vcf \
-ST JexlExpression --select_names "snp" --select_exps "vc.isSNP()" \
-ST JexlExpression --select_names "indel" --select_exps "vc.isIndel()" \
-L /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-nt 12 \
-dt NONE

#Genotype Concordance
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T GenotypeConcordance \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_genotype_evaluation.txt \
--eval:"$seqId" "$seqId"_filtered_meta.vcf \
--comp:omni2.5 /state/partition1/db/human/gatk/2.8/b37/1000G_omni2.5.b37.vcf \
--comp:hapmap3.3 /state/partition1/db/human/gatk/2.8/b37/hapmap_3.3.b37.vcf \
--comp:cosmic78 /state/partition1/db/human/cosmic/b37/cosmic_78.b37.vcf \
-ST JexlExpression --select_names "snp" --select_exps "vc.isSNP()" \
-ST JexlExpression --select_names "indel" --select_exps "vc.isIndel()" \
-L /data/diagnostics/pipelines/SomaticAmplicon/SomaticAmplicon-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-nt 12 \
-dt NONE

### Clean up ###

#delete unused files
rm "$seqId"_variants.lcr.vcf "$seqId"_variants.lcr.vcf.idx "$seqId"_filtered.vcf "$seqId"_filtered.vcf.idx

#log with Trello
phoneTrello "$seqId" "Analysis complete"