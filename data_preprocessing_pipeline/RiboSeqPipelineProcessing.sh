#!/bin/bash
#!/usr/bin/perl
# Perl (v5.18.2)

# Ribosome Profiling Pipeline For Processing
# ----------------------------------------------------

export species=$1
export tissue=$2
export library=$3
export type=$4
export replicate=$5
export nproc=$6


# USER NAME
USER=zhongyi

# GTF_EXON
export GTF_EXON=${species}_${tissue}.exons.of.dominant.isoforms.gtf

# GTF_REFINED
export GTF_REFINED=${species}_${tissue}.refined.gtf

# DOMINANT_ISOFORMS
export DOMINANT_ISOFORMS=${species}_${tissue}.tissue.specific.isoforms.txt

# TRANSCRIPT_LENGTHS
export TRANSCRIPT_LENGTHS=${species}_${tissue}.5UTR.CDS.3UTR.length.txt

## Check If $USER Is Set Otherwise Exit
[[ -z "$USER" ]] && exit 1

## Create The Aliases
export WORK=your/own/path/...
export HOME=your/own/path/...
export pathResultDir=${WORK}/Results
export pathWorkDir=${pathResultDir}/${species}/${tissue}/${replicate}_${library}/Analysis
export pathMappingData=${pathResultDir}/${species}/${tissue}/${replicate}_${library}/Tophat/${library}.unique.bam
export pathGTFExon=${HOME}/Genome_data/${species}/${GTF_EXON}
export pathGTFRefined=${HOME}/Genome_data/${species}/${GTF_REFINED}
export pathDominantIsoforms=${HOME}/Genome_data/${species}/${DOMINANT_ISOFORMS}
export pathTranscriptLengths=${HOME}/Genome_data/${species}/${TRANSCRIPT_LENGTHS}
export pathScripts=${HOME}/scripts

########################################

## Create Work Directory ##
if [ -d ${pathWorkDir} ]; then  # -d: Return true value if exists and is a directory
echo "Path Exists"
else
mkdir -p ${pathWorkDir}
fi

## ... and test if it exists and if it is writable
[[ ! -d "${pathWorkDir}" ]] && exit 1
[[ ! -w "${pathWorkDir}" ]] && exit 1

cd ${pathWorkDir} || exit 1

#########################################
echo "Job started at $(date)"

## 1. Aggregate reads that Mapping to the Same Genomic Position
echo "### 1 ### Aggregate reads that Mapping to the Same Genomic Position ### 1 ###"
# output file: "unique_mappers_frame_analysis_index.sam"

samtools view -h -o ${pathWorkDir}/temp.sam ${pathMappingData}

if [ "$type" == "RP" ]; then

echo "Ribo-Seq!!!"

perl ${pathScripts}/Processing1_Refine_Index_Sam_File.pl ${pathWorkDir}/temp.sam 26 34 ${pathWorkDir} 1> ${pathWorkDir}/Processing1_Refine_Index_Sam_File.out 2>&1

else

echo "RNA-Seq!!!"

perl ${pathScripts}/Processing1_Refine_Index_Sam_File.pl ${pathWorkDir}/temp.sam 20 29 ${pathWorkDir} 1> ${pathWorkDir}/Processing1_Refine_Index_Sam_File.out 2>&1

fi

rm -f ${pathWorkDir}/temp.sam

echo "Finished Step 1 at $(date)"


## 2. Get Ribo-seq or RNA-seq Mapping Profiles for Each Gene
echo "### 2 ### Get Ribo-seq or RNA-seq Mapping Profiles for Each Gene ### 2 ###"
#  output dir: "distribution_of_ribosome_footprints" OR "distribution_of_mRNA_reads"
#  output file: "read_distribution_for_individual_genes.txt"
if [ "$type" == "RP" ]; then

echo "Ribo-Seq!!!"

rm -rf ${pathWorkDir}/distribution_of_ribosome_footprints

perl ${pathScripts}/Processing2_Obtain_Read_Distributions_for_Each_Gene.pl ${pathDominantIsoforms} ${pathGTFExon} ${pathGTFRefined} ${pathWorkDir}/unique_mappers_frame_analysis_index.sam ${species} ${type} ${pathWorkDir} 1> ${pathWorkDir}/Processing2_Obtain_Read_Distributions_for_Each_Gene.out 2>&1

else

echo "RNA-Seq!!!"

rm -rf ${pathWorkDir}/distribution_of_mRNA_reads

perl ${pathScripts}/Processing2_Obtain_Read_Distributions_for_Each_Gene.pl ${pathDominantIsoforms} ${pathGTFExon} ${pathGTFRefined} ${pathWorkDir}/unique_mappers_frame_analysis_index.sam ${species} ${type} ${pathWorkDir} 1> ${pathWorkDir}/Processing2_Obtain_Read_Distributions_for_Each_Gene.out 2>&1

fi

echo "Finished Step 2 at $(date)"


## 3. Obtain for Each Gene the Read Counts for 5'UTR, CDS and 3'UTR
echo "### 3 ### Obtain for Each Gene the Read Counts for 5'UTR, CDS and 3'UTR ### 3 ###"
# output file: "distribution_of_ribosome_footprints_for_5UTR_CDS_3UTR.txt" OR "distribution_of_mRNA_reads_for_5UTR_CDS_3UTR.txt"
if [ "$type" == "RP" ]; then

echo "Ribo-Seq!!!"

perl ${pathScripts}/Processing3_Obtain_Read_Distributions_for_5UTR_CDS_3UTR.pl ${pathDominantIsoforms} ${pathWorkDir}/distribution_of_ribosome_footprints ${pathWorkDir} 1> ${pathWorkDir}/Processing3_Obtain_Read_Distributions_for_5UTR_CDS_3UTR.out 2>&1

else

echo "RNA-Seq!!!"

perl ${pathScripts}/Processing3_Obtain_Read_Distributions_for_5UTR_CDS_3UTR.pl ${pathDominantIsoforms} ${pathWorkDir}/distribution_of_mRNA_reads ${pathWorkDir} 1> ${pathWorkDir}/Processing3_Obtain_Read_Distributions_for_5UTR_CDS_3UTR.out 2>&1

fi

## Get The Number of Reads for Each Gene for FPKM calculation
if [ "$type" == "RP" ]; then
echo "Ribo-Seq!!!"
export number_of_alignments="$(grep "CDS Reads: " ${pathWorkDir}/how_many_reads_are_mapped_to_each_transcript.txt | awk '{ print $3 }')"
else
echo "total RNA-Seq!!!"
export number_of_alignments="$(grep "CDS Reads: " ${pathWorkDir}/how_many_reads_are_mapped_to_each_transcript.txt | awk '{ print $3 }')"
fi

million=1000000
export millions_of_alignments=$(echo "scale = 6; ${number_of_alignments} / ${million}" | bc)
echo "millions of mapped reads: ${millions_of_alignments}"


# Copy Numbers to The Final Report
cat ${pathWorkDir}/how_many_reads_are_mapped_to_each_transcript.txt >> ${pathResultDir}/${species}/${tissue}/${replicate}_${library}/${type}\.${species}\.${tissue}\.${replicate}\.${library}\.report.txt

echo "Finished Step 3 at $(date)"


## 4. FPKM Calculation
# output file: "FPKM_and_density.txt"

echo "### 4 ### FPKM Calculation ### 4 ###"
if [ "$type" == "RP" ]; then

echo "Ribo-Seq!!!"

perl ${pathScripts}/Processing4_FPKM_Calculation.pl ${pathTranscriptLengths} ${pathWorkDir}/distribution_of_ribosome_footprints_for_5UTR_CDS_3UTR.txt ${type} ${millions_of_alignments} ${pathWorkDir} 1> ${pathWorkDir}/Processing4_FPKM_Calculation.out 2>&1

else

echo "RNA-Seq!!!"

perl ${pathScripts}/Processing4_FPKM_Calculation.pl ${pathTranscriptLengths} ${pathWorkDir}/distribution_of_mRNA_reads_for_5UTR_CDS_3UTR.txt ${type} ${millions_of_alignments} ${pathWorkDir} 1> ${pathWorkDir}/Processing4_FPKM_Calculation.out 2>&1

fi

echo "Finished Step 4 at $(date)"

echo "###################"
echo "Everything is done!"
exit 0; ## exit(0) indicates successful termination.

