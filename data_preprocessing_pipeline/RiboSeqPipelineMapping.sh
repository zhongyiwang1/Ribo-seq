#!/bin/bash
#!/usr/bin/perl
# Perl (v5.18.2)

# Ribosome Profiling Pipeline For Mapping
# ----------------------------------------------------

export species=$1
export tissue=$2
export library=$3
export type=$4
export replicate=$5
export nproc=$6


# USER NAME
USER=zhongyi

# SAMPLE NUMBER
SAMPLE_NUMBER=$(echo ${library} | sed 's/[^0-9]*//g')

# ADAPTOR TRIMMING CONFIGURATION
ADAPTOR=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

# FASTQC CONFIGURATION
FASTQC_PARAMS="--threads ${nproc} --format fastq"

# CUTADAPT CONFIGURATION
CUTADAPT_PARAMS="--adapter=${ADAPTOR} --minimum-length=6 --maximum-length=60 -q 20"

# BOWTIE2 CONFIGURATION
BOWTIE2_PARAMS="--phred33 -p ${nproc} -L 20 -N 1 -t --no-unal"

# TOPHAT CONFIGURATION
TOPHAT_PARAMS="--no-novel-juncs --library-type fr-firststrand --read-realign-edit-dist 0 --segment-length 20 --min-anchor-length 5 --min-intron-length 50 --num-threads ${nproc}"

# TOPHAT_GTF
export TOPHAT_GTF=${species}_${tissue}.gtf

## Check If $USER Is Set Otherwise Exit
[[ -z "$USER" ]] && exit 1 ## exit(1) indicates unsucessful termination.

## Create The Aliases
export WORK=your/own/path/...
export HOME=your/own/path/...
export pathRawData=${WORK}/Raw_data_NoFruitfly/${species}/${tissue}/${library}.fastq.gz
export pathBowtie2Index=${HOME}/Bowtie2_index/${species}
export pathGenomeData=${HOME}/Genome_data/${species}
export pathWorkDir=${WORK}/Results/${species}/${tissue}/${replicate}_${library}
export pathFastqc=${pathWorkDir}/Fastqc
export pathTrimmed=${pathWorkDir}/Trimmed
export pathNonrRNA=${pathWorkDir}/NonrRNA
export pathNontRNA=${pathWorkDir}/NontRNA
export pathNonsnoRNA=${pathWorkDir}/NonsnoRNA
export pathTophat=${pathWorkDir}/Tophat
export pathReport=${pathWorkDir}/${type}\.${species}\.${tissue}\.${replicate}\.${library}\.report.txt

########################################

## Create Work Directory ##
if [ -d ${pathWorkDir} ]; then  # -d: Return true if the directory exists
echo "Path Exists"
else
mkdir -p ${pathWorkDir}
fi

## ... and test if it exists and it's writable
[[ ! -d "${pathWorkDir}" ]] && exit 1
[[ ! -w "${pathWorkDir}" ]] && exit 1

cd ${pathWorkDir} || exit 1

#########################################
echo "Job started at $(date)"

### 1. Quality Control
echo "### 1 ### Quality Control ### 1 ###"
mkdir -p ${pathFastqc}
mkdir -p ${pathFastqc}/RawData

fastqc ${FASTQC_PARAMS} ${pathRawData} --outdir ${pathFastqc}/RawData &


### 2. Adapter Trimming
echo "### 2 ### Adapter Trimming ### 2 ###"
mkdir -p ${pathTrimmed}

cutadapt ${CUTADAPT_PARAMS} ${pathRawData} \
--output=${pathTrimmed}/${library}.trimmed.fastq \
--info-file=${pathTrimmed}/${library}.trimming.info.txt \
--too-short-output=${pathTrimmed}/${library}.trimmed.short.reads.fastq \
--too-long-output=${pathTrimmed}/${library}.trimmed.long.reads.fastq \
1>${pathTrimmed}/report.txt 2>&1


### 3. Remove rRNA Reads
echo "### 3 ### Remove rRNA Reads ### 3 ###"
mkdir -p ${pathNonrRNA}

bowtie2 ${BOWTIE2_PARAMS} \
-x ${pathBowtie2Index}/${species}_all_rrna/${species}_all_rrna \
-U ${pathTrimmed}/${library}.trimmed.fastq \
-S ${pathNonrRNA}/${library}.rrna.sam \
--un ${pathNonrRNA}/${library}.nonrrna.fastq \
--al ${pathNonrRNA}/${library}.rrna.fastq \
1>${pathNonrRNA}/report.txt 2>&1


### 4. Remove Human, Mouse and Rat rRNA Reads
echo "### 4 ### Remove Human, Mouse and Rat rRNA Reads ### 4 ###"

bowtie2 ${BOWTIE2_PARAMS} \
-x /beegfs/home/hd/hd_hd/hd_fd141/Bowtie2_index/Human/Human_Mouse_Rat_all_rrna/Human_Mouse_Rat_all_rrna \
-U ${pathNonrRNA}/${library}.nonrrna.fastq \
-S ${pathNonrRNA}/${library}.rrna.human.mouse.rat.sam \
--un ${pathNonrRNA}/${library}.nonrrna.nonhuman.nonmouse.nonrat.fastq \
--al ${pathNonrRNA}/${library}.rrna.human.mouse.rat.fastq \
1>${pathNonrRNA}/report.human.mouse.rat.rrna.txt 2>&1


### 5. Remove tRNA Reads
echo "### 5 ### Remove tRNA Reads ### 5 ###"
mkdir -p ${pathNontRNA}

bowtie2 ${BOWTIE2_PARAMS} \
-x ${pathBowtie2Index}/${species}_all_trna/${species}_all_trna \
-U ${pathNonrRNA}/${library}.nonrrna.nonhuman.nonmouse.nonrat.fastq \
-S ${pathNontRNA}/${library}.nonrrna.trna.sam \
--un ${pathNontRNA}/${library}.nonrrna.nontrna.fastq \
--al ${pathNontRNA}/${library}.nonrrna.trna.fastq \
1>${pathNontRNA}/report.txt 2>&1


### 6. Remove snoRNA Reads
echo "### 6 ### Remove snoRNA Reads ### 6 ###"
mkdir -p ${pathNonsnoRNA}

bowtie2 ${BOWTIE2_PARAMS} \
-x ${pathBowtie2Index}/${species}_all_snorna/${species}_all_snorna \
-U ${pathNontRNA}/${library}.nonrrna.nontrna.fastq \
-S ${pathNonsnoRNA}/${library}.nonrrna.nontrna.snorna.sam \
--un ${pathNonsnoRNA}/${library}.nonrrna.nontrna.nonsnorna.fastq \
--al ${pathNonsnoRNA}/${library}.nonrrna.nontrna.snorna.fastq \
1>${pathNonsnoRNA}/report.txt 2>&1


### 7. Trim RNA-Seq Reads to 29 nt When Longer Than 29 nt
echo "### 7 ### Trim RNA-Seq Reads to 29 nt when longer than 29 nt ### 7 ###"

if [ "${type}" == "TR" ]; then

echo "Library type is RNA-Seq."

awk ' { if (length($1)<3 || $2~/:/) print $0 ; else print substr($0, 1, 29) ;} ' ${pathNonsnoRNA}/${library}.nonrrna.nontrna.nonsnorna.fastq > ${pathNonsnoRNA}/${library}.nonrrna.nontrna.nonsnorna.trimmed.fastq

# only keep TR reads between 20 and 29 nt
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 20 && length(seq) <= 29) {print header, seq, qheader, qseq}}' < ${pathNonsnoRNA}/${library}.nonrrna.nontrna.nonsnorna.trimmed.fastq > ${pathNonsnoRNA}/${library}.20_29.fastq

export InputRPTRFastq=${pathNonsnoRNA}/${library}.20_29.fastq

else

# only keep RP reads between 26 and 34 nt
echo "Library type is Ribo-Seq."

awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 26 && length(seq) <= 34) {print header, seq, qheader, qseq}}' < ${pathNonsnoRNA}/${library}.nonrrna.nontrna.nonsnorna.fastq > ${pathNonsnoRNA}/${library}.26_34.fastq

export InputRPTRFastq=${pathNonsnoRNA}/${library}.26_34.fastq

fi

###################### mapping statistics ############################
# how many raw reads?
gzip -d -c ${pathRawData} | wc -l | \
awk '{print "Raw Reads: ", $1 / 4;}' >> ${pathReport}

# how many reads left after trimming?
wc -l < ${pathTrimmed}/${library}.trimmed.fastq | \
awk '{print "Trimmed Reads: ", $1 / 4;}' >> ${pathReport}

# how many reads left after removing rRNA reads?
wc -l < ${pathNonrRNA}/${library}.nonrrna.fastq | \
awk '{print "NonrRNA Reads: ", $1 / 4;}' >> ${pathReport}

# how many reads left after removing human-mouse-rat rRNA reads?
wc -l < ${pathNonrRNA}/${library}.nonrrna.nonhuman.nonmouse.nonrat.fastq | \
awk '{print "NonHumanMouseRatrRNA Reads: ", $1 / 4;}' >> ${pathReport}

# how many reads left after removing tRNA reads?
wc -l < ${pathNontRNA}/${library}.nonrrna.nontrna.fastq | \
awk '{print "NontRNA Reads: ", $1 / 4;}' >> ${pathReport}

# how many reads left after removing snoRNA reads?
wc -l < ${pathNonsnoRNA}/${library}.nonrrna.nontrna.nonsnorna.fastq | \
awk '{print "NonsnoRNA Reads: ", $1 / 4;}' >> ${pathReport}

# how many useable reads left (RP: 26-34 nt and TR: 20-29 nt)?
wc -l < ${InputRPTRFastq} | \
awk '{print "Useable Reads: ", $1 / 4;}' >> ${pathReport}
######################################################################


### 8. Align Remaining Reads to the Genome and Tissue-Specific Transcriptome
echo "### 8 ### Align Remaining Reads to the Genome and Tissue-Specific Transcriptome ### 8 ###"
mkdir -p ${pathTophat}

tophat2 ${TOPHAT_PARAMS} \
--output-dir ${pathTophat} \
--GTF ${pathGenomeData}/${TOPHAT_GTF} \
${pathBowtie2Index}/${species}_toplevel_dna/${species}_toplevel_dna \
${InputRPTRFastq} \
1>${pathTophat}/report.txt 2>&1

# copy Tophat2 mapping summary to the report
cat ${pathTophat}/align_summary.txt >> ${pathReport}

### 9. Extract Unique Alignments and Count Reads for Each Gene (up to one mismatch)
echo "### 9 ### Extract Unique Alignments and Count Reads for Each Gene ### 9 ###"
samtools view -h ${pathTophat}/accepted_hits.bam | grep -E -i "(NM:i:[01]\b)|(^@)" | \
grep -E -i "(NH:i:1\b)|(^@)" -> ${pathTophat}/${library}.unique.sam

# how many uniquely mapped reads?
grep -v "^@" ${pathTophat}/${library}.unique.sam | wc -l | \
awk '{print "Uniquely Mapped Reads: ", $1;}' >> ${pathReport}


### 10. Convert Sam to Bam
echo "### 10 ### Convert Sam to Bam ### 10 ###"
samtools view -b -S ${pathTophat}/${library}.unique.sam -o ${pathTophat}/${library}.unique.bam


### 11. Sort Bam
echo "### 11 ### Sort Bam ### 11 ###"
samtools sort ${pathTophat}/${library}.unique.bam ${pathTophat}/${library}.unique.sorted


### 12. Index Bam
echo "### 12 ### Index Bam ### 12 ###"
samtools index ${pathTophat}/${library}.unique.sorted.bam ${pathTophat}/${library}.unique.sorted.bai

echo "###################"
echo "Everything is done!"
exit 0; ## exit(0) indicates successful termination.
