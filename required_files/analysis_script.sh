#!/bin/bash

# Preparation and setup of required files

FILES=$(pwd)
WKDIR=$(echo $FILES | sed 's:/required_files::g')

read -p 'Do you want to retrieve genomic data from the CGD? (yes or no): ' GENEDATA

# Ask for raw data file format

echo 'Are the data in bam or fastq format?'
read -p 'Specify file format (bam or fastq): ' FORMAT

read -p 'Do you want to do a quality control of the raw data (yes or no): ' QCRAW

read -p 'Are the libraries prepared in a strand-specific way? (yes or no): ' STRANDED

read -p 'Are the data paired-end? (yes or no): ' PAIRED

read -p 'How many threads (cores) should be used for the analysis (use 1, if you are not sure): ' THREAD

if [ $GENEDATA == 'yes' ]
then
	echo 'Retrieving genomic data from CGD.'

	wget http://www.candidagenome.org/download/sequence/C_auris_B8441/current/C_auris_B8441_current_chromosomes.fasta.gz     ## include WKDIR/required_files folder in wget!!!
	wget http://www.candidagenome.org/download/gff/C_auris_B8441/C_auris_B8441_current_features.gff
	gunzip C_auris_B8441_current_chromosomes.fasta
else
	echo 'No genomic data are retrieved.'
fi


GENOME=$WKDIR/required_files/C_auris_B8441_current_chromosomes.fasta
FEATURES=$WKDIR/required_files/C_auris_B8441_current_features_rm_rRNA.gff
ADAPT1=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
ADAPT2=GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG
rRNA=$WKDIR/required_files/rRNAloci.bed
mkdir $WKDIR/QC
PICARD=$WKDIR/required_files/picard.jar


# QC of raw data

if [ $QCRAW == 'yes' ]
then
	mkdir $WKDIR/QC_raw
	echo 'Quality control of raw data:'
	if [ $FORMAT == 'bam' ]
	then
		for i in $WKDIR/*.bam
		do
			fastqc -o $WKDIR/QC_raw $i
		done
	else
		for SNAME in $(ls $WKDIR | egrep '(\.f.*q$)|(q\.gz$)')
		do
			i=$WKDIR/$SNAME
			fastqc -o $WKDIR/QC_raw $i
		done
	fi
	multiqc -o $WKDIR/QC_raw $WKDIR/QC_raw
else
	echo 'No QC of raw data done.'
fi


# Convert .bam to .fastq format

if [ $FORMAT == 'bam' ]
then
	echo 'File format is bam.'
	for i in $WKDIR/*.bam
	do
		bamToFastq -i $i -fq $i.fq
	done
elif [ $FORMAT == 'fastq' ]
then
	echo 'File format is fastq.'
else
	echo 'Invalid file format! Options are "bam" or "fastq".'
	exit
fi


# Adapter removal with cutadapt and mapping of all files with NGM

for SNAME in $(ls $WKDIR | egrep '(\_1.f.*q$)|(_1.f.*q\.gz$)')

do
	i=$WKDIR/$SNAME # sets sample name and file path

	SAMPLE=$(echo ${SNAME} | sed "s/_1\.fq.gz//") 
	
	if [ $PAIRED == 'yes' ]
	then
		echo ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz

		cutadapt -q 30 -O 1 -a $ADAPT1 -A $ADAPT2 -o ${SAMPLE}_trimmed_1.fq.gz -p ${SAMPLE}_trimmed_2.fq.gz ${SAMPLE}_1.fq.gz  ${SAMPLE}_2.fq.gz >$WKDIR/QC/Cutadapt_$SAMPLE.txt  

	else
		cutadapt -q 30 -O 1 -a $ADAPT1 $i > $i.trimmed.fq 2>$WKDIR/QC/Cutadapt_$SNAME.txt
	fi
	rm ${SAMPLE}_1.fq.gz  
	rm ${SAMPLE}_2.fq.gz

	if [ $PAIRED == 'yes' ]
	then
		ngm -r $GENOME -1 ${SAMPLE}_trimmed_1.fq.gz -2 ${SAMPLE}_trimmed_2.fq.gz -o $SAMPLE.trimmed.fq.bam -p -b -Q 30 -t $THREAD -g

	else
		ngm -q $i.trimmed.fq -r $GENOME -o $i.trimmed.fq.bam -b -Q 30 -t $THREAD
	fi
	rm ${SAMPLE}_trimmed_1.fq.gz
	rm ${SAMPLE}_trimmed_2.fq.gz

	samtools sort -@ $THREAD $SAMPLE.trimmed.fq.bam -o $SAMPLE.trimmed.fq.bam.sort.bam   # sort .bam files using samtools
	rm $SAMPLE.trimmed.fq.bam

	# Labelling of duplicated reads and removal of optical duplicates
	java -jar $PICARD MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT  I=$SAMPLE.trimmed.fq.bam.sort.bam O=$SAMPLE.final.bam M=$WKDIR/QC/$SAMPLE.markdup.metrics.txt
	rm $SAMPLE.trimmed.fq.bam.sort.bam

	#Quality control and statistics about mapped samples
	samtools flagstat $SAMPLE.final.bam >> $WKDIR/QC/$SAMPLE.final.flagstat_analysis.txt   # flagstat analysis

	fastqc -o $WKDIR/QC $SAMPLE.final.bam

done

multiqc -s -o $WKDIR/QC $WKDIR/QC

# Preparation of coverage files for visualization in IGV

mkdir $WKDIR/IGV_files

for i in $WKDIR/*.final.bam
do
	samtools index $i
	SNAME=$(echo $i | sed 's:/.*/::g')
	if [ $PAIRED == 'yes' ]
	then
		bamCoverage -b $i -o $WKDIR/IGV_files/$SNAME.bw --normalizeUsing CPM -bs 5 -p $THREAD -e
	else
		bamCoverage -b $i -o $WKDIR/IGV_files/$SNAME.bw --normalizeUsing CPM -bs 5 -p $THREAD
	fi
done



mkdir $WKDIR/count
mkdir $WKDIR/diff_expr_analysis

for i in $WKDIR/*.final.bam
do
	if [ $PAIRED == 'yes' ]
	then
		htseq-count -f bam -s $STRANDED -r pos -t gene -i ID $i $FEATURES > $i.count.txt  # read count  for each gene with htseq-count
	else
		htseq-count -f bam -s $STRANDED -t gene -i ID $i $FEATURES > $i.count.txt
	fi
	mv $i.count.txt $WKDIR/count
done

for i in $WKDIR/count/*.count.txt
do
	head -n -5 $i > $i.crop.txt  # clear count files for flags
done

cp $WKDIR/count/*.crop.txt $WKDIR/diff_expr_analysis
#cp $FILES/edgeR_analysis.R $WKDIR/diff_expr_analysis
#cp $FILES/Targets.txt $WKDIR/diff_expr_analysis
