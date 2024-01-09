# RNAseq_analysis_Cauris

RNAseq analysis pipeline for analysis of Illumina short-read paired-end sequencing data from *Candida auris* modified from **Jenull mBio (2022)** which was established by Michael Tscherner at https://github.com/kakulab/RNAseq_analysis_Cauris. This script is used in the Kuchler lab (http://cdl.univie.ac.at/) at Max Perutz Labs Vienna ([https://www.mfpl.ac.at/de.html](https://www.maxperutzlabs.ac.at/)).

*C. auris* genome sequence and annotation file (obtained from Candida Genome Database) are provided. For the current version, see CGD.

# Tools required for analysis:

samtools (http://www.htslib.org/)

FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

MultiQC (https://multiqc.info/)

cutadapt (https://cutadapt.readthedocs.io/en/stable/)

NextGenMap (https://github.com/Cibiv/NextGenMap/wiki)

DeepTools (https://deeptools.readthedocs.io/en/develop/)

HTSeq (https://htseq.readthedocs.io/en/)

All the above-mentioned tools have to be included in your PATH environment.

# Usage:

Clone the repository by typing "git clone https://github.com/kakulab/CSP2024.git" and copy the raw data into the RNAseq_analysis directory.

Clone the repository and add the read files in the base directory.

Change the adapter sequence for read trimming in the `analysis_script.sh` file if necessary. By default, it contains the Illumina TrueSeq adapter.

Change into the required_files directory and run the analysis script (by typing: `bash analysis_script.sh`).

After the pipeline has finished, change into the diff_expr_analysis directory and use the `edgeR_analysis.R` script as a basis for differential expression analysis in R.
