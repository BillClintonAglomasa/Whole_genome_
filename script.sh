#!/bin/bash

#Exit when there is an error or an undefined variable
set -eu

#Set our working directory
data="/home1/amr28/FDCO/Fastq/G4_Project"

#Check working directory
if [ ! -d "${data}" ]
then
    echo "Error: '${data}' working directory not found"
    echo "Kindly check the path you gave the data variable"
    exit 1
else
    echo "The path to your data directory is corrrect"
    echo "I will proceed with the analysis"
fi

#Creating all the directories necessary for analysis


#Function for running fastqc and multiqc
fastqc_multiqc_run () {
    local data="$1"
    local base="$2"

    fastqc "${data}/*.fastq.gz" -o "${data}/fastqc_dir"
    multiqc "${data}/fastqc_dir/*fastqc*"
}

#Function for running fastp for trimming reads
fastp_run () {
    local data="$1"
    local base="$2"

    fastp -i "${base}"_R1_001.fastq.gz -I "${base}"_R2_001.fastq.gz -l 30 --cut_front --cut_tail -W 4 -M 20 -o \
	  "${data}/fastp_dir/${base}"_1.trimmed.fastq.gz -O "${data}/fastp_dir/${base}"_2.trimmed.fastq.gz --html \
	  "${data}/fastp_dir/${base}".html --json "${data}/fastp_dir/${base}".json
}

#Function for running spades to assembly reads
spades_run () {
    local data="$1"
    local base="$2"

    spades.py -1 "${data}/${base}"_1.trimmed.fastq.gz -2 "${data}/${base}"_2.trimmed.fastq.gz --careful -t 50 -o "${data}/spades_dir/${base}"
}

#Function for checking the quality of assembled reads
quast_run () {
    local data="$1"
    local base="$2"

    quast.py "${data}/${base}/contigs.fasta" -o "${data}/quast_dir/${base}"
}

prokka_run () {
    local data="$1"
    local base="$2"

    prokka --kingdom Bacteria --cpus 50 --addgenes --prefix "${base}" --force "${data}/${base}/contigs.fasta" --outdir "${data}/prokka_dir/${base}"
}

    
#Create working fastqc if not found
if [ ! -d "${data}/fastqc_dir" ]
then
    echo "fastqc_dir not found. Creating it"
    mkdir -p  fastqc_dir
else
    echo "fastqc_dir already exist"
fi

#Running fastqc
for file in ./*.fastq.gz
do
    echo
    fastqc $file -o "${data}/fastqc_dir"
done


#Generate Multiqc report
#multiqc "${data}/fastqc_dir"/*fastqc*

#Create fastp_dir directory if not found
if [ ! -d "${data}/fastp_dir" ]
then
    echo "fastp_dir not found. Creating one"
    mkdir -p  fastp_dir
else
    echo "fastp_dir already directory"
fi

#Run the fastp command
#Selecting paired reads for analysis
for read1 in ./*_R1_001.fastq.gz
do
    if [ ! -f "${read1}" ]
    then
	echo "Error: '${read1}' not found. Kindly check the path"
	exit 1
    fi
    base=$(basename "${read1%_R1_001.fastq.gz}")
    read2="${base}_R2_001.fastq.gz"

    if [ ! -f "${read2}" ]
    then
	echo "Error '${read2}' not found"
	echo "Kindly check path to read2 or if there are errors in the base name"
	exit 1
    fi
    echo
    echo "Matching reads found for '$read1' and '$read2'"
    fastp -i "${base}"_R1_001.fastq.gz -I "${base}"_R2_001.fastq.gz -l 30 --cut_front --cut_tail -W 4 -M 20 -o \
	"${data}/fastp_dir/${base}"_1.trimmed.fastq.gz -O "${data}/fastp_dir/${base}"_2.trimmed.fastq.gz --html \
	"${data}/fastp_dir/${base}".html --json "${data}/fastp_dir/${base}".json
done

#Creating directory and running the spade command
#Creating output directory for spades analysis if absent
if [ ! -d "${data}/spades_dir" ]
then
    echo "spades_dir directory not found. Creating one"
    mkdir -p spades_dir
fi

#Selecting paired reads for analysis
for trim_read1 in "${data}/fastp_dir"/*_1.trimmed.fastq.gz
do
    if [ ! -f "${trim_read1}" ]
    then
	echo
	echo "Error: '$trim_read1' not found"
	echo "Kindly check the path to your trimmed reads"
	exit 1
    fi
    base=$(basename "${trim_read1%_1.trimmed.fastq.gz}")
    trim_read2="${data}/fastp_dir/${base}_2.trimmed.fastq.gz"
    if [ ! -f "${trim_read2}" ]
    then
	echo
        echo "Error: '${trim_read2}' not found"
	echo "Kindly check if your basename and/or the path to trim_read2 is correct"
	exit 1
    fi
    echo
    echo "Matching reads found for '${trim_read1}' and '${trim_read2}'"
    spades.py -1 "${base}"_1.trimmed.fastq.gz -2 "${base}"_2.trimmed.fastq.gz --careful -t 50 -o "${data}/spades_dir/${base}"
done

echo "Spades has been successful"

#Running quast
#Checking output_dir and creating if absent
if [ ! -d "${data}/quast_dir" ]
then
    echo "quast_dir directory not found. Creating one"
    mkdir -p quast_dir
fi

#Looping over assembled genomes to evaluate genome assembly
#run the quast
for folder in "${data}/spades_dir"/*
do
    if [ ! -d "$folder" ]
    then
	echo "Error: Folder not found. Kindly check the path to make sure you are in the right directory"
	exit 1
    fi
    echo "'$folder' found"
    base=$(basename "${folder}")
    echo "${base}"
    quast.py "${base}/contigs.fasta" -o "${data}/quast_dir/${base}"
done

#Running Prokka
#Checking output_dir and creating if absent
if [ ! -d "${data}/prokka_dir" ]
then
    echo "prokka_dir directory not found. Creating one"
    mkdir -p prokka_dir
fi

#Looping over files to run prokka
for folder in "${data}/spades_dir"/*
do
    if [ ! -d "$folder" ]
    then
	echo "Error: Folder not found. Kindly check the path to make sure you are in the right directory"
	exit 1
    fi
    echo "'$folder' found"
    base=$(basename "${folder}")
    echo "${base}"
    prokka --kingdom Bacteria --cpus 50 --addgenes --prefix "${base}" --force "${base}/contigs.fasta" --outdir "${data}/prokka_dir/${base}"
done

#Checking output_dir and creating if absent
if [ ! -d "${data}/blast_dir" ]
then
    echo "blast_dir directory nto found. Creating one"
    mkdir -p blast_dir/megablast_db
fi

#Copying fasta file contaiining resistant genes that
#will be used to create database
cp /home1/amr28/BLAST/TARGETGENES/amrgenes.fasta "${data}/blast_dir/megablast_db"

#Creating a local database
#Enter the appropriate directory to create database
cd "${data}/blast_dir/megablast_db"
#makeblastdb -dbtype nucl -in "${data}/blast_dir/megablast_db"/amrgenes.fasta

#Getting back to my working directory
cd "${data}"

#Looping over files to run blast
for folder in "${data}/spades_dir"/*
do
    if [ ! -d "$folder" ]
    then
       echo "Folder not found. Kindly check the path to the spades_dir"
    fi
    echo "'$folder' found"
    base=$(basename "${folder}")
    echo "${base}"
    blastn -query "${data}/spades_dir/${base}/contigs.fasta" -db "${data}/blast_dir/megablast_db"/amrgenes.fasta -outfmt 10 -out "${data}/blast_dir${base}".csv 
done

#Running staramr                                                                                                                                  
#Checking output_dir and creating if absent                                         
if [ ! -d "${data}/staramr_dir" ]
then
    echo "staramr_dir directory not found. Creating one"
    mkdir -p staramr_dir 
fi


#Looping over contigs to identify amr genes using staramr
#run the quast 
for folder in "${data}/spades_dir"/*
do
    if [ ! -d "$folder" ]
    then
        echo "Error: Folder not found. Kindly check the path to make sure you are in the right directory"
        exit 1
    fi
    echo "'$folder' found"
    base=$(basename "${folder}")
    echo "${base}"
    staramr search "${data}/spades_dir/${base}/contigs.fasta" -o "${data}/staramr_dir/${base}"
done


#Running abricate 
#Checking output_dir and creating if absent
if [ ! -d "${data}/abricate_dir" ]
then
    echo "abricate_dir directory not found. Creating one and its subdirectories"
    mkdir -p abricate_dir/ncbi abricate_dir/card abricate_dir/plasmidfinder abricate_dir/argannot abricate_dir/resfinder abricate_dir/vfdb
fi

#Looping over contigs to identify amr genes using abricate
#run the abricate
for folder in "${data}/spades_dir"/*
do
    if [ ! -d "$folder" ]
    then
        echo "Error: Folder not found. Kindly check the path to make sure you are in the right directory"
	exit 1
    fi
    echo "'$folder' found"
    base=$(basename "${folder}")
    echo "${base}"
    abricate "${data}/spades_dir/${base}/contigs.fasta" --db ncbi --nopath > "${data}/abricate_dir/ncbi/${base}.csv"
    abricate "${data}/spades_dir/${base}/contigs.fasta" --db card --csv > "${data}/abricate_dir/card/${base}.csv"
    abricate "${data}/spades_dir/${base}/contigs.fasta" --db plasmidfinder --csv > "${data}/abricate_dir/plasmidfinder/${base}.csv"
    abricate "${data}/spades_dir/${base}/contigs.fasta" --db argannot --csv > "${data}/abricate_dir/argannot/${base}.csv"
    abricate "${data}/spades_dir/${base}/contigs.fasta" --db resfinder --csv > "${data}/abricate_dir/resfinder/${base}.csv"
    abricate "${data}/spades_dir/${base}/contigs.fasta" --db vfdb --csv > "${data}/abricate_dir/vfdb/${base}.csv"
done
