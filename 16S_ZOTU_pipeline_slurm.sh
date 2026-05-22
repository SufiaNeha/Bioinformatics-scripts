#!/bin/bash
#SBATCH -J zotu_pipeline
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH -p compute
#SBATCH --export=ALL

###########################################################
# 16S ZOTU Pipeline
# Author: Sufia Akter Neha
#
# Workflow:
# 1. Merge paired-end reads with PEAR
# 2. Quality filtering with USEARCH
# 3. Relabel reads by sample
# 4. Pool reads
# 5. Dereplicate sequences
# 6. Generate ZOTUs using UNOISE3
# 7. Build OTU table
# 8. Assign taxonomy with SINTAX
# 9. Multiple sequence alignment
# 10. Build phylogenetic tree
###########################################################

############################
# SOFTWARE PATHS
############################

pear="/path/to/pear"
usearch="/path/to/usearch"
sintax_db="/path/to/16S_database.fa"
fasttree="/path/to/FastTree"

############################
# CREATE OUTPUT DIRECTORIES
############################

mkdir -p stitched
mkdir -p filtered

############################
# GENERATE SAMPLE LIST
############################

ls raw | sed 's/_.*//' | sort | uniq > samples

############################
# MERGE PAIRED-END READS
############################

for s in $(cat samples)
do
    echo "Merging reads for ${s}"

    $pear \
        -f raw/${s}_R1.fastq \
        -r raw/${s}_R2.fastq \
        -o stitched/${s}
done

############################
# QUALITY FILTERING
############################

for s in $(cat samples)
do
    echo "Filtering reads for ${s}"

    $usearch \
        -fastq_filter stitched/${s}.assembled.fastq \
        -fastaout filtered/${s}_filtered.fa \
        -fastq_maxee 1.0 \
        -threads 12
done

############################
# RELABEL READS
############################

for s in $(cat samples)
do
    echo "Relabeling reads for ${s}"

    sed "s/^>\(.*\)/>\1;barcodelabel=${s};/" \
        filtered/${s}_filtered.fa > filtered/${s}.fa
done

############################
# POOL ALL READS
############################

cat filtered/*.fa > combined.fa

############################
# DEREPLICATION
############################

$usearch \
    -fastx_uniques combined.fa \
    -fastaout uniques.fa \
    -sizeout \
    -relabel Uniq \
    -threads 12

############################
# GENERATE ZOTUS
############################

$usearch \
    -unoise3 uniques.fa \
    -zotus zotus.fa

############################
# BUILD OTU TABLE
############################

$usearch \
    -otutab combined.fa \
    -zotus zotus.fa \
    -otutabout zotutab.txt \
    -mapout zmap.txt

############################
# TAXONOMY ASSIGNMENT
############################

$usearch \
    -sintax zotus.fa \
    -db $sintax_db \
    -tabbedout zotus.tax \
    -strand both \
    -threads 12

############################
# ALIGNMENT
############################

# Optional alignment step
# Uncomment if SSU-ALIGN is installed

# export PATH="$PATH:/path/to/ssu-align/bin"
# export SSUALIGNDIR="/path/to/ssu-align"

# ssu-align zotus.fa zotus.aln
# ssu-mask --pf 0.95 --pt 0.95 --afa zotus.aln
# ssu-draw zotus.aln

############################
# PHYLOGENETIC TREE
############################

# Example FastTree command
# Uncomment if alignment generated

# $fasttree \
#     -nt zotus.aln/zotus.aln.bacteria.mask.afa > tree.tre

echo "Pipeline completed successfully!"