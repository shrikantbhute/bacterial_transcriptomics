# Bacterial Transcriptomics Pipeline with Salmon

This pipeline provides a high-performance workflow for quantifying bacterial gene expression. It is optimized for prokaryotic genomes, which are gene-dense and often contain plasmids.

## Sequencing Compatibility

This pipeline is designed and optimized for the following data types:

- Paired-End Reads: Required for maximum mapping accuracy and fragment length estimation.
- Illumina Platforms: Compatible with HiSeq, NextSeq, and MiSeq data.
- Stranded or Unstranded: Automatically detects library orientation (Sense/Antisense), which is critical for bacterial operon analysis.
- Ribosomal RNA (rRNA) Depleted: Suitable for total RNA-seq where rRNA has been removed (typically via Ribo-Zero).

# 1. Reference Preparation (Indexing)
Bacterial genomes consist mostly of coding regions (ORFs), but reads can still map to intergenic spaces. We use a "Gentrome" approach with Decoy-Aware Indexing to prevent these reads from being incorrectly counted as gene expression.

## Step 1: Prepare the "Gentrome"
The Gentrome contains all target transcripts (CDS) followed by the full genome sequence.
```
# Example: Combine Coding Sequences and Whole Genome. CDS must come FIRST so Salmon prioritizes gene mapping

cat cds_from_genomic.fna genome.fna > gentrome.fna
```

## Step 2: Extract Decoy Names
Salmon needs to know which sequences in the file are the "background" (the chromosomes/plasmids).

```
# Extract headers from the genome file to define what is a "decoy"
grep "^>" genome.fna | cut -d " " -f 1 | sed 's/>//g' > decoys.txt
```

## Step 3: Build the Index

```
# Import the Salmon path before running the command below

salmon index -t gentrome.fna -d decoys.txt -i bacterial_index -p 8

-t (Transcriptome): The combined Gentrome file.
-d (Decoys): The list of IDs Salmon should treat as genomic background.
-i (Index): The output directory for the built index.

```

# 2. Cluster Quantification Script (salmon_quant.sh)
This script is designed to run on UCLA Hoffman. It uses positional arguments to process samples in parallel across the cluster. Changes may be needed to make it work on your HPC system.

```
#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
#$ -l h_rt=24:00:00,h_data=8G
#$ -pe shared 8

# Load environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3

# --- PATH CONFIGURATION ---
# Use ABSOLUTE paths to ensure compute nodes find the files
SALMON="/u/home/j/jpjacobs/project-jpjacobs/software_rna_seq/salmon/salmon-latest_linux_x86_64/bin/salmon"
INDEX="/u/home/j/jpjacobs/project-jpjacobs/yang/shewanella_ref_genome/data/shewanella_index"

# --- EXECUTION ---
# $1 = Read 1 file, $2 = Read 2 file
$SALMON quant -i "$INDEX" \
             -l A \
             -1 "$1" \
             -2 "$2" \
             -p 8 \
             --gcBias \
             --validateMappings \
             -o "${1%.fastq.gz}_quant"


-l A (Library Type): Automatically determines if the library is stranded (ISR/ISF) or unstranded (U).
--gcBias: Bacterial genomes have high variability in GC content. This corrects for coverage biases introduced during PCR/sequencing.
--validateMappings: Uses a more rigorous alignment check to ensure reads are assigned to the correct gene, especially useful for similar genes in different operons.
-o ${1%.fastq.gz}_quant: Dynamically creates a results folder based on the R1 filename.
```

# 3. Batch Submission Loop
Run this from your data folder to submit every sample as an independent job.
```
for f in *R1_001.fastq.gz; do 
    # Extract the base name (e.g., Ag1_)
    name=$(basename "$f" R1_001.fastq.gz)
    
    # Identify the matching R2 mate
    r2="${name}R2_001.fastq.gz"
    
    echo "Submitting Sample: $name"
    
    # Pass R1 and R2 filenames to the script as $1 and $2
    qsub ../rna_scripts/salmon_quant.sh "$f" "$r2"
done
```

# 4. Expected Outputs
- quant.sf: The primary results file containing TPM (normalized) and raw counts.
- logs/: Check salmon_quant.log to see the Mapping Rate (aim for >70% in well-sequenced samples).
- lib_format_counts.json: Confirms whether your library was stranded or not.

