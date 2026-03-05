# Bacterial Transcriptomics Pipeline with Salmon

This repository contains a high-performance workflow for quantifying bacterial gene expression. It is optimized for prokaryotic genomes (e.g., Shewanella oneidensis MR-1) and specifically addresses removal of rRNA content from the libraries.

## Sequencing Compatibility

This pipeline is designed and optimized for the following data types:

- Paired-End Reads: Required for maximum mapping accuracy and fragment length estimation.
- Illumina Platforms: Compatible with HiSeq, NextSeq, and MiSeq data.
- Stranded or Unstranded: Automatically detects library orientation (Sense/Antisense), which is critical for bacterial operon analysis.
- Ribosomal RNA (rRNA) Depleted: Suitable for total RNA-seq where rRNA has been removed (typically via Ribo-Zero).

# 1. Reference Preparation (rRNA-aware Indexing)
In many bacterial samples, rRNA can make up >90% of the reads. Instead of filtering them out with general databases, we extract the isolate-specific rRNA and include them in the Salmon index as targets. This "sinks" the rRNA reads, providing an accurate mapping rate and cleaner mRNA data.


## Step 1: Install gffread
If not available on your cluster, download the standalone binary to your project space:


```
cd /u/home/j/jpjacobs/project-jpjacobs/software_rna_seq/

wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.7.Linux_x86_64.tar.gz

tar -xvzf gffread-0.12.7.Linux_x86_64.tar.gz

```

## Step 2: Extract rRNA Sequences from the genome file.

```
# Path: /u/home/j/jpjacobs/project-jpjacobs/yang/shewanella_ref_genome/data/GCF_000146165.2/
/u/home/j/jpjacobs/project-jpjacobs/software_rna_seq/gffread-0.12.7.Linux_x86_64/gffread \
-w rrna.fa \
-g GCF_000146165.2_ASM14616v2_genomic.fna \
-r rRNA \
genomic.gff
```


## Step 3: Create the "Gentrome" & Decoys
The Gentrome must contain all transcripts (CDS + rRNA) followed by the whole genome.


```
# Combine sequences (Order: CDS -> rRNA -> Genome)
cat cds_from_genomic.fna rrna.fa GCF_000146165.2_ASM14616v2_genomic.fna > gentrome_plus_rrna.fna

# Extract headers for the decoy list (Genome only)
grep "^>" GCF_000146165.2_ASM14616v2_genomic.fna | cut -d " " -f 1 | sed 's/>//g' > decoys.txt

```

## Step 4: Build the Index
```
/u/home/j/jpjacobs/project-jpjacobs/software_rna_seq/salmon/salmon-latest_linux_x86_64/bin/salmon index \
-t gentrome_plus_rrna.fna \
-d decoys.txt \
-i shewanella_rrna_index \
-p 8
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
# --- ABSOLUTE PATHS ---
SALMON="/u/home/j/jpjacobs/project-jpjacobs/software_rna_seq/salmon/salmon-latest_linux_x86_64/bin/salmon"
INDEX="/u/home/j/jpjacobs/project-jpjacobs/yang/shewanella_ref_genome/data/GCF_000146165.2/shewanella_rrna_index"

# --- EXECUTION ---
# $1 = R1 file, $2 = R2 file
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

