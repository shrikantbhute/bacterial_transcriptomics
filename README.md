# Bacterial Transcriptomics Pipeline with Salmon

This repository contains a high-performance workflow for quantifying bacterial gene expression. It is optimized for prokaryotic genomes (e.g., Shewanella oneidensis MR-1) and features an rRNA-aware indexing strategy to cleanly account for and isolate residual ribosomal RNA contamination from true coding mRNA sequences.

# Sequencing Compatibility

This pipeline is designed and optimized for the following data types:

- Paired-End Reads: Required for maximum structural mapping accuracy and discordant/dovetail fragment filtering.

- Illumina Platforms: Natively compatible with HiSeq, NextSeq, and NovaSeq output files.

- Automated Strandedness Detection: Salmon automatically infers library orientation (e.g., ISR, ISF), protecting downstream bacterial operon architectures.

- Total RNA-Seq (Ribosomal Depleted): Tailored for datasets treated with wet-lab rRNA depletion kits (e.g., Ribo-Zero) where variable extraction efficiency occurs.

 # 1. Reference Preparation (rRNA-Aware Indexing)
In bacterial transcriptomics, residual rRNA can dominate sequencing space if wet-lab depletion varies. Rather than discarding these or allowing them to distort coding alignments, we bundle standard coding sequences (CDS) and ribosomal features (rna-) alongside whole-genome decoys. This provides a competitive "sink" that absorbs technical noise while yielding non-ambiguous mRNA profiles.

## Step 1: Establish Your Core Genomes and Transcripts
Acquire your reference assembly (.fna) and annotated transcript target sequences from NCBI or RefSeq.

## Step 2: Build the Decoy Index
Save the following blueprint as 1_build_salmon_index.sh to generate a decoy-aware Shewanella index.

```
#!/bin/bash
#$ -cwd
#$ -o build_index.log
#$ -j y
#$ -l h_rt=01:00:00,h_data=4G,slots=4

set -e

# --- Configuration ---
GENOME="Shewanella_oneidensis_MR-1_genome.fna"
TRANSCRIPTOME="Shewanella_oneidensis_MR-1_rna.fna"
DECOY_HEADERS="decoys.txt"
INDEX_NAME="shewanella_decoy_index"

echo "Building structural decoys..."
grep "^>" "$GENOME" | cut -d ' ' -f 1 | sed 's/^>//' > "$DECOY_HEADERS"

echo "Merging transcript targets and genomic scaffolds..."
cat "$TRANSCRIPTOME" "$GENOME" > combined_ref.fna

echo "Executing Salmon Indexer..."
salmon index \
    -t combined_ref.fna \
    -d "$DECOY_HEADERS" \
    -i "$INDEX_NAME" \
    -p 4 \
    --gencode

echo "Salmon index built successfully."
```

# 2. Cluster Quantification Script
Rather than manually managing parallel execution loops over your cluster terminal, use this integrated array handler. It isolates files cleanly inside a designated output directory (quant_results/), catches missing mate pairs, and tracks stream corruption.

Save this file as 2_run_salmon_quant.sh:
```
#!/bin/bash
#$ -cwd
#$ -o salmon_batch_quant.log
#$ -j y
#$ -l h_rt=08:00:00,h_data=4G,slots=4

# --- Configuration ---
FASTQ_DIR="./raw_reads"
OUTPUT_DIR="./quant_results"
INDEX="shewanella_decoy_index"

mkdir -p "$OUTPUT_DIR"

echo "Beginning mass quantification..."

for r1_file in "${FASTQ_DIR}"/*_R1_001.fastq.gz; do
    [ -e "$r1_file" ] || continue

    # Derive matching reverse mate cleanly
    r2_file="${r1_file/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    base_name=$(basename "$r1_file" _R1_001.fastq.gz)
    sample_out_dir="${OUTPUT_DIR}/${base_name}_quant"

    if [ ! -f "$r2_file" ]; then
        echo "[ERROR] Missing reverse pair file for: ${base_name}. Skipping sample."
        continue
    fi

    echo "Quantifying Sample: ${base_name}"
    
    salmon quant \
        -i "$INDEX" \
        -l A \
        -1 "$r1_file" \
        -2 "$r2_file" \
        -p 4 \
        --validateMappings \
        -o "$sample_out_dir"
done

echo "Batch processing completed successfully."
```

# 3. Results Summary and Extraction
Because bacterial total RNA libraries exhibit variable Ribo-Zero depletion efficiency across samples, we parse the mapping metrics to verify coding depth before initiating downstream calculations in R (DESeq2 / edgeR).

This script uses a fast awk pass targeting specific structural annotation flags (^rna-) to isolate mRNA reads from residual rRNA signals.

Save this file as 3_summarize_quant_results.sh:

```
#!/bin/bash
#$ -cwd
#$ -o summary_job.log
#$ -j y
#$ -l h_rt=00:15:00,h_data=2G

RESULTS_DIR="./quant_results"
OUTPUT_FILE="mapping_summary.txt"

# Formulate column alignments
printf "%-30s %-15s %-15s %-15s %-10s\n" "Sample" "Total_Mapped" "rRNA_Reads" "mRNA_Reads" "rRNA_%" > "$OUTPUT_FILE"
printf "%-30s %-15s %-15s %-15s %-10s\n" "------------------------------" "---------------" "---------------" "---------------" "----------" >> "$OUTPUT_FILE"

if [ ! -d "$RESULTS_DIR" ]; then
    echo "[ERROR] Target output directory ${RESULTS_DIR} not found!" >> summary_job.log
    exit 1
```text
for dir in "${RESULTS_DIR}"/*_quant/; do
    [ -e "$dir" ] || continue
    quant_file="${dir}quant.sf"

    if [ -f "$quant_file" ]; then
        sample=$(basename "$dir" _quant)

        # Evaluate mapped target distributions
        stats=$(awk '
            BEGIN {total=0; rrna=0}
            NR > 1 {
                total += $5;
                if ($1 ~ /^rna-/) { rrna += $5 }
            }
            END {
                mrna = total - rrna;
                perc = (total > 0) ? (rrna / total) * 100 : 0;
                printf "%.0f %.0f %.0f %.2f", total, rrna, mrna, perc
            }' "$quant_file")

        read total rrna mrna perc <<< "$stats"
        printf "%-30s %-15s %-15s %-15s %-10s%%\n" "$sample" "$total" "$rrna" "$mrna" "$perc" >> "$OUTPUT_FILE"
    fi
done
```
# 4. Expected Outputs
Running the pipeline yields the following files within your output folders:

- quant.sf: The primary quantification matrix containing estimated transcript lengths, effective lengths, TPM (Transcripts Per Million), and raw sequence counts. Use this directly for tximport into R.

- logs/salmon_quant.log: Contains internal statistics. Review this file if your overall mapping efficiency shifts unexpectedly.

- mapping_summary.txt: A cleanly tabularized summary file providing an informative bird's-eye view of your coding depth vs. ribosomal contamination.
- 
