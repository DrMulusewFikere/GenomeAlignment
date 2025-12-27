# GenomeAlignment
RNA-seq Alignment and Variant Calling Pipeline (SLURM)
# Overview
GenomeAlignment is a SLURM-based RNA-seq alignment and variant calling pipeline designed for high-performance computing (HPC) environments. The pipeline processes pre-trimmed RNA-seq reads, aligns them to a reference genome using STAR, and performs RNA-seq–aware variant calling using GATK SplitNCigarReads and bcftools.
The workflow is automated to process multiple sequencing runs and samples with minimal manual intervention.
# Pipeline Workflow
For each RNA-seq sample, the pipeline performs:
## Genome Alignment
Aligns reads to a reference genome using STAR and produces coordinate-sorted BAM files with gene-level read counts.
## RNA-seq BAM Processing
Uses GATK SplitNCigarReads to properly handle spliced alignments prior to variant calling.
## Variant Calling
Calls SNPs and indels using bcftools mpileup and bcftools call, producing compressed and indexed VCF files.
# Requirements
## System
Linux HPC environment
SLURM workload manager
## Software
STAR 2.7.11b
samtools 1.22.1
bcftools 1.17
GATK 4.6.0.0

Each Run* directory is processed automatically
FASTQ files must be pre-trimmed (e.g., using fastp)
Single-end reads are assumed
Reference Genome and Index
The STAR genome index must be generated before running the pipeline.
# Running the Pipeline
## Submit the job using SLURM:
sbatch rnaseq_variant_calling.sh
The pipeline will:
Loop through all Run* directories
Process all *_trimmed.fastq.gz files
Skip samples where a .vcf.gz already exists
## Output Files
For each sample, the pipeline generates:
*_Aligned.sortedByCoord.out.bam — aligned reads
*.split.bam — RNA-seq processed BAM
*.vcf.gz — variant calls
*.vcf.gz.tbi — VCF index
## STAR gene count files
SLURM log files:
pipeline.out
pipeline.err
## Notes
Input FASTQ files must be trimmed before running the pipeline
Variant calls are unfiltered; downstream filtering is recommended
Designed for RNA-seq variant discovery, not DNA-seq
Modify STAR parameters for paired-end data if required
