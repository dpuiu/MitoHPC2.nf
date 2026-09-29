# MitoHPC2 Short Read Nextflow Pipeline

`MitoHPC2.sr.nf` is a **Nextflow DSL2 workflow for mitochondrial DNA (mtDNA) analysis from paired-end short-read sequencing data**.

The workflow takes paired-end FASTQ files, aligns them to a whole-genome reference, identifies mitochondrial reads while distinguishing mitochondrial sequences from NUMTs, generates a mitochondrial BAM, calls variants using multiple callers, normalizes and selects variants, and produces downstream mitochondrial annotations and a sample-specific mitochondrial sequence.

## Workflow Overview

```text
Paired-end FASTQ
       │
       ▼
ALIGN_REFERENCE
       │
       ▼
INDEX_ALIGNMENT
       │
       ├──────────────► COMPUTE_ALIGNMENT_STATS
       │                         │
       │                         ▼
       │                 COMPUTE_MTDNA_COPY_NUMBER
       │                         │
       │                         ▼
       │                 CALCULATE_SUBSAMPLING_RATE
       │                         │
       ▼                         │
SUBSAMPLE_AND_TRIM ◄─────────────┘
       │
       ├──────────────► REALIGN_CIRCULARIZED_MT ──┐
       │                                           │
       └──────────────► REALIGN_NUMTS ─────────────┤
                                                   ▼
                                           COMPARE_SCORES
                                                   │
                                                   ▼
                                           SELECT_MT_READS
                                                   │
                                                   ▼
                                        FILTER_MT_ALIGNMENTS
                                                   │
                              ┌────────────────────┼───────────────────┐
                              │                    │                   │
                              ▼                    ▼                   ▼
                    COMPUTE_MT_COVERAGE   IDENTIFY_SPLIT_ALIGNMENTS  SNV calling
                                                                          │
                                         ┌────────────────────────────────┤
                                         │        │       │       │       │
                                         ▼        ▼       ▼       ▼       ▼
                                      Mutect2  Mutserve FreeBayes VarScan BCFtools
                                         │        │       │       │       │
                                         └────────┴───────┴───────┴───────┘
                                                           │
                                                           ▼
                                                   NORMALIZE_SNVS
                                                           │
                                                           ▼
                                               IDENTIFY_DOMINANT_SNVS
                                                           │
                                                           ▼
                                                     ANNOTATE_SNVS
                                                           │
                              ┌────────────────────────────┼─────────────────────┐
                              ▼                            ▼                     ▼
                       IDENTIFY_HAPLOGROUP          CHECK_CONTAMINATION   COMPUTE_SAMPLE_MT
                                                                                │
                                                                                ▼
                                                                         INDEX_SAMPLE_MT
```

## Input Data

The workflow expects paired-end FASTQ files in the `data/` directory:

```text
data/
├── SAMPLE1_1.fastq.gz
├── SAMPLE1_2.fastq.gz
├── SAMPLE2_1.fastq.gz
└── SAMPLE2_2.fastq.gz
```

Files are discovered using:

```nextflow
channel.fromFilePairs("${projectDir}/data/*_{1,2}.fastq.gz")
```

The filename prefix becomes the sample ID.

For example:

```text
data/NA12878_1.fastq.gz
data/NA12878_2.fastq.gz
```

produces the sample ID:

```text
NA12878
```

## Reference Files

The workflow uses several reference files:

| Parameter               | Purpose                                                                        |
| ----------------------- | ------------------------------------------------------------------------------ |
| `params.reference`      | Whole-genome reference used for initial alignment                              |
| `params.mt_reference`   | Standard mitochondrial reference used for variant calling                      |
| `params.mtc_reference`  | Circularized mitochondrial reference used to improve mtDNA read identification |
| `params.numt_reference` | NUMT reference used to identify nuclear mitochondrial sequences                |
| `params.mt_fai`         | mtDNA reference index/length information                                       |

The whole-genome alignment is performed first, followed by a targeted analysis of reads potentially originating from mtDNA.

## Main Processing Steps

### 1. Whole-genome alignment

`ALIGN_REFERENCE`

Paired-end reads are aligned to the whole-genome reference using `bwa mem`.

The output is coordinate-sorted BAM:

```text
SAMPLE.bam
```

Read groups are added using the sample ID.

### 2. BAM indexing

`INDEX_ALIGNMENT`

The BAM is indexed using `samtools index`:

```text
SAMPLE.bam
SAMPLE.bai
```

### 3. Alignment statistics

`COMPUTE_ALIGNMENT_STATS`

`samtools idxstats` is used to count reads aligned to each reference sequence.

Output:

```text
SAMPLE.idxstats
```

### 4. Estimate mitochondrial copy number

`COMPUTE_MTDNA_COPY_NUMBER`

The alignment statistics are processed to estimate mitochondrial DNA copy number.

Output:

```text
SAMPLE.idxstats.count
```

### 5. Calculate subsampling rate

`CALCULATE_SUBSAMPLING_RATE`

The estimated mtDNA copy number is used to determine whether the reads should be subsampled before the targeted mtDNA analysis.

The calculation uses the configured `HP_L` parameter.

### 6. Subsample and trim reads

`SUBSAMPLE_AND_TRIM`

Reads are extracted from the whole-genome BAM and processed through:

* `samtools view`
* `samtools sort`
* mate/read processing
* `samblaster`
* `bedtools bamtofastq`
* `fastp`

The resulting interleaved FASTQ file is:

```text
SAMPLE.fq
```

### 7. Realign against circularized mtDNA

`REALIGN_CIRCULARIZED_MT`

The extracted reads are aligned against a circularized mitochondrial reference.

The process generates:

```text
SAMPLE.sam
SAMPLE.score
```

The alignment score information is used to distinguish mtDNA-derived reads from reads that may originate from NUMTs.

### 8. Realign against NUMTs

`REALIGN_NUMTS`

The same extracted reads are aligned against a NUMT reference.

Output:

```text
SAMPLE.numt.score
```

### 9. Compare mtDNA and NUMT scores

`COMPARE_SCORES`

The mitochondrial and NUMT alignment scores are compared to identify reads more consistent with mtDNA.

Output:

```text
SAMPLE.MT.ids
```

### 10. Select mitochondrial reads

`SELECT_MT_READS`

The selected read IDs are intersected with the mitochondrial alignments.

Output:

```text
SAMPLE.sam
```

### 11. Generate filtered mitochondrial BAM

`FILTER_MT_ALIGNMENTS`

The selected mitochondrial alignments are converted and sorted into a BAM file.

The process also handles circular mtDNA coordinates.

Outputs:

```text
SAMPLE.bam
SAMPLE.bam.bai
```

This BAM is the primary input for downstream mtDNA variant analysis.

## Mitochondrial Coverage

`COMPUTE_MT_COVERAGE`

Mitochondrial alignment coverage is calculated using `bedtools`.

Outputs:

```text
SAMPLE.cvg
SAMPLE.cvg.stat
```

The coverage statistics can be used to assess sequencing depth and coverage distribution across the mitochondrial genome.

## Split Alignments

`IDENTIFY_SPLIT_ALIGNMENTS`

Split and supplementary alignments are extracted from the mitochondrial BAM.

Output:

```text
SAMPLE.sa.bed
```

This provides information useful for identifying reads with split alignment patterns.

## SNV Calling

The workflow uses multiple variant callers:

* **GATK Mutect2**
* **Mutserve**
* **FreeBayes**
* **VarScan**
* **BCFtools**

All callers operate on the filtered mitochondrial BAM.

```text
                    mtDNA BAM
                       │
        ┌──────────────┼──────────────┐
        ▼              ▼              ▼
     Mutect2        Mutserve       FreeBayes
        │              │              │
        ├──────────────┼──────────────┤
        ▼              ▼
      VarScan       BCFtools
        │              │
        └──────────────┴──────────────►
                       │
                       ▼
                  NORMALIZE_SNVS
```

### Mutect2

`CALL_SNVS_MUTECT2`

Calls mitochondrial variants with GATK Mutect2 and applies `FilterMutectCalls`.

### Mutserve

`CALL_SNVS_MUTSERVE`

Calls mitochondrial variants using Mutserve, including deletion and insertion detection.

### FreeBayes

`CALL_SNVS_FREEBAYES`

Calls variants using FreeBayes in pooled-continuous mode.

### VarScan

`CALL_SNVS_VARSCAN`

Uses `samtools mpileup` followed by VarScan for SNV and indel detection.

### BCFtools

`CALL_SNVS_BCFTOOLS`

Uses the BCFtools `mpileup`/`call` workflow.

## Variant Normalization

`NORMALIZE_SNVS`

Variant representations from the different callers are normalized using:

```text
bcftools norm
```

followed by the pipeline's custom VCF normalization utility.

Output:

```text
SAMPLE.fix.vcf
```

The goal is to bring variant representations into a consistent format before downstream comparison and selection.

## Dominant Variant Selection

`IDENTIFY_DOMINANT_SNVS`

The normalized VCF is processed to identify the dominant SNVs.

Outputs include:

```text
SAMPLE.fix.max.vcf
SAMPLE.fix.max.vcf.gz
SAMPLE.fix.max.vcf.gz.tbi
```

The compressed VCF is indexed with `tabix`.

## Variant Annotation

`ANNOTATE_SNVS`

The selected variants are annotated using:

```text
annotateVcf.sh
```

Output:

```text
SAMPLE.fix.max.annotated.vcf
```

## Haplogroup Assignment

`IDENTIFY_HAPLOGROUP`

Mitochondrial haplogroup classification is performed with Haplogrep.

Output:

```text
SAMPLE.fix.max.annotated.haplogroup
```

## Contamination Check

`CHECK_CONTAMINATION`

Mitochondrial contamination is assessed using Haplocheck.

Output:

```text
SAMPLE.fix.max.annotated.haplocheck
```

## Sample-Specific Mitochondrial Sequence

`COMPUTE_SAMPLE_MT`

A sample-specific mitochondrial consensus sequence is generated using:

```text
bcftools consensus
```

Output:

```text
SAMPLE.fix.max.annotated.fa
```

## Index the Sample-Specific Reference

`INDEX_SAMPLE_MT`

The generated FASTA is indexed with:

```text
samtools faidx
```

and a GATK sequence dictionary is generated.

Outputs:

```text
SAMPLE.fix.max.annotated.fa.fai
SAMPLE.fix.max.annotated.dict
```

## External Tools

The workflow invokes the following command-line tools and custom utilities:

### Alignment and BAM processing

* BWA
* samtools
* bedtools
* samblaster
* fastp

### Variant calling

* GATK
* Mutserve
* FreeBayes
* VarScan
* BCFtools

### Mitochondrial analysis utilities

* `idxstats2count.pl`
* `bed2bed.pl`
* `count.pl`
* `intersectSam.pl`
* `circSam.pl`
* `st.pl`
* `sam2bedSA.pl`
* `uniq.pl`
* `fix*Vcf.pl`
* `maxVcf.pl`
* `annotateVcf.sh`

### Mitochondrial interpretation

* Haplogrep
* Haplocheck

The custom Perl and shell utilities must be available in the execution environment.

## Configuration Parameters

The workflow references parameters controlling references, thresholds, memory, Java options, and tool-specific options.

Important parameters include:

```text
reference
mt_reference
mtc_reference
numt_reference
mt_fai

HP_MT
HP_L
HP_RMT
HP_RNUMT
HP_MM
HP_DOPT
HP_FOPT
HP_JOPT
HP_GOPT
HP_DP

MINAF
MAXDP
M
java_options
sort_memory
```

These parameters should be defined through the project's Nextflow configuration or command-line parameters.

For example:

```bash
nextflow run MitoHPC2.part.nf \
    --reference reference.fa \
    --mt_reference rCRS.fa \
    --mtc_reference circular_mt.fa \
    --numt_reference numt.fa \
    --mt_fai rCRS.fa.fai
```

The exact parameter values and paths depend on the reference datasets and execution environment.

## Running the Workflow

From the repository directory:

```bash
nextflow run MitoHPC2.part.nf
```

For a configuration file:

```bash
nextflow run MitoHPC2.part.nf -c nextflow.config
```

For a specific work directory:

```bash
nextflow run MitoHPC2.part.nf \
    -c nextflow.config \
    -work-dir work/
```

Before running, make sure that:

1. Paired FASTQ files are present under `data/`.
2. All reference files are available.
3. Reference indexes have been generated where required.
4. The external command-line tools are available.
5. The custom MitoHPC utilities are in `PATH`.
6. Required Nextflow parameters are defined.

## Workflow Outputs

The major output products are:

| Output            | Description                                |
| ----------------- | ------------------------------------------ |
| `*.bam`           | Whole-genome and filtered mtDNA alignments |
| `*.bai`           | BAM indexes                                |
| `*.idxstats`      | Alignment statistics                       |
| `*.count`         | mtDNA copy-number/count information        |
| `*.fq`            | Extracted/interleaved reads                |
| `*.score`         | Circularized mtDNA alignment scores        |
| `*.numt.score`    | NUMT alignment scores                      |
| `*.MT.ids`        | Selected mitochondrial read IDs            |
| `*.cvg`           | mtDNA coverage                             |
| `*.cvg.stat`      | Coverage statistics                        |
| `*.sa.bed`        | Split-alignment information                |
| `*.vcf`           | Variant calls                              |
| `*.fix.vcf`       | Normalized variants                        |
| `*.max.vcf`       | Selected/dominant variants                 |
| `*.annotated.vcf` | Annotated variants                         |
| `*.haplogroup`    | Haplogroup assignment                      |
| `*.haplocheck`    | Contamination analysis                     |
| `*.fa`            | Sample-specific mitochondrial FASTA        |
| `*.fa.fai`        | FASTA index                                |
| `*.dict`          | GATK sequence dictionary                   |

## Pipeline Structure

The workflow is organized into the following major stages:

```text
1. Whole-genome alignment
2. Alignment statistics
3. mtDNA copy-number estimation
4. Read subsampling
5. mtDNA / NUMT realignment
6. mtDNA read selection
7. mtDNA BAM generation
8. Coverage and split-alignment analysis
9. Multi-caller SNV detection
10. Variant normalization
11. Dominant variant selection
12. Variant annotation
13. Haplogroup classification
14. Contamination assessment
15. Sample-specific mtDNA consensus generation
```

## Notes

`MitoHPC2.part.nf` is a **pipeline component rather than a completely self-contained software distribution**. In particular, it depends on project-specific scripts and externally installed bioinformatics tools.

The workflow also assumes that the required `params.*` values are supplied through the surrounding project configuration.

For the current implementation, see:

[MitoHPC3 — MitoHPC2.part.nf](https://github.com/dpuiu/MitoHPC3/blob/main/MitoHPC2.part.nf?utm_source=chatgpt.com)

