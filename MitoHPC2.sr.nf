workflow {

    reference_ch = channel.value(file(params.reference))
    mt_reference_ch = channel.value(file(params.mt_reference))
    mtc_reference_ch = channel.value(file(params.mtc_reference))
    numt_reference_ch = channel.value(file(params.numt_reference))
    mt_fai_ch = channel.value(file(params.mt_fai))

    reads_ch = channel
        .fromFilePairs("${projectDir}/data/*_{1,2}.fastq.gz", checkIfExists: true)

    bam_ch = ALIGN_REFERENCE(reads_ch, reference_ch)
    indexed_bam_ch = INDEX_ALIGNMENT(bam_ch)

    stats_ch = COMPUTE_ALIGNMENT_STATS(
        indexed_bam_ch.map { id, bam, bai -> tuple(id, bam) }
    )

    mt_count_ch = COMPUTE_MTDNA_COPY_NUMBER(stats_ch, params.HP_MT)
    rate_ch = CALCULATE_SUBSAMPLING_RATE(mt_count_ch)

    trimmed_ch = SUBSAMPLE_AND_TRIM(
        indexed_bam_ch
            .map { id, bam, bai -> tuple(id, bam) }
            .join(rate_ch, by: 0),
        reference_ch
    )

    mito_ch = REALIGN_CIRCULARIZED_MT(trimmed_ch, mtc_reference_ch)
    numt_ch = REALIGN_NUMTS(trimmed_ch, numt_reference_ch)

    ids_ch = COMPARE_SCORES(
        mito_ch.map { id, sam, score -> tuple(id, score) }
            .join(numt_ch, by: 0)
    )

    selected_ch = SELECT_MT_READS(
        ids_ch.join(
            mito_ch.map { id, sam, score -> tuple(id, sam) },
            by: 0
        )
    )

    mt_bam_ch = FILTER_MT_ALIGNMENTS(selected_ch, mt_fai_ch)

    COMPUTE_MT_COVERAGE(
        mt_bam_ch.map { id, bam, bai -> tuple(id, bam) },
        mt_fai_ch
    )

    IDENTIFY_SPLIT_ALIGNMENTS(
        mt_bam_ch.map { id, bam, bai -> tuple(id, bam) }
    )

    /*
     * SNV calling
     */
    snvs_ch = CALL_SNVS_MUTECT2(mt_bam_ch, mt_reference_ch)
        .mix(CALL_SNVS_MUTSERVE(mt_bam_ch, mt_reference_ch))
        .mix(CALL_SNVS_FREEBAYES(mt_bam_ch, mt_reference_ch))
        .mix(CALL_SNVS_VARSCAN(mt_bam_ch, mt_reference_ch))
        .mix(CALL_SNVS_BCFTOOLS(mt_bam_ch, mt_reference_ch))

    normalized_ch = NORMALIZE_SNVS(snvs_ch, mt_reference_ch)

    dominant_ch = IDENTIFY_DOMINANT_SNVS(normalized_ch)

    annotated_ch = ANNOTATE_SNVS(dominant_ch)

    IDENTIFY_HAPLOGROUP(annotated_ch)
    CHECK_CONTAMINATION(annotated_ch)

    sample_mt_ch = COMPUTE_SAMPLE_MT(annotated_ch, mt_reference_ch)
    INDEX_SAMPLE_MT(sample_mt_ch)
}

//////////////////////////////
// ALIGN READS

process ALIGN_REFERENCE {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path reference

    output:
    tuple val(sample_id), path("output.bam")

    script:
    """
    bwa mem \
        -v 1 \
        -t ${task.cpus} \
        -Y \
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' \
        ${reference} \
        ${reads} \
    | samtools view -bu \
    | samtools sort \
        -m ${params.sort_memory} \
        -@ ${task.cpus} \
        -o output.bam
    """
}

//////////////////////////////
// INDEX

process INDEX_ALIGNMENT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam.simpleName}.bai")

    script:
    """
    samtools index \
        -@ ${task.cpus} \
        ${bam}
    """
}

//////////////////////////////////////////
// Count reads alihned to each chr
process COMPUTE_ALIGNMENT_STATS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${bam.simpleName}.idxstats")

    script:
    """
    samtools idxstats \
        ${bam} \
        > ${bam.simpleName}.idxstats
    """
}

--- 
//////////////////////////////////
// GET MTDNA-CN; error messagel exist if very loaw counts
process COMPUTE_MTDNA_COPY_NUMBER {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(idxstats)
    val chrM

    output:
    tuple val(sample_id), path("${idxstats.simpleName}.count")

    script:
    """
    idxstats2count.pl \
        --sample ${sample_id} \
        --chrM ${chrM} \
        < ${idxstats} \
        > ${idxstats.simpleName}.count
    """
}

/////////////////////////////////
// CALCULATE_SUBSAMPLING_RATE

process CALCULATE_SUBSAMPLING_RATE {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(count)

    output:
    tuple val(sample_id), val(subsampling_rate)

    script:
    """
    subsampling_rate=""

    if [ -n "${params.HP_L}" ]; then
        subsampling_rate=\$(tail -1 ${count} | \
            perl -ane '\$rate=${params.HP_L}/(\$F[-1]+1); print \$rate if (\$rate < 1)')
    fi

    echo "subsampling_rate=\$subsampling_rate"
    """
}

/////////////////////////////////
// SUBSAMMPLE READS

process SUBSAMPLE_AND_TRIM {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), val(subsampling_rate)
    path reference

    output:
    path "${sample_id}.fq"

    script:
    def subsample = subsampling_rate ? "-s ${subsampling_rate}" : ""

    """
    samtools view \
        ${subsample} \
        ${bam} \
        ${params.HP_RMT} \
        ${params.HP_RNUMT} \
        -bu \
        -F 0x900 \
        -T ${reference} \
        -@ ${task.cpus} \
    | samtools sort \
        -n \
        -O SAM \
        -m ${params.HP_MM} \
        -@ ${task.cpus} \
    | perl -ane '
        if (/^@/) {
            print
        }
        elsif (\$P[0] eq \$F[0]) {
            print \$p, \$_
        }
        @P = @F;
        \$p = \$_;
    ' \
    | samblaster \
        ${params.HP_DOPT} \
        --addMateTags \
    | samtools view \
        -bu \
    | bedtools bamtofastq \
        -i /dev/stdin \
        -fq /dev/stdout \
        -fq2 /dev/stdout \
    | fastp \
        --stdin \
        --interleaved_in \
        --stdout \
        ${params.HP_FOPT} \
        > ${sample_id}.fq
    """
}

process REALIGN_CIRCULARIZED_MT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(fastq)
    path mtc_reference

    output:
    tuple val(sample_id),
          path("${sample_id}.sam"),
          path("${sample_id}.score")

    script:
    """
    bwa mem \
        ${mtc_reference} \
        - \
        -p \
        -v 1 \
        -t ${task.cpus} \
        -Y \
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' \
    < ${fastq} \
    | samtools view \
        -F 0x10C \
        -h \
    | tee ${sample_id}.sam \
    | samtools view \
        -bu \
    | bedtools bamtobed \
        -i /dev/stdin \
        -tag AS \
    | bed2bed.pl \
        --rmsuffix \
    | count.pl \
        -i 3 \
        -j 4 \
    | sort \
    > ${sample_id}.score
    """
}

process REALIGN_NUMTS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(fastq)
    path numt_reference

    output:
    tuple val(sample_id), path("${sample_id}.numt.score")

    script:
    """
    bwa mem \
        ${numt_reference} \
        - \
        -p \
        -v 1 \
        -t ${task.cpus} \
        -Y \
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' \
    < ${fastq} \
    | samtools view \
        -bu \
        -F 0x10C \
    | bedtools bamtobed \
        -i /dev/stdin \
        -tag AS \
    | bed2bed.pl \
        --rmsuffix \
    | count.pl \
        -i 3 \
        -j 4 \
    | sort \
    > ${sample_id}.numt.score
    """
}

process COMPARE_SCORES {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(mito_score), path(numt_score)

    output:
    tuple val(sample_id), path("${sample_id}.MT.ids")

    script:
    """
    join \
        ${mito_score} \
        ${numt_score} \
        -a 1 \
        --nocheck-order \
    | perl -ane '
        next if (@F == 3 && \$F[2] > \$F[1]);
        print join "\\t", @F;
        print "\\n";
    ' \
    > ${sample_id}.MT.ids
    """
}

process SELECT_MT_READS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(MT_ids), path(mito_sam)

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    intersectSam.pl \
        ${mito_sam} \
        ${MT_ids} \
    > ${sample_id}.sam
    """
}

process FILTER_MT_ALIGNMENTS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(sam)
    path mt_fai

    output:
    tuple val(sample_id),
          path("${sample_id}.bam"),
          path("${sample_id}.bam.bai")

    script:
    """
    cat ${sam} \
    | circSam.pl \
        --ref_len ${mt_fai} \
        --offset 0 \
    | samtools view \
        -bu \
    | samtools sort \
        -m ${params.HP_MM} \
        -@ ${task.cpus} \
        -o ${sample_id}.bam

    samtools index \
        -@ ${task.cpus} \
        ${sample_id}.bam
    """
}

process COMPUTE_MT_COVERAGE {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)
    path mt_fai

    output:
    tuple val(sample_id),
          path("${sample_id}.cvg"),
          path("${sample_id}.cvg.stat")

    script:
    """
    bedtools bamtobed \
        -cigar \
        -i ${bam} \
    | grep '^${params.HP_MT}' \
    | bedtools genomecov \
        -i - \
        -g ${mt_fai} \
    | tee ${sample_id}.cvg \
    | cut -f3 \
    | st.pl \
        --sample ${sample_id} \
    > ${sample_id}.cvg.stat
    """
}

process IDENTIFY_SPLIT_ALIGNMENTS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sa.bed")

    script:
    """
    samtools view \
        -h \
        -@ ${task.cpus} \
        ${bam} \
    | sam2bedSA.pl \
    | uniq.pl \
        -i 3 \
    | sort \
        -k2,2n \
        -k3,3n \
    > ${sample_id}.sa.bed
    """
}

process CALL_SNVS_MUTECT2 {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)
    path mt_reference
    path mtr_reference

    output:
    path "${sample_id}.vcf"

    script:
    """
    # Call SNVs against the standard mitochondrial reference
    gatk --java-options "${params.HP_JOPT}" Mutect2 \
        -R ${mt_reference} \
        -I ${bam} \
        -O ${sample_id}.orig.vcf \
        ${params.HP_GOPT} \
        --native-pair-hmm-threads ${task.cpus} \
        --callable-depth 6 \
        --max-reads-per-alignment-start 0 \
        -min-AF ${params.MINAF}

    gatk --java-options "${params.HP_JOPT}" FilterMutectCalls \
        -R ${mt_reference} \
        -V ${sample_id}.orig.vcf \
        -O ${sample_id}.filtered.vcf \
        --min-reads-per-strand 2

    mv ${sample_id}.filtered.vcf ${sample_id}.vcf
    """
}

process CALL_SNVS_MUTSERVE {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)
    path mt_reference

    output:
    path "${sample_id}.vcf"

    script:
    """
    mutserve call \
        --deletions \
        --insertions \
        --level ${params.MINAF} \
        --output ${sample_id}.vcf \
        --reference ${mt_reference} \
        ${bam}
    """
}

process CALL_SNVS_FREEBAYES {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)
    path mt_reference

    output:
    path "${sample_id}.vcf"

    script:
    """
    freebayes \
        -p 1 \
        --pooled-continuous \
        --min-alternate-fraction ${params.MINAF} \
        ${bam} \
        -f ${mt_reference} \
        > ${sample_id}.vcf
    """
}

process CALL_SNVS_VARSCAN {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)
    path mt_reference
    path mt_fai
    path varscan_vcf

    output:
    path "${sample_id}.vcf"

    script:
    """
    samtools mpileup \
        -f ${mt_reference} \
        ${bam} \
        -r ${params.HP_MT} \
        -B \
        -d ${params.MAXDP} \
    | varscan mpileup2snp \
        --min-coverage ${params.HP_DP} \
        -B \
        --variants \
        --min-var-freq ${params.MINAF} \
        --output-vcf 1 \
    > ${sample_id}.orig.vcf

    samtools mpileup \
        -f ${mt_reference} \
        ${bam} \
        -r ${params.HP_MT} \
        -B \
        -d ${params.MAXDP} \
    | varscan mpileup2indel \
        --min-coverage ${params.HP_DP} \
        -B \
        --variants \
        --min-var-freq ${params.MINAF} \
        --output-vcf 1 \
    | grep -v '^#' \
    >> ${sample_id}.orig.vcf

    cat ${varscan_vcf} > ${sample_id}.vcf

    cat ${mt_fai} \
    | perl -ane '
        print "##contig=<ID=\$F[0],length=\$F[1]>\\n"
    ' \
    >> ${sample_id}.vcf

    printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t%s\\n' \
        '${sample_id}' \
    >> ${sample_id}.vcf

    bcftools query \
        -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t.\\tGT:DP:AD:AF\\t[%GT:%DP:%AD:%FREQ]\\n' \
        ${sample_id}.orig.vcf \
    | perl -lane '
        print "$1:", int($2 * 100 + .5) / 10000 if (/(.+):(.+)%$/);
    ' \
    | sort -k2,2n \
    >> ${sample_id}.vcf
    """
}


process CALL_SNVS_BCFTOOLS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)
    path mt_reference

    output:
    path "${sample_id}.vcf"

    script:
    """
    bcftools mpileup \
        -f ${mt_reference} \
        ${bam} \
        -d ${params.MAXDP} \
    | bcftools call \
        --ploidy 2 \
        -mv \
        -Ov \
    > ${sample_id}.vcf
    """
}

process NORMALIZE_SNVS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf)
    path mt_reference

    output:
    path "${sample_id}.fix.vcf"

    script:
    """
    bcftools norm \
        -m-any \
        -f ${mt_reference} \
        ${vcf} \
    | fix${params.M}Vcf.pl \
        --file ${mt_reference} \
    | bedtools sort \
        -header \
    > ${sample_id}.fix.vcf
    """
}

process IDENTIFY_DOMINANT_SNVS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${vcf.simpleName}.max.vcf")
    tuple val(sample_id), path("${vcf.simpleName}.max.vcf.gz")
    tuple val(sample_id), path("${vcf.simpleName}.max.vcf.gz.tbi")

    script:
    """
    maxVcf.pl ${vcf} \
        | bedtools sort -header \
        | tee ${vcf.simpleName}.max.vcf \
        | bgzip -f -c > ${vcf.simpleName}.max.vcf.gz

    tabix -f ${vcf.simpleName}.max.vcf.gz
    """
}

process ANNOTATE_SNVS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf)

    output:
    path "${vcf.simpleName}.annotated.vcf"

    script:
    """
    annotateVcf.sh ${vcf}
    """
}

process IDENTIFY_HAPLOGROUP {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf)

    output:
    path "${vcf.simpleName}.haplogroup"

    script:
    """
    haplogrep classify \
        --in ${vcf} \
        --format vcf \
        --out ${vcf.simpleName}.haplogroup
    """
}

process CHECK_CONTAMINATION {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf)

    output:
    path "${vcf.simpleName}.haplocheck"

    script:
    """
    haplocheck \
        --out ${vcf.simpleName}.haplocheck \
        ${vcf}
    """
}

process COMPUTE_SAMPLE_MT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf)
    path reference

    output:
    path "${vcf.simpleName}.fa"

    script:
    """
    bcftools consensus \
        -f ${reference} \
        ${vcf} \
        -H A \
        | perl -ane '
            chomp;
            if (\$. == 1) {
                print ">\$ENV{sample_id}\\n"
            } else {
                s/N//g;
                print
            }
            END {
                print "\\n"
            }
        ' > ${vcf.simpleName}.fa
    """
}

process INDEX_SAMPLE_MT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(fasta)

    output:
    path "${fasta}.fai"
    path "${fasta.baseName}.dict"

    script:
    """
    samtools faidx ${fasta}

    gatk --java-options "${params.java_options}" \
        CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT ${fasta.baseName}.dict
    """
}

