#!/usr/bin/env nextflow
// This file defines individual processes (separated for portability)

// taking input bam files for sorting and indexing
process "preprocessing" {

    label "high"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    input:
    tuple sample, path(bam)
    // eg. [sample, /path/to/sample.bam]
    path fasta

    output:
    tuple sample, path("sorted.bam")
    // eg. [sample, /path/to/sorted.bam]

    script:
    """
    samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
    -o sorted.bam ${bam} || exit \$?
    """
}


// eg. perform samtools calmd + index
process "process1" {

    label "low"
    label "finish"
    tag "$sample - $type"

    maxForks "${params.fork}".toInteger()

    input:
    tuple type, sample, path(bam)
    // eg. [process1, sample, /path/to/sorted.bam]
    path fasta
    path fai

    output:
    tuple type, sample, path("${sample}.${type}.bam"), path("${sample}.${type}.bam.bai")
    // eg. [process1, sample, /path/to/process1.bam, /path/to/process1.bam.bai]

    script:
    """
    samtools calmd -b sorted.bam ${fasta} 1> ${sample}.${type}.bam 2> /dev/null
    samtools index ${sample}.${type}.bam
    """
}


// eg. perform samtools index + make an additional file
process "process2" {

    label "low"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    input:
    tuple type, sample, path(bam)
    // eg. [process2, sample, /path/to/sorted.bam]

    output:
    tuple type, sample, path("${sample}.${type}.bam"), path("${sample}.${type}.bam.bai")
    // eg. [process2, sample, /path/to/process2.bam, /path/to/process2.bam.bai]
    path "${sample}.txt"

    // a directive when you only want the process to execute under certain conditions
    when:
    params.process2

    script:
    """
    cp sorted.bam ${sample}.${type}.bam
    samtools index ${sample}.${type}.bam
    touch ${sample}.txt
    """
}

// eg. perform samtools fastq
process "process3" {

    label "low"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    input:
    tuple type, sample, path(bam), path(bai)
    // eg. [process2, sample, /path/to/process2.bam, /path/to/process2.bam.bai]

    output:
    tuple sample, path("*.fastq.gz")
    // eg. [sample, /path/to/type.fastq.gz]
    // eg. [sample, [/path/to/type_1.fastq.gz, /path/to/type_2.fastq.gz]]

    script:
    if (params.PE)
        """
        samtools fastq -c 6 -1 ${sample}.${type}_1.fastq.gz -2 ${sample}.${type}_2.fastq.gz -0 /dev/null -s /dev/null -n ${bam}
        """
    else
        """
        samtools fastq -c 6 -0 /dev/null ${bam} > ${sample}.${type}.fastq.gz
        """
}
