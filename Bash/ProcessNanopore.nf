nextflow.enable.dsl=2

workflow {
    // Define the initial input channel
    def fastqFolders = Channel.fromPath("${params.inputDir}/${params.folderPattern}*", type: 'dir')

    // Process to find the largest directory size for normalization
    def largestDirSize = findLargestDirectorySize(fastqFolders)
        .map { it.toString().toInteger() }

    // Process to merge FASTQ files and run FastQC
    def mergedFastqs = mergeFastq(fastqFolders.combine(largestDirSize))

    // Process to perform SARS-CoV-2 mapping with minimap2
    def mappedBams = mapReads(mergedFastqs)

    // Process to trim primers using TilingAmpliconParser
    def trimmedBams = trimPrimers(mappedBams)

    // Process to sort and index trimmed BAM files
    def sortedBams = sortAndIndex(trimmedBams)

    // Process to generate sequencing depth and variant frequency files
    def (depthFiles, variantTables) = generateDepthAndVariants(sortedBams)

    // Optional process to delete intermediate BAM files
    if (params.deleteBams) {
        deleteBamFiles(sortedBams)
    }

    // Optional summary report
    generateReport(depthFiles, variantTables)
}

process findLargestDirectorySize {
    input:
    path fastqFolders

    output:
    stdout

    script:
    """
    du -sb ${fastqFolders.join(' ')} | sort -nr | head -1 | cut -f1
    """
}

process mergeFastq {
    input:
    tuple path(folder), val(largestSize)

    output:
    path "${folder}/merged.fastq.gz"

    script:
    def folderName = folder.getName()
    def fastpOutput = "${folder}/merged.fastq.gz"
    
    """
    if [ "${params.normalize}" = "true" ]; then
        echo "Normalizing reads for folder ${folderName}"
        pigz -dc ${folder}/*.fastq.gz | head -n ${largestSize} | ${params.fastp} --thread 50 --stdin --stdout -q 18 --length_limit 5000 | pigz > ${fastpOutput}
    else
        pigz -dc ${folder}/*.fastq.gz | ${params.fastp} --thread 50 --stdin --stdout -q 18 --length_limit 5000 | pigz > ${fastpOutput}
    fi
    
    ${params.fastqc} -o ${folder} -t 70 ${fastpOutput}
    """
}

process mapReads {
    input:
    path mergedFastq

    output:
    path "${mergedFastq}.bam"

    script:
    """
    ${params.minimap2} --secondary-seq -a -t 90 ${params.referenceFasta} ${mergedFastq} | \
    ${params.samtools} view -bS -F 4 - | \
    ${params.samtools} sort -@ 30 -m 6G -o ${mergedFastq}.bam -
    """
}

process trimPrimers {
    input:
    path bamFile

    output:
    path "${bamFile}_Trimmed.bam"

    script:
    """
    ${params.javabin} -Xmx100G -jar ${params.TilingAmpliconParserJar} trimprimers -i ${bamFile} -b ${params.bedFiles} -f 35
    """
}

process sortAndIndex {
    input:
    path bamFile

    output:
    path "${bamFile}_Sorted.bam"

    script:
    """
    ${params.samtools} sort ${bamFile} -@ 90 -o ${bamFile}_Sorted.bam
    ${params.samtools} index ${bamFile}_Sorted.bam
    """
}

process generateDepthAndVariants {
    input:
    path sortedBam

    output:
    path "${sortedBam}.depth.tsv"
    path "${sortedBam}.ivar.tsv"

    script:
    """
    ${params.samtools} depth ${sortedBam} > ${sortedBam}.depth.tsv
    ${params.samtools} mpileup -A -d 600000 -B -Q 0 --reverse-del -f ${params.referenceFasta} ${sortedBam} > ${sortedBam}.mpileup
    ${params.javabin} -jar ${params.MpileupParserJar} -i ${sortedBam}.mpileup -o ${sortedBam}.ivar.tsv
    rm ${sortedBam}.mpileup
    """
}

process deleteBamFiles {
    input:
    path bamFile

    script:
    """
    rm ${bamFile}
    rm ${bamFile}.bai
    """
}

process generateReport {
    input:
    path depthFiles
    path variantTables

    output:
    path "mapping_summary.tsv"

    script:
    """
    echo -e "Sample\tReads\tMappedReads" > mapping_summary.tsv
    for file in ${depthFiles}; do
        sample_name=\$(basename \$file .depth.tsv)
        read_count=\$(cat \$file | wc -l)
        mapped_count=\$(cat \${variantTables} | wc -l)
        echo -e "\${sample_name}\t\${read_count}\t\${mapped_count}" >> mapping_summary.tsv
    done
    """
}