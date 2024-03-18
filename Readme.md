# Workflow for analyzing Nanopore SARS-Cov-2 amplicon sequencing data from wastewater samples

$indir: directory with input fastq files

$outdir: directory for output

$sarsCovRef: SARS-Cov-2 reference genome fasta file

## 1. Concatenate fastq files

cat \$indir/*.fastq.gz >  $outdir/merged.fastq.gz
 
## 2. Map to SARS-Cov-2 genome

    minimap2 -a  -t 90 $sarsCovRef $outdir/merged.fastq.gz | samtools view -bS -F 4 - | samtools sort -@ 30 -m 6G -o $outDir/tmp.bam -	&& $samtools index $outDir/tmp.bam 

## 3. Trim amplification primers

    java -jar  <path to jar>/CovidVariantFilter-1.0.jar  trimprimers -i $outDir/tmp.bam -b <path to folder with amplicification primer bed files> -f 20  -o $outDir/outbam.bam

### Parameters
**--inbam, -i** input bam file
**--outbam, -o** output bam file
**--fuzzyness, -f** fuzzyness for amplicon ends allowed. e.g. -f 20 will accept bam records that start matching +/- 20 nt.  from expected pos
**--bedfiles, -b** directory with bed files defining the amplification primers. Typically should contain one bed file per primer panel
