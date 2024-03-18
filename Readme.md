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
<<<<<<< HEAD
**--bedfiles, -b** directory with bed files defining the amplification primers. Typically should contain one bed file per primer panel
MN908947.3	25	50	SARS-CoV-2_1_LEFT	1	+	AACAAACCAACCAACTTTCGATCTC
MN908947.3	408	431	SARS-CoV-2_1_RIGHT	1	-	CTTCTACTAAGCCACAAGTGCCA
MN908947.3	324	344	SARS-CoV-2_2_LEFT	2	+	TTTACAGGTTCGCGACGTGC
MN908947.3	705	727	SARS-CoV-2_2_RIGHT	2	-	ATAAGGATCAGTGCCAAGCTCG
1st Column: Name of reference genome
2nd Column: Position of the most 5' nucleotide of the primer on the reference genome
3rd Column: Name of the primer. must end with _\<amplicon number>_<LEFT if forward primer RIGHT if reverse primer>. There must be no other "_" in the primer name. E.g.  SARS-CoV-2_2_LEFT means : Amplicon 2, forward primer
A row with a forward primer is  followed by a  row with a reverse primer. A row with a reverse primer is  followed by a  row with a forward primer.
Only files with the extension .bed are considered.
Example bed files are in the AmpliconPanels directory

