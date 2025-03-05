
# Workflow for analyzing Nanopore SARS-Cov-2 amplicon sequencing data from wastewater samples

## Authors
Pauline Steichen <steichen@ipmc.cnrs.fr>

Rainer Waldmann <waldmann@ipmc.cnrs.fr>

# Mapping, generation of variation frequency tables (each sample)
$indir: directory with input fastq files

$outdir: directory for output

$sarsCovRef: SARS-Cov-2 reference genome fasta file

## 1. Concatenate fastq files 

cat \$indir/*.fastq.gz >  $outdir/merged.fastq.gz
 
## 2. Map to SARS-Cov-2 genome

    minimap2 -a  -t 90 $sarsCovRef $outdir/merged.fastq.gz | samtools view -bS -F 4 - | samtools sort -@ 30 -m 6G -o $outDir/tmp.bam -	&& $samtools index $outDir/tmp.bam 

## 3. Trim amplification primers
custom java softwarer [TilingAmpliconParser](Java/TilingAmpliconParser/TilingAmpliconParser-1.0.jar) 

    java -jar  <path to jar>/TilingAmpliconParser-1.0.jar  trimprimers -i $outDir/tmp.bam -b <[path to folder with amplicification primer bed files](AmpliconPanels/)> -f 20  -o $outDir/outbam.bam

generates an output bam file where primers are soft-clipped in the SAM records. Only SAM records consistent with an amplicon defined in the bed files are written to the output bam file.
Also generates a tsv file *.TrimStats.tsv with the number of SAM record matching each amplicon.


### Parameters 
(required parameters are in bold)

**--inbam, -i** input bam file

**--outbam, -o** output bam file

**--fuzzyness, -f** fuzzyness for amplicon ends allowed. e.g. -f 20 will accept sam records that start and end matching +/- 20 nt.  from expected pos

**optional: --writesplitbams, -s : ** writes one additional output bam for each amplicon panel, optional flag

**optional; --writeNonMatching, -n : ** writes an additional output bam for non assigned records, Use this flag to find errors in amplicon bed files. optional flag

**--bedfiles, -b** directory with bed files defining the amplification primers. Typically should contain one bed file per primer panel. Columns are tab delimited

Only files with the extension .bed are considered.

Bed files for Artic, Varskip and Spike primer panels are in the [AmpliconPanels](AmpliconPanels) directory


Example for two amplicons:

    MN908947.3	25	50	SARS-CoV-2_1_LEFT	1	+	AACAAACCAACCAACTTTCGATCTC
    MN908947.3	408	431	SARS-CoV-2_1_RIGHT	1	-	CTTCTACTAAGCCACAAGTGCCA
    MN908947.3	324	344	SARS-CoV-2_2_LEFT	2	+	TTTACAGGTTCGCGACGTGC
    MN908947.3	705	727	SARS-CoV-2_2_RIGHT	2	-	ATAAGGATCAGTGCCAAGCTCG

    1st Column: Name of reference genome

    2nd Column: Most 5' pos on the reference genome

    3rd Column: Most 3' pos on the reference genome

    4th Column: Name of the primer. must end with _<amplicon number>_<LEFT if forward primer RIGHT if reverse primer>. There must be no other "_" in the primer name. E.g.  SARS-CoV-2_2_LEFT means : Amplicon 2, forward primer

    5th Column: primer pool (typically 1 or 2)

    6th column: strand , + for forward primer, - for reverse primer

    7th column: primer sequence

    A row with a forward primer is  followed by a  row with a reverse primer. A row with a reverse primer is  followed by a  row with a forward primer.** If there are e.g. two optional reverse primers for an amplicon (e.g. older versions of the Artic panel) the forward primer has to be repeated with another name.



## 4. Generate variation frequency data

generate a TSV file with frequencies for each variation

**mpileup** generates a text pileup file

    samtools mpileup -A -d 600000 -B -Q 0 --reverse-del -f $sarsCovRef $outDir/outbam.bam > $outDir/out.mpileup


**mpileupparser** parses the pileup file and generates a variation frequency table

 [MpileupParser](Java/MpileupParser/MpileupParser-1.0.jar)

    java -jar /home/rainer/apps/wastewater/MpileupParser-1.0.jar -i $outDir/out.mpileup -o $outDir/ivar.tsv

an alternate option is the iVar software https://github.com/andersen-lab/ivar. However, we use our custom software since we filter out most of the background variations due to low quality reads by filtering out variations which are only backed by low quality reads or have a skewed forward reverse balance.  iVar does not provide those informations for indels.

## 4. Generate sequencing depth data

    samtools depth $outDir/outbam.bam > $outDir/depth.tsv

# Variant Analysis 

**Python script parses folders in --inputDir , The folder name is used as the sample name** 

## Usage:
[main.py](Python/CagableaParser/main.py)

    > python main.py --inputDir [path to folder that contains one folder per sample]


**Each folder within the specified folder corresponds to one sample and must contain:**

 - a variation frequency file with a name ending with 'ivar.tsv'      
 - a sequencing depth data file with a name ending with 'depth.tsv'
 
 ## Comments:
 - ignores folders that do not contain a variation frequency file with a name ending with 'ivar.tsv'
- ignores folders with names starting with '_'.
- if removeFirstCharsFromSampleName = True in [settings.py](Python/CagableaParser/settings.py), the first characters from the sample name until the first underscore at max pos 5 in the name will be removed. 
- Used to define sorting order in outputs
e.g.  a sample (folder name) AAA_samplexxxx  AAB_sampleaaaa . AAA_ and AAB_ will be removed in sample names in graphs, reports ... samplexxxx will appear before sampleaaaa in the output  

## Output
**Output will be written to the folder provided as --inputDir**

Output includes:
- a html report with heatmaps, coverage plots, variant frequencies
- a tsv file with variation frequencies, quality values .... : VarFrequencies.tsv . Column naming is similar to that used by the ivar software
- a tsv file with just variation frequencies: VarFrequenciesLight.tsv
- a tsv file with sequencing depth data for all positions "CoverageAllPositions.tsv"
 - a tsv file with sequencing depth data for every tenth position "CoveragePositionsReduced.tsv"
- a tsv file with mean coverage values: coverageMeans.tsv
- optionally pdf files with heat maps, pie charts .... See [settings.py](Python/CagableaParser/settings.py) for options 


uses multithreading as much as this is possible with Python

## Data filtering

First, the files from the different folders (samples) are read and the variation frequency data are filtered
The following variations are filtered out:
- Variation frequencies are set to 0 if the the mutation is both supported by a read QV in the given position of below 17 and the QV of the basecall supporting the variation is at least 7 below the QV of the reads supporting the wild-type sequence. Can be modified by changing conditionALT_QUAL in [settings.py](CagableaParserV0.3/settings.py)
- Variations with an imbalanced forward/reverse read ratio of at least 3.5 . Optional, active if settings.do_filter_FWD_REV_balance=True in [settings.py](Python/CagableaParser/settings.py). Max Imbalance factor is set with settings.maxALT_FwdRevImbalance.

The data are merged into one dataframe.

## Variant Frequencies

The Mutation frequencies for each lineage's characteristic mutations [Variants.tsv.txt](CagableaParserV0.3/Data/Variants.tsv.txt) are extracted from the variation frequency tables for each sample. The frequency of each lineage is set to the mean of its variant-defining mutation frequencies. Outliers were identified and excluded using an interquartile range (IQR)-based filter, ignoring values exceeding 0.1 times the IQR below the first quartile or above the third quartile as outliers. Can be adjusted with interq_range_outliers. To prevent filtering values near the mean when the IQR was small, values deviating less than 5% from the mean are preserved. 

### Variant definition file  
[Variants.tsv.txt](Python/CagableaParser/Data/Variants.tsv.txt), a tsv file containing variant defining mutations , hierarchy and instructions for the graphs (colors ...) 

**Lines of Variants.tsv.txt:**

- Variant: Name of the variant  
- variant_altname: Alternate name, can be empty. Currently not used by script.  
- comment: Optional, comment which is not used by script.  
- parent: The parent lineage (must be a lineage in this table) or empty if no parent lineage or if parent not in table.  
- excempt from child sum: Defaults to false. Means the variant is ignored when sums of all child lineages of a parent lineage are calculated. Is typically set true for 'helper' lineages such as BA.4_BA.5 in the table. BA.4 and BA.5 don't have a parent child lineage relationship, however all mutations of BA.5 are found in BA.4. We are using BA.4_BA.5 as a virtual helper lineage, which is not printed in the output, that contains mutations common in BA.4 and BA.5.  
- virtual parent: BA.4_BA.5 has the mutations of BA.2. For frequency calculations BA.2 is used as a virtual parent of BA.4_BA.5  
- remove counts from virtual parent: True for BA.4_BA.5 since they share all BA.2 mutations and BA4_BA5 frequency needs to be substracted from BA.2 frequency  
- min muts for pass: Variant is only considered present if the minimal number of mutations are present after outlier filtering.  
- min * muts for pass: Mutations can be marked as important by preceeding the mutation with a '*'. The number here defines the number of those * - flagged mutatants required for presence of the variant 
- calcstrategy: BA.5 freq is calculated by substracting the BA.4 frequency from the frequency of the virtual BA.4_BA.5 lineage. The instruction is provided in the following format: BA.4_BA.5:-:BA.4. Means substract BA.4 from BA.4_BA.5; ':' is just a delimiter 
- print: Whether variant should be printed in pie charts, histograms ... Is e.g. set false for the virtual lineage BA.4_BA.5  
- histogram group: Instructions on how histogram bars with variant frequencies are organized. E.g. delta/0 (name/[index]) means group named delta is first bar in histogram, omicron/1 -> second histogram bar (the omicron sublineages that have the omicron/1 group will be shown in this histogram bar).  The order how the variants are shown within a bar can be defined by the second index. E.G. XBB/3-1 -- XBB/3-2 -- XBB/3-3 ... 
- color: A CSS color name or hex color, the color used for the variant in histograms and pie charts  
- mutations: Mutations characterizing a variant.  Substitutions: e.g. C16466T , deletions: e.g. 28362del9 means nine nucleotides are deleted , first deleted position is 28362. Optionally a '/' followed by the amino acid mutation can be added for better lisibility of the table - won't be used by the script. Insertions are currently not treated

# Identification of multiple mutations in one read

[TilingAmpliconParser-1.0.jar](Java/TilingAmpliconParser/TilingAmpliconParser-1.0.jar)

    java -jar  TilingAmpliconParser-1.0.jar  splitvariants -inbam \<input bam file\> --outDir \<output directory\> --variantData  \<TSV file with info on mutations\> --minmuts 5 --fractionmutsrequired 0.8

    
will parse the -inbam and write results into --outDir
Output:
Stats.tsv : number of hits for each variant
One bam file per variant if the --writeBams option is provided

### Parameters

- --variantData, -v :
    tsv file that contains:

    a) first row: names of variants starting in the second column

    b) any row below in the first coulumn the word "mutations"

    c) for each variant the mutations that should be searched for. Same format as in [Variant definition file](Python/CagableaParser/Data/Variants.tsv.txt) used in the Python script for lineage deconvolution. The same file can be used. Just delete column for variants that should not be searched for. Only the row with the variant names and the rows starting with the row "mutations" are read by the program

- --inbam, -i : the input bam file to parse
- --outDir, -o : the output directory
- --minmuts, -m : the minimal number of variant specific mutations that need to be potentially present in one read
- --fractionmutsrequired, -f : The fraction of the mutations defined with --minmuts that have to be found in a read.
- --writeBams, -w (Optional) No parameters: If present, will write a bam file with matching records for each variant. 

