#!/bin/bash

readonly javabin=/opt/java/java/bin/java
readonly samtools=/opt/bioinfo/samtools/samtools-1.16.1/bin/samtools
readonly sarsCovRef=/data/REFERENCES/minimap/SarsCov/SARS-CoV-2.reference.fasta
readonly fastpProg=fastp
readonly minimapProg="minimap2 --secondary-seq -a  -t 90"
readonly CovidVariantFilterApp=/home/rainer/apps/wastewater/TilingAmpliconParser-1.0.jar
readonly FastpQV=18
readonly variantsTSV="/home/rainer/apps/wastewater/variant_data/Variants.tsv.txt"
readonly ampliconbedfiles=/home/rainer/apps/wastewater/Bedfiles
analysisOutBase="analysis"
mappingStatsFile=./$analysisOutBase/MappingStats.tsv
MinDirSizePctMax=12 # minimal size of a directory as fraction of the biggest directory

maxFastqLines=999999999999
maxDirSizeFound=0

cd $1
foundFastqFolders=($(find  $1 -maxdepth 1 -name $2"*" -type d ))
#---------------------------------------------------------
#find biggest dir
for folder in "${foundFastqFolders[@]}"
  do
    dirSize=$(du -sb $folder | cut -f1)
    if [[ $dirSize -gt $maxDirSizeFound ]]
      then
        maxDirSizeFound=$dirSize
      fi
  done
  
minDirSize=$(($maxDirSizeFound / $MinDirSizePctMax))

# if --normalize is set (last parameter), find minimal number of fastq lines  
#find minimal number of fastq lines
#----------------------------------------------------
if [[ "$*" == *"--normalize"* ]];then
for folder in "${foundFastqFolders[@]}"
  do
  echo "Found folder "$folder
	dirSize=$(du -sb $folder | cut -f1)	
  echo -e "counting reads in $folder \t Dirsize: $dirSize"
	if [[ $dirSize -gt $minDirSize &&  $(du -s  $folder |awk ' {print $1}') -gt 0 ]]
	  then 
        
        x=$(pigz -dc $folder/*.fastq.gz | wc -l)
        echo "found $(($x / 4)) reads"
        #c=$(($x / 4))
        if [[ "$x" -lt "$maxFastqLines" ]]
        then
          maxFastqLines=$x
        fi
  fi    
 done
 c=$(($maxFastqLines / 4))
 echo "!!!!!!!!!!!!!   Normalizing to $c reads !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
fi

if [ ! -d "$analysisOutBase" ]; then
    mkdir ./$analysisOutBase
  fi

#---------------------------------------------------------------------------------------
echo -e "Sample\tReads\tMappedReads" > $mappingStatsFile
#------------------------------------------------------------------------

for folder in "${foundFastqFolders[@]}" 
do
	dirSize=$(du -sb $folder | cut -f1)
	
	if [[ $dirSize -gt $minDirSize &&  $(du -s  $folder |awk ' {print $1}') -gt 0 ]]
	then
		folderBase=$(basename $folder)
		outDir=./$analysisOutBase/$folderBase
		echo "******************** PROCESSING $folder ---> $outDir ********************************************"
	  if [ ! -d "$outDir" ]; then
      mkdir $outDir
    fi

		echo "*****  MERGING FASTQ ***********************************************"
		pigz -dc $folder/*.fastq.gz | head -n $maxFastqLines | $fastpProg --thread 50 --stdin --stdout -q $FastpQV --length_limit 5000 | pigz > $outDir/${folderBase}.fastq.gz
    x=$(pigz -dc $outDir/${folderBase}.fastq.gz | wc -l)
    echo reads after fastp: $x
    reads=$(($x / 4))
    fastqc -o $outDir -t 70 $outDir/${folderBase}.fastq.gz

		echo " ******SARSCOV2 MAPPING ***********************************************"
			$minimapProg $sarsCovRef $outDir/${folderBase}.fastq.gz | $samtools view -bS -F 4 - | $samtools sort -@ 30 -m 6G -o $outDir/${folderBase}.bam -	&& $samtools index $outDir/${folderBase}.bam # $samtools sort -o $outDir/${folderBase}.bam - sort should not be necessary since it is sorted later
      # trim primers
      $javabin -Xmx100G -jar  $CovidVariantFilterApp  trimprimers -i $outDir/${folderBase}.bam -b $ampliconbedfiles -f 35  # -n -s # generates a <inbam>Trimmed.bam output file
		
		#remove non trimmed bam
		rm $outDir/${folderBase}.ba*
		#sort and index all bam files 
		foundbams=($(find  $outDir -name "*.bam" -type f ))
		for bam in "${foundbams[@]}" 
		do
			n=$(basename $bam .bam)
			samtools sort $bam -@ 90 -o $outDir/${n}Sorted.bam 
			rm  $bam
			#rm $bam.bai
			#samtools index $outDir/${n}Sorted.bam			
		done
    # --- create sequencing depth tsv --------------------
		$samtools depth $outDir/${folderBase}TrimmedSorted.bam > $outDir/${folderBase}.depth.tsv
    
    c=$($samtools view -c  $outDir/${folderBase}TrimmedSorted.bam)
    echo -e "$folderBase\t$reads\t$c" >> $mappingStatsFile
   
echo " ******** Generating Variation Frequency Table ***********************************************"
    $samtools mpileup  -A -d 600000 -B -Q 0   --reverse-del -f $sarsCovRef $outDir/${folderBase}TrimmedSorted.bam > $outDir/${folderBase}.mpileup
    $javabin -jar /home/rainer/apps/wastewater/MpileupParser-1.0.jar -i $outDir/${folderBase}.mpileup -o $outDir/${folderBase}.ivar.tsv
    rm $outDir/${folderBase}.mpileup

		rm $outDir/${folderBase}.fastq.gz
         
     if [[ "$*" == *"--deletebams"* ]];then
       rm $outDir/*.bam
       rm $outDir/*.bai
     fi   
  
	fi
done


 
