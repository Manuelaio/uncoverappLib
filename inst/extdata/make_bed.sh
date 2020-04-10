#!bin/bash


 . $1

Rscript --vanilla $pathscript/preprocessing.R $geneList $genome $output $notation_bam


mkdir $output/$(date +%Y_%m_%d_%H:%M:%S)

samtools depth -b $output/teMpFoldeR/*.bed -f $bamList -Q -q | awk 'BEGIN { FS=" "; OFS=" " } { $2=$2 "\t" $2 } 1' | sed -e 's/\(^[0-9XY]\)/chr\1/' -e 's/^MT/chrM/'| bgzip  > $output/$(date +%Y_%m_%d_%H:%M:%S)/multisample.bed.gz


rm -rf teMpFoldeR
