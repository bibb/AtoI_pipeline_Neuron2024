## WGS processing

#### Read quality control -> FASTQC

#### Loading 16 fastq files, one pair per sample (2 samples, GM and GC)
```
module load fastqc

find . -name "*.fq.gz" | xargs fastqc -t 4```
```

#### Generate a single HTML report with multiqc
```
module load multiqc

multiqc .
```

#### Alingment to the reference genome (example for GC3 sample)
````
module load bwa/0.7.15 samtools spark/2.3.0
                        
bwa mem -M -t 22 -R $'@RG\tID:GC3\tSM:GC3\tLB:GC3\tPL:ILLUMINA' Homo_sapiens_assembly38.fasta <(zcat GC3_CSFP200004795-1a_H5MVWDSXY_L2_1.fq.gz) <(zcat GC3_CSFP200004795-1a_H5MVWDSXY_L2_2.fq.gz) > GC3_1.sam
 
samtools view -bS -t Homo_sapiens_assembly38.fasta.fai -@22 -o GC3_1.bam GC3_1.sam
 
samtools sort -@22 GC3_1.bam -o GC3_1_sorted.bam
 
gatk MarkDuplicates --ASSUME_SORTED TRUE --REMOVE_DUPLICATES FALSE --VALIDATION_STRINGENCY LENIENT --INPUT GC3_1_sorted.bam --OUTPUT GC3_1_sorted_unique.bam --METRICS_FILE GC3_1_picard_metrics.out
 
gatk BQSRPipelineSpark -R Homo_sapiens_assembly38.fasta -I GC3_1_sorted_unique.bam --known-sites dbsnp_146.hg38.vcf.gz -O GC3_1_sorted_unique_recal.bam -- --spark-runner LOCAL --spark-master 'local[22]'
 

        
#Remove Intermediate BAM files leaving only final processed_GATK_bam file:

                if [ -e GC3_1_sorted_unique_recal.bam ] 
                        then
                                rm GC3_1_sorted.bam
                                rm GC3_1.sam
                                rm GC3_1.bam
                                rm GC3_1_sorted_unique.bam

                                echo 'old BAM and intermediate bam.bai all removed'
                        else
                                echo '! something has gone wrong '
                        fi
````
#### Haplotype caller (example for chromosome 1)
````
gatk HaplotypeCaller -R Homo_sapiens_assembly38.fasta -L chr1 -ERC GVCF -I GC3_1_sorted_unique_recal.bam -O GC3.raw.snps.indels.chr1.updated.g.vcf
````
#### Gatger gVCF into a single gVCF
````
gatk GatherVcfs -I GC3.raw.snps.indels.chr1.updated.g.vcf \
-I GC3.raw.snps.indels.chr2.updated.g.vcf \
-I GC3.raw.snps.indels.chr3.updated.g.vcf \
-I GC3.raw.snps.indels.chr4.updated.g.vcf \
-I GC3.raw.snps.indels.chr5.updated.g.vcf \
-I GC3.raw.snps.indels.chr6.updated.g.vcf \
-I GC3.raw.snps.indels.chr7.updated.g.vcf \
-I GC3.raw.snps.indels.chr8.updated.g.vcf \
-I GC3.raw.snps.indels.chr9.updated.g.vcf \
-I GC3.raw.snps.indels.chr10.updated.g.vcf \
-I GC3.raw.snps.indels.chr11.updated.g.vcf \
-I GC3.raw.snps.indels.chr12.updated.g.vcf \
-I GC3.raw.snps.indels.chr13.updated.g.vcf \
-I GC3.raw.snps.indels.chr14.updated.g.vcf \
-I GC3.raw.snps.indels.chr15.updated.g.vcf \
-I GC3.raw.snps.indels.chr16.updated.g.vcf \
-I GC3.raw.snps.indels.chr17.updated.g.vcf \
-I GC3.raw.snps.indels.chr18.updated.g.vcf \
-I GC3.raw.snps.indels.chr19.updated.g.vcf \
-I GC3.raw.snps.indels.chr20.updated.g.vcf \
-I GC3.raw.snps.indels.chr21.updated.g.vcf \
-I GC3.raw.snps.indels.chr22.updated.g.vcf \
-I GC3.raw.snps.indels.chrX.updated.g.vcf \
-I GC3.raw.snps.indels.chrY.updated.g.vcf \
-I GC3.raw.snps.indels.chrM.updated.g.vcf  \
-O GC3.raw.snps.indels.g.v2.vcf
````
#### Merge gVCFs from the two samples with genomics DB Import. Here it is recommended to use more than 15 samples. Here is an example for chromosome 1. The mapfils needs to be a tab separated file with the first column showing the sample ID (e.g. GC3) and the second column the sample full path (e.g. /path/to/GC3.raw.snps.indels.g.v2.vcf). You need to also provide a temporary directory (tmp_dir)
```
for i in chr {1..22} X Y 
do

gatk GenomicsDBImport \
   --genomicsdb-workspace-path chr${chr}_database \
   --batch-size 15 \
   -L chr${chr} \
   --sample-name-map mapfile.txt \
   --tmp-dir ${tmp_dir} \
   --reader-threads 7
   
 done
  
````

#### Genotype GVCF from database. Here you use each chromosome database from the previous step and call raw variants in a .vcf.gz output, one for chromosome.

```
for i in chr {1..22} X Y 
do

gatk GenotypeGVCFs \
   -R Homo_sapiens_assembly38.fasta \
   -V gendb:///chr${chr}_database \
   -O chr${chr}_output.vcf.gz
done
```   

#### Concatenate and filter

```
for chr in {1..22} X Y
do

echo "chr${chr}_output.vcf.gz" >> all_chrs.txt

done


bcftools concat -n `cat all_chrs.txt` \
-Oz \
-o genomes.filtered_vars.recal.snps.indels.vcf_alleleReduction_VQRS.recode.vcf.gz \
--threads 10


bcftools view -f PASS genomes.filtered_vars.recal.snps.indels.vcf_alleleReduction_VQRS.recode.vcf -Ou \| 
bcftools filter -S . -e 'FMT/DP<10 \
| FMT/GQ<20' -Ou \
| bcftools filter -e 'F_MISSING > 0.1' -Ou \
| bcftools norm -m-both -Ou \
| bcftools norm -f hg38.fa -Ou \
| bcftools filter -e 'AC==0 \
| AC==AN' -Ou \
| bcftools view -e 'ALT="*"' \
-v 'snps,indels' -Ou \
| bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ou \
-Oz \
-o genomes.GM3_GC3.DP10_GC20_CR90_allVars.vcf.gz
```

## RNA-Seq processing

### Read quality controls -> FastQC


#### Loading 16 fastq files, one pair per sample (8 samples, M and C samples from 1 to 4)
```
module load fastqc

find . -name "*.fq.gz" | xargs fastqc -t 16```
```

#### Generate a single HTML report with multiqc
```
module load multiqc

multiqc .
```
### Alignment to the reference genome

#### Align to the hg38 reference genome using STAR

```
module load STAR/2.7.5

STAR \
--runThreadN 22 \
--genomeDir STAR_genome_hg38/ \
--sjdbGTFfile gencode.v37.annotation.gtf \
--sjdbOverhang 149 \
--readFilesIn $inputFiles \
--readFilesCommand zcat \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${nID}.
```
#### Bam filtering and recalibration

```
module load bwa/0.7.15 samtools spark/2.3.0


gatk MarkDuplicates --ASSUME_SORTED TRUE --REMOVE_DUPLICATES FALSE --VALIDATION_STRINGENCY LENIENT --INPUT Aligned.sortedByCoord.out.bam --OUTPUT Aligned_Sorted_MarkDup.bam --METRICS_FILE Aligned_MarkDup_picard_metrics.out

gatk SplitNCigarReads -R Homo_sapiens_assembly38.fasta -I Aligned_Sorted_MarkDup.bam -O Aligned_Sorted_MarkDup_NCigar.bam

gatk AddOrReplaceReadGroups -I Aligned_Sorted_MarkDup_NCigar.bam -O Aligned_Sorted_MarkDup_NCigar_wRG.bam --RGID C1 --RGLB C1 --RGPL ILLUMINA --RGSM C1 --RGPU unit1

gatk BQSRPipelineSpark -R Homo_sapiens_assembly38.fasta -I Aligned_Sorted_MarkDup_NCigar_wRG.bam --known-sites dbsnp_146.hg38.vcf.gz -O Aligned_Sorted_MarkDup_NCigar_wRG_recal.bam -- --spark-runner LOCAL --spark-master 'local[22]'

if [ -e Aligned_Sorted_MarkDup_NCigar_wRG_recal.bam ] 
                        then
                                rm Aligned.sortedByCoord.out.bam
                                rm Aligned_Sorted_MarkDup.bam
                                rm Aligned_Sorted_MarkDup_NCigar.bam
                                rm Aligned_Sorted_MarkDup_NCigar_wRG.bam

                                echo 'old BAM and intermediate bam.bai all removed'
                        else
                                echo '! something has gone wrong '
                        fi
```

#### Haplotype caller, example for chromosome 1

```
gatk HaplotypeCaller -R Homo_sapiens_assembly38.fasta \
--dont-use-soft-clipped-bases \
--standard-min-confidence-threshold-for-calling 20 \
--dbsnp dbsnp_146.hg38.vcf.gz \
-L chr1 \
-I Aligned_Sorted_MarkDup_NCigar_wRG_recal.bam \
-O raw.snps.indels.chr1.updated.allConfidentSites.vcf.gz \
--output-mode EMIT_ALL_CONFIDENT_SITES
```

#### Variant filtration

```
gatk IndexFeatureFile -I raw.snps.indels.chr1-M.updated.allConfidentSites.vcf.gz

gatk VariantFiltration --R Homo_sapiens_assembly38.fasta \
--V raw.snps.indels.chr1-M.updated.allConfidentSites.vcf.gz \
--window 35 \
--cluster 3 \
--filter-name "FS" \
--filter "FS > 30.0" \
--filter-name "QD" \
--filter "QD < 2.0" \
-O raw.snps.indels.chr1-MT.updated.allConfidentSites.HardFiltered.vcf.gz
```

#### Merge all chromosomes and merge all samples into a single .vcf.gz file

```
for i in chr{1..22} chrX chrY chrM
do
echo "C1.raw.snps.indels.${i}.updated.allConfidentSites.vcf.gz" >> vcf_input.AllSites.txt
done

bcftools concat -n `cat vcf_input.AllSites.txt` -Ou | bcftools sort -Oz -o C1.raw.snps.indels.chr1-M.updated.allConfidentSites.vcf.gz
```

#### QC per sample RNASeq VCF files

```
bcftools view -f PASS C1.raw.snps.indels.chr1-MT.updated.allConfidentSites.HardFiltered.vcf.gz \
-Ou \| 
bcftools filter -S . -e 'FMT/DP<10 \| 
FMT/GQ<20' \
-Ou \| 
bcftools norm -m-both -Ou \| 
bcftools norm -f hg38/hg38.fa -Ou \|
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' \
-Oz \
-o C1.OnlyPASS_DP10_GQ20.MultiSplit.SNPsInDels.RawID.bcftools.vcf.gz \
--threads 10

tabix -p vcf C1.OnlyPASS_DP10_GQ20.MultiSplit.SNPsInDels.RawID.bcftools.vcf.gz
```

#### Merge all QCed RNASeq VCF files into a single VCF file
```
bcftools merge -m both  M1/M1.OnlyPASS_DP10_GQ20.MultiSplit.SNPsInDels.RawID.bcftools.vcf.gz M2/M2.OnlyPASS_DP10_GQ20.MultiSplit.SNPsInDels.RawID.bcftools.vcf.gz M3/M3.OnlyPASS_DP10_GQ20.MultiSplit.SNPsInDels.RawID.bcftools.vcf.gz M4/M4.OnlyPASS_DP10_GQ20.MultiSplit.SNPsInDels.RawID.bcftools.vcf.gz C1/C1.OnlyPASS_DP10_GQ20.MultiSplit.SNPsInDels.RawID.bcftools.vcf.gz C2/C2.OnlyPASS_DP10_GQ20.MultiSplit.SNPsInDels.RawID.bcftools.vcf.gz C3/C3.OnlyPASS_DP10_GQ20.MultiSplit.SNPsInDels.RawID.bcftools.vcf.gz C4/C4.OnlyPASS_DP10_GQ20.MultiSplit.SNPsInDels.RawID.bcftools.vcf.gz  \
--threads 10 \
-Oz \
-o Merged.renamedID.DP10_GQ20.RNASeq_D60_wREfs.vcf.gz

```

#### Annotate variant IDs and filter

```
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' Merged.RNASeq_D90v2_wREfs.vcf.gz -Oz -o Merged.renamedID.RNASeq_D90v2_wREfs.vcf.gz

bcftools filter -S . -e 'FMT/DP<10 | FMT/GQ<20' Merged.renamedID.RNASeq_D90v2_wREfs.vcf.gz -Ou | bcftools filter -e 'F_MISSING > 0.1'-Oz -o RNASeq.renamedID.RNASeq_D90v2_wREfs.DP10_GQ20_CR90.vcf.gz --threads 10

```

#### Merge RNASeq with WGS data

```
bcftools merge -m all genomes.GM3_GC3.DP10_GC20_CR90_allVars.vcf.gz RNASeq.renamedID.DP10_GQ20.RNASeq_D60_wREfs.vcf.gz \
-Ou \| 
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' \
-Oz \
-o  D60.Merged.WGS_NoREf_RNASeQ_wREFs.DP10.vcf.gz \
--threads 10
```

### A-to-I analysis

#### Parse WGS+RNASeq merged VCf to plain text.
```
bcftools query -f '%ID\t%REF\t%ALT [\t%GT]\n' D60.Merged.WGS_NoREf_RNASeQ_wREFs.DP10.vcf.gz | sed 's/|/\//g' > D60.Merged.WGS_NoREf_RNASeQ_wREFs.DP10.formatted.txt
```

#### Parse WGS+RNASeq merged VCf to plain text. The following lines go through the genotype lines that start at column 4 onwards, where column 4 and 5 are the two whole-genomes, and the rest the M and C samples respestively.

#### This part checks that the reference allele ($2) is an A, and the alternative ($3) is a G, therefore an A to I site. Then it makes sure that the WGS sample at $4 is reference (0/0, or A/A), then it chacks that that the RNASeq M samples contains at least one alternative allele (/1/). If it checks that the WGS sample ($5) is an A, then it will check RNASeq samples C from 10 onwards for at least an alternative allele.

#### It also checks for A to I sites in the reverse strand, so $3 must be a T and $3 a C, and then if the M and C samples have an alternative allele. 

## Finally, it checks for G to A sites, where the WGS is now alternative (G/G) and the RNASeq samples are references (A), and also the reverse, where WGS are alternative (C/C) and the RNAseq samples are references (T/T).

```
FILENAME="D60.Merged.WGS_NoREf_RNASeQ_wREFs.DP10.formatted.txt"

awk '$2 == "A" && $3 == "G" && $4 == "0/0" {print}' $FILENAME | awk '$6 ~ /1/ || $7 ~ /1/  || $8 ~ /1/ || $9 ~ /1/ {print}' > Ms_AtoG.allPosSites.txt

awk '$2 == "A" && $3 == "G" && $5 == "0/0" {print}' $FILENAME | awk '$10 ~ /1/ || $11 ~ /1/  || $12 ~ /1/ || $13 ~ /1/ {print}' > Cs_AtoG.allPosSites.txt

awk '$2 == "T" && $3 == "C" && $4 == "0/0" {print}' $FILENAME | awk '$6 ~ /1/ || $7 ~ /1/  || $8 ~ /1/ || $9 ~ /1/ {print}' > Ms_TtoC.allPosSites.txt

awk '$2 == "T" && $3 == "C" && $5 == "0/0" {print}' $FILENAME | awk '$10 ~ /1/ || $11 ~ /1/  || $12 ~ /1/ || $13 ~ /1/ {print}' > Cs_TtoC.allPosSites.txt

awk '$2 == "G" && $3 == "A" && $4 == "1/1" {print}' $FILENAME | awk '$6 ~ /0/ || $7 ~ /0/  || $8 ~ /0/ || $9 ~ /0/ {print}' > Ms_GtoA.allPosSites.txt

awk '$2 == "G" && $3 == "A" && $5 == "1/1" {print}' $FILENAME | awk '$10 ~ /0/ || $11 ~ /0/  || $12 ~ /0/ || $13 ~ /0/ {print}' > Cs_GtoA.allPosSites.txt

awk '$2 == "C" && $3 == "T" && $4 == "1/1" {print}' $FILENAME | awk '$6 ~ /0/ || $7 ~ /0/  || $8 ~ /0/ || $9 ~ /0/ {print}' > Ms_CtoT.allPosSites.txt

awk '$2 == "C" && $3 == "T" && $5 == "1/1" {print}' $FILENAME | awk '$10 ~ /0/ || $11 ~ /0/  || $12 ~ /0/ || $13 ~ /0/ {print}' > Cs_CtoT.allPosSites.txt

# Include RefMissings

awk '$2 == "A" && $3 == "G" && $4 == "./." {print}' $FILENAME | awk '$6 ~ /1/ || $7 ~ /1/  || $8 ~ /1/ || $9 ~ /1/ {print}' > Ms_AtoG.RefMissing.allPosSites.txt

awk '$2 == "A" && $3 == "G" && $5 == "./." {print}' $FILENAME | awk '$10 ~ /1/ || $11 ~ /1/  || $12 ~ /1/ || $13 ~ /1/ {print}' > Cs_AtoG.RefMissing.allPosSites.txt

awk '$2 == "T" && $3 == "C" && $4 == "./." {print}' $FILENAME | awk '$6 ~ /1/ || $7 ~ /1/  || $8 ~ /1/ || $9 ~ /1/ {print}' > Ms_TtoC.RefMissing.allPosSites.txt

awk '$2 == "T" && $3 == "C" && $5 == "./." {print}' $FILENAME | awk '$10 ~ /1/ || $11 ~ /1/  || $12 ~ /1/ || $13 ~ /1/ {print}' > Cs_TtoC.RefMissing.allPosSites.txt
```

#### Merge files

```

cat Ms_AtoG.allPosSites.txt Cs_AtoG.allPosSites.txt Ms_TtoC.allPosSites.txt Cs_TtoC.allPosSites.txt Ms_GtoA.allPosSites.txt Cs_GtoA.allPosSites.txt Ms_CtoT.allPosSites.txt Cs_CtoT.allPosSites.txt Ms_AtoG.allPosSites.wWGS_REFS.txt Ms_TtoC.allPosSites.wWGS_REFS.txt Cs_AtoG.allPosSites.wWGS_REFS.txt Cs_TtoC.allPosSites.wWGS_REFS.txt | sort -V -k1,1 -u | grep -F -v "./." > AllPossibleSites.NoMissingness.txt


cat Ms_AtoG.RefMissing.allPosSites.txt Ms_TtoC.RefMissing.allPosSites.txt Cs_AtoG.RefMissing.allPosSites.txt Cs_TtoC.RefMissing.allPosSites.txt | sort -V -k1,1 -u | cut -f1 > RefMissing.allPosSites.IDs.txt

sed 's/_/\t/g' RefMissing.allPosSites.IDs.txt | awk '{print $1 "_" $2 "_" $3 "_" "."}' > RefMissing.allPosSites.IDs.wRefFormat.txt
```

#### Counting sites by missingness degree. Here we Selected 1 missing site.
```
D60

5053 AllPossibleSites.NoMissingness.txt
10118 AllPossibleSites.1Missing.txt
17428 AllPossibleSites.2Missing.txt
28180 AllPossibleSites.3Missing.txt
44258 AllPossibleSites.4Missing.txt
67218 AllPossibleSites.5Missing.txt
103429 AllPossibleSites.6Missing.txt
185554 AllPossibleSites.7Missing.txt
```

### Post processing
#### Extract all possible A to I sites allowing 1 missing site from the merged WGS+RNASeq VCF file and parse the output to a txt file.

```
bcftools view -i 'I=@AllPossibleSites.1Missing.txt' D60.Merged.WGS_OnlyREf_RNASeQ_wREFs.DP10.wCorrectID.vcf.gz -Ou | bcftools query -f '%ID\t%REF\t%ALT [\t%GT\t%AD\t%DP]\n' | awk '{OFS="\t"} {for(i=5;i<=NF;i=i+3) if ($i=="\.") $i="\.,\."; print}' | sed 's/,/\t/g' > AllPossibleSites.1Missing.WGS_NoRef_RNASeqYesRef.GTADDP.separated.txt
```

#### Extract selected columns
````
cut -f1,14,15,18,19,22,23,26,27,30,31,34,35,38,39,42,43 AllPossibleSites.1Missing.WGS_NoRef_RNASeqYesRef.GTADDP.separated.txt > AllPossibleSites.1Missing.wID.txt
````

#### Prepare file for Excel
````
cut -f2- AllPossibleSites.1Missing.wID.txt | sed 's/\./NA/g' > AllPossibleSites.1Missing.wID.ForCalculation.txt

````

#### Processed in Excel

````
File = AllPossibleSites.1Missing.wID.Processed.txt 9983 vars
````
#### Add ID
````
cut -f1 AllPossibleSites.1Missing.wID.txt | paste - AllPossibleSites.1Missing.wID.Processed.txt > AllPossibleSites.1Missing.wID.Processed.wID.txt

````

#### Generate allele counts
````
cut -f2-5 AllPossibleSites.1Missing.wID.Processed.wID.txt | awk '{for(i=1;i<=NF;i++) t+= $i == "NA" ; print $0 "\t" t; t=0}' | cut -f5 > counts1.txt

cut -f6- AllPossibleSites.1Missing.wID.Processed.wID.txt | awk '{for(i=1;i<=NF;i++) t+= $i == "NA" ; print $0 "\t" t; t=0}' | cut -f5 > counts2.txt

cut -f1 AllPossibleSites.1Missing.wID.Processed.wID.txt | paste - counts1.txt counts2.txt | awk '$2 < 3 && $3 < 3{print}' | cut -f1 > with2orLessNAs.txt

grep -F -wf with2orLessNAs.txt AllPossibleSites.1Missing.wID.Processed.wID.txt > AllPossibleSites.1Missing.with2orLessNAs.txt

awk '{for(i=2;i<=NF;i++) t+= $i == 1 ; print $0 "\t" t; t=0}' AllPossibleSites.1Missing.with2orLessNAs.txt | cut -f10 > ones.counts.txt

awk '{for(i=2;i<=NF;i++) t+= $i == "NA" ; print $0 "\t" t; t=0}' AllPossibleSites.1Missing.with2orLessNAs.txt | cut -f10 > NAS.counts.txt

cut -f1 AllPossibleSites.1Missing.with2orLessNAs.txt | paste - NAS.counts.txt ones.counts.txt | awk '$2 == 0 && $3 < 8{print $1}' > CerosSietes.txt

cut -f1 AllPossibleSites.1Missing.with2orLessNAs.txt | paste - NAS.counts.txt ones.counts.txt | awk '$2 == 1 && $3 < 7{print $1}' > OneSiete.txt

cut -f1 AllPossibleSites.1Missing.with2orLessNAs.txt | paste - NAS.counts.txt ones.counts.txt | awk '$2 == 2 && $3 < 6{print $1}' > TwoSeis.txt

cut -f1 AllPossibleSites.1Missing.with2orLessNAs.txt | paste - NAS.counts.txt ones.counts.txt | awk '$2 == 3 && $3 < 5{print $1}' > ThreeCinco.txt

cut -f1 AllPossibleSites.1Missing.with2orLessNAs.txt | paste - NAS.counts.txt ones.counts.txt | awk '$2 == 4 && $3 < 4{print $1}' > CuatroCuatro.txt

cat CerosSietes.txt OneSiete.txt TwoSeis.txt ThreeCinco.txt CuatroCuatro.txt | sort -V -k1,1 -u > NAs_NoConstant.txt
````

#### File for R
````
grep -F -wf NAs_NoConstant.txt AllPossibleSites.1Missing.with2orLessNAs.txt | cut -f2- > AllPossibleSites.1Missing.with2orLessNAs.fractions.NoConstant.txt
````

#### R t-test
````
data <- read.table("AllPossibleSites.1Missing.with2orLessNAs.fractions.NoConstant.txt", header = FALSE)
for (i in 1:nrow(data)) {
row <- data[i,]
dataA <- matrix(c(row$V1, row$V2, row$V3, row$V4),nrow = 1)
dataB <- matrix(c(row$V5, row$V6, row$V7, row$V8),nrow = 1)
Pvalues<-t.test(dataA, dataB, alternative="greater")$p.value
write.table(Pvalues, file="editing_sites_percentages_M-C_Pvalues.1tail-greaterM.August2021.txt", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
}
````

#### After R
````
grep -F -wf NAs_NoConstant.txt AllPossibleSites.1Missing.wID.Processed.wID.txt

grep -F -wf NAs_NoConstant.txt AllPossibleSites.1Missing.WGS_NoRef_RNASeqYesRef.GTADDP.separated.txt > AllPossibleSites.1Missing.WGS_NoRef_RNASeqYesRef.GTADDP.separated.NAs_NoConstant.txt

paste AllPossibleSites.1Missing.WGS_NoRef_RNASeqYesRef.GTADDP.separated.NAs_NoConstant.txt AllPossibleSites.1Missing.with2orLessNAs.fractions.NoConstant.txt editing_sites_percentages_M-C_Pvalues.1tail-greaterM.August2021.txt > August2021_results.D60.DP10_CR90.wFractions.wPvalue.txt
````

#### Annotate with ANNOVAR and generate final master table for analysis
````
FILE="D60.Merged.WGS_OnlyREf_RNASeQ_wREFs.DP10.wCorrectID"

/projects/b1049/genetics_programs/annovar_2017/annovar/table_annovar.pl $FILE.vcf.gz /projects/b1049/genetics_programs/annovar_2017/annovar/humandb/ -buildver hg38 \
--thread 10 \
-out $FILE.anno \
-remove \
-protocol refGene,genomicSuperDups,gnomad211_exome,gnomad30_genome \
-operation g,r,f,f \
-nastring . \
-otherinfo \
-vcfinput


cut -f1 August2021_results.D60.DP10_CR90.wFractions.wPvalue.txt | grep -F -wf - D60.Merged.WGS_OnlyREf_RNASeQ_wREFs.DP10.wCorrectID.anno.hg38_multianno.txt | cut -f1-9,11,12,29,47 > August2021_results.D60.DP10_CR90.wFractions.wPvalue.annotated.txt

paste August2021_results.D60.DP10_CR90.wFractions.wPvalue.annotated.txt August2021_results.D60.DP10_CR90.wFractions.wPvalue.txt > August2021_results.D60.DP10_CR90.wFractions.wPvalue.MasterTable.tsv
````