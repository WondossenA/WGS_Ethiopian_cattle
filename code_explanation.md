1.	FASTQC v0.11.8

```bash
/usr/local/bin/fastqc -q -t 12 *fastq.gz -o ${OUTPUTdirectory}

```

-q is for quiet mode

 -t is the number of thread

2. Trimmomatic 0.39
```bash
java -jar /usr/local/trimmomatic/0.39/trimmomatic-0.39.jar PE -threads 8 -summary ${Fastq}. summary ${Fastq_R1}.fastq ${Fastq_R2}.fastq ${Fastq_R1}-trimmed.fastq ${Fastq_R1}-singleton.fastq ${Fastq_R2}-trimmed.fastq ${Fastq_R2}-singleton.fastq MINLEN:${MINLEN} ILLUMINACLIP:${ADAPTER.fasta}:2:30:3:1:true LEADING:20 TRAILING:20 SLIDINGWINDOW:3:15 AVGQUAL:20 MINLEN:${MINLEN}

```
PE: Specifies input data is paired-end sequencing data

-threads:  Specifies the number of threads 

-summary  is the filename for the summary report

MINLEN:  Specifies the minimum length a read must be to be retained after trimming

ILLUMINACLIP: Cut adapter and other Illumina-specific sequences from the read

LEADING:20: Removes low-quality bases from the beginning of reads until a base with a quality of at least 20 is encountered.

TRAILING:20: Removes low-quality bases from the end of reads until a base with a quality of at least 20 is encountered.

SLIDINGWINDOW:3:15: Performs a sliding window trimming approach. Trims the end of a read when the average quality within a 3-base window drops below 15.

AVGQUAL:20: Removes reads with an average quality below 20.

MINLEN:  Specifies the minimum length a read must be to be retained after all trimming steps.

3. BWA 0.7.17-r1188

```bash
   
bwa mem -M -t 12 -R '@RG\tID: ${Sample}\tPL: ILLUMINA\tLB: ${Sample}\tSM:
${Sample}' ${Reference} ${Read_1. fq.gz} ${Read_2.fq.gz} | samtools view -bS - >
${Sample}. Bam

```
-M: This option marks shorter split hits as secondary (for Picard compatibility).

-t: Specifies the number of threads or CPU cores to be used during the alignment

-R '@RG\tID:${Sample}\tPL:ILLUMINA\tLB:${Sample}\tSM:${Sample}': This option adds read group information to the alignment.

ID: Sample ID

PL: Platform used

LB: Library ID

SM: Sample name

samtools view -bS -: Uses samtools to convert the SAM format output from BWA into the binary BAM  format.

4. Samtools 1.6
   
```bash

samtools sort ${Sample}.sam > ${Sample}_sorted.bam

```

```bash

samtools index ${Sample}.}_sorted.bam -@ 12

```

```bash
samtools flagstat ${Sample}_sorted.bam > ${Sample}_sorted.bam. mappingstat.txt

```
samtools sort: This command sorts a SAM/BAM file

samtools index: This command creates an index file for a BAM file

-@: Specifies the number of threads or CPU cores to be used during indexing

samtools flagstat: This command generates statistics about the flag fields of SAM/BAM records

5. Picard 2.27.4
   
```bash
picard MarkDuplicates I= ${Sample}.sorted.bam O= ${Sample}_dedup.bam
CREATE_INDEX=true M=${Sample}_dedup_metrics.txt
OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=true
VALIDATION_STRINGENCY=LENIENT

```

I: Specifies the input BAM file to be processed

O: Specifies the output BAM file where duplicate-marked reads will be written

CREATE_INDEX=true: Creates an index for the output BAM file

M: Specifies the file where Picard will write metrics (statistics) related to the duplicate marking process

OPTICAL_DUPLICATE_PIXEL_DISTANCE=100: Specifies the maximum offset between two duplicate clusters to be considered optical duplicates. 

VALIDATION_STRINGENCY=LENIENT: Sets the stringency of validation during processing. "LENIENT" means that Picard will attempt to continue processing even if it encounters potential issues.

6. GATK 4.3.0.0  ### BQSR

```bash

 gatk BaseRecalibrator \
 -I ${Sample}. dedup.bam\
 -R ${Reference}\
 --known-sites ${dbSNP}\
 -o ${Sample}_recal_data. table1

```

gatk BaseRecalibrator: Initiates the GATK tool for Base Quality Score Recalibration (BQSR)

-I: Specifies the input BAM file, which should be the BAM file after duplicate marking (deduplication). 

-R: Specifies the path to the reference genome used for the alignment and subsequent analysis.

--known-sites : Specifies known variant sites to be used in the recalibration process. 

-o: Specifies the output file where the recalibration data will be written

BQSR 2

```bash

 gatk ApplyBQSR -R ${Reference}\
 -I ${Sample}.dedup.bam \
 --bqsr-recal-file ${Sample}_recal_data. table1\
 -O ${Sample}_recal.bam

 ```

gatk ApplyBQSR: Initiates the GATK tool for applying Base Quality Score Recalibration (BQSR).

-R: Specifies the path to the reference genome used during the alignment and the initial BQSR process.

-I: Specifies the input BAM file, which should be the deduplicated BAM file.

--bqsr-recal-file: Specifies the path to the table file containing the recalibration data generated in the previous BaseRecalibrator step.

-o: Specifies the output file where the recalibration data will be written

7. HaplotypeCaller Variant calling in GVCF mode by HaplotypeCaller

```bash

gatk HaplotypeCaller
 -R ${Reference}\
 -I ${Sample}_recal.bam\
 -O ${Sample}.g.vcf.gz -ERC GVCF

 ```
gatk HaplotypeCaller: Initiates the GATK tool for variant calling, specifically for calling germline variants (SNPs and indels) in individual samples.

-R: Specifies the path to the reference genome used for alignment and variant calling.

-I: Specifies the input BAM file, which should be the recalibrated BAM file containing information about base quality score recalibration.

-O: Specifies the output file in variant call format (VCF), which contains information about the detected variants. 

-ERC GVCF: Specifies the mode for emitting reference confidence scores in the VCF. 

8. CombineGVCFs

```bash
gatk CombineGVCFs
 -R ${Reference}
 --variant  sample1.g.vcf.gz \
 --variant  sample2.g.vcf.gz \
 --variant  sample3.g.vcf.gz \
 -O All-combined.vcf.gz

 ```
gatk CombineGVCFs: Initiates the GATK tool for combining multiple genomic VCF (GVCF) files into a single GVCF.

-R: Specifies the path to the reference genome used for the variant calling.

--variant: Specifies the input GVCF files. 

-O: Specifies the output file where the combined GVCF will be written. 

9. GenotypeGVCFs

```bash

 gatk GenotypeGVCFs \
 -R Homo_sapiens_assembly38.fasta \
 -V All-combined.vcf.gz \
 -O All-genotyped.vcf.gz

 ```
gatk GenotypeGVCFs: Initiates the GATK tool for genotyping genomic VCF (GVCF) files. This step involves determining the genotype of each variant, including homozygous or heterozygous variants.

-R: Specifies the path to the reference genome used for the variant calling and genotyping.

-V: Specifies the input combined GVCF file.

-O: Specifies the output file where the genotyped VCF will be written.

 10. GATK 4.3.0 VQSR steps

```bash

java -d64 -Xmx48g -jar ${GenomeAnalysisTK.jar} VariantRecalibrator -R
${Reference}\
-V ${All-genotyped}. vcf.gz 
--resource: dbSNP, known=false, training=true, truth=true, prior=15.0${TRUEVAR}
--resource: dbSNP, known=true, training=false, truth=false, prior=2.0${KNOWNVAR}
-an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP
--tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 
-mode SNP \
-O ${All-genotyped} _recalibrate_SNP.recal 
--tranches-file ${All-genotyped}_recalibrate_SNP.tranches 
--rscript-file ${All-genotyped}_recalibrate_SNP_plots.R

 ```
-R: Specifies the path to the reference genome used for variant calling and recalibration.

-V: Specifies the input genotyped VCF file that needs to be recalibrated

--resource:dbSNP,known=false,training=true,truth=true,prior=15.0${TRUEVAR} and --resource:dbSNP,known=true,training=false,truth=false,prior=2.0${KNOWNVAR}: Specifies resources for recalibration, in this case, a resource from dbSNP. These are used for truth sensitivity and training sensitivity with associated priors.

-an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR: Specifies annotations to use during variant recalibration, such as depth, quality by depth, mapping quality rank sum test, read position rank sum test, FisherStrand, and StrandOddsRatio.

-mode SNP: Specifies the mode for variant recalibration, in this case, it's set to SNP
--tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 90.0: Specifies the tranche sensitivity thresholds for filtering variants during recalibration.

-O: Specifies the output file where the recalibrated SNP variants will be written

--tranches-file: Specifies the output file where the sensitivity tranche information will be written

--rscript-file: Specifies the output R script file for generating recalibration plots


10. Variant Selection 

```bash

 gatk SelectVariants \
 -R ${Reference}
 -V ${All-genotyped} _snp_VQSR_ApplyRecal_filtered.vcf.gz \
 -O ${SNP_filtered}. vcf.gz \
 --restrict-alleles-to BIALLELIC \
 --select-type-to-include SNP \
 --exclude-non-variants \
 --exclude-filtered \

 ```
-R: Specifies the reference genome or sequence against which the variants are being called

-V: Indicates the input variant call set. 

-O: Specifies the output VCF file where the selected variants will be written

--restrict-alleles-to BIALLELIC: This option restricts the output to contain only biallelic variants

--select-type-to-include SNP: Specifies that only SNPs (Single Nucleotide Polymorphisms) should be included in the output.

--exclude-non-variants: Excludes non-variant sites from the output

--exclude-filtered: Excludes variants that are marked as filtered










