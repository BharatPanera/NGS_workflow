# NGS Workflow
**Links to download the data:**
* https://www.ncbi.nlm.nih.gov/sra/?term=SRR24927859
* https://www.ncbi.nlm.nih.gov/sra/?term=SRR24927860

**1. Download SRA file using SRA-toolkit the command line:**
```shell
prefetch SRR24927859
```
**2. Split SRA file into R1 and R2 fastq files**
```shell
fastq-dump --split-files SRR24927859.sra
```
**3. Get the count of reads from fastq files**
```shell
grep -c "^@" SRR24927859_R1.fastq > SRR24927859_total_reads_R1.txt
grep -c "^@" SRR24927859_R2.fastq > SRR24927859_total_reads_R2.txt
```
**4. Perform basic QC checks on both files**
```shell
mkdir qc_results
fastqc SRR24927859_R1.fastq SRR24927859_R2.fastq --outdir=qc_results
```
**5. Download human genome fasta file from UCSC**
```shell
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```
**6. Index the reference file**
```shell
bwa index hg38.fa.gz
```
**7. Align both sets of fastq files onto the genome. Use BWA for alignment.**
```shell
bwa mem -t 8 hg38.fa SRR24927859_R1.fastq SRR24927859_R2.fastq | samtools view -h -bS - > SRR24927859_aligned_reads.bam
```
**8. Calculate alignment statistics and store them in a file.**
```shell
samtools sort SRR24927859_aligned_reads.bam -o SRR24927859_aligned_reads_sorted.bam
samtools flagstat SRR24927859_aligned_reads_sorted.bam > SRR24927859_alignment_statistics.txt

```
**9. Perform mark duplicates and store Mark duplicate statistics.**
```shell
picard MarkDuplicates INPUT=SRR24927859_aligned_reads_sorted.bam OUTPUT=SRR24927859_marked_duplicates.bam METRICS_FILE=mark_duplicate_stats.txt
```
**10. Variant calling.**
```shell
bcftools mpileup -Ou -f hg38.fa.gz SRR24927859_marked_duplicates.bam | bcftools call -mv -Ob -o SRR24927859.bcf
bcftools view SRR24927859.bcf | vcfutils.pl varFilter - > SRR24927859.vcf
```
**11. Separate out SNPs and Indels into individual vcf files.**
```shell
#For SNPs
samtools mpileup --skip-indels -d 250 -m 1 -E --BCF --output-tags DP,DV,DP4,SP -f hg38.fa.gz -o SRR24927859_snps.bcf SRR24927859_marked_duplicates.bam
bcftools index SRR24927859_snps.bcf ${file_name}_indexed_snps.bcf
bcftools call --skip-variants indels --multiallelic-caller --variants-only -O v SRR24927859_snps.bcf -o SRR24927859_snps.vcf

#For INDELs
samtools mpileup -d 250 -m 1 -E --BCF --output-tags DP,DV,DP4,SP -f hg38.fa.gz -o SRR24927859_indels.bcf SRR24927859_marked_duplicates.bam
bcftools index SRR24927859_indels.bcf SRR24927859_indexed_indels.bcf
bcftools call --skip-variants snps --multiallelic-caller --variants-only  -O v SRR24927859_indels.bcf -o SRR24927859_indels.vcf
```
