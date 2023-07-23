#Usage: sh ngs_pipeline.sh


# https://www.ncbi.nlm.nih.gov/sra/?term=SRR24927859
# https://www.ncbi.nlm.nih.gov/sra/?term=SRR24927860


file_name="Sample2"
fastq_files_path="/path/to/fastq_files/"
ref_file_path="/path/to/ref/hg38/hg38.fa"

mkdir ${file_name}
cd ${file_name}

#Download SRA files from NCBI
#prefetch SRR24927859

#fastq-dump --split-files SRR24927859.sra

get total reads
grep -c "^@" ${fastq_files_path}/${file_name}_R1.fastq > ${file_name}_total_reads_R1.txt
grep -c "^@" ${fastq_files_path}/${file_name}_R2.fastq > ${file_name}_total_reads_R2.txt
echo -e "\n################################ Read count is completed ################################\n"

# get total bases
total_bases_r1=$(awk '{s+=$NF} END {print s}' "${fastq_files_path}/${file_name}_R1.fastq") > ${file_name}_total_base_R1.txt
total_bases_r2=$(awk '{s+=$NF} END {print s}' "${fastq_files_path}/${file_name}_R2.fastq") > ${file_name}_total_base_R2.txt
echo -e "\n################################ Base count is completed ################################\n"

#fastqc
mkdir qc_results
fastqc ${fastq_files_path}/${file_name}_R1.fastq ${fastq_files_path}/${file_name}_R2.fastq --outdir=qc_results
echo -e "\n################################ Fastqc is completed ################################\n"

#adapter trimming:
cutadapt -q 20 -m 36 -o ${file_name}_trimmed_output_R1.fastq -p ${file_name}_trimmed_output_R2.fastq ${fastq_files_path}/${file_name}_R1.fastq ${fastq_files_path}/${file_name}_R2.fastq
echo -e "\n################################ Adapter trimming is completed ################################\n"

#Download genome fasta file from UCSC
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip ${ref_file_path}
echo -e"\n################################ Reference genome is downloaded ################################\n"

#Reference index
bwa index ${ref_file_path}

#Alignment
# bwa mem -t 8 ${ref_file_path} R1.fastq R2.fastq > ${file_name}_aligned_reads.sam
# bwa mem -t 8 ${ref_file_path} ${fastq_files_path}/${file_name}_R1.fastq ${fastq_files_path}/${file_name}_R2.fastq > ${file_name}_aligned_reads.sam
bwa mem -t 8 ${ref_file_path} ${fastq_files_path}/${file_name}_R1.fastq ${fastq_files_path}/${file_name}_R2.fastq | samtools view -h -bS - > ${file_name}_aligned_reads.bam
echo -e"\n################################ Alignment is completed ################################\n"

#sort bam file:
samtools sort ${file_name}_aligned_reads.bam -o ${file_name}_aligned_reads_sorted.bam
echo -e"\n################################ .bam file is sorted ################################\n"

#get alignment stats
samtools flagstat ${file_name}_aligned_reads_sorted.bam > ${file_name}_alignment_statistics.txt
echo -e"\n################################ Alignment stats are generated ################################\n"

#mark duplicates
. /home/bharat/anaconda3/bin/activate
picard MarkDuplicates INPUT=${file_name}_aligned_reads_sorted.bam OUTPUT=${file_name}_marked_duplicates.bam METRICS_FILE=mark_duplicate_stats.txt
conda deactivate
echo -e"\n################################ Duplicate marking is completed ################################\n"

#variant calling

#command for GATk
# gatk HaplotypeCaller -R genome_fasta/corona.fa -I ${file_name}_marked_duplicates.bam -O ${file_name}_raw_variants.vcf --native-pair-hmm-threads 8

#command for bcftools
bcftools mpileup -Ou -f ${ref_file_path} ${file_name}_marked_duplicates.bam | bcftools call -mv -Ob -o ${file_name}.bcf
bcftools view ${file_name}.bcf | vcfutils.pl varFilter - > ${file_name}.vcf
echo -e"\n################################ vcf file is generated ################################\n"

#fetching SNPs from vcf
samtools mpileup --skip-indels -d 250 -m 1 -E --BCF --output-tags DP,DV,DP4,SP -f ${ref_file_path} -o ${file_name}_snps.bcf ${file_name}_marked_duplicates.bam
bcftools index ${file_name}_snps.bcf ${file_name}_indexed_snps.bcf
bcftools call --skip-variants indels --multiallelic-caller --variants-only  -O v ${file_name}_snps.bcf -o ${file_name}_snps.vcf
echo -e"\n################################ SNPs are saved in ${file_name}_snps.vcf ################################\n"

#fetching INDELS from vcf
samtools mpileup -d 250 -m 1 -E --BCF --output-tags DP,DV,DP4,SP -f ${ref_file_path} -o ${file_name}_indels.bcf ${file_name}_marked_duplicates.bam
bcftools index ${file_name}_indels.bcf ${file_name}_indexed_indels.bcf
bcftools call --skip-variants snps --multiallelic-caller --variants-only  -O v ${file_name}_indels.bcf -o ${file_name}_indels.vcf
echo -e"\n################################ Indels are saved in ${file_name}_indels.vcf ################################\n"

