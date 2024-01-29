# Transposable-Elements-Expression-Analysis-Pipeline
In this project, we have utilized the RNA-Seq data to detect the Transposable Elements (TEs) expression using TEtranscript package (https://academic.oup.com/bioinformatics/article/31/22/3593/240793?login=true). Here, in this project, we have analyzed the Morc3 alleles (Morc3b and Morc3ab) mutants post fertilization time points (i:e 6 hpf and 24hpf) bulk RNA-Seq data. 
### Required Tools 
1) &nbsp; Pre-processing <br />
2) &nbsp; Genome alignment with STAR (with Multimappers) <br />
3) &nbsp; Post-Alignment formatting and filtering with SAMtools <br />
4) &nbsp; TE expression identification with TEtranscript (pre-installed at Sapelo2 GACRC)
5) &nbsp; DESeq2 R-Package to detect the differentially expressed TEs (Integrated installed in TEtranscript )

# General steps for data processing 

### 1) Pre-processing: 
Trim_Glore was used for the pre-processing of raw fastq files.
```
mkdir trim
```
```
trim_galore --fastqc -j 24 --output_dir ./trim --paired read_1.fastq read_2.fastq 
```
### 2) Genome alignment with STAR
--> Reference Genome fasta seq retrieval
```
mkdir ref
```

```
 curl -s ftp://ftp.ensembl.org/pub/release-98/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz | gunzip -c > danio_refseq.fa
```
--> Genome Index building 
```
  STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./ref --genomeFastaFiles ./ref/genome.fa
```
--> Genome Mapping
```
mkdir bams
```
Repeat this step for wild and mutant samples
```
STAR --runThreadN 24 --genomeDir ./ref/ --outFileNamePrefix ./bams/ \
      --readFilesCommand zcat --readFilesIn wild_read_val_1.fastq wild_read_val_2.fastq  --outSAMtype BAM SortedByCoordinate \
      --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000

STAR --runThreadN 24 --genomeDir ./ref/ --outFileNamePrefix ./bams/ \
      --readFilesCommand zcat --readFilesIn mutant_read_val_1.fastq mutant_read_val_2.fastq  --outSAMtype BAM SortedByCoordinate \
      --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
```      
### 3) Post alignment formatting and filtering with SAMtools

```
samtools view -bS -h -bq1 wild_Aligned.sortedByCoord.out.bam | samtools sort - > wild_filter_sorted.bam

samtools view -bS -h -bq1 mutant_Aligned.sortedByCoord.out.bam | samtools sort - > mutant_filter_sorted.bam
```
### 4) File conversion from bam to bigwig for visualization 
```
samtools index -M wild_filter_sorted.bam -> wild_filter_sorted.bam.bai
samtools index -M mutant_filter_sorted.bam -> mutant_filter_sorted.bam.bai
```
#### IGV visualization 
```
bamCoverage -b wild_filter_sorted.bam -o wild_filter.bw
bamCoverage -b mutant_filter_sorted.bam -o mutant_filter.bw


