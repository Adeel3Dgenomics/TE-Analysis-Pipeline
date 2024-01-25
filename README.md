# Transposable-Elements-Expression-Analysis-Pipeline
In this project, we have utilized the RNA-Seq data to detect the Transposable Elements (TEs) expression using TEtranscript package (The Gale Hammell Lab). Here, in this project, we have analyzed the Morc3 alleles (Morc3b and Morc3ab) mutants post fertilization time points (i:e 6 hpf and 24hpf) bulk RNA-Seq data. 
### Required Tools 
1) &nbsp; Pre-processing <br />
2) &nbsp; Genome alignment with STAR  <br />
3) &nbsp; Post-Alignment formatting and filtering with SAMtools <br />
4) &nbsp; TE expression identification with TEtranscript (pre-installed at Sapelo2 GACRC)
5) &nbsp; DESeq2 R-Package to detect the differentially expressed TEs (Integrated installed in TEtranscript )


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
Repeat this step for control and treated samples
```
STAR --runThreadN 24 --genomeDir ./ref/ --outFileNamePrefix ./bams/ \
      --readFilesCommand zcat --readFilesIn read_val_1.fastq read_val_2.fastq  --outSAMtype BAM SortedByCoordinate \
      --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
```
### 3) Post alignment formatting and filtering with SAMtools
