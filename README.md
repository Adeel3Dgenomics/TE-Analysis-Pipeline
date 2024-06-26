# Transposable Elements Expression Analysis
Bulk RNA-Seq data was used to detect the Transposable Elements (TEs) expression using TEtranscript package [Publication](https://academic.oup.com/bioinformatics/article/31/22/3593/240793?login=true). This software package is originally build by [Hammell Lab](https://hammelllab.labsites.cshl.edu/software/#TEtranscripts) and available at their lab [Github](https://github.com/mhammell-laboratory/TEtranscripts).This software calculates the TE and gene expression using RNA-Seq data, with the integrated DESeq2 R-Package it can calculates the differential expression of TE/Genes and give output in refined tabular form. Here, we have used test RNA-Seq data of Morc3 mutant and wild alleles.  
## General Steps and tools 
1) &nbsp; Pre-processing <br />
2) &nbsp; Genome alignment with _STAR_ (with Multimappers) <br />
3) &nbsp; Post-Alignment formatting and filtering with _SAMtools_ <br />
4) &nbsp; TE expression identification with _TEtranscript_ (pre-installed at Sapelo2 GACRC)
5) &nbsp; _DESeq2_ R-Package to detect the differentially expressed TEs (Integrated installed in TEtranscript )

# General steps for data processing 

### 1) Pre processing: 
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
#### batch script for Pre-processing and Mapping at HPC (Parallel Mode)  (Alternative method) 

```
sbatch mapping.sh
```
### 4) File conversion from bam to bigwig for visualization 
```
samtools index -M wild_filter_sorted.bam -> wild_filter_sorted.bam.bai
samtools index -M mutant_filter_sorted.bam -> mutant_filter_sorted.bam.bai
```
--> IGV visualization 
```
bamCoverage -b wild_filter_sorted.bam -o wild_filter.bw
bamCoverage -b mutant_filter_sorted.bam -o mutant_filter.bw
```
### 5) TE counts using TEtranscipt (Method 1) 
```
TEtranscripts --format BAM -t mutant_rep1_filter_sorted.bam \
mutant_rep2_filter_sorted.bam mutant_rep3_filter_sorted.bam \
-c wild_rep1_filter_sorted.bam \
wild_rep2_filter_sorted.bam wild_rep3_filter_sorted.bam \
--GTF refann.gtf \
--TE TE_annotations_DEclusters_new.gtf \
--mode multi --project TE_out --minread 1 -i 10 --padj 0.05 --sortByPos
```
 -t  &nbsp;   treated/mutant sample bam file <br />
-c   &nbsp;   control/wild sample bam file <br />
--GTF&nbsp;   gene annotation file corresponding to the reference genome assembly (refann.gtf) <br />
--TE &nbsp;   TE annotatation (download from [here](https://www.dropbox.com/sh/1ppg2e0fbc64bqw/AACUXf-TA1rnBIjvykMH2Lcia?dl=0) )<br />

#### batch script for TE counts at HPC (Parallel Mode) (Alternative method)

```
sbatch TE_analysis.sh
```

#### Test output 
The output contains read counts table (*.cnTable), Expressed Genes (*_gene_TE_analysis.txt), and Significant Differentially expressed Genes table (*sigdiff_gene_TE_Only.txt).   
```
"ENSDARG00000000001"    54      35      48      43      41      31
"ENSDARG00000000002"    0       2       1       0       0       0
"ENSDARG00000000018"    776     793     490     410     490     621
"ENSDARG00000000019"    624     726     740     1074    992     719
"ENSDARG00000000068"    608     683     662     503     613     631
"ENSDARG00000000069"    556     598     581     532     575     500
"ENSDARG00000000086"    279     226     256     220     308     134
"ENSDARG00000000103"    1306    1270    1366    1230    1522    1031
"ENSDARG00000000142"    137     116     140     119     123     96
"ENSDARG00000000151"    6       8       3       2       21      4
"ENSDARG00000000161"    2       3       3       2       3       6
"ENSDARG00000000175"    0       2       0       0       0       0
piggyBac-N2_DR:PiggyBac:DNA	167	264	135	147	95	75
piggyBac-N3_DR:PiggyBac:DNA	240	289	217	123	74	84
piggyBac-N4_DR:PiggyBac:DNA	121	159	125	88	50	45
piggyBac-N5B_DR:PiggyBac:DNA	130	140	124	103	69	62
piggyBac-N5C_DR:PiggyBac:DNA	56	41	49	56	22	28
piggyBac-N5_DR:PiggyBac:DNA	69	74	71	51	32	35
piggyBac-N6_DR:PiggyBac:DNA	17	25	18	7	5	5
piggyBac-N7_DR:PiggyBac:DNA	343	444	319	249	175	171
piggyBac-N8_DR:PiggyBac:DNA	231	264	188	142	84	92
piggyBac-N9_DR:PiggyBac:DNA	207	195	185	128	94	75
```
#### Columns names
1)  Gene/TE
2)  mutant_rep1_filter_sorted.bam.T
3)  mutant_rep2_filter_sorted.bam.T
4)  mutant_rep3_filter_sorted.bam.T
5)  wild_rep1_filter_sorted.bam.C
6)  wild_rep2_filter_sorted.bam.C
7)  wild_rep3_filter_sorted.bam.C

#### Note: 
This tool is designed to handle both gene counts and TE counts, leading to the inclusion of read counts from both genes and TEs in the output table "cntTable." However, our primary objective is to exclusively process TEs. Therefore, we will exclude the gene counts and retain only the read counts associated with TEs.
### 6) TE counts filtering 
```
grep -v "ENS*" out.cntTable > TEs_only_count.cntTable
```
### 7) TE counts using FeatureCount (Method 2)

Another very handy tool 'FeatureCount' can be use to quantify the read counts of TEs. It uses the TE Annotations and calcultaes the read counts. 
```
featureCounts -s 0 -p -M -a TE_ann.gtf  -o count.txt mutant_rep1_filter_sorted.bam mutant_rep2_filter_sorted.bam mutant_rep3_filter_sorted.bam wild_rep1_filter_sorted.bam wild_rep2_filter_sorted.bam wild_rep3_filter_sorted.bam
```
Next step,  is to filter out the irrelavent extraAttributes as follows. 

```
awk '{print $1"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' count.txt  > final_count.txt
```
Use the 'final_count.txt' in down stream processing. 

### 8) Differentially Expressed TEs detection 

```
Rscript Diff_TE_DESeq2.R
```
Arrange the tables, merge the count informationa and already available 'Reference Table' data. Filter out the 'Up and down-regulated TEs' based on your criteria. We used FC > 2 and P-value < 0.05.  






