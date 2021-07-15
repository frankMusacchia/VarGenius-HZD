# VarGenius-HZD

Rare homozygous/hemizygous deletions caller within a set of BAM files

----------------------------------------

Please, clone or download the entire folder of VarGenius-HZD from this repository. 
Installation is not required, but you need some dependencies...

**Dependencies**

R

PERL (Needed libraries: Parallel::ForkManager; Getopt::Long

bedtools

samtools 

once those are installed, prepare your input:


**Input to prepare**

- genome in FASTA format that you used for the alignment
- the genome file for your genome. You can produce it with:
	1. samtools faidx ucsc.hg19.fa (that produces ucsc.hg19.fa.fai)
	2.  cut -f1,2 ucsc.hg19.fasta.fai | sort -k1,1 -k2,2n  > ucsc.hg19.genomefile

	VarGenius-HZD uses bedtools coverage with the -sorted option to accelerate the process. This assumes that the records are sorted lexicographically (e.g., chr1, chr10, etc.). 
	Thus, you should take a look at your BAM files to check the chromosomes order. You could do this with the following command:
	samtools view -H sample.bam | grep SQ | cut -f 2 | awk '{ sub(/^SN:/, ""); print;}'

	For example 1KGP samples are sorted "numerically" (e.g., 1, 2, etc.).  If so, you need to provide the tool with a "genome file" (-g option) defining the expected order. 
	You can find more information here:
	http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html#g-define-an-alternate-chromosome-sort-order-via-a-genome-file
	We did include a genomefile within the data folder which was used for the 1KGP dataset.
	 
- the target BED file for your sequenced samples
- a file with a list of sample names and paths to the BAM files tab separated (bamlist.txt). 
  E.g.:

```

samplename	path
sampleA	/path/to/SampleA.bam
sampleB	/path/to/SampleB.bam

```

Now you should be able to run VarGenius-HZD pre processing and algorith:

**Step 1: Get the exons on-target (to be run only once per-target)**

perl HZD_launch.pl -f GET\_EXONS\_ON\_TARGET -p /path/where/you/downloaded/thistool/ -o /path/to/outfolder/ -ref /path/to/genome.fasta -b //path/to/bedtools -l /path/to/bamlist.txt  -t /path/to/target.bed

**Step 2: Execute the breadth/depth of coverage step for each sample**

We suggest to run Step 2 and Step 3 within a screen or in background as they will require a while to finish.

perl HZD_launch.pl -f HZD\_PREPROCESSING -p /path/where/you/downloaded/thistool/ -o /path/to/outfolder/ -ref /path/to/genome.fasta -b //path/to/bedtools -l /path/to/bamlist.txt -g /path/to/genomefile  -t /path/to/target.bed

**Step 3. Launch HZD algorithm**

Here the -m parameter (max_lowBoC) allows to define the number of samples which might have a specific exon HD.
Default is 2, Please change accordingly to your case (e.g. siblings in the dataset) 

perl HZD_launch.pl -f DETECT\_HDs -p /folder/where/you/downloaded/thistool/ -o /path/to/outfolder/ -ref /path/to/genome.fasta -b //path/to/bedtools -l /path/to/bamlist.txt  -t /path/to/target.bed -m 2


**OUTPUT Description**

*VarGenius-HZD* outputs a TAB separated file with the following fields:

compid: an identifier of the HD

samplename: (self-explained)

pBoC: sample Breadth-of-coverage

pDoC: sample Depth-of-coverage
fDoC: father Depth-of-coverage (only in VarGenius)

mDoC: mather Depth-of-coverage (only in VarGenius)

avgDoC: average Depth-of-coverage across samples

cnv_type: type of CNV

gene: UCSC gene symbol




**TEST WITH SYNTHETIC HDs in 1KGP DATA**

This is the test that we conducted on 1KGP data to detect synthetic deletions using *VarGenius-HZD*. 
It requires the use of bedtools and samtools.

1. Download UCSC Hg19 FASTA genome from: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/
1. Download both BAM and BAI files from 1KGP repository (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/ - exome\_alignment folder) as in the bam_list.txt file in the **VarGenius-HZD/data** folder: 
	
2. The bam_list.txt file will be used by *VarGenius-HZD* please change paths accordingly. The samples names containing *\_deleted* indicate those where the deletions will be introduced.
3. Download the 1KGP target BED file here: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/exome_pull_down/20120518.analysis_exome_targets.consensus.annotation.bed
4. Once downloaded the BAM file you will need to insert the deletions. Samples names and corresponding region deleted are the following:

```	
NA06989 21 48063447 48063551
NA07347 21 27326904 27327003
NA12058 21 35091133 35091161
NA12748 21 10906904 10907040
NA12830 21 40188932 40189015
```

Please create 5 new BED files containing the 5 intervals above and use the following commands to generate the new BAM files
	
```
	
bedtools intersect -a sample.bam -b deletion_i.bed -v > sample_deleted.bam; 
samtools sort sample_deleted.bam > sample_deleted_sort.bam; 
samtools index sample_deleted_sort.bam.

```

5. Now that you have the samples with the synthetic deletion you can run *VarGenius-HZD*. 

	5a. Generate the exons on target file:
	
	perl /path/where/you/downloaded/VarGenius-HZD/HZD\_launch.pl -f GET\_EXONS\_ON\_TARGET -p /folder/where/you/downloaded/thistool/ -o /path/to/outfolder/ -ref /path/to/hg19.fa -b /path/to/bedtools -l /path/to/bamlist.txt  -t /path/to/20120518.analysis_exome_targets.consensus.annotation.bed --nochr

	For this test we will remove the *chr* prefix of chromosome names as in the 1KGP BAM files this is not present.

	5b. For the pre-processing step we suggest to use *screen* command or to launch it in background:
	
	perl /path/where/you/downloaded/VarGenius-HZD/HZD\_launch.pl -f HZD\_PREPROCESSING -p /folder/where/you/downloaded/thistool/ -o /path/to/outfolder/ -ref /path/to/hg19.fa -b /path/to/bedtools -s /path/to/samtools -l /path/to/bamlist.txt  -t /path/to/20120518.analysis_exome_targets.consensus.annotation.bed
	
	A genomefile will be automatically generated for 1KGP format file. Then *bedtools coverage* will be run for all samples and needed files will be generated within the same folder.
	
	5c. Detect putative HDs:
	
	perl /path/where/you/downloaded/VarGenius-HZD/HZD_launch.pl -f DETECT\_HDs -p /folder/where/you/downloaded/thistool/ -o /path/to/outfolder/ -ref /path/to/hg19.fa -b //path/to/bedtools -l /path/to/bamlist.txt  -t /path/to/20120518.analysis_exome_targets.consensus.annotation.bed -m 2


	5d. Filter the HDs from the output selecting only the "deleted" samples and average DoC>50:
	
	grep 'deleted' /path/to/outfolder/20120518.analysis_exome_targets.consensus_20120518.analysis_exome_targets.consensus_VG-HZD_final_suspect_table_gene.txt | awk '$7>50'
	

