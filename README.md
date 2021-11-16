# VarGenius-HZD

Rare homozygous/hemizygous deletions caller within a set of BAM files

----------------------------------------

Please, clone or download the entire folder of VarGenius-HZD from this repository. 
Installation is not required, but you need some dependencies...


**Installation with Singularity**

- Download entire folder of VarGenius-HZD from this repository
- Install singularity and PERL
- Create a directory to work in

mkdir /home/frank/hzd_analyses
cd VarGenius-HZD
perl install.pl 

- When asked provide 
	- 1. your work directory (e.g. /path/to/hzd_analyses).
	- 2. the path of the BAM files 
	- 3. the folder of the reference genome

VarGenius-HZD automatically detects your singularity installation and sets the paths to samtools, bedtools and R.
If you won't plan to install singularity then you have to manually install them.



**Input to prepare**

- the reference genome in FASTA format that you used for the alignment
- the genome file for your genome. It is generated automatically, otherwise see FAQs. 
- the target BED file for your sequenced samples
- a file with a list of sample names and paths to the BAM files tab separated (bamlist.txt). 
  E.g.:

```

samplename	path
sampleA	/path/to/SampleA.bam
sampleB	/path/to/SampleB.bam

```

**BE CAREFUL! :  The BAM file, the genome fasta and the target BED file should have the same naming convention for chromosomes! 
Read the FAQs for more details on this and how to check prior to run the analysis.**

If the naming convention is not with "chr1" but with "1" to enumerate the chromosomes you must always add the parameter -nochr to the following commands.


Now you should be able to run VarGenius-HZD pre processing and algorithm:

**Step 1: Get the exons on-target (to be run only once per-target)**

Go to the work folder

cd /path/to/hzd_analyses

perl /path/to/VarGenius-HZD/HZD_launch.pl -f GET\_EXONS\_ON\_TARGET  -t /path/to/target.bed [-nochr]


**Step 2: Execute the breadth/depth of coverage step for each sample**

We suggest to run Step 2 and Step 3 within a screen or in background as they will require a while to finish.

perl /path/to/VarGenius-HZD/HZD_launch.pl -f HZD\_PREPROCESSING  -l /path/to/bamlist.txt  -t /path/to/target.bed  [-nochr]


**Step 3. Launch HZD algorithm**

Here the -m parameter (max_lowBoC) allows to define the number of samples which might have a specific exon HD.
Default is 2, Please change accordingly to your case (e.g. siblings in the dataset) 

perl /path/to/VarGenius-HZD/HZD_launch.pl -f DETECT\_HDs -l /path/to/bamlist.txt  -t /path/to/target.bed -m 2


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
	
2. The bam_list.txt file will be used by *VarGenius-HZD* please change paths accordingly. The samples names containing *\_deleted* indicate those where the deletions were introduced.
3. VarGenius-HZD folder contains the 1KGP target BED file with changed naming convention as it is within the BAM files (VarGenius-HZD/data folder).
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

if you planned to use our singularity container, it includes bedtools and samtools. Hence you can run the previous commands with command like this:

```
singularity run VarGenius-HZD/LIB/SOFTWARE/VarGenius-HZD_container.sif bedtools

```

5. Now that you have the samples with the synthetic deletion you can run *VarGenius-HZD*. 

	5a. Generate the exons on target file:
	
	perl /path/to/VarGenius-HZD/HZD\_launch.pl -f GET\_EXONS\_ON\_TARGET  -l /path/to/bamlist.txt  -t /path/to/VarGenius-HZD/data/20120518.analysis_exome_targets.consensus.annotation.nochr.bed -nochr

	For this test we will remove the *chr* prefix of chromosome names as in the 1KGP BAM files this is not present.

	5b. For the pre-processing step we suggest to use *screen* command or to launch it in background:
	
	perl /path/to/VarGenius-HZD/HZD\_launch.pl -f HZD\_PREPROCESSING -l /path/to/bamlist.txt  -t /path/to/VarGenius-HZD/data/20120518.analysis_exome_targets.consensus.annotation.nochr.bed -nochr
	
	A genomefile will be automatically generated for 1KGP format file. Then *bedtools coverage* will be run for all samples and needed files will be generated within the same folder.
	
	5c. Detect putative HDs:
	
	perl /path/to/VarGenius-HZD/HZD_launch.pl -f DETECT\_HDs -l /path/to/bamlist.txt  -t /path/to/VarGenius-HZD/data/20120518.analysis_exome_targets.consensus.annotation.nochr.bed -m 2


	5d. Filter the HDs from the output selecting only the "deleted" samples and average DoC>50:
	
	grep 'deleted' /path/to/workfolder/20120518.analysis_exome_targets.consensus_20120518.analysis_exome_targets.consensus_VG-HZD_final_suspect_table_gene.txt | awk '$7>50'
	


Frequently Asked Questions

**Issues with the genome file for bedtools*

Please read carefully here. 
VarGenius-HZD uses bedtools coverage with the -sorted option to accelerate the process. This assumes that the records are sorted lexicographically (e.g., chr1, chr10, etc.). bedtools requires the genome file that must be produced using the exact reference genome file that was used for the alignment (http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html#g-define-an-alternate-chromosome-sort-order-via-a-genome-file).

It is better that prior to launch VarGenius-HZD you check that the reference fasta file, the BAM file and your target BED keep the same nomenclature for the chromosomes i.e. chromosome 1 can be "chr1" or "1". The BAM file keeps the nomenclature of the fasta file, but the target BED file might not. Please check also this latter.

You can open with an editor both the fasta and the target BED while you will require samtools to read the BAM:

e.g.	samtools view -H sample.bam | grep SQ | cut -f 2 | awk '{ sub(/^SN:/, ""); print;}'


If the reference folder is not writable you must produce manually the genome file. Please execute these commands:
	1. samtools faidx ucsc.hg19.fa (that produces ucsc.hg19.fa.fai)
	2. awk -v OFS='\t' {'print $1,$2'} ucsc.hg19.fasta.fai > ucsc.hg19.genomefile



For example 1KGP samples are sorted "numerically" (e.g., 1, 2, etc.).  If so, you need to provide the tool with a "genome file" (-g option) defining the expected order. 
	


**Issues with naming convention in bedtools**

If you get an error like the followiing:

***** WARNING: File SampleA_ucsc_exons_nc.bed has inconsistent naming convention for record:
1       896073  896180  0       0       107     0.0000000

means that bedtools is comparing something with "chrX" and "X" and doesn't get it. 
