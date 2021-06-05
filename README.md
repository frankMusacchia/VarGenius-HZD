# VarGenius-HZD

Rare homozygous/hemizygous deletions caller from targeted sequencing

----------------------------------------

**Dependencies**

R

PERL (Needed libraries: Parallel::ForkManager; Getopt::Long)

bedtools


**Input to prepare**

- genome in FASTA format that you used for the alignment
- the genome-file for bedtools using your genome. E.g. for UCSC Hg19:
	samtools faidx ucsc.hg19.fa (that produces ucsc.hg19.fa.fai)
	cut -f1,2  ucsc.hg19.fa.fai > ucsc.hg19.genomefile
- the target BED file for your sequenced samples
- a file with a list of sample names and paths to the BAM files tab separated (bamlist.txt). 
  E.g.:

```

samplename	path
sampleA	/path/to/bam/
sampleB	/path/to/bam/

```
**Execution**

**Step 1: Get the exons on-target (to be run only once per-target)**

perl HZD_launch.pl -p /path/where/you/downloaded/VarGenius-HZD -o /path/to/outfolder/ -ref /path/to/genome.fasta -b //path/to/bedtools -l /path/to/bamlist.txt -f GET\_EXONS\_ON\_TARGET -t /path/to/target.bed

**Step 2: Execute the breadth/depth of coverage step for each sample**

perl HZD_launch.pl -p /path/where/you/downloaded/VarGenius-HZD -o /path/to/outfolder/ -ref /path/to/genome.fasta -b //path/to/bedtools -l /path/to/bamlist.txt -f HZD\_PREPROCESSING -t /path/to/target.bed

**Step 3. Launch HZD algorithm**

perl HZD_launch.pl -p /path/where/you/downloaded/VarGenius-HZD -o /path/to/outfolder/ -ref /path/to/genome.fasta -b //path/to/bedtools -l /path/to/bamlist.txt -f DETECT\_HDs -t /path/to/target.bed


