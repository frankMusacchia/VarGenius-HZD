#!/usr/bin/perl

#VarGenius-HZD path
use lib '/home/francesco/bin/VarGenius-HZD/';
use lib '/home/francesco/bin/VarGenius-HZD/LIB/';
use lib '/home/francesco/bin/VarGenius-HZD/LIB/PERL/';

####PLATFORM_SPECIFIC_SETTINGS_TERMINATED

#!/usr/bin/perl

#VarGenius-HZD searches rare homozygous and hemizygous deletions in targeted sequencing
#Copyright (C) 2021 Francesco Musacchia

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.



use strict;
use warnings;

#Using a library to manage files
use LIB::files_management qw( save_hash load_hash file_not_present
				download_file dl_and_extract delete_file extract_columns_from_file
				append_str_2_file check_presence file_num_rows 
				invert_cols_position get_col_index insert_col_in_file_table delete_columns
				extract_col_from_file file_name_wostrange_chars append_hash_to_file
				delete_rows_containing file_list_to_array extract_colnum_from_file_linux
				shift_column compress_folder extract_name append_str_2_file_if_path_notexist
				join_files_with_cat separate_elements_on_col join_vcf_files_with_cat
				vcf_insert_string_into_field vcf_filter_identical_calls
				overwrite_str_in_file_if_exists tsv_2_xls list_to_array
				count_lines_file compid_intervals_2_bed insert_line_into_file
				is_folder_empty append_file_2_file merge_rows_with_same_col
				print_str_2_file extract_file_folder delete_file_list
				append_list_2_file);
				
				
#Using a library to manage configuration files and their global variables
use LIB::programs_management qw(hash2ConfigFile configFile2Hash try_exec_command 
				try_exec_job initialize_folders print_and_log log_and_exit
				separate_input_ids remove_sample_sheet show_analyses_from_sample_sheet
				get_ped_file remove_analyses clean_folders save_analyses checkConfigVariables
				separate_bed_perchr show_analysesids_from_analysesnames show_samples_from_group
				check_config_parameters check_genedb_links get_flowcell_and_lane
				get_ped_file_for_enlarged_analysis store_and_compress_analyses
				exec_store_and_compress_analyses check_datasets_links
				show_samplesids_from_samplesnames get_extended_regions
				JOB_get_status_field JOB_is_active_queue);
				
##Get the file with the parameters for the module
my $main_data_folder = $ARGV[0];
#my $target_bed = $ARGV[1];
my $ucsc_genes = $ARGV[1];
my $ucsc_exons = $ARGV[2];
my $samplename = $ARGV[3];
my $bam_path = $ARGV[4];
my $target_exons = $ARGV[5];
my $config_file = $ARGV[6];
my $log_file = $ARGV[7];
my $out_folder = $ARGV[8]; #extract_name($bam_path,"path");


#Get the parameters from input
my $cfg_hash = load_hash($config_file);
	
my $hum_ref = $cfg_hash->{'reference_f'};
my $genomefile = $cfg_hash->{'genome_file'};
my $nc_threshold = $cfg_hash->{'min_BoC_HDs'};#"0.2";#percent of interval covered

my $bedtools_path = $cfg_hash->{'bedtools_path'};

my $remove_temp = "NO";
my $bed_ext = "bed";

my $coverbed_prog = "coverage";
my $sortbed_prog = "sort";
my $intersectbed_prog = "intersect";

my $noncov_ext = "_exons.noncov";
#Exonic raw coverage outname:
#samplename_$target_exons.bed
#Non covered exons output
#
		
#I give in input a genome file that has been obtained with 
# samtools faidx ucsc.hg19.fa that produces ucsc.hg19.fa.fai
#and then cut -f1,2  ucsc.hg19.fa.fai >  ucsc.hg19.genomefile
#and with the -sorted parameter the computation will be very fast
#
#Target exons were previously sorted with
#sortBed -i Agilent_SureSelect_CCP_v1_exons_UCSC.bed -faidx names.txt  > Agilent_SureSelect_CCP_v1_exons_UCSC.sort.bed
#where names.txt is a file:
#chrM
#chr1
#chr2
#..
#chrY		
if ( -e $target_exons and !(-z $target_exons) ){
	
	#my $target_name = extract_name($target_bed,1);
	#2: Get exon coverage
	#my $target_exons_cov = $out_folder."/$samplename\_".extract_name($target_exons,1).".$bed_ext";
	my $target_exons_cov = $out_folder."/$samplename.$bed_ext";

	#my $covbedparams = " -g ".extract_name($hum_ref,"noext").".genomefile -sorted ";
	my $covbedparams = " -g $genomefile -sorted ";
	run_BEDTOOLS_coverageBed_gen($target_exons,$bam_path,$target_exons_cov,$covbedparams,$log_file);
	
	
	#3: Restrict to ipotetically low coverage exons only
	my $target_exons_nc = extract_name($target_exons_cov,"noext")."_nc.".$bed_ext;
	print_and_log( "Selecting putative homozygous deletions with BoC < $nc_threshold...",$log_file);#DEBUGCODE
	#my $command = "sed '".'s/\./,/'."' $target_exons_cov | awk "."'".'$7<'.$nc_threshold."'  > $target_exons_nc";
	my $command = " awk "."'".'$7<'.$nc_threshold."' $target_exons_cov  > $target_exons_nc";
	print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
	try_exec_command($command) or die "Unable to execute command: $command\n";	


	#4: intersection with ucsc genes (intersectBed)
	#I use -wo because I want the dimension of the intersection and the Gene name.
	#But I don't want the regions that do not have intersection with Genes because I will not use them
	my $noncovexons_out = extract_name($target_exons_cov,"noext").$noncov_ext;
	my $inters_params = " -wo ";
	print "...intersecting UCSC exons with Genes...";#DEBUGCODE
	run_BEDTOOLS_intersectBed($target_exons_nc,$ucsc_genes,$noncovexons_out,$inters_params,$log_file);#######	

	print_and_log( "...DONE!",$log_file);#DEBUGCODE
	#if ( $remove_temp eq 'YES' ){
	delete_file($target_exons_nc);
	#}
	
}else{
	print "ERROR: $target_exons file does not exist...";#DEBUGCODE
}



#######################UTILITIES



=head2 run_BEDTOOLS_coverageBed_gen

 Title   : run_BEDTOOLS_coverageBed_gen
 Usage   : run_BEDTOOLS_coverageBed_gen(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs bedtools intersect
								
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BEDTOOLS_coverageBed_gen{
	my $bed1 = shift;
	my $bed2 = shift;
	my $outFile = shift;
	my $params = shift;
	my $log_file = shift;
	
	#When comparing alignments in BAM format (-abam) to features in BED format (-b),
	# bedtools wants you to write the output in BED format. Otherwise the error
	# "ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option."
	#will be thrown
	my $files_str = "";
	$files_str = " -a $bed1 -b $bed2 ";
		
	#Execute the command
	my $command = $bedtools_path." ".$coverbed_prog." $files_str $params > $outFile";
	if (defined $log_file){	print_and_log( " [ Executing command: $command ] ",$log_file);}
	else {print " [ Executing command: $command ] "; }
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}



=head2 run_BEDTOOLS_sortBed

 Title   : run_BEDTOOLS_sortBed
 Usage   : run_BEDTOOLS_sortBed(   );

 Function: Runs sortBed taking in input the bam file
								
			INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BEDTOOLS_sortBed{
	my $file = shift;
	my $params = shift;
	my $outFile = shift;
	my $log_file = shift;
	
	#Execute the command
	my $command = "$bedtools_path $sortbed_prog -i $file $params  > $outFile";
	print_and_log( "[ Executing command: $command ]",$log_file);
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}


=head2 run_BEDTOOLS_intersectBed

 Title   : run_BEDTOOLS_intersectBed
 Usage   : run_BEDTOOLS_intersectBed(   );

 Function: Takes parameters from the configuration hash and references to 
					filenames and programs runned before from the database and
					Runs bedtools intersect
								
					INPUT: needs bam output from the alignment against the reference
					
 Returns : nothing
=cut
sub run_BEDTOOLS_intersectBed{
	my $bed1 = shift;
	my $bed2 = shift;
	my $outFile = shift;
	my $params = shift;
	my $log_file = shift;
	
	#When comparing alignments in BAM format (-abam) to features in BED format (-b),
	# bedtools wants you to write the output in BED format. Otherwise the error
	# "ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option."
	#will be thrown
	my $files_str = "";
	if ( $bed1 =~ /\.bam$/ and $bed2 =~ /\.bed$/ ){
		#$params .= " -bed ";
		$files_str = " -abam $bed1 -b $bed2 ";
	}elsif ($bed2 =~ /\.bam$/ and $bed1 =~ /\.bed$/){
		$files_str = " -abam $bed2 -b $bed1 ";
	}else{
		$files_str = " -a $bed1 -b $bed2 ";
	}
		
	#Execute the command
	my $command = $bedtools_path." ".$intersectbed_prog." $files_str $params > $outFile";
	if (defined $log_file){	print_and_log( " [ Executing command: $command ] ",$log_file);}
	else {print " [ Executing command: $command ] "; }
	try_exec_command($command) or die "Unable to execute command: $command\n";		
}



