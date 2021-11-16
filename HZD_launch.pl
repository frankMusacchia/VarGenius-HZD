#!/usr/bin/perl

#VarGenius-HZD path
use lib '/home/francesco/bin/VarGenius-HZD/';
use lib '/home/francesco/bin/VarGenius-HZD/LIB/';
use lib '/home/francesco/bin/VarGenius-HZD/LIB/PERL/';

####PLATFORM_SPECIFIC_SETTINGS_TERMINATED

#VarGenius-HZD searches rare homozygous and hemizygous deletions in targeted sequencing
#Copyright (C) 2022 Francesco Musacchia

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
use Getopt::Long;
use Parallel::ForkManager;#For parallel executions
use File::Copy;#To manage files

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
			
#my $progDir;
#BEGIN { $progDir = getcwd; }
#use lib "$progDir/LIB/PERL/";

my $program_name = "VarGenius-HZD";
my $folders_file = "folders.txt";

#Get working and program folders
my $workDir = ""; 
my $progDir = "";

($workDir,$progDir) = initialize_folders($folders_file);

my $parallelism = "threads";

#Get the confighash and pre-defined variables
my $cfg_hash;
my $program_config = "$progDir/CONFIGURATION/program_config.txt";
#checkConfigVariables($program_config,$variables_file,1);
configFile2Hash($program_config,\$cfg_hash);

#Set needed folders
my $lib_folder = $cfg_hash->{'lib_folder'};#"LIB";
my $config_folder = $cfg_hash->{'config_folder'};#"CONFIGURATION";
my $data_folder = $workDir."/".$cfg_hash->{'data_folder'};
my $log_folder = $workDir."/".$cfg_hash->{'log_folder'};

my $perl_libs_paths_f = $cfg_hash->{'perl_libs_paths_f'};# "perl_libs_paths.txt";#Into config folder

my @main_pls = split(",",$cfg_hash->{'vg_pl_scripts_2_update'});


my $vargenius_sing_cont = $cfg_hash->{'vargenius_sing_cont'};
my $perl_mod_folder = "PERL";
my $log_file = "";


my $program_folder = "";
my $out_folder = "";
my $function = "";#function to execute						
my $bam_list = "";
my $target_bed = "";
#human reference
my $hum_ref = ""; 
#Path to bedtools
my $bedtools_path = $cfg_hash->{'bedtools_path'};# 
my $samtools_path = $cfg_hash->{'samtools_path'};

my $HZD_algorithm = "$progDir/HZD_algorithm.R";

my $coverbed_prog = "coverage";
my $sortbed_prog = "sort";
my $intersectbed_prog = "intersect";

#max_samples_with_lowcov
my $max_lowBoC = "2";

#If the BAM files do not contain the "chr"
my $nochr ="";

#The genome file provided in input
my $genomefile = "";


#Get input parameters
parse_command_line_args();

#Folder containing needed data
my $main_data_folder = "$progDir/data/";
my $program_call = $program_folder."/$program_name/get_nc_exons.pl";


#UCSC exons and genes bed file
my $ucsc_genes = $main_data_folder."/ucsc_genes_hg19";
my $ucsc_exons = $main_data_folder."/ucsc_exons_hg19";
if ($nochr){
	$ucsc_genes .= "_nochr";
	$ucsc_exons .= "_nochr";
}
$ucsc_genes .= ".bed";
$ucsc_exons .= ".bed";

my $target_exons = "$data_folder/".extract_name($target_bed,"1")."_ucsc_exons.bed";#MODIFED
print "Setting the datasets to use:\n";
print "ucsc_genes: $ucsc_genes\n";
print "ucsc_exons: $ucsc_exons\n";
print "target_exons: $target_exons\n";
#my $target_exons = "coverage.bed";


$log_file = $log_folder."/".$function.".log";



##STEP 1
##The target BED file will be used with UCSC exons to select only exons on target.
if ($function eq 'GET_EXONS_ON_TARGET'){
	
	if (file_not_present($target_bed)>0 ){ die "ERROR. The target BED is required in input. Please use --target.\n";}

	if ( !(-e $target_exons) or (-z $target_exons) ){
		print_and_log( "Getting the exons on the target...",$log_file);#DEBUGCODE

		
		#1. Intersect the exons with the target so that we have the exons on-target
		#to see how much is the overlap and also when overlap is missing	
		#Parameters used are:
		 # -a the exons
		 # -b the target file
		 # -wa which allows to write the original entry in A . 
		my $inters_params = " -wa ";
		print_and_log( "...intersecting UCSC exons with target (intersectBed)...",$log_file);#DEBUGCODE
		run_BEDTOOLS_intersectBed($ucsc_exons,$target_bed,$target_exons.".temp",$inters_params,$log_file);#######	
		
		#Get only the chrom and interval, sort and get uniq
		print_and_log( "Get only the chrom and interval, sort and get unique lines...",$log_file);#DEBUGCODE
		my $command = " cut -f1-3 $target_exons.temp | sort -k1,1V -k2,2n | uniq > $target_exons.temp2";
		print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
		try_exec_command($command) or die "Unable to execute command: $command\n";		

		#Removing "chr" from target exons since the BAM files contain 1 instead of chr1
		if ($nochr){
			print_and_log( "Removing chr...",$log_file);#DEBUGCODE
			#$command = "  sed 's/chr//'  $target_exons.temp3 > $target_exons";
			$command = "  sed 's/chr//'  $target_exons.temp2 > $target_exons";	
			print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
			try_exec_command($command) or die "Unable to execute command: $command\n";	
		}else{
			move("$target_exons.temp2",$target_exons);
		}
		print_and_log( "\nExons-on-target file generation is done!\n",$log_file);#DEBUGCODE
	
			
		#if (  $remove_temp eq 'YES' ){
			delete_file("$target_exons.temp");
			delete_file("$target_exons.temp2");
			delete_file("$target_exons.temp3");
		#}			
	}else{print "Cool! $target_exons already present! Go ahead with VarGenius-HZD!\n";}	
}


#STEP 2
#A Genome file is generated if not provided for bedtools and for each sample in bam_list.txt
#will be executed bedtools coverage to obtain BoC and DoC for exons_on_target
if ($function eq 'HZD_PREPROCESSING'){

my $config_file = $log_folder."/".$function.".hash";

save_hash(\$cfg_hash,$config_file);


	my $log_file = $log_folder."/$program_name.log";

	#This analysis is performed for each sample, hence the first operation is to get the samples list 
	open (FILE,"<$bam_list") or die "Cannot open $bam_list\n";
	my @lines = <FILE>;
	close(FILE);
	
	#When using threads instead of jobs declare commands
	my @commands = ();
	
	my $num_row = 0;
	foreach my $row (@lines){
		if ($num_row > 0){
			chop($row);
			my @pieces = split("\t",$row);
			my $samplename = $pieces[0];
			my $bam_path = $pieces[1];
			
			my 	$program_call = "perl $progDir/get_nc_exons.pl".
				" $main_data_folder $ucsc_genes $ucsc_exons $samplename $bam_path $target_exons ".
				" $config_file $log_file $workDir";
				
			push(@commands,$program_call);

		}
		$num_row++;
	}
	execute_threads(\@commands);
	print "Commands executed and running...\n";	
}

#STEP 3
#suspect HDs are detected and annotation is provided
#Launch the algorithm
if ( $function eq 'DETECT_HDs'){
	
	#Check for samtools
	if (! (-e $bam_list) ){
		die "list of samples to analyze was not provided. Please use --bam_list\n";
	}
	#print "Target exons: $target_bed\n";
	my $target_suff = extract_name($target_exons,1);
	print "targetsuff $target_suff\n";
	my $out_suff = extract_name($target_exons,1)."_VG-HZD";
	print "out_suff $out_suff\n";

	
	#my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $bam_list $target_suff ucsc_exons _exons.noncov $out_folder $out_suff $max_lowBoC' $HZD_algorithm";
	my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $bam_list $target_suff exons.noncov $workDir $out_suff $max_lowBoC ".$cfg_hash->{'analtype'}."' $HZD_algorithm";

	print  "Launching $command\n";
	try_exec_command($command);
}



=head2 try_exec_slurm_job

 Title   : try_exec_slurm_job
 Usage   : try_exec_slurm_job(    
	$qsub_account => account to use for execution;
	$qsub_queue => queue to use;
	$env_vars =s> variables to pass to subject program
	$program_to_run => program to be run;
	$lnodes => number of nodes needed;
	$ppn => number of cpus needed;
                               );

 Function:  Given in input a command it will try to execute it with system function more than one times.
						
						If the variable qsub_param_for_all=YES then it takes the parameters for the resources
						request from the program_config. This means that all the tasks will ask for the same
						number of resources. Otherwise the suffix of the variable will change with the task name 
						(qc,align, refine, etc..) and the resources will be ask depending by the task 
						from the user_config.
						The following are the parameters for qsub that can be changed
							qsub_mem = 120GB
							qsub_ncpus = 12
							qsub_nodes =
							qsub_select = 1
							
							
 Returns : 1 if succeed, -1 otherwise

=cut
sub try_exec_slurm_job{
	my $cfg_hash = shift;
	my $env_vars =shift;
	my $task = shift;
	my $program_path = shift;
	my $job_name = shift;
	my $dependencies = shift;
	my $std_out = shift;
	my $std_err = shift;
	my $log_file = shift;
	
	my $qsub_command;
	#Account to be used to launch job on the cluster
	if (defined $cfg_hash->{'qsub_cmd'} ){
		$qsub_command = $cfg_hash->{'qsub_cmd'};
	}	
	#A string with parameters for the job taken from the program config file
	my $qsub_params = "";
	#Account to be used to launch job on the cluster
	if (defined $cfg_hash->{'qsub_account'} ){
		$qsub_params .= " -p ".$cfg_hash->{'qsub_account'};
	}


	##Active the restartable job
	#if ( $cfg_hash->{'qsub_restartable'} eq 'YES'){
		#$qsub_params .= " -r y";
	#}
	######Given in input
	#Define the jobname
	if ($job_name ne 'null' ){
		$qsub_params .= " --job-name=".$job_name;
	}

	#if (defined $cfg_hash->{'qsub_account'} ){
		#$qsub_params .= " -W group_list=".$cfg_hash->{'qsub_account'}." ";
	#}
	##Define the standard output
	#if ($std_out ne 'null' ){
		#$qsub_params .= " -o ".$std_out;
	#}
	##Define the standard error
	#if ($std_err ne 'null' ){
		#$qsub_params .= " -e ".$std_err;
	#}	
	#From the user configuration file (dependent by program)
	my $qsub_resources = "";
	#Define the resources needed
	my $suff ="";
	if ( $cfg_hash->{'qsub_param_for_all'} eq 'YES' ){
		$suff = $qsub_command;
	}else{
		$suff = $task;
	}
	
	#The queue to use. First check the one for the tast, then the default one
	if (defined $cfg_hash->{$suff.'_queue'}){
		$qsub_params .= " -q ".$cfg_hash->{$suff.'_queue'};
	}else{	
		if (defined $cfg_hash->{'qsub_queue'} ){
			$qsub_params .= " -q ".$cfg_hash->{'qsub_queue'};
		}
	}
		
	#QSUB RESOURCES
	if ( defined $cfg_hash->{$suff.'_cpus'} 
				or defined $cfg_hash->{$suff.'_nodes'} or defined $cfg_hash->{$suff.'_job_memory'}
				or defined $cfg_hash->{$suff.'_ncpus'}){
		$qsub_resources .= " -l ";
	}
	if (defined $cfg_hash->{$suff.'_ncpus'}){
		$qsub_resources .= "--ntasks-per-node=".$cfg_hash->{$suff.'_ncpus'}.":";
	}
	if (defined $cfg_hash->{$suff.'_mem'}){
		$qsub_resources .= "--mem=".$cfg_hash->{$suff.'_mem'}.":";
	}

	chop($qsub_resources);	
			
	#The walltime to use
	if (defined $cfg_hash->{$suff.'_walltime'} ){
		$qsub_params .= " --time=".$cfg_hash->{$suff.'_walltime'};
	}
	#Use dependencies
	if ($dependencies ne 'no_depend'){
		#$qsub_params .= " -W depend=".$cfg_hash->{$suff.'_depend'}.":$dependencies ";
	}		
	
	#Add environmental variables
	$qsub_params .= " $qsub_resources ";
	if ( $env_vars ne 'NONE'){
		 $qsub_params .= " --export=$env_vars ";
	}
	my $command = $qsub_command.$qsub_params.$program_path;
	print_and_log( "Executing command: $command\n",$log_file);
	
	my $qsub_id = `$command`;
	chomp($qsub_id);
	if ($qsub_id ne ''){
		#print "Output from qsub: $qsub_id\n";#DEBUGCODE
		#The qsub id should come out with number.host. Hence I remove the host
		$qsub_id =~ s/\..*//;
		#print "qsub_id: $qsub_id\n";#DEBUGCODE
	}else{
		print "ERROR: Unable to start $program_path. Please check error from qsub..\n ";
	}
	return $qsub_id;
	die "Unable to execute: $command\n";
}



=head2 grep_file
 Title  : grep_file
 Usage  : grep_file(  -fileToCheck => 'path to the file to check');

 Function: 	Checks if a file  contains the string given in input

 Returns : 1 if the string is present, otherwise 0

=cut
sub grep_file{
	my $fileToCheck = shift;
	my $string = shift;
	
	my $retval = 0;
	open( FILE, "<$fileToCheck" ) or die "$!";
	while (my $line = <FILE>) {
		if ( $line =~ /$string/ ) {
			print "found one!\n";#DEBUGCODE
			$retval = 1;
			last;
		}
		#else{
		 #print "$line does not match $string\n";
		#}
	}
	close(FILE);
	return $retval;
}


=head2 execute_threads

 Title  : execute_threads
 Usage  : execute_threads( - commands => 'the commands to be executed',
                              );

 Function: for each command it will launch a thread performing that operation
			
 Returns : nothing

=cut	 
sub execute_threads{
  my $commands = shift;#reference to array of commands
  
  my $processes = scalar(@$commands);
  print "Executing $processes commands \n";# I'll write a log in $log.\n"; 
  #Creates an object to fork processes in parallel
  my $pm = new Parallel::ForkManager($processes);

  #You can define a subroutine which is called when a child is terminated. It is called in the parent process.
  #  - pid of the process, which is terminated
  # - exit code of the program
  $pm->run_on_finish(
    sub {
      my($pid,$exit_code) = @_;
      #print "** Just got out of the pool with PID $pid and exit code: $exit_code\n";#DEBUCODE
    }
  );

	foreach my $cmd (@$commands){
    # Forks and returns the pid for the child:
    my $pid = $pm->start and next;
   
    # Here is the parallelized block
    # -----------
    #print "$pid ---> running \n";
    print  "Launching $cmd\n";
    try_exec_command($cmd);
    
    # Terminates the child process
    $pm->finish;
    print "Thread $pid finished. \n";
	}

  $pm->wait_all_children; 
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

#

=head2 substitute_start_of_file

 Title   : substitute_start_of_file
 Usage   : substitute_start_of_file( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: substitutes the start of the file given in input.
					It searches for the string_to_find, and until that string
					writes only the string_to_ins in the new file. Once that the string
					is found it copies all the rest of the input file into the output
					
 Returns : nothing

=cut
sub substitute_start_of_file{
  my $file = shift;
  my $string_to_find = shift;
  my $string_to_ins = shift;
 
  my $temp_f = $file.".temp";
  
  open(TEMP, ">$temp_f") or die "Couldn't open file $temp_f";
  open(FILE, "<$file") or die "Couldn't open file $file";
	
	#Print the input string in the new file
	print TEMP $string_to_ins."\n";
	
	# print the lines before the change
	while( my $line = <FILE> ) {
    chop($line);
    last if $line eq $string_to_find; # chosen string to find found
   }

	print TEMP $string_to_find."\n";
	# print the rest of the lines
	while( my $line = <FILE> )  {
     print TEMP $line;
   }
  close(TEMP);  
  close(FILE);
  
  #Overwrite input file
  move($temp_f,$file) or print "ERROR: unable to move $temp_f into $file\n";
	#Remove the temp file
	unlink($temp_f);
}

=head2 substitute_text_in_file

 Title   : substitute_text_in_file
 Usage   : substitute_text_in_file( - path1 -> the file to be copied,
                                - path2 => the file to be written)

 Function: substitutes a text within a file given a start and an end strings
					
 Returns : nothing

=cut
sub substitute_text_in_file{
  my $file = shift;
  my $start = shift;
  my $end = shift;
  my $string_to_ins = shift;
 
  my $temp_f = $file.".temp";
  
  open(TEMP, ">$temp_f") or die "Couldn't open file $temp_f";
  open(FILE, "<$file") or die "Couldn't open file $file";
	

	# print the lines before the change
	while( my $line = <FILE> ) {
 	#Print the input string in the new file
	print TEMP $line;
	#Remove the newline for the check
	chop($line);
    last if $line eq $start; # chosen string to find found
   }
	#Print the input string in the new file
	print TEMP $string_to_ins;

	
	# jump all that is between start and end
	while( my $line = <FILE> )  {
	chop($line);
     last if $line eq $end; # chosen string to find found
   }
   
   #Now print the end and the rest of the file
   print TEMP $end."\n";
   	while( my $line = <FILE> )  {
	print TEMP $line;
   }
   
  close(TEMP);  
  close(FILE);
  
  #Overwrite input file
  move($temp_f,$file) or print "ERROR: unable to move $temp_f into $file\n";
	#Remove the temp file
	delete_file($temp_f);
}


=head2 parse_command_line_args

 Title   : parse_command_line_args
 Usage   : parse_command_line_args(   );

 Function:  Parses the arguments specified upon the command line.
						
 Returns : nothing

=cut
sub parse_command_line_args{
  my $HELP  = 0;# Shows help overview.
  my $VERSION = 0;# Shows version number and exit.

	
	my $howToUse = "Use with: perl HZD_launch.pl \n\n";
	
	#  Parse options
  GetOptions(
           "help" =>        \$HELP,
           "version" =>      \$VERSION,
           "p|program_folder=s" => \$program_folder,
           "o|out_folder=s" => \$out_folder,
           "ref|human_ref=s" => \$hum_ref,
           "b|bedtools=s" => \$bedtools_path,
           "l|bam_list=s" => \$bam_list,
           "f|function=s" => \$function,
           "t|target=s" => \$target_bed,
			"m|max_lowBoC=s" => \$max_lowBoC, 
			"nc|nochr" => \$nochr, 
			"g|genomefile=s" => \$genomefile,
			"s|samtools=s" => \$samtools_path
	          );
           
           #Print a little help
  if ( $HELP ){
    print $howToUse;
    pod2usage(1);
    exit;
  }

  #Print version
  if ( $VERSION ){
    print "Version: 1.0 \n";
    exit;
  }
  
}
