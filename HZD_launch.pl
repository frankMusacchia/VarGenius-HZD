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
use Getopt::Long;
use Parallel::ForkManager;#For parallel executions

my $parallelism = "threads";

my $cfg_hash;
$cfg_hash->{'qsub_cmd'} = "sbatch";
$cfg_hash->{'R_path'} = "R";

my $program_name = "VarGenius-HZD";
my $program_folder = "";
my $out_folder = "";
my $function = "";#function to execute						
my $bam_list = "";
my $target_bed = "";
#human reference
my $hum_ref = ""; 
#Path to bedtools
my $bedtools_path = "";# 
my $coverbed_prog = "coverage";
my $sortbed_prog = "sort";
my $intersectbed_prog = "intersect";
	
#Get input parameters
parse_command_line_args();

#Folder containing needed data
my $main_data_folder = $program_folder."/$program_name/data/";
my $program_call = $program_folder."/$program_name/get_nc_exons.pl";
#UCSC genes bed file
my $ucsc_genes = $main_data_folder."/ucsc_genes.bed";
#Prepare Target Exons:
my $ucsc_exons = $main_data_folder."/UCSC_coding_exons.bed";
my $target_exons = extract_name($target_bed,"noext")."_ucsc_exons_bed";

if (file_not_present($target_bed) ){ die "Cannot proceed Check: $target_bed.\n";}

	
if (! (-d $out_folder) ){
	die "ERROR: $out_folder does not exist. Please use --out_folder\n";
}
		
if ($function eq 'GET_TARGET_EXONS'){
	my $log_file = $out_folder."/".$function.".log";

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
		my $command =	"  cut -f1-3 $target_exons.temp | sort -k1,1 -k2,2n | uniq > $target_exons.temp2";
		print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
		try_exec_command($command) or die "Unable to execute command: $command\n";		
		
		#Sort using bedtools and a reference chromosome order
		my $chr_order_f = "chr_order.txt";
		my $sortparams = " -faidx $main_data_folder/$chr_order_f ";
		run_BEDTOOLS_sortBed("$target_exons.temp2",$sortparams,"$target_exons.temp3",$log_file);
		
		#Removing "chr" from target exons since the BAM files at the CNAG contain 1 instead of chr1
		print_and_log( "Removing chr...",$log_file);#DEBUGCODE
		$command =	"  sed 's/chr//'  $target_exons.temp3 > $target_exons";
		print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
		try_exec_command($command) or die "Unable to execute command: $command\n";		
			
		#if (  $remove_temp eq 'YES' ){
			delete_file("$target_exons.temp");
			delete_file("$target_exons.temp2");
			delete_file("$target_exons.temp3");
		#}			
	}else{print "Cool! $target_exons already present!\n";}	
}
		
if ($function eq 'VGII_COV_STEP'){


	my $log_file = $main_data_folder."/$program_name.log";
	my $task = "align";

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
			
			my 	$program_call = "perl $program_folder/VarGeniusII-HZD/get_nc_exons.pl".
				" $main_data_folder $target_bed $ucsc_genes $ucsc_exons $samplename $bam_path $target_exons $hum_ref ".
				" $bedtools_path $log_file";
				
			push(@commands,$program_call);

		}
		$num_row++;
	}
	execute_threads(\@commands);
	print "Commands executed and running...\n";	
}


#Launch the algorithm
if ( $function eq 'VGII-HZD'){
	
	my $HZD_algorithm = "$program_folder/$program_name/HZD_algorithm.R";
	my $log_file = $HZD_algorithm.".log";

	my $suff_from_target = extract_name($target_exons,1);
	my $out_suff = extract_name($target_exons,1)."_VGII-HZD";

	my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $bam_list $suff_from_target ucsc_exons exons.noncov $out_folder $out_suff' $HZD_algorithm";			
	#my $command = $cfg_hash->{'R_path'}." CMD BATCH --no-save --no-restore '--args $bam_list $suff_from_target exons.noncov $out_folder $out_suff' $HZD_algorithm";			
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




=head2 extract_name

 Title   : extract_name
 Usage   : extract_name( -filePath => 'complete path of the file',
                        -type => 'a number saying what you want to extract'
			       );

 Function: extract the name from a complete path of a file. Even the file name only with the extension
              0: the complete name with the extension
              1: the name only (remove the last extension: string separated by dot)
              2: the first two words joined by a dot
              noext: just remove the last extension from the path
              no2ext: remove the last 2 extensions from the path
              gz: the name from the .gz
              targz: the name from the .tar.gz
              zip: the name from the .zip
              tar: the name from the .tar
              fqgz: the name before two extensions

 Returns : the name only

=cut
sub extract_name {
  my $filePath = shift;#Path to the file
  my $type = shift;#The type of file
  
  #Separate the path in pieces using slashes
  my @list = split("/",$filePath);
  my $complName = pop(@list);
  
  
  my @nameElements = split (/\./, $complName);
  my $name;
  if ($type eq "0"){ $name = $complName;}
  elsif ($type eq "1"){ 
	 pop(@nameElements);
	 my $name = join(".",@sams);
	 }
  elsif ($type eq "2"){ $name = $nameElements[0].'.'.$nameElements[1];}
  elsif ($type eq 'noext'){    
		my @parts = split(/\./,$filePath);
    pop @parts;
    $name = join '.', @parts;
	}
  elsif ($type eq 'no2ext'){    
		my @parts = split(/\./,$filePath);
    pop @parts;
    pop @parts;
    $name = join '.', @parts;
	}
  elsif ($type eq "gz"){ $complName =~ /(\S+).gz/;
                 $name= $1;##Name for uncompressed file
                 }
  elsif ($type eq "targz"){$complName =~ /(\S+).tar.gz/;
                 $name= $1;##Name for uncompressed file
                 }
   elsif ($type eq "zip"){$complName =~ /(\S+).zip/;
                 $name= $1;##Name for uncompressed file
                 }
   elsif ($type eq "tar"){$complName =~ /(\S+).tar/;
                 $name= $1;##Name for uncompressed file
                 }
   elsif ($type eq "fqgz"){
		$complName =~ /(\S+)\.\S+\.\S+/;
		
                 $name= $1;##Name for uncompressed file
                 }
   else { die	"ERROR [$?]: $type is not a valid input extracting a name: ?\n";}
  return $name;
  
}


=head2 print_and_log

 Title   : print_and_log
 Usage   : print_and_log( - string -> the sentence that have to be print_and_log 
											- onlyLog -> a number);

 Function: will print (the string in input always in the log file and on the STDOUT
					if onlyLog is used then the print will be only in the log file
 Returns : nothing

=cut
sub print_and_log{
  my $string = shift;    
  my $logFile = shift; 
  my $onlyLog = shift;
  
  open(LOG, ">>$logFile") or die "ERROR [$!]: Cannot open $logFile! Check permissions.\n";
  if ( defined $onlyLog){
    print LOG $string;
  }else{
    my $STDOUT = *STDOUT;
    my $LOG = *LOG;
    #Prints on both OUT
    for ($LOG, $STDOUT) { print $_ $string; }
  }
  #Close the log
	close(LOG)
}




=head2 try_exec_command

 Title   : try_exec_command
 Usage   : try_exec_command( -sysCall => is the string that should be repeated
                               );

 Function:  Given in input a command it will try to execute it with system function more than one times.

 Returns : 1 if succeed, -1 otherwise

=cut
sub try_exec_command{
    my $command = shift;

    my $maxTimes = 5;
    my $success = -1;
    my $timesCount = 0;

    while ($success == -1 and $timesCount < $maxTimes){
        if ( (system $command) == 0) {
          $success = 1;
        }
        else{
         if ($? == -1) {
              print "failed to execute: $!\n";
          }
          elsif ($? & 127) {
              printf "child died with signal %d, %s coredump\n",
                  ($? & 127),  ($? & 128) ? 'with' : 'without';
          }
          elsif ($? == 35 ) {
              printf "child exited normally\n",
              $success = 1;   
          }          
          else {
              printf "child exited with value %d\n", $? >> 8;
          }
         $timesCount++;
        }
    }
    return $success;
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

=head2 file_not_present
 Title  : file_not_present
 Usage  : file_not_present(  -fileToCheck => 'path to the file to check');

 Function: 	Checks if a file is present and its dimensions are more than zero

 Returns : 1 if the filename is an empty string
					 2 if the file is not present
					 3 if the file is present but the size is 0

=cut
sub file_not_present{
 my $fileToCheck = shift;
 
 my $retVal = 0;
 #print $fileToCheck." \n";
 if ( $fileToCheck ne ''){
   if(-z $fileToCheck){
    $retVal = 1;
    print "$fileToCheck is present but empty...!\n";
   }
   #-e could return false if the directory is not readable by the program giving a misleading behaviour
	
   if(! (-e $fileToCheck) ){
    $retVal = 2;
		#print "$fileToCheck does not exists..\n";
   }
 }else{
	 print "The file to be checked is an empty string. Please check the code...\n";
	 $retVal = 3;
	}
 
 return $retVal;
}

=head2 delete_file

 Title   : delete_file
 Usage   : delete_file( -filePath => 'a path to a file to delete'
			       );

 Function:  Delete a file in a given location. It can erase also more 
						than one file starting with the same name. This
            can be done by using the usual '*'. When it is found the 
            'glob' keyword will be used.
 Returns : Error code: 1 -> correctly deleted; -1 -> error with perl function; -2 file does not exist

=cut
sub delete_file{
  my $filePath = shift;
 
  my $retVal = 0;
  
	#With this first IF we can delete a set of file with the same name
	if ($filePath =~ /\*/){
		$retVal = 1;
		unlink glob $filePath or $retVal = -1;
	 }elsif ( -e $filePath ){
			if  ( unlink($filePath) == 1) {
				$retVal = 1;
				#deleted successfully
			}else{
			# not deleted for some problems with unlink subroutine or the file does not exist anymore
				$retVal = -1;
			}
	}else{
	# does not exist
		$retVal = -2;
	}
    
  return $retVal;
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
           "t|target=s" => \$target_bed 
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
