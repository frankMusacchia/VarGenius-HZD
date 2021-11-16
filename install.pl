#!/usr/bin/perl

#VarGenius-HZD searches rare homozygous and hemizygous deletions in targeted sequencing
#Copyright (C) <2022>  <Francesco Musacchia>

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
    
#Getting the folder path of VarGenius
use Cwd;
my $progDir;
BEGIN { $progDir = getcwd; chmod 0755,$progDir."/".$program_name;} 

print "$progDir\n";
#LIBRARIES NEEDED HERE
use lib "$progDir";
use lib "$progDir/LIB/PERL/";

use utf8;
use JSON;
  
use strict;
use warnings;

use File::Copy;
use Getopt::Long;


my $update = "";
my $program_config = "CONFIGURATION/program_config.txt";
my $variables_file = "CONFIGURATION/variables.txt";
				
my $program_name = "VarGenius-HZD";
my $configUserFile = "user_config.txt";

print "This script will prepare $program_name to work in your directory.\n";

#Get the confighash and pre-defined variables
my $cfg_hash;
checkConfigVariables($program_config,$variables_file,1);
configFile2Hash($program_config,\$cfg_hash);
my $data_folder = $cfg_hash->{'data_folder'};
my $lib_folder = $cfg_hash->{'lib_folder'};#"LIB";
my $config_folder = $cfg_hash->{'config_folder'};#"CONFIGURATION";
my $log_folder = $cfg_hash->{'log_folder'};#"LOG";
my $repository_folder = $cfg_hash->{'repository_folder'};#"REPOSITORY";

#Files included ito be inserted during the installation
my $perl_libs_paths_f = $cfg_hash->{'perl_libs_paths_f'};# "perl_libs_paths.txt";#Into config folder
my $folders_file = $cfg_hash->{'folders_file'};
my $perl_mod_folder = "PERL";
my $vg_version = $cfg_hash->{'VarGeniusVer'};
my $vargenius_singularity_cont = $lib_folder."/SOFTWARE/VarGenius-HZD_container.sif";

my $bam_folder = "";
my $reference_folder = "";
my $reference_file = "";

my $sure= getInput("Are you sure you want to enjoy this software?(y or n)",'^[YynN]$');

my @main_pls = split(",",$cfg_hash->{'vg_pl_scripts_2_update'});


parse_command_line_args();

if ( ($sure eq "y") or ($sure eq "Y")){ 
  
  #Asks to the user what is the directory where he wants to work
  my $workDir = '';    
  while ( !(-e $workDir) ){
    print "\nWrite an existing complete path where you want to play with $program_name (/home/username/vargenius_analyses): ";
    $workDir= <STDIN>;
    chomp $workDir; #to cancel the return 
  }

	#Change the configuration according to user
	initial_config_vargenius($workDir);
  
    
  #Writing a file with the paths to the current folder and the one
  open (FOLD, ">$workDir/$folders_file") or die "Cannot create file $folders_file. Exiting.. \n";
  print FOLD "$workDir $progDir";
  close(FOLD);
	
 #print "UPDATE: $update\n";
	#Set the LIB folder complete path
	$lib_folder = "$progDir/$lib_folder";
  #Set the file with libraries if it exists
	my $perl_libs_paths_f_p = "$progDir/$config_folder/$perl_libs_paths_f";
	
		
  #Creating the LOG folder into the working directory
  my $main_log_folder = $workDir."/$log_folder";
	#Check if directory exists, otherwise it creates it
	unless(-d $main_log_folder){
		print "Creating folder $main_log_folder...\n";
		mkdir $main_log_folder or die "ERROR: can't create folder $main_log_folder. Check permissions. \n";
	}

  #Creating the DATA folder into the working directory
  my $main_data_folder = $workDir."/$data_folder";
	#Check if directory exists, otherwise it creates it
	unless(-d $main_data_folder){
		print "Creating folder $main_data_folder...\n";
		mkdir $main_data_folder or die "ERROR: can't create folder $main_data_folder. Check permissions. \n";
	}
		
	#Now change the starting part of the .pl scripts in VarGenius according to
	#system settings
	print "Updating the perl scripts with specific platform settings...\n";
	my $start_string = "####PLATFORM_SPECIFIC_SETTINGS_TERMINATED";
	my $perl_location = `which perl`;
	my $vargenius_lib = "use lib '".$progDir."/';";
	chop($perl_location);
	my $new_start_string = "#!".$perl_location."\n\n#$program_name path\n".$vargenius_lib."\n";
	$new_start_string .= "use lib '$lib_folder/';\nuse lib '$lib_folder/$perl_mod_folder/';\n";
  	
	my @perl_libs_paths = ();
	if ( -e $perl_libs_paths_f_p and ! (-z $perl_libs_paths_f_p) ){
		@perl_libs_paths = list_to_array($perl_libs_paths_f_p,"NO_NEW_LINE");
		foreach  my $perl_libs_path (@perl_libs_paths){ 
			$new_start_string .= "use lib '$lib_folder/$perl_mod_folder/$perl_libs_path';\n";
		}
	}
	
	#Each .pl which should change its "use lib" will be modified
	foreach my $main_pl (@main_pls){
		print "$main_pl...";
		substitute_start_of_file($progDir."/".$main_pl,$start_string,$new_start_string);			
	}

	if ( $update > 0 ){ 
		print "Update Completed! You will now use VarGenius $vg_version\n";
	}else{
		print "Done! Now you can work with VarGenius-HZD from this folder!\n";
	}
}else {die "Change directory and start again the install script.\n";}
  

=head2 initial_config_vargenius
 Title   : initial_config_vargenius
 Usage   : initial_config_vargenius( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  sets the very important parameters for VarGenius asking them to the 
			user which install
 
 
 Returns : nothing
=cut
sub initial_config_vargenius{
	my $workDir = shift;
	
	####SECTION 1: DB AND USER ACCOUNT CONFIGURATION
	
	my $string = "#DATABASE configuration\n";
	my $input = "";

	
	#Folders: input, scratch and storage
	$input = "";
	print "\nNow you need to set the reference genome fasta path.\n";
	

	#References area
	$input = "";
	while ( length($input) == 0 ){
		print "\nPlease provide the path to the reference genome: ";
		$input= <STDIN>;
		chomp $input; #to cancel the return 
	}
	$string .= "\n#Reference genome\n";
	if ( $input eq 'NA'){
		$string .= "reference_f =\n";		
	}else{
		while ( !(-e $input) ){
			print "\n$input file does not exist. Please indicate an existing reference: ";
			$input= <STDIN>;
			chomp $input; #to cancel the return 
		}
		$string .= "reference_f = $input\n";
		#Get the reference folder path
		$reference_folder = extract_name($input,"path");
		$reference_file = $input;
	}
	
	
    #all the other tools are within the container
    #Getting singularity path
    my $singularity_path = `which singularity`;
    chomp($singularity_path);
    
    my $samtools_cmd = "";
    if ($singularity_path =~ /singularity/){

		#Input area
		$input = "";
		while ( length($input) == 0 ){
			print "\nSince you are using singularity, it needs to bind the area were the BAM files are localized. ".
			"Please type the folder where all the BAM files are located: ";
			$input= <STDIN>;
			chomp $input; #to cancel the return 
		}
		#$string .= "\n#Input area (for BAM files\n";
		#if ( $input eq 'NA'){
		#	$string .= "input_f =\n";		
		#}else{
			while ( ! (-d $input) ){
				print "\n$input folder does not exist. Please indicate an existing folder: ";
				$input= <STDIN>;
				chomp $input; #to cancel the return 
			}
			#$string .= "input_f = $input\n";
			$bam_folder = $input;
		#}
			
		my $pwd = `pwd`;
		chomp($pwd);
		my $singularity_call = " --bind ".extract_name($workDir,"path")." --bind $bam_folder --bind $reference_folder ";

		$singularity_call .= " $pwd/$lib_folder/SOFTWARE";
		
		print "\nSingularity call will be $singularity_call:\n";
		 
		$string .= "\n##PROGRAM PATHS\n";
		$string .= "R_path = $singularity_path run $singularity_call/".$cfg_hash->{'vargenius_sing_cont'}." R\n";
		$string .= "samtools_path = $singularity_path run $singularity_call/".$cfg_hash->{'vargenius_sing_cont'}." samtools\n";
		$string .= "bedtools_path = $singularity_path run $singularity_call/".$cfg_hash->{'vargenius_sing_cont'}." bedtools\n";	
		$samtools_cmd = "$singularity_path run $singularity_call/".$cfg_hash->{'vargenius_sing_cont'}." samtools";
		
	}else{
		print "WARNING: you don't have singularity tool as reccomended. Looking if you have installed R, samtools and bedtools...\n";
		my $samtools_path = `which samtools`;
		chomp($samtools_path);
		if ($samtools_path !~ /samtools/){
			die "samtools is not present. Please install it.."; 
		}
		my $bedtools_path = `which bedtools`;
		chomp($bedtools_path);
		if ($bedtools_path !~ /bedtools/){
			die "bedtools is not present. Please install it.."; 

		}
		my $R_path = `which R`;
		chomp($R_path);
		if ($R_path !~ /R/){
			die "R is not present. Please install it.."; 
		}
		
		#If everything is ok
		$string .= "\n##PROGRAM PATHS\n";
		$string .= "R_path = R\n";
		$string .= "samtools_path = samtools\n";
		$string .= "bedtools_path = bedtools\n";
		$samtools_cmd = "samtools";	
		
	}
	
	#Generate the genome file
	if( not -w $reference_folder){
		print "VarGenius-HZD tried to produce the genome file but $reference_folder is not writeable. You must do it manually ".
		" and add into the configuration file ($program_config)\n";
	}else{
		my $reference_f_ind =  $reference_file.".fai";
		my $genome_file = "$reference_file.genomefile";
		if (! -e $reference_f_ind){
			print "Indexing the genome...\n";
			my $command = "$samtools_cmd faidx $reference_file";
			try_exec_command($command) >0 or die "Error: Unable to execute: $command.\n" ;
		}
		print "Producing the genome file...\n";
		my $command = "awk -v OFS='".'\t'."' {'print ".'$1,$2'."'} $reference_f_ind > $genome_file ";
		try_exec_command($command) >0 or die "Error: Unable to execute: $command.\n" ;
		
		if ( -e $genome_file){
			$string .= "genome_file = $genome_file\n";
		}
	}
		
	print "Printing this text: $string\n";
	 #Now change these parameters within the program_config.txt
    substitute_text_in_file($program_config,"*******GROUP 1: START","*******GROUP 1: END",$string);
    
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



=head2 PBS_get_queues_info

 Title   : PBS_get_queues_info
 Usage   : PBS_get_queues_info(   );

 Function:  uses qstat to info about queue
 
 Returns : a string with status which can be
		
		Obtains the output of qstat -Q in JSON format for parsing
		and 
		
=cut
sub PBS_get_maximum_queue_resources{
	my $queue = shift;

	my $qstat_out_json = `qstat -Q -f -Fjson`;

	my $data = decode_json($qstat_out_json);
	#print Dumper $data;

	my $walltime = "";
	my $ncpus = "";
	my $mem = "";

	my $queues = $data->{'Queue'};

	if ( defined $data->{'Queue'}->{$queue}){
		my $q_data = $data->{'Queue'}->{$queue};
		if (defined $q_data->{'queue_type'}){
				#print "Queue: $queue\n";
				my $resources = $q_data->{'resources_max'};
				if (defined $resources->{'walltime'}){
					$walltime = $resources->{'walltime'};
					print "Maximum walltime: ".$resources->{'walltime'}."\n";
				}else{
					$walltime = $cfg_hash->{'qsub_walltime_max'};
				}
				if (defined $resources->{'ncpus'}){
					$ncpus = $resources->{'ncpus'};
					print "Maximum ncpus: ".$resources->{'ncpus'}."\n";	
				}else{
					$ncpus = $cfg_hash->{'qsub_ncpus'};
				}
				if (defined $resources->{'mem'}){
					$mem = uc($resources->{'mem'});
					print "Maximum mem: ".$resources->{'mem'}."\n";
				}else{
					$mem = $cfg_hash->{'qsub_mem'};
				}
												
		}		
	}else{
		die "ERROR: the queue $queue does not exists\n";
	} 
	return ($walltime, $ncpus, $mem);
}


=head2 PBS_get_queues_info

 Title   : PBS_get_queues_info
 Usage   : PBS_get_queues_info(   );

 Function:  uses qstat to info about queue
 
 Returns : a string with status which can be
		
		Obtains the output of qstat -Q in JSON format for parsing
		and 
		
=cut
sub PBS_get_queues_info{

	my $qstat_out_json = `qstat -Q -f -Fjson`;

	my $data = decode_json($qstat_out_json);
	#print Dumper $data;

   my $queues = $data->{'Queue'};
   foreach my $queue (keys %$queues){
        my $q_data = $data->{'Queue'}->{$queue};
        if (defined $q_data->{'queue_type'}){
			if ($q_data->{'queue_type'} eq "Execution"){
				print $queue."\n";
			}
        }
   }
}

=head2 my_extract_file_sys
 Title   : my_extract_file_sys
 Usage   : my_extract_file_sys( -input = the file compressed
                            - output = the file uncompressed
                               );

 Function:  takes in input the path of the compressed file and that of the file to create.
                  Then insert the file uncompressed where choosen
                It extracts by using system calls.

 Returns : nothing
=cut
sub my_extract_file_sys{
  my $input = shift;
  my $outDir = shift;
  
  my $command = '';
  my $outName = '';
  
  my $removeText = 'Please go in the folder and remove all the files coming out from this erroneous execution and try again.'.
  'If you get again an error please consider to unzip manually the files and leave them inside the folder. Then restart Annocript.';
  if ($input =~ /\.tar.gz$/){
      $command = "tar -zxf $input -C $outDir";
      try_exec_command($command) >0 or die "Error: Unable to execute: $command. $removeText \n" ;
    }elsif ($input =~ /\.gz$/) {
        $outName = extract_name($input,"gz");
        $outName = $outDir."/".$outName;
        $command = "gunzip -c $input > $outName";
        try_exec_command($command) >0 or die "Error: Unable to execute: $command. $removeText\n" ;
      }elsif ($input =~ /\.zip$/){
        $outName = extract_name($input,"zip");
        $outName = $outDir."/".$outName;
        $command = "gunzip -c $input > $outName";
        try_exec_command($command) >0 or die "Error: Unable to execute: $command. $removeText\n" ;
      }else{
        #extract_archive($input,$outDir);
      } 
      
  printf "Decompression of ".$input." finished\n";
}



=head2 extract_name

 Title   : extract_name
 Usage   : extract_name( -filePath => 'complete path of the file',
                        -type => 'a number saying what you want to extract'
			       );

 Function: extract the name from a complete path of a file. Even the file name only with the extension
              0: the complete name with the extension
              1: the name only
              2: the first two words joined by a dot
              noext: just remove the last extension from the path
              no2ext: remove the last 2 extensions from the path
              gz: the name from the .gz
              targz: the name from the .tar.gz
              zip: the name from the .zip
              tar: the name from the .tar
              fqgz: the name before two extensions
			path: the last element after the last slash is removed 
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
  elsif ($type eq "1"){ $name = $nameElements[0];}
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
      elsif ($type eq "path"){
		#The last element is removed
        $name= join '/', @list;
       }              
   else { die	"ERROR [$?]: $type is not a valid input extracting a name: ?\n";}
  return $name;
  
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

=head2 getInput

 Title   : getInput
 Usage   : getInput( - sentence: a sentence that will printed in input to ask something to the user;
					 - regex: what the answer of the user have to respect
                               );

 Function: Takes in input a sentence and a regex. Asks to the user the sentence and controls its input with regex
 
 Returns : input given by the user

=cut
sub getInput{
  my $sentence = shift;
  my $regex = shift;
  my $input='';
		
  while (!($input =~ /$regex/) or ($input eq "")){  
    print $sentence." ";
    $input = <STDIN>;
    chomp $input;
    print "\n";
  }	
  return $input;
}


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

=head2 file_list_to_array

 Title   : file_list_to_array
 Usage   : file_list_to_array(    );

 Function:  puts lines of a file inside an array and returns the array
 Returns : nothing

=cut
sub file_list_to_array{ 
	my $listFile = shift;
	
	my @list = ();
	
	open(LIST,"<$listFile") or die "Cannot open $listFile..\n";
	while ( my $row = <LIST>){
			chomp($row);
			if (-e $row ){
				push(@list,$row);
			}else{
				die "ERROR in $listFile: $row does not exist.\n ";
			}
	}
	close(LIST);
	return @list;
}


=head2 parse_command_line_args

 Title   : parse_command_line_args
 Usage   : parse_command_line_args(   );

 Function:  Parses the arguments specified upon the command line.
 Returns : nothing

=cut
sub parse_command_line_args{
  my $HELP  = 0;# Shows help overview.

	my $howToUse = "Use with: \nperl complete_install.pl \n\n".
	"\t-u|--update: Use this command if you just want to update VarGenius with the same dependencies you installed before\n";
 
  #  Parse options
  GetOptions(
           "help" => \$HELP,
           "u|update" => \$update );
	#print "UPDATE: $update\n";
	if ($update){$update=1;}else{$update=-1;}
	#print "UPDATE: $update\n";
  #Print a little help
  if ( $HELP ){
    #pod2usage(1);
    print $howToUse;
    exit;
  }

}



=head2 append_str_2_file_if_path_notexist

 Title   : append_str_2_file_if_path_notexist
 Usage   : append_str_2_file_if_path_notexist( - path -> the file where to write,
                                - string => the string to be written)

 Function: 
						will append the string at the file path and put a newline
						Both string and file list must be this format: 
							AnyInformationAlsoMultipleColumns\tpathToFile
						but only if the file does not exist or there is not the title
						
						NOTICE: The path MUST be always the last element in the row!

 Returns : nothing

=cut
sub append_str_2_file_if_path_notexist{
  my $file = shift;
  my $string = shift;
  
  #print "String: $string will be appended in $file..\n";
  #Check if the file exists
  if (-e $file){
		#Get all the paths of the files in an array
		open(DATA,"<$file") or die "Could not open file $file";
		my @pres = ();
		while (my $row = <DATA>){
			chomp($row);
			my @fields = split("\t",$row);
			#The path must be always the last element in the row!
			my $last_el_ind = scalar(@fields)-1;
			push(@pres,$fields[$last_el_ind]);
		}
		close(DATA);
		
		#print_array(\@pres);#DEBUGCODE
		#Insert the element if it does not exists
		#Get the file name from the string in input
		my @in_str_flds = split("\t",$string);
		my $last_el_ind = scalar(@in_str_flds)-1;
		my $filename = $in_str_flds[$last_el_ind];
		
		#Append to file only if not present
		if ( !(grep {/\Q$filename\E/} @pres) ){
			open(DATA,">>$file") or die "Couldn't open file $file";
			#print DATA "index: $last_el_ind $filename \n";#DEBUGCODE
			print DATA $string."\n";
			close(DATA);			
		}
	}else{
		open(DATA,">$file") or die "Couldn't open file $file";
		print DATA $string."\n";
		close(DATA);
	}
}



=head2 configFile2Hash

 Title   : configFile2Hash
 Usage   : configFile2Hash( - configFilePath = path of the config file
                             - configHash = the pointer to the hash to be filled
                               );

 Function:  gets the hash table with all the path and names in input from the config file in input
 Returns : nothing

=cut
sub configFile2Hash{
  my $configFilePath=shift;
  my ($configHash) = shift;

  my $start = 0;
  #Here we open config file and read all its line to find elements belonging to each of the executers
  open (configFile,$configFilePath) or die "ERROR: The file $configFilePath doesn't exists. The program will exit..\n";
  while (my $line = <configFile>){
		#This IF is useful if you want to put a description above in the text file. Delimit it with a set of hashtags
    if ($line =~ /#########/){
      $start = 1;
    }
   # if( ($line =~ /(\S+)\s*=\s*(.+)/) and ($start == 1) and !($line =~ /#/) ){
   if( ($line =~ /(\w+)\s*=\s*(.+)/) and ($start == 1) and !($line =~ /^#/) ){
      if ( $2 ne ''){
			  $$configHash->{$1} = $2;
				#annoPrint ("$1 = $2\n") ;#DEBUGCODE
			}
    #}elsif ( $line =~ /(\S+)\s*=/ ){
    }elsif ( $line =~ /(\w+)\s*=$/ ){
			delete($$configHash->{$1});
		}
  }
  close(configFile);
  #print Dumper\$configHash; #DEBUGCODE
}


=head2 checkConfigVariables

 Title   : checkConfigVariables
 Usage   : checkConfigVariables( - configFile -> file with the user configuration
                              - variablesFile -> the path to a file with all variables written
                              - lineToCheck -> line in the variables file to be used
          );

 Function: this subroutine reads the config files and check if all variables are there and are well written.
            The variables.txt file is needed fot this operation.

 Returns : nothing

=cut
sub checkConfigVariables {
  my $configFile = shift;
  my $variablesFile = shift;
  my $lineToCheck = shift;

  my $hashToCheck;

  if (! open(VARF,"<$variablesFile")){ die "ERROR: Failure opening '$variablesFile'. Your program version is corrupted - $!";}
  if (! open(CFILE,"<$configFile")){ die "ERROR: Cannot find '$configFile' - Your program version is corrupted - $!";}

  #Stores the variables in the config user file inside the hash
  my $start = 0;
  while (my $line = <CFILE>){
    #This code is to jump to the line where there are some ###
    if ($line =~ /#########/){
      $start = 1;
    }
    #Put variables in a hash
    if( ($line =~ /(\S+)\s*=/) and ($start == 1) and !($line =~ /#/) ){
     #annoPrint ($line."\n");#DEBUGCODE
     $hashToCheck->{$1} = "OK";
    }
  }
  close(CFILE);

  my @confVars = ();

  #Variables that are in the configuration file must be also in the variables.txt file
  my $errors=0;
  my $lines=0;
  #For each line of the variables file
  while (my $line = <VARF>){
    $line =~ s/\n//;#Remove \n in the end of the line

    #get the variables in the line
    my @variables = split (/;/,$line);

    $lines++;
    #For each of the variables in the variables file
    foreach my $var (@variables){
      if ($lines == $lineToCheck){
        #print "program Variable: $var - - value: ".$hashToCheck->{$var}."\n";#DEBUGCODE
        #put the variable inside an array for program config
        push (@confVars, $var);
        if( !(defined($hashToCheck->{$var})) ){

          die "ERROR: in $configFile variable $var is missing. Please check the file. Closing...\n ";
          $errors=1;
        }#else{ annoPrint ("From the hash: ".$hashCheck->{$var}."\n";)}
      }
    }
  }

  #print_array(\@allVars);
  #print Dumper\$hashCheck;
  #Now check if all the elements in the hash are also in the array
  foreach my $key (keys %$hashToCheck){
     # print "Search $key in array...\n";#DEBUGCODE
     if (!(grep {/$key/} @confVars )){
       die "ERROR: Variable $key is in the config file $configFile and not in $variablesFile file. This is completely wrong. Program will not work...\n ";
     }
  }
  #if ($errors == 0){annoPrint ("ok";}
  close(VARF);
}



=head2 list_to_array

 Title   : list_to_array
 Usage   : list_to_array( listFile = file path
													new_line = take or not the new line
										);

 Function:  puts lines of a file inside an array and returns the array.
					You can use it taking or not the new line character '\n' at the end of eac line
					by using the parameter new_line (NO_NEW_LINE, to take or nothing or NO to
					not take.
					DEFAULT: will take
 Returns : an array with the lines

=cut
sub list_to_array{ 
	my $listFile = shift;
	my $new_line = shift;
	
	my @list = ();
	if ( -e $listFile ){
		#you can do it using the newline or not
		if ($new_line eq 'NO_NEW_LINE'){
			open(LIST,"<$listFile") or die "Cannot open $listFile..\n";
			while ( my $row = <LIST>){
					chomp($row);
					push(@list,$row);
			}
			close(LIST);
		}else{
			#Reading all in one with <FILE> command
			open (FILE,"<$listFile") or die "Cannot open $listFile\n";
			@list = <FILE>;
			close(FILE);
		}		
	}

	return @list;
}
