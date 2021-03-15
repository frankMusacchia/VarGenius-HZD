#!/usr/bin/perl

use strict;
use warnings;


##Get the file with the parameters for the module
my $main_data_folder = $ARGV[0];
my $target_bed = $ARGV[1];
my $ucsc_genes = $ARGV[2];
my $ucsc_exons = $ARGV[3];
my $samplename = $ARGV[4];
my $bam_path = $ARGV[5];
my $target_exons = $ARGV[6];
my $hum_ref = $ARGV[7];
my $bedtools_path = $ARGV[8];
my $log_file = $ARGV[9];
my $out_folder = extract_name($bam_path,"path");


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
		
					
if ( -e $target_exons and !(-z $target_exons) ){
	
	my $target_name = extract_name($target_bed,1);
	#2: Get exon coverage
	my $target_exons_cov = $out_folder."/$samplename\_".extract_name($target_exons,1).".$bed_ext";
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
	my $covbedparams = " -g ".extract_name($hum_ref,"noext").".genomefile -sorted ";
	run_BEDTOOLS_coverageBed_gen($target_exons,$bam_path,$target_exons_cov,$covbedparams,$log_file);
	
	
	#3: Restrict to ipotetically low coverage exons only
	my $nc_threshold = "0.2";#percent of interval covered
	my $target_exons_nc = extract_name($target_exons_cov,"noext")."_nc.".$bed_ext;
	print_and_log( "Get only the chrom and interval, sort and get unique lines...",$log_file);#DEBUGCODE
	my $command = "sed '".'s/\./,/'."' $target_exons_cov | awk "."'".'$7<'.$nc_threshold."'  > $target_exons_nc";
	print_and_log( "\nExecuting command: $command\n",$log_file);#DEBUGCODE
	try_exec_command($command) or die "Unable to execute command: $command\n";	

	my $noncovexons_out = extract_name($target_exons_cov,"noext").$noncov_ext;

	#4: intersection with ucsc genes (intersectBed)
	#I use -wo because I want the dimension of the intersection and the Gene name.
	#But I don't want the regions that do not have intersection with Genes because I will not use them
	my $inters_params = " -wo ";
	print "...intersecting UCSC exons with Genes...";#DEBUGCODE
	run_BEDTOOLS_intersectBed($target_exons_nc,$ucsc_genes,$noncovexons_out,$inters_params,$log_file);#######	

	print_and_log( "...DONE!",$log_file);#DEBUGCODE
	if ( $remove_temp eq 'YES' ){
		delete_file($target_exons_cov);
		delete_file($target_exons_nc);
	}
	
}else{
	print "ERROR: $target_exons file does not exist...";#DEBUGCODE
}



#######################UTILITIES



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
  elsif ($type eq "path"){ 
                 $name = join("/",@list);
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
