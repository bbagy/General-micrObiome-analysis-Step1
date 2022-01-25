#!/usr/bin/perl
# 날짜를 잘 적어야 하는데...
# Gokbk 에서 braken 속도를 위해 뺐다. 
# run_kraken_macproV3.pl
# kraken2를 사용할때가 되었다.
# ./Gokraken2.pl -thread 4 -o kraken_out -b bracken_out -c classified -u unclassified -r kraken_report -log kraken2/kraken2_log.txt --db ${minikraken2}  fastq/*.fastq
# Apr 2nd 2020
# Added gzip for out put files



use warnings;
use strict;

use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $help; 
my $version = "0.1";
my $version_marker;

my $parallel = 1;
my $db ;
my $out_dir = "kraken_output";
my $log = "kraken_log";
my $contam = "contam_list";
my $filter1 = "filter1_files";
my $filter2 = "filter2_files";
#my $out_dir_classified = "classified_out";
#my $out_dir_unclassified = "unclassified_out" ;

#my $report = "report";

#added
#my $out_dir_report = "report";
#my $out_map_dir_report = "map_report";
#my $out_dir_bracken = "bracken_out" ;




my $keep;  

### delimiter used to identify sample name (either "_" or ".")
my $delimiter = "_";

### kraken2 parameters:
my $preload;
my $quick;
my $min_hits;
my $paired;
my $check_names;
my $classified_out;
my $unclassified_out;
my $threads;
my $mpa;


### initialized variables not user-specified
my %s2f = (); ### samples as keys and fastqs as values (2 for PE and 1 for SE)
my @log_files = (); ### array to contain all log files, which will be looped over at the end
my @out_files = (); ### out files needed for kraken-mpa-report
my @contam_list = (); ### out files  for contam list
my @filter1_files = (); ### out files  for filter1_files
my @filter2_files = (); ### out files  for filter2_files

my $kraken_options = ""; ### options to be passed to command-line

my $res = GetOptions("database|db=s"=>\$db,
		     "out_dir|o=s" => \$out_dir,
                     
             #"classified_out|c=s" => \$out_dir_classified,
             #"unclassified_out|u=s" => \$out_dir_unclassified,
             #"report|r=s" => \$out_dir_report,
             #"bracken|b=s" => \$out_dir_bracken,

		     "log=s"=>\$log,
		     "threads:i"=>\$parallel,
		     "preload"=>\$preload,
		     "quick"=>\$quick,
		     "min-hits:i"=>\$min_hits,
		     "paired"=>\$paired,
		     "check-names"=>\$check_names,
		     "keep"=>\$keep,
             #"report=s"=>\$report,
		     "delimiter=s"=>\$delimiter,
		     "help|h"=>\$help,
		     "version|v"=>\$version_marker,

	  )	or pod2usage(2);



pod2usage(-verbose=>2) if $help;

if ( $version_marker )	{	print "version $version\n";	exit	}

### the "2>&1" syntax means output STDERR to same output as STDOUT (i.e. same output as file descripter 1).
if ( index( `kraken2 2>&1` , "Need to specify input filenames!" ) == -1 )      { die "Stopping job, because \"kraken2\" is not in your path.\n";	}

if ( ! defined $db )	{	die "need database directory\n";	}

if ( ( $delimiter ne "_" ) and ( $delimiter ne "." ) )	{	die "delimiter $delimiter needs to be \"_\" or \".\"\n";	}

my @files = @ARGV;

### link sample names to fastqs
&fastqs2sample( @files );

pod2usage($0.': You must provide a list of fastq files.') unless @files;

# check that all files exist: 
foreach my $f ( @files )	{	 
	if ( ! -e $f )	{ die "file $f doesn't exist\n";	}	
}

system( "mkdir -p kraken2/$out_dir");
system( "mkdir -p kraken2/filtered_fastq");


# 확인된 option
if ( $threads )    {    $kraken_options = $kraken_options . " --threads ";    }
#if ( $classified_out )    {    $kraken_options = $kraken_options . " --classified_out";    }
#if ( $unclassified_out )    {    $kraken_options = $kraken_options . " --unclassified_out";    }
#if ( $report )    {    $kraken_options = $kraken_options . " --report";    }
if ( $mpa )    {    $kraken_options = $kraken_options . " --use-mpa-style";    }
if (  defined $min_hits )    {    $kraken_options = $kraken_options . " --min-hits $min_hits";    }
if ( $paired )    {    $kraken_options = $kraken_options . " --paired";    }


# 확인 안된 option
if (  $preload )	{	$kraken_options = $kraken_options . " --preload";	}
if (  $quick )	{	$kraken_options = $kraken_options . " --quick";	}
if ( $check_names )	{	$kraken_options = $kraken_options . " --check-names";	}

# test


foreach my $sample ( keys %s2f )	{
	
	my @fastqs = @{ $s2f{$sample} };
    
    my $cpu_count=1;
    #if the option is set
    if(defined($parallel)){
        #option is set but with no value then use the max number of proccessors
        if($parallel ==0){
            #load this module dynamically
            eval("use Sys::CPU;");
            $cpu_count=Sys::CPU::cpu_count();
        }else{
            $cpu_count=$parallel;
        }
    };
	my $out = "kraken2/$out_dir" . "/" . $sample . "_out.txt";
    my $log = "kraken2/$out_dir" . "/" . $sample . "_out.log";
    my $contam = "kraken2/$out_dir" . "/" . $sample . "_contam.txt";
    my $filter1 = "kraken2/filtered_fastq" . "/" . $sample . "_R1_filtered.fastq.gz";
    my $filter2 = "kraken2/filtered_fastq" . "/" . $sample . "_R2_filtered.fastq.gz";
    
 
	push( @out_files , $out );
	push( @log_files , $log );
    push( @contam_list , $contam );
    push( @filter1_files , $filter1 );
    push( @filter2_files , $filter2 );
    
	print STDERR "\nRunning:\nkraken2 --db $db --threads $cpu_count --output $out @fastqs 2> $log\n";
    system( "kraken2 --db $db --threads $cpu_count --output $out @fastqs 2> $log" );
    print STDERR "\nRunning:\n  $_ \n";
    
    print STDERR "\nRunning:\ngrep -v '\@' $out | cut -f2 > $contam\n";
    system( "grep -v '\@' $out | cut -f2 > $contam" ); #"grep -v "^@" $out | cut -f2 > $contam"
    
    
    print STDERR "\nRunning:\nseqkit grep -f $contam -v @fastqs -o $filter1 $filter2\n";
    #system( "seqkit grep -f $contam -v @fastqs -o $filter1 $filter2" );
    
    
    #system( "gzip $filter1" );
    #system( "gzip $filter2" );
	
}





### now combine all kraken output files into mpa report:
#print STDERR "\nRunning:\nmerge_metaphlan_tables.py kraken2/$out_dir_report/*_mpa.txt > kraken2/kraken2_mpa.txt\n";

#system( "merge_metaphlan_tables.py kraken2/$out_dir_report/*_mpa.txt > kraken2/kraken2_mpa.txt" );





### create single log file for all samples by looping through individual log files
open( 'MASTERLOG' , '>' , $log ) or die "cant create MASTERLOG file $log\n";
print MASTERLOG "sample	total_reads	classified	unclassified	classified_percent	unclassified_percent\n";

my $i = 0;

foreach my $err ( @log_files )	{

	### read in log file for kraken jobs:
	open( 'LOG' , '<' , $err ) or die "cant open LOG $err\n";
	chomp( my @lines = <LOG> );
	close( 'LOG' );
	
	my ($name,$path,$suffix) = fileparse( $err , (".fq" , ".fastq") );
	
	### scan 2nd last line of file for info on classified sequences:
	$lines[$#lines-1] =~ m/(\d+) sequences classified \((\S+)%\)/;
	my $classified = $1;
	my $classified_percent = $2;
	
	### scan last line of file for info on unclassified sequences:
	$lines[$#lines] =~ m/(\d+) sequences unclassified \((\S+)%\)/;
	my $unclassified = $1;
	my $unclassified_percent = $2;
	
	my $total = $classified + $unclassified;
	
	print MASTERLOG "$name	$total	$classified	$unclassified	$classified_percent	$unclassified_percent\n";

	### remove individual log files for each vsearch job:
	if ( ! $keep )	{	system( "rm $err" );	}
	++$i;
}

close( 'MASTERLOG' );


sub fastqs2sample {
	
	foreach my $f ( @_ )	{
		my @f_split = ();
		if ( $delimiter eq "_" )	{
			@f_split = split( '_' , basename($f) );
		} elsif ( $delimiter eq "." )	{
			@f_split = split( '\.' , basename($f) );
		}
		my $s = $f_split[0]; #sample
		my $base = basename($f);

		if ( ! exists $s2f{$s} )	{
			my @dummy = ();
			$s2f{$s} = \@dummy;
		}
		push( @{$s2f{$s}} , $f );
	}

	### keep track of #s of SE and PE samples, shouldn't be both types of samples
	my $se_marker = 0;
	my $pe_marker = 0;
	
	### check that if there are 2 fastqs per sample that they are R1 and R2. Also, no more than 2 fastqs per sample.
	foreach my $s ( keys %s2f )	{
		my @f = @{$s2f{$s}};
		my $num_fastq = scalar @f;
		if ( $num_fastq == 1 )	{
			++$se_marker;
			### only 1 fastq for this sample, assuming it's SE and moving on.
		} elsif ( $num_fastq == 2 )	{
			++$pe_marker;
			### assuming this sample has 2 fastqs because they are PE, check that either "_R1." or "_R1_" or "forward" is in one filename and "_R2." or "_R2_" or "reverse" is in the other's filename
			if ( ( basename($f[0]) =~ m/_R1\.|_R1_|forward/ ) and ( basename($f[1]) =~ m/_R2\.|_R2_|reverse/ ) ) {
				### the filenames seem to be PE reads and already ordered as forward and reverse, so continue
			} elsif ( ( basename($f[1]) =~ m/_R1\.|_R1_|forward/ ) and ( basename($f[0]) =~ m/_R2\.|_R2_|reverse/ ) ) 	{
				### filenames are PE, but reorder so they are forward and reverse
				my @tmp = ($f[1] , $f[0] );
				@{$s2f{$s}} = \@tmp;
			} else {
				die "the 2 fastqs for sample $s don't have PE read ids in their names, are you sure these are PE reads? They need _R1. and _R2. (or either _R1_ and _R2_ or forward and reverse) in their names. The filenames are: @f\n";
			}
		} elsif ( ($num_fastq < 1 ) or ( $num_fastq > 2 ) )	{
			die "$num_fastq fastqs for sample $s. This script assumed that fastq prefixes (when delimited by \"$delimiter\") are sample IDs, check your filenames\n";
		}
	}

	if ( ( $se_marker > 0 ) and ( $pe_marker > 0 ) )	{	
		die "Stopping job - $se_marker SE sample(s) and $pe_marker PE sample(s) input. You should run SE and PE samples through this script separately\n";	
	} elsif ( ( $se_marker > 0 ) and ( $paired ) )	{
		die "Stopping job - $se_marker SE sample(s), but '--paired' option used. You should run SE and PE samples through this script separately\n";
	} elsif ( ( $pe_marker > 0 ) and ( ! $paired ) ) 	{
		die "Stopping job - $pe_marker PE sample(s), but '--paired' option was not used.\n";
	}
}
