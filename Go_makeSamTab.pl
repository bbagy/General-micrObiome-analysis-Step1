#!/usr/bin/perl -w
# 20211222
# creating count table from the sam files
# Go_makeSamTab.pl samDIR outfile
# 이네스가 보내준 코드를 바탕으로 완성.
use strict;
use Data::Dumper qw(Dumper);

use warnings 'all';

use File::Basename qw( basename );
use Getopt::Long   qw( GetOptions );

# add argument
sub usage {
   if (@_) {
      my ($msg) = @_;
      chomp($msg);
      print(STDERR "$msg\n");
   }

   my $prog = basename($0);
   print(STDERR "$prog --help for usage\n");
   exit(1);
}

sub help {
   my $prog = basename($0);
   print(STDERR "$prog [options] samdir outfile\n");
   print(STDERR "$prog --help\n");
   exit(0);
}

Getopt::Long::Configure(qw( posix_default ));  # Optional, but makes the argument-handling consistent with other programs.
GetOptions(
    'help|h|?' => \&help,
)
    or usage();

@ARGV == 2
    or usage("Incorrect number of arguments\n");

my ($in_dir) = $ARGV[0];
my ($out) = $ARGV[1];




# add dir infor below
system( "mkdir $in_dir/mapped");
system( "mkdir $in_dir/counts");


##################
# step1
##################
#Select seqs that have been mapped to the contigs

opendir (my $DIR, "$in_dir") or die "Cannot open directory Count_mapped_r\n";


my @files = grep {/^[\w\d]+/} readdir ($DIR);
closedir $DIR;

foreach my $files (@files) {
  open (my $infile, "$in_dir/$files") or die "Cannot open the file $files\n";
  open (my $outfile1, ">$in_dir/mapped/mapped_$files") or die "Cannot open $files\n";

  foreach my $line (<$infile>) {
    if ($line =~ /^\s*$/) {
                next;
   }
    	else {
              	#my @a;
                chomp ($line);
                if ($line =~ /^@/) {
                        next;
                }
                else {
                      	my @split = split ('\t',$line);
                        my $contig = $split[2];
                        my $name = $split[0];
                                if ($contig ne "*") {
                                        print $outfile1 "$name\t$contig\n";
                                }
                }
        }
  }
close ($infile) or die "Cannot close the infile $infile because $! \n";
close ($outfile1) or die "Cannot close the outfile $outfile1 because $! \n";
}

#Count how many times each contig is represented in the sample

opendir (my $DIR2, "$in_dir/mapped") or die "Cannot open directory count_Badol/m\n";
my @files2 = grep {/^[\w\d]+/} readdir ($DIR2);
closedir $DIR2;

foreach my $files2 (@files2) {
        open (my $infile2, "$in_dir/mapped/$files2") or die "Cannot open files a\n";
        open (my $outfile2, ">$in_dir/counts/counts_$files2") or die "Cannot ope\n";

        my %main_cnt;

        foreach my $line2 (<$infile2>) {

                #print "processing $files2\n";

                my @split2 = split ('\t',$line2);

                $main_cnt{$split2[1]} += 1;

        }
	while ((my $keys, my $values) = each(%main_cnt)){
                chomp ($keys);
                print $outfile2 "$keys\t$values\n";
        }


close ($infile2) or die "Cannot close the infile $infile2 because $! \n";
close ($outfile2) or die "Cannot close the outfile $outfile2 because $! \n";
}

print "\nStep1...done\n";

##################
# step2
##################
#Get unique elements from contig names array
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


#Get contig names into an array
opendir (my $DIR, "$in_dir/counts/") or die "Cannot open directory Make_1_table/infiles because $! \n";
my @files = grep {/^[\w\d]+/} readdir ($DIR);


closedir $DIR;

my @samples = "Contig\t";
my @matches;

foreach my $files (@files) {
    #print STDERR "\nRunning:\n $in_dir/counts/$files \n";
    
  open (my $infile, "$in_dir/counts/$files") or die "Cannot open the 1file $files because $! \n";
  push (@samples, "$files\t");
  foreach my $line (<$infile>) {
    if ($line =~ /^\s*$/) {
        next;
    }
    else {
        #my @a;
        chomp ($line);
        my @split = split ('\t',$line);
           my $name = $split[0];
#           print ("$name\n");
           push (@matches, "$name");
    }
  }
close ($infile) or die "Cannot close the infile $infile because $! \n";
}




my @unique_names = uniq (@matches);
#foreach my $id (@unique_names) {
#    chomp ($id);
#    print "$id\n";
#}

################################
### Print header on output file
open (my $outfile, ">$out") or die "Cannot open files at Count_mapped_reads/counts directory because $! \n";
print $outfile "@samples\n";


# Make an array containing contig names as keys and number of observations in the samples as values

foreach my $name (@unique_names) {
      chomp ($name);
      my @array;
      foreach my $files (@files) {
        my $counter=0;
        
        # print STDERR "\nRunning:\n $in_dir/counts/$files \n";
        
          open (my $infile, "$in_dir/counts/$files") or die "Cannot open the 2file $files because $! \n";
          foreach my $line (<$infile>) {
            if ($line =~ /^\s*$/) {
                next;
            }
            else {
                chomp ($line);
                my @split = split ('\t',$line);
#                print ("$split[0]\n");
                if ($name eq $split[0]) {
#                    print "match $name was found\n";
                    push (@array, "$split[1]\t");
                    $counter ++;
                }
            }
    }
    if ($counter == 0) {
        push (@array, "0\t");
    }
    close ($infile) or die "Cannot close the infile $infile because $! \n";
    
  }
print $outfile "$name\t@array\n";
}

close ($outfile) or die "Cannot close the outfile $outfile because $! \n";


print "\nStep2...done\n";
