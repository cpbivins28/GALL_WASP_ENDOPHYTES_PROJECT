#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use File::Basename;

sub reverse_complement($);
sub trim($);
sub usage;
sub isint($);

my @File_array;

my $opt_string='hd:p:n:l:t:s:m:c:';
my $prg = basename $0;
my %opt;
my $Seq_directory;
my $Trim_path;
my $Threads;
my $Lead;
my $Trail;
my $Sliding_window;
my $Min_length;
my $Headcrop;

getopts("$opt_string", \%opt) or usage();
$opt{h} and usage();
usage() unless ($opt{d} && $opt{p});
print "Could not read directory: $opt{d}\n" and die(usage()) unless (-e $opt{d});
print "File does not exist: $opt{p}\n" and die(usage()) unless (-e $opt{p});
if($opt{n}){print "Number of threads must be an integer\n" and die(usage()) unless(isint($opt{n}));}
if($opt{c}){print "Headclip argument must be an integer\n" and die(usage()) unless(isint($opt{c}));}
if($opt{l}){print "Leading trim quality score minimum must be an integer\n" and die(usage()) unless(isint($opt{l}));}
if($opt{t}){print "Trailing trim quality score minimum must be an integer\n" and die(usage()) unless(isint($opt{t}));}
if($opt{s}){print "Sliding window parameter must follow the format (Integer):(Integer)" and die(usage()) unless($opt{s} =~ m/\d+:\d+/);}
if($opt{m}){print "Minimum sequence length must be an integer\n" and die(usage()) unless(isint($opt{m}));}

$Seq_directory = $opt{d};
$Trim_path= $opt{p};
$Threads = $opt{n} ? $opt{n} : 1;
$Lead = $opt{l} ? $opt{l} : 20;
$Trail = $opt{t} ? $opt{t} : 20;
$Sliding_window = $opt{s} ? $opt{s} : "4:20";
$Min_length = $opt{m} ? $opt{m} : 200;
$Headcrop = $opt{c} ? $opt{c} : 10;

print STDERR "Creating trimmomatic command script for all .fastq files stored in the directory: $opt{d}\n";
print STDERR "Path to trimmomatic .jar file (including filename): $Trim_path\n";
print STDERR "Threads for each analysis: $Threads\n";
print STDERR "Number of bases at the head of the sequence: $Headcrop\n";
print STDERR "Minimum quality score for bases at the leading edge of the sequence: $Lead\n";
print STDERR "Minimum quality score for bases at the trailing edge of the sequence: $Trail\n";
print STDERR "Sliding_window parameter set to: $Sliding_window\n";
print STDERR "Minimum sequence length to be retained after trimming: $Min_length\n";

opendir(SEQDIR, $Seq_directory) or die "can't open directory: $Seq_directory $!\n";

if ($Seq_directory !~ /\/$/){
    $Seq_directory .= "/";
}

my %pairs;

foreach my $f (readdir (SEQDIR)){
    next if ($f =~ /^\./);
    next if ($f !~ /\.fastq$/);
    
    if ($f =~ /_R1_/) {
        my $pair = $f;
        $pair =~ s/_R1_/_R2_/;
        if (-e "$Seq_directory$pair") {
            push(@File_array, [$f, $pair]);
        } else {
            print STDERR "Warning: Missing pair for $f\n";
        }
    }
}

@File_array = sort { $a->[0] cmp $b->[0] } @File_array;

my $Output_script = "my_Make_Trimmomatic_forAMPtk.sh";
open(OUTPUT, ">$Output_script");
print OUTPUT "#!/bin/bash\n";
print OUTPUT "mkdir -p my_Trimmomatic_OUTPUT\n";
print OUTPUT "PID=\"\"\n";

foreach my $pair (@File_array) {
    my ($f1, $f2) = @$pair;
    my $sample_id = $f1;
    $sample_id =~ s/_R1_.*//;
    my $infile1 = $Seq_directory . $f1;
    my $infile2 = $Seq_directory . $f2;
    print OUTPUT "java -jar $Trim_path PE -threads $Threads $infile1 $infile2 my_Trimmomatic_OUTPUT/${sample_id}_R1_paired.fastq my_Trimmomatic_OUTPUT/${sample_id}_R1_unpaired.fastq my_Trimmomatic_OUTPUT/${sample_id}_R2_paired.fastq my_Trimmomatic_OUTPUT/${sample_id}_R2_unpaired.fastq HEADCROP:$Headcrop LEADING:$Lead TRAILING:$Trail SLIDINGWINDOW:$Sliding_window MINLEN:$Min_length &\n";
    print OUTPUT "PID=\$!\n";
    print OUTPUT "wait \$PID\n";
    print OUTPUT "cat my_Trimmomatic_OUTPUT/${sample_id}_R1_paired.fastq my_Trimmomatic_OUTPUT/${sample_id}_R1_unpaired.fastq | gzip > my_Trimmomatic_OUTPUT/${sample_id}_R1_ALL.fastq.gz\n";
    print OUTPUT "cat my_Trimmomatic_OUTPUT/${sample_id}_R2_paired.fastq my_Trimmomatic_OUTPUT/${sample_id}_R2_unpaired.fastq | gzip > my_Trimmomatic_OUTPUT/${sample_id}_R2_ALL.fastq.gz\n";
}

close(OUTPUT);

sub isint($){
  my $val = shift;
  return ($val =~ m/^\d+$/);
}

sub reverse_complement($) {
    my $dna = shift;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/NACGTacgt/NTGCAtgca/;
    return $revcomp;
}

sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

sub usage
{
    print STDERR << "EOF";

Name $prg - This Perl script will print a .sh script that can be executed from the command line
to process a folder full of paired-end .fastq read files 
with Trimmomatic. As argument this script takes the directory of fasta 
files, the path to the Trimmomatic software and several parameters to use for 
running Trimmomatic. Naming conventions are expected in the fastq files 
in order to process the folder correctly.
Files must be named as follows:
Unique-sample-ID_any-text_R(1|2)_any-text.fastq
The placement of "_" characters is significant. Everything before the first "_" character
is considered as the unique sample ID. the rest of the filename can follow any convention
but the string "_R1_" or "_R2_" must be present. This string identifies members of paired
read files where R1 = (forward reads) and R2 = (reverse reads). Read files should be placed
as pairs into a single directory and read files with the same ID, and the R1 and R2 strings
will be paired in commands generated for the .sh script
Output of the shell script produced by this script creates files that are named using the
conventions expected by AMPtk software. The forward_paired and forward_unpaired files for
Trimmomatic command will also be concatenated into a single file by this script.

NOTE: file names of fastq files must end in .fastq

usage: $prg -d directory -p path [-c clip_bases -n threads -l leading_min_score -t trailing_min_score -s sliding-window -m min-length]

-h      : print this help message
-d directory : Path to a directory which contains paired end fastq sequence files in .fastq format
-p path : Full path to the Trimmomatic software
-n threads : The number of threads to execute under trimmomatic commands (default: 1)
-c clip_bases : The number of bases to trim off the front of the read before proceeding with any other 
trimming command (default: 10).
-l min_score : The minimum quality score to keep when cutting from the leading edge of 
a sequence (default: 20)
-t min_score : The minimum quality score to keep when cutting from the trailing edge of 
a sequence (default: 20)
-s sliding_window : The parameters for sliding window trimming (default: 4:20)
-m min_length : The minimum length allowable when keeping a trimmed sequence. (default: 200)

ex: $prg -d Some_directory/ -p path_to_trimmomatic
This will pair all files in "Some_directory" into Trimmomatic commands using the default setting above.

EOF
    exit;
}

