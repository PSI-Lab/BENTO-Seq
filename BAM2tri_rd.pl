#!/usr/bin/perl -w
# this is the master perl script to generate triplet junction read distributions
# from a BAM/SAM file or previosly processed all junction rd file
#
# Written by Leo J Lee, University of Toronto, 2013


use strict;
use Getopt::Std;

my $start = time();
my $Usage = "Usage: $0 [-o min_overhang -d edit_distance -m multi_hit] TRI_file [BAM/SAM_file]\n";
@ARGV >= 1 or die $Usage;
my %options=();
getopts("l:o:d:m:", \%options);

# parse input arguments
my ($TRI_file, $BAM_file) = @ARGV[-2..-1]; #print "$TRI_file\t$BAM_file\n";
my ($length, $overhang, $dist, $multi) = (75, 8, 2, 1);
$overhang = $options{o} if (defined $options{o});
$dist = $options{d} if (defined $options{d});
$multi = $options{m} if (defined $options{m});
#die "$BAM_file, $TRI_file, $length, $overhang, $dist, $multi\n";

# run the read processing pipeline
my $cmd; my $file_out = "all_junc_rd.tsv";
if ($BAM_file =~ /bam$/)
{
    my $tmp_SAM = "temp$$.sam";
    $cmd = "samtools view $BAM_file | head -n1 > $tmp_SAM"; system("$cmd");
    open (IN, $tmp_SAM) or die ("cannot open $tmp_SAM to read: $!");
    while (<IN>) { my @line = split /\t/; $length = length($line[9]); last; }
    close(IN); 
    unlink $tmp_SAM or warn "Cannot remove $tmp_SAM: $!";
    $cmd = "samtools view $BAM_file | ./filter_junc_read.pl | ./build_junc_rd.pl -l $length > $file_out";
    print "$cmd\n"; system("$cmd");
} elsif ($BAM_file =~ /sam$/) {
    open (IN, $BAM_file) or die ("cannot open $BAM_file to read: $!");
    while (<IN>) {
	next if (/^\@/);
	my @line = split /\t/; $length = length($line[9]);
	last;
    }
    close(IN);
    $cmd = "./filter_junc_read.pl $BAM_file | ./build_junc_rd.pl -l $length > $file_out";
    print "$cmd\n"; system("$cmd");
} else {
    $TRI_file = $BAM_file;
    if (-e $file_out) { 
	print "Using previoulsy processed junction read density file \"$file_out\" ...\n"; 
    } else { 
	print "Error: junction read density file \"$file_out\" not found!\n";
	print "Please provide a BAM/SAM file in addition to the TRI_File.\n"; 
	print $Usage; exit;
    }
}

unless ($TRI_file) { 
    print "Error: TRI_file not specified!\n"; print $Usage; exit; 
}
$cmd = "./extract_tri_junc_rd.pl $TRI_file $file_out";
print "$cmd\n"; system("$cmd");   

my $end = time();
my $diff = ($end-$start)/60;
#print "$start\t$end\n";
print "Finished processing $BAM_file in $diff minutes.\n";
