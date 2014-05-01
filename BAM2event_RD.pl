#!/usr/bin/perl -w
# this is the master perl script to generate inclusion/exclusion read 
# distributions from a BAM/SAM file or previosly processed junction RD file
#
# Written by Leo J Lee, University of Toronto, 2014


use strict;
use Getopt::Std;

my $start = time();
my $Usage = "Usage: $0 [-o -d -m] event_file [BAM/SAM_file]\n";
$Usage = $Usage."\t-o: minimum overhang on each side of the junction, default 5\n";
$Usage = $Usage."\t-d: maximum edit distance when mapping a read, default 2\n";
$Usage = $Usage."\t-m: maximum number of mapped positions, default 1 (uniquely mapped)\n";

@ARGV >= 1 or die $Usage;
my %options=();
getopts("o:d:m:", \%options);

# parse input arguments
my ($event_file, $BAM_file) = @ARGV[-2..-1]; #print "$event_file\t$BAM_file\n";
my ($length, $overhang, $dist, $multi) = (0, 5, 2, 1);
$overhang = $options{o} if (defined $options{o});
$dist = $options{d} if (defined $options{d});
$multi = $options{m} if (defined $options{m});
#die "$BAM_file, $event_file, $length, $overhang, $dist, $multi\n";

# run the read processing pipeline
my $cmd; my $file_out = "all_juncs.rd";
if ($BAM_file =~ /bam$/)
{
    #obtain read length from the BAM file
    my $tmp_SAM = "temp$$.sam";
    $cmd = "samtools view $BAM_file | head -n1 > $tmp_SAM"; system("$cmd");
    open (IN, $tmp_SAM) or die ("cannot open $tmp_SAM to read: $!");
    while (<IN>) { my @line = split /\t/; $length = length($line[9]); last; }
    close(IN); 
    unlink $tmp_SAM or warn "Cannot remove $tmp_SAM: $!";
    $cmd = "samtools view $BAM_file | ./filter_junc_read.pl | ./build_junc_RD.pl -l $length > $file_out";
    print "$cmd\n"; system("$cmd");
} elsif ($BAM_file =~ /sam$/) {
    #obtain read length from the SAM file
    open (IN, $BAM_file) or die ("cannot open $BAM_file to read: $!");
    while (<IN>) {
	next if (/^\@/);
	my @line = split /\t/; $length = length($line[9]);
	last;
    }
    close(IN);
    $cmd = "./filter_junc_read.pl $BAM_file | ./build_junc_RD.pl -l $length > $file_out";
    print "$cmd\n"; system("$cmd");
} else {
    $event_file = $BAM_file;
    if (-e $file_out) { 
	print "Using previoulsy processed junction read distribution file \"$file_out\" ...\n"; 
    } else { 
	print "Error: junction read distribution file \"$file_out\" not found!\n";
	print "Please provide a BAM/SAM file in addition to the event_file.\n"; 
	print $Usage; exit;
    }
}

unless ($event_file) { 
    print "Error: event_file not specified!\n"; print $Usage; exit; 
}
$cmd = "./extract_event_RD.pl $event_file $file_out";
print "$cmd\n"; system("$cmd");   

my $end = time();
my $diff = ($end-$start)/60;
#print "$start\t$end\n";
print "Finished processing $BAM_file in $diff minutes.\n";
