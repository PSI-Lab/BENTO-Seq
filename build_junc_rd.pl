#!/usr/bin/perl -w
# build the read distribution matrix for each junction
# Options: -l: read length (default 75); -o: minimum overhang (default 8)
#
# Now index reads by left start position instead of the balance factor so that
# indels (as well as slightly different read lengths) can be handled better
# (although not perfectly!)
#
# Written by Leo J Lee, University of Toronto, 2013

use strict;
use Getopt::Std;

my %options=();
getopts("l:o:", \%options);

my ($r_len, $min_overhang) = (75, 8);
$r_len = $options{l} if (defined $options{l});
$min_overhang = $options{o} if (defined $options{o});

# Only include junctions on the reference genome at the moment
my @chrs = qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM);
my %ref_chr;
foreach my $chr (@chrs) { $ref_chr{$chr} = 1; }
 
my %junc_rd;
my @x_rd;
# my $x_max = $r_len - $min_overhang - $min_overhang;
# my $x_min = -$x_max;
# my ($i, $x) = (0, $x_min);
# while ($x <= $x_max) {
#     $x_rd[$i] = $x;
#     $i++; $x += 2;
# }
my ($i, $x_min, $x_max) = (0, $min_overhang - $r_len, -$min_overhang);
foreach my $x ($x_min..$x_max) {
    $x_rd[$i] = $x; $i++;
}
#die "@x_rd\n";

# make this work like a standard unix utility
while (<>) 
{
    next if /^#/; chomp;
    my @line = split /\t/;
    my ($chr, $n_junc) = ($line[3], $line[5]);
    next unless (exists $ref_chr{$chr});
    #print "$_\n" if ($n_junc==2);
    my @pos = split /,/, $line[6]; my @seg;
    foreach my $i (1..$n_junc+1) {
	$seg[$i] = $pos[2*$i-1] - $pos[2*$i-2] + 1;
    }
    foreach my $i (1..$n_junc) {
	my $junc_pos = "$pos[2*$i-1]:$pos[2*$i]";
	my $junc_key = "$chr:$junc_pos";
	unless (exists $junc_rd{$junc_key}) {
	    foreach my $x (@x_rd) {
		$junc_rd{$junc_key}{$x} = 0;
	    }
	}
	my ($left, $right) = (0, 0);
	foreach my $j (1..$i) {
	    $left += $seg[$j];
	}
	foreach my $j ($i+1..$n_junc+1) {
	    $right += $seg[$j];
	}
	# my $balance = $right - $left;
	# $junc_rd{$junc_key}{$balance}++;
	my $l_start = -$left; $junc_rd{$junc_key}{$l_start}++;
	#print "$pos_junc\t$left\t$right\t$balance\n" if ($n_junc==2);
    }
    #last if ($n_junc==2);
}

print "#junction\tpositions\treads\n"; 
# sort the junctions by chromosome first, then start position and end position
my @junc_keys = sort {
    my ($chr1, $start1, $end1) = split(/:/, $a, 3);
    my ($chr2, $start2, $end2) = split(/:/, $b, 3);
    $chr1 =~ s/chr//; $chr2 =~ s/chr//;
    if ($chr1 eq $chr2)
    {
	if ($start1 == $start2) {
	    return ($end1 <=> $end2);
	} else {
	    return ($start1 <=> $start2);
	}
    }
    elsif ($chr1 eq "X")
    {
	if ($chr2 eq "Y" or $chr2 eq "M") { return -1; }
	else { return 1; }
    }
    elsif ($chr1 eq "Y")
    {
	if ($chr2 eq "M") { return -1; }
	else { return 1; }
    }
    elsif ($chr1 eq "M")
    {
	return 1;
    }
    elsif ($chr2 eq "X" or $chr2 eq "Y" or $chr2 eq "M")
    {
	return -1;
    }
    else
    {
	return ($chr1 <=> $chr2);
    }
} keys %junc_rd;
foreach my $junc_key (@junc_keys) 
{
    print "$junc_key\t";
    foreach my $x (@x_rd) {
	print "$x,";
    }
    print "\t";
    foreach my $x (@x_rd) {
	print "$junc_rd{$junc_key}{$x},";
    }
    print "\n";
}
