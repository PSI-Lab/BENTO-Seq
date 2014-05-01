#!/usr/bin/perl -w
# build the read distribution matrix for each junction
# Options: -l: read length (default 75); -o: minimum overhang (default 5)
#
# Reads are indexed by left start position (on the plus strand), indels 
# (and reads with slightly different length?) are handled approximately
# (all positions here are 1-based!)
#
# Written by Leo J Lee, University of Toronto, 2014

use strict;
use Getopt::Std;

my %options=();
getopts("l:o:", \%options);

my ($r_len, $min_overhang) = (75, 5);
$r_len = $options{l} if (defined $options{l});
$min_overhang = $options{o} if (defined $options{o});
 
my %junc_rd;
my @x_rd; #a list of possible left start positions (all negative)
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
	my $l_start = -$left; $junc_rd{$junc_key}{$l_start}++;
	#print "$pos_junc\t$left\t$right\t$balance\n" if ($n_junc==2);
    }
    #last if ($n_junc==2);
}

print "#junction_ID\tpositions\tread_distribution\n"; 
# sort the junctions by chromosome first, then start position and end position
my @junc_keys = sort {
    my ($chr1, $start1, $end1) = split(/:/, $a, 3);
    my ($chr2, $start2, $end2) = split(/:/, $b, 3);
    if ($chr1 eq $chr2)
    {
	if ($start1 == $start2) {
	    return ($end1 <=> $end2);
	} else {
	    return ($start1 <=> $start2);
	}
    }
    else
    {
	return ($chr1 cmp $chr2);
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
