#!/usr/bin/perl -w
# Extract junction read coverage and distribution for exon triplets
# The current version only works with a very strict exon triplet file format,
# will relax that later plus adding more error checks!
#
# Written by Leo J Lee, University of Toronto, 2013

use strict;
use Cwd;
use List::Util qw(sum);

my $Usage = "Usage: $0 TRI_file RD_file\n";
@ARGV == 2 or die $Usage;

my ($TRI_file, $RD_file) = @ARGV;
#die "$TRI_file, $RD_file";

my %junc_rd; my $pos = -1;
open (IN, $RD_file) or die ("cannot open $RD_file to read: $!");
while (<IN>)
{
    next if /^#/; chomp;
    my @line = split /\t/;
    my $junc_id = $line[0];
    my $rd_matrix = $line[2];
    $pos = $line[1] if ($pos eq "-1");
    my @rd = split /,/, $rd_matrix;
    $junc_rd{$junc_id} = [ @rd ];
}
close(IN);
my @tmp = split /,/, $pos;
my $n_pos = $#tmp;

my $file_out = "triplet_junc_rd.tsv";
my ($n_event, $n_detect) = (0, 0);
open (IN, "$TRI_file") or die ("cannot open $TRI_file to read: $!");
open (OUT, ">$file_out") or die ("cannot open $file_out to write: $!");
print OUT "#ID\tn_C1A\tn_AC2\tn_C1C2\tpositions\tRD_C1A\tRD_AC2\tRD_C1C2\n";
while (<IN>)
{
    next if (/^#/);
    $n_event++;
    chomp;
    my @line = split /\t/;
    @line >= 9 or die "Invalid exon triplet format: not enough tab delimited fields!";
    my ($event_id, $chr, $strand, $C1_s, $C1_e, $A_s, $A_e, $C2_s, $C2_e) = @line[0..8];
    my ($C1A, $AC2, $C1C2);
    if ($strand eq "+") {
	$C1A = "$C1_e:$A_s";
	$AC2 = "$A_e:$C2_s";
	$C1C2 = "$C1_e:$C2_s";
    } else {
	$C1A = "$A_e:$C1_s";
	$AC2 = "$C2_e:$A_s";
	$C1C2 = "$C2_e:$C1_s";
    }
    my ($n_C1A, $n_AC2, $n_C1C2) = (0, 0, 0);
    my (@rd_C1A, @rd_AC2, @rd_C1C2);
    foreach my $k (0..$n_pos) {
	$rd_C1A[$k] = 0; $rd_AC2[$k] = 0; $rd_C1C2[$k] = 0;
    }
    my $junc_id = "$chr:$C1A";
    if (exists $junc_rd{$junc_id}) {
	@rd_C1A = @{ $junc_rd{$junc_id} };
	@rd_C1A = reverse @rd_C1A if ($strand eq "-");
	$n_C1A = sum @rd_C1A;
    }
    $junc_id = "$chr:$AC2";
    if (exists $junc_rd{$junc_id}) {
	@rd_AC2 = @{ $junc_rd{$junc_id} };
	@rd_AC2 = reverse @rd_AC2 if ($strand eq "-");
	$n_AC2 = sum @rd_AC2;
    }
    $junc_id = "$chr:$C1C2";
    if (exists $junc_rd{$junc_id}) {
	@rd_C1C2 = @{ $junc_rd{$junc_id} };
	@rd_C1C2 = reverse @rd_C1C2 if ($strand eq "-");
	$n_C1C2 = sum @rd_C1C2;
    }
    print OUT "$event_id\t$n_C1A\t$n_AC2\t$n_C1C2\t$pos\t";
    foreach my $n (@rd_C1A) { print OUT "$n,"; }
    print OUT "\t";
    foreach my $n (@rd_AC2) { print OUT "$n,"; }
    print OUT "\t";
    foreach my $n (@rd_C1C2) { print OUT "$n,"; }
    print OUT "\n";
    $n_detect++ unless ($n_C1A+$n_AC2+$n_C1C2 == 0);
}
close(OUT); close(IN);
print "$n_detect out of $n_event exon triplets are detected.\n";
