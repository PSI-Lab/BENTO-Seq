#!/usr/bin/perl -w
# Build inclusion and exclusion read distributions for different AS events.
# The supported event types are:
#     CAS: cassette or exon skipping events
#     A5SS: alternative 5' splice site
#     A3SS: alternative 3' splice site
#     MXE: mutually exclusive exons
#     AFE: alternative first exon
#     ALE: alternative last exon
#     SPR: splicing ratio (an upstream exon can splice with one of two 
#          downstream exons)
# For all types of events, exons must be specified in the order of 5' to 3'
#     (4 exons for MXE while 3 for all others)
#
# Written by Leo J Lee, University of Toronto, 2014

use strict;
use Cwd;
use List::Util qw(max min sum);

my $Usage = "Usage: $0 event_file RD_file\n";
my $FORMAT = "The event file must be tab-delimited with the following fields:\n";
$FORMAT = $FORMAT."\ttype: one of CAS, A5SS, A3SS, MXE, AFE, ALE and SPR\n";
$FORMAT = $FORMAT."\tID: a unique ID for the event\n";
$FORMAT = $FORMAT."\tchr: chromosome name matching the BAM file\n";
$FORMAT = $FORMAT."\tstrand: + or -\n";
$FORMAT = $FORMAT."\texon coordinates: in the form of start:end for each exon\n";
$FORMAT = $FORMAT."   Exons must be listed from 5' to 3' and tab-delimited\n";
$FORMAT = $FORMAT."   (4 for MXE and 3 for all other event types).\n";
@ARGV == 2 or die $Usage;
my ($event_file, $RD_file) = @ARGV;

my @types = qw(CAS A5SS A3SS MXE AFE ALE SPR);
my %support;
foreach my $type (@types) { $support{$type} = 1; }

my %junc_rd; 
my $pos = -1; #possible left start positions, same for all junctions
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
my @l_pos = split /,/, $pos;
# deduce min_overhang and read length from the left start positions
my $o = -max(@l_pos); my $l_r = $o-min(@l_pos);
my $n_pos = $#l_pos;

my $file_out = "event_inc_exc.rd";
my ($n_event, $n_detect) = (0, 0);
open (IN, "$event_file") or die ("cannot open $event_file to read: $!");
open (OUT, ">$file_out") or die ("cannot open $file_out to write: $!");
print OUT "#ID\tRD_inc\tRD_exc\n";
while (<IN>)
{
    next if (/^#/);
    next if (/^\s+$/); #skip blank lines
    chomp;
    my @line = split /\t/;
    @line >= 7 or die "Too few tab delimited fields in $event_file!\n$FORMAT";
    my ($type, $id, $chr, $strand, $e1, $e2, $e3) = @line[0..6];
    if (exists $support{uc($type)}) {
	$type = uc($type); $n_event++;
    } else {
	print "Event type $type is not supported for $id, skipping...\n";
	next;
    }

    my ($e1_s, $e1_e) = split /:/, $e1;
    my ($e2_s, $e2_e) = split /:/, $e2;
    my ($e3_s, $e3_e) = split /:/, $e3;
    # use negative coordinates on the minus strand for internal processing 
    if ($strand eq "-") { 
	$e1_s = -$e1_s; $e1_e = -$e1_e;
	$e2_s = -$e2_s; $e2_e = -$e2_e;
	$e3_s = -$e3_s; $e3_e = -$e3_e;
    }
    ($e1_s, $e1_e) = ($e1_e, $e1_s) if ($e1_s > $e1_e);
    ($e2_s, $e2_e) = ($e2_e, $e2_s) if ($e2_s > $e2_e);
    ($e3_s, $e3_e) = ($e3_e, $e3_s) if ($e3_s > $e3_e);
    unless (($e2_s > $e1_s or $e2_e > $e1_e) and ($e3_s > $e2_s or $e3_e > $e2_e)) {
	print "Exons must be listed from 5' to 3' on the transcribed strand!\n";
	print "Skipping $id due to invalid exon list...\n";
	$n_event--; next;
    }
    my ($e1e2, $e1e3, $e2e3);
    if ($strand eq "+") {
	$e1e2 = "$e1_e:$e2_s";
	$e2e3 = "$e2_e:$e3_s";
	$e1e3 = "$e1_e:$e3_s";
    } else {
	$e1e2 = "$e2_s:$e1_e"; $e1e2 =~ s/\-//g;
	$e2e3 = "$e3_s:$e2_e"; $e2e3 =~ s/\-//g;
	$e1e3 = "$e3_s:$e1_e"; $e1e3 =~ s/\-//g;
    }
    my $l_e1 = $e1_e - $e1_s + 1;
    my $l_e2 = $e2_e - $e2_s + 1;
    my $l_e3 = $e3_e - $e3_s + 1;

    my ($e4, $e4_s, $e4_e, $e2e4, $e3e4, $l_e4);
    if ($type eq "MXE") {
	$e4 = $line[7]; ($e4_s, $e4_e) = split /:/, $e4;
	if ($strand eq "-") { $e4_s = -$e4_s; $e4_e = -$e4_e; }
	($e4_s, $e4_e) = ($e4_e, $e4_s) if ($e4_s > $e4_e);
	unless ($e4_s > $e3_s or $e4_e > $e3_e) {
	    print "Exons must be listed from 5' to 3' on the transcribed strand!\n";
	    print "Skipping $id due to invalid exon list...\n";
	    $n_event--; next;
	}
	if ($strand eq "+") { 
	    $e2e4 = "$e2_e:$e4_s"; $e3e4 = "$e3_e:$e4_s";
	} else {
	    $e2e4 = "$e4_s:$e2_e"; $e2e4 =~ s/\-//g;
	    $e3e4 = "$e4_s:$e3_e"; $e3e4 =~ s/\-//g;
	}
	$l_e4 = $e4_e - $e4_s + 1;
    }
    #print "$id\n"; print "e1e2: $e1e2\ne2e3: $e2e3\ne1e3: $e1e3\n";
    my (@rd_e1e2, @rd_e2e3, @rd_e1e3);
    foreach my $k (0..$n_pos) {
	$rd_e1e2[$k] = 0; $rd_e2e3[$k] = 0; $rd_e1e3[$k] = 0;
    }
    my $junc_id = "$chr:$e1e2";
    if (exists $junc_rd{$junc_id}) {
	@rd_e1e2 = @{ $junc_rd{$junc_id} };
	@rd_e1e2 = reverse @rd_e1e2 if ($strand eq "-");
    }
    $junc_id = "$chr:$e2e3";
    if (exists $junc_rd{$junc_id}) {
	@rd_e2e3 = @{ $junc_rd{$junc_id} };
	@rd_e2e3 = reverse @rd_e2e3 if ($strand eq "-");
    }
    $junc_id = "$chr:$e1e3";
    if (exists $junc_rd{$junc_id}) {
	@rd_e1e3 = @{ $junc_rd{$junc_id} };
	@rd_e1e3 = reverse @rd_e1e3 if ($strand eq "-");
    }

    # Process different types of events
    my (@rd_inc, @rd_exc);
    if ($type eq "CAS") {
	unless ($e2_s > $e1_e and $e3_s > $e2_e) {
	    print "$id is not a valid CAS event, skipping...\n";
	    $n_event--; next;
	}
	#avoid double counting reads spanning two splice junctions
	my $l = $l_r - 2*$o - $l_e2;
	if ($l>0) { splice @rd_e1e2, -$l; }
	@rd_inc = (@rd_e1e2, @rd_e2e3);
	@rd_exc = @rd_e1e3;
	
    } elsif ($type eq "A5SS") {
	unless ($e2_s == $e1_s and $e3_s > $e2_e) {
	    print "$id is not a valid A5SS event, skipping...\n";
	    $n_event--; next;
	}
	@rd_inc = @rd_e2e3;
	@rd_exc = @rd_e1e3;
    } elsif ($type eq "A3SS") {
	unless ($e2_s > $e1_e and $e3_e == $e2_e) {
	    print "$id is not a valid A3SS event, skipping...\n";
	    $n_event--; next;
	}
	@rd_inc = @rd_e1e2;
	@rd_exc = @rd_e1e3;
    } elsif ($type eq "MXE") {
	unless ($e2_s > $e1_e and $e3_s > $e2_e and $e4_s > $e3_e) {
	    print "$id is not a valid MXE event, skipping...\n";
	    $n_event--; next;
	}
	my (@rd_e2e4, @rd_e3e4);
	foreach my $k (0..$n_pos) { $rd_e2e4[$k] = 0; $rd_e3e4[$k] = 0; }
	$junc_id = "$chr:$e2e4";
	if (exists $junc_rd{$junc_id}) {
	    @rd_e2e4 = @{ $junc_rd{$junc_id} };
	    @rd_e2e4 = reverse @rd_e2e4 if ($strand eq "-");
	}
	$junc_id = "$chr:$e3e4";
	if (exists $junc_rd{$junc_id}) {
	    @rd_e3e4 = @{ $junc_rd{$junc_id} };
	    @rd_e3e4 = reverse @rd_e3e4 if ($strand eq "-");
	}
	#avoid double counting reads spanning two splice junctions
	my $l = $l_r - 2*$o - $l_e2;
	if ($l>0) { splice @rd_e1e2, -$l; }
	$l = $l_r - 2*$o - $l_e3;
	if ($l>0) { splice @rd_e1e3, -$l; }
	@rd_inc = (@rd_e1e2, @rd_e2e4);
	@rd_exc = (@rd_e1e3, @rd_e3e4);
    } elsif ($type eq "AFE") {
	unless ($e3_s > $e1_e and $e3_s > $e2_e) {
	    print "$id is not a valid AFE event, skipping...\n";
	    $n_event--; next;
	}
	#reads shouldn't extend beyond the first exon start position
	my $l = $l_r - $o - $l_e1;
	if ($l>0) { splice @rd_e1e3, 0, $l; }
	$l = $l_r - $o - $l_e2;
	if ($l>0) { splice @rd_e2e3, 0, $l; }
	@rd_inc = @rd_e2e3;
	@rd_exc = @rd_e1e3;
    } elsif ($type eq "ALE") {
	unless ($e2_s > $e1_e and $e3_s > $e1_e) {
	    print "$id is not a valid ALE event, skipping...\n";
	    $n_event--; next;
	}
	#reads shouldn't extend beyond the last exon end position
	my $l = $l_r - $o - $l_e2;
	if ($l>0) { splice @rd_e1e2, -$l; }
	$l = $l_r - $o - $l_e3;
	if ($l>0) { splice @rd_e1e3, -$l; }
	@rd_inc = @rd_e1e2;
	@rd_exc = @rd_e1e3;
    } elsif ($type eq "SPR") {
	unless ($e2_s > $e1_e and $e3_s > $e1_e) {
	    print "$id is not a valid SPR event, skipping...\n";
	    $n_event--; next;
	}
	@rd_inc = @rd_e1e2;
	@rd_exc = @rd_e1e3;
    } else {
	die "This shouldn't happen! Please debug $0";
    }
	
    print OUT "$id\t";
    foreach my $n (@rd_inc) { print OUT "$n,"; }
    print OUT "\t";
    foreach my $n (@rd_exc) { print OUT "$n,"; }
    print OUT "\n";
    $n_detect++ unless (sum(@rd_inc)+sum(@rd_exc) == 0);
}
close(OUT); close(IN);
print "$n_detect out of $n_event events are expressed.\n";
