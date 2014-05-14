#!/usr/bin/perl -w
# filter and reformat junction reads
# Options: 
#    -d: maximum edit distance (default 2, include mismatch, indel and clipping)
#    -m: maximum multiple hits (default 1)
# (could use the optional MD field if that becomes widely available later on)
# The current version works for TopHat and STAR BAM files
#
# Written by Leo J Lee, with help from Alice Gao
# PSI lab, University of Toronto, 2013

use strict;
use Getopt::Std;

my %options=();
getopts("d:m:", \%options);

my ($max_NM, $max_NH) = (2, 1);
# maximum edit distance allowed
$max_NM = $options{d} if (defined $options{d});
# maximum mapping positions allowed
$max_NH = $options{m} if (defined $options{m});
my $header = "#Read\tNM\tNH\tchr\tstrand\tn_junc\tpos\tCIGAR\n";
print $header;

# make this work like a standard unix utility
while (<>) 
{
    chomp;
    next if (/^\@/);
    my @line = split /\t/;
    next unless ($line[5] =~ /(\d+)N/);
    
    # join all optional tags and do a preliminary check
    my $optional_tags = join(",", @line[11..$#line]);
    #use [Nn] to accomondate both STAR and Tophat outputs
    $optional_tags =~ /([Nn]M):i:(\d+)/ or die "Incompatible BAM/SAM format: optional TAG [Nn]M is not present!";
    my ($type, $n_NM) = ($1, $2);
    next if ($n_NM > $max_NM); 
    $optional_tags =~ /NH:i:(\d+)/ or die "Incompatible BAM/SAM format: optional TAG NH is not present!";
    my $n_NH = $1;
    next if ($n_NH > $max_NH);

    my ($read, $flag, $chr, $start1, $CIGAR) = ($line[0], $line[1], $line[2], $line[3], $line[5]);
    # count soft clipping towards the edit distance
    my $CIGAR_old = $CIGAR; my $n_tmp;
    if ($CIGAR =~ /^(\d+)S(\d+)M/) {
	#print "$read\t$CIGAR_old\t$n_NM\t";
	$n_NM += $1; $n_tmp = $1 + $2; $CIGAR =~ s/^(\d+)S(\d+)M/${n_tmp}M/;
	#print "$CIGAR\t$n_NM\n";
    }
    if ($CIGAR =~ /(\d+)M(\d+)S$/) {
	#print "$read\t$CIGAR_old\t$n_NM\t";
	$n_NM += $2; $n_tmp = $1 + $2; $CIGAR =~ s/(\d+)M(\d+)S$/${n_tmp}M/;
	#print "$CIGAR\t$n_NM\n";
    }
    # also count indels with STAR input
    if ($type eq "nM") {
	while ($CIGAR =~ /(\d+)[ID]/g) { 
	    #print "$read\t$CIGAR\t$n_NM\t";
	    $n_NM += $1;
	    #print "$n_NM\n";
	}
    }
    next if ($n_NM > $max_NM);

    # indels right at the splice junction are not allowed
    next if ($CIGAR =~ /(\d+)[ID](\d+)N/ or $CIGAR =~ /(\d+)N(\d+)[ID]/);
    # skip reads with indels, for debugging
    #next if ($CIGAR =~ /(\d+)[ID]/);

    # pre-process indels to properly find read positions relative to junctions
    if ($CIGAR =~ /(\d+)[ID]/)
    {
    	while ($CIGAR =~ /(\d+)M(\d+)([ID])/g) {
	    #print "$read\t$CIGAR_old:\t$CIGAR -> ";
    	    if ($3 eq "I") { $n_tmp = $1; } #Alice used "$1 - $2"
    	    else { $n_tmp = $1 + $2; }
     	    $CIGAR =~ s/(\d+)M(\d+)([ID])/${n_tmp}M/;
	    #print "$CIGAR\n";
    	}
	# This was used to accomodate some bugs in STAR
    	# while ($CIGAR =~ /(\d+)N(\d+)([ID])(\d+)M/g) {
	#     #print "$read\t$CIGAR_old:\t$CIGAR -> ";
    	#     if ($3 eq "I") { $n_tmp = $4; } #Alice used "$4 - $2"
    	#     else { $n_tmp = $4 + $2; }
    	#     $CIGAR =~ s/(\d+N)(\d+)([ID])(\d+)M/$1${n_tmp}M/;
	#     #print "$CIGAR\n";
    	# }
	# merge consecutive matches resulted from previous steps
	while ($CIGAR =~ /(\d+)M(\d+)M/g) {
	    #print "$read\t$CIGAR_old:\t$CIGAR -> ";
	    $n_tmp = $1 + $2;
	    $CIGAR =~ s/(\d+)M(\d+)M/${n_tmp}M/;
	    #print "$CIGAR\n";
	}
    }

    # compute and ouput junction read positions
    my $strand = "+"; my $pos = "$start1";
    $strand = "-" if ($flag & 16);
    my $n_junc = 0;
    while ($CIGAR =~ /(\d+)M(\d+)N/g) {
    	$n_junc++;
    	my $end1 = $start1 + $1 - 1;
    	my $start2 = $start1 + $1 + $2;
    	$pos = "$pos,$end1,$start2";
    	$start1 = $start2;
    }
    $CIGAR =~ /(\d+)M$/;
    my $end1 = $start1 + $1 - 1;
    $pos = "$pos,$end1";
    # all coordinates in $pos are 1-based!
    #print modified $CIGAR in comment line as debug reference
    print "## $read:\t$CIGAR_old -> $CIGAR\n" unless ($CIGAR_old eq $CIGAR);
    print "$read\t$n_NM\t$n_NH\t$chr\t$strand\t$n_junc\t$pos\t$CIGAR\n";

    # for debugging multi-junction reads
    # if ($n_junc == 2) {
    # 	print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[5]\t$line[11]\t$line[13]\n";
    # 	print "$read\t$n_NM\t$n_NH\t$chr\t$strand\t$n_junc\t$pos\n";
    # 	last;
    # }
    
}
