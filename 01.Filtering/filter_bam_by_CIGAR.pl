#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($bam, $ide,$mlen,$mapq,$out,$rev,$sim);
die "
	$0
	
	-b  bam file.
	-i  identity [mapped/tot_len].
	-s  similarity [matched/mapped].
	-m  minimum insert length [mapped_len].
	-q  minimum mapping quality.
	-o  out bam file name.
	-r  reverse operation.
	
" unless @ARGV >= 5;

GetOptions (
    "bam=s"        => \$bam,
    "identity=f"   => \$ide,
    "similarity=f" => \$sim,
	"minlen=i"     => \$mlen,
	"quality=i"    => \$mapq,
	"outfile=s"    => \$out,
	"reverse:1"    => \$rev
) or die ("Error in command line arguments\n");

open(IN, "samtools view --threads 8 -h $bam |") or die "$!";
open(OUT, "|samtools view --threads 8 -b >$out") or die "$!";

while(<IN>) {
	my (@hclips, @sclips, @nMs, @nIs, @nDs);
	my ($hclip, $sclip, $nM, $nI, $nD) = (0,0,0,0,0);
	if (/^@/) {print OUT $_; next}
	my @p = split(/\t/);
	my ($q,$cigar,$i, $seq, $nm) = ($p[4], $p[5], $p[8], $p[9], $p[11]);
    
    if (!$rev) {
        next if $q < $mapq;
        next if abs($i) < $mlen;
    } else {
        
        my $len = length($seq);
        # number of mismatches.
        if ($nm =~ /NM:i:(\d+)/) {$nm = $1} else {$nm=0}
        # number of insertions.
        if (@nIs = $cigar =~ /(\d+)I/g) {$nI = sum(\@nIs)} else {$nI=0}
        # number of deletions.
        if (@nDs = $cigar =~ /(\d+)D/g) {$nD = sum(\@nDs)} else {$nD=0}
        # number of matches.
        if (@nMs = $cigar =~ /(\d+)M/g) {$nM = sum(\@nMs)} else {$nM=0}
        $nm -= ($nI + $nD);
		if ($q < $mapq || abs($i) < $mlen || (1-$nm/$nM) < $sim) {
			print OUT $_;
			next;
		}
		if ($cigar =~ /[SH]/) {
			if (@hclips = $cigar =~ /(\d+)H/g) {
				$hclip = sum(\@hclips); 
				$len += $hclip;
			}
	
			if (@sclips = $cigar =~ /(\d+)S/g) {
				$sclip = sum(\@sclips);
			}
			
			print OUT $_ if 1-($hclip+$sclip)/$len < $ide;
			
		}
		
		next;
	}
    
    my $len = length($seq);
    # number of mismatches.
    if ($nm =~ /NM:i:(\d+)/) {$nm = $1} else {$nm=0}
    # number of insertions.
    if (@nIs = $cigar =~ /(\d+)I/g) {$nI = sum(\@nIs)} else {$nI=0}
    # number of deletions.
    if (@nDs = $cigar =~ /(\d+)D/g) {$nD = sum(\@nDs)} else {$nD=0}
    # number of matches.
    if (@nMs = $cigar =~ /(\d+)M/g) {$nM = sum(\@nMs)} else {$nM=0}
    $nm -= ($nI + $nD);
    next if ($nM != 0 && $nm != 0 && (1-$nm/$nM) < $sim);
	if ($cigar !~ /[SH]/) {print OUT $_; next}
	
	if (@hclips = $cigar =~ /(\d+)H/g) {
		$hclip = sum(\@hclips); 
		$len += $hclip;
	}
	
	if (@sclips = $cigar =~ /(\d+)S/g) {
		$sclip = sum(\@sclips);
	}
	next if 1-($hclip+$sclip)/$len < $ide;
	print OUT $_; 
}

sub sum {
	my $tmp = shift;
	my $s   = 0;
	map {$s+=$_} @$tmp;
	return $s;
}