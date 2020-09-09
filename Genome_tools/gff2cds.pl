#!/usr/bin/perl
use strict;
use Getopt::Std;
use vars qw($opt_g $opt_f $opt_o);
my $usage = "usage: perl $0 -g gff file -f genome fasta file -o output name\n\
		notice : your gff file must to have CDS region !\n\
		your system need to have bedtools software (yum install bedtools)\n";
getopts('g:f:o:') or die $usage;
my $name = $opt_o;
open(IN,"$opt_g") or die $usage;
open(OUT,">$opt_o\.bed") or die "open $opt_o file failed\n";
my %hash;
while (<IN>) {
	chomp;
	next if(/^#/);
	my @l=split /\t/;
	$l[3]=$l[3]-1;
	if($l[2] eq "CDS")
	{
		my ($tran_id)=$_=~/Parent=([\w\.\-]+)/;
		$hash{$tran_id} = $l[6];
		print OUT "$l[0]\t$l[3]\t$l[4]\t$tran_id\t.\t$l[6]\n";
	}
}
close IN;
`bedtools getfasta -fi $opt_f -fo $name\.tab -bed $name\.bed -name+ -tab`;
my %hash_seq;
my %hash_cdna;
open(IN,"$name\.tab") or die "open $name\.tab failed\n";
while (<IN>) {
	chomp;
	my @l=split /\t/;
	my @m=split /\:|\-/,$l[0];
	$hash_seq{$m[0]}{$m[3]} =$l[1];
}
close IN;

foreach my $k(keys %hash_seq)
{
	foreach my $kk(sort{$a <=> $b} keys %{$hash_seq{$k}})
	{
		$hash_cdna{$k} .= "$hash_seq{$k}{$kk}"
	}
}
open(OUT,">$name\.cds.fa") or die $!;
foreach my $k(keys %hash_cdna)
{
	if($hash{$k} eq "+")
	{
		print OUT ">$k\n$hash_cdna{$k}\n";
	}
	else
	{
		my $str = &reverse_complementary($hash_cdna{$k});
		print OUT ">$k\n$str\n";
	}
}
sub reverse_complementary{
	my ($base) =@_;
	$base = reverse($base);
	$base =~tr/ATCGatcg/TAGCtagc/;
	return $base;
}

