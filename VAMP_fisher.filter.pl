# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
	'a=f' => \(my $alpha = 0.05),
	'D' => \(my $noDriveCount = ''),
	'O' => \(my $noOddsratio = ''),
	'm' => \(my $includeMutants = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl VAMP_fisher.filter.pl VAMP_fisher.txt > VAMP_fisher.filter.txt

Options: -h       display this help message
         -a FLOAT alpha, p-value cutoff
         -D       do not consider drive count
         -O       do not consider odds ratio
         -m       include mutants of selected clusters

EOF
}
my ($fisherFile) = @ARGV;
open(my $reader, $fisherFile);
my @columnList = ();
my @tokenHashList = ();
while(my $line = <$reader>) {
	chomp($line);
	if($line =~ s/^#//) {
		@columnList = split(/\t/, $line);
		print '#', join("\t", @columnList), "\n";
	} else {
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line);
		if($noDriveCount || $tokenHash{'driveCount'}) {
			if($tokenHash{'pvalue'} <= $alpha) {
				$tokenHash{'selected'} = 1;
			} elsif($noOddsratio eq '' && ($tokenHash{'oddsratio'} eq "Inf" || $tokenHash{'oddsratio'} == 0)) {
				$tokenHash{'selected'} = 1;
			}
		}
		push(@tokenHashList, \%tokenHash);
	}
}
close($reader);

if($includeMutants) {
	my %clusterPhenotypeHash = ();
	foreach(grep {$_->{'selected'}} @tokenHashList) {
		foreach my $genotype (split(/,/, $_->{'genotypes'})) {
			(my $cluster = $genotype) =~ s/\|.*$//;
			if($cluster eq $genotype) {
				$clusterPhenotypeHash{$cluster} = $_->{'phenotype'};
			}
		}
	}
	foreach(grep {!$_->{'selected'}} @tokenHashList) {
		my @genotypeList = ();
		foreach my $genotype (split(/,/, $_->{'genotypes'})) {
			(my $cluster = $genotype) =~ s/\|.*$//;
			if($cluster ne $genotype) {
				if(defined(my $clusterPhenotype = $clusterPhenotypeHash{$cluster})) {
					push(@genotypeList, $genotype) if($_->{'phenotype'} ne $clusterPhenotype && $_->{'phenotype'} ne '');
				}
			}
		}
		if(@genotypeList) {
			$_->{'genotypes'} = join(',', @genotypeList);
			$_->{'selected'} = 1;
		}
	}
}

print join("\t", @$_{@columnList}), "\n" foreach(grep {$_->{'selected'}} @tokenHashList);
