use strict;
use warnings;
use Text::NSP::Measures::2D::Fisher2::twotailed;

unless ($ARGV[2]) {
	die "USAGE: perl e-Driver.pl FILE_WITH_MUTATIONS FILE_WITH_ANNOTATIONS FILE_WITH_PROTEIN_SEQUENCES OUTPUT_FILE FLAG_TISSUE_MODE\n";
}

my %features;

open IN, $ARGV[1] or die "CAN'T OPEN INPUT FILE WITH FEATURE \"$ARGV[1]\"\n";

while (<IN>) {
	
	chomp;
	my @array = split ("\t", $_);
	
	my $ensp = shift @array;
	
	foreach my $feature (@array) {
		$features{$ensp}{$feature} = 1;
	}
}

close IN or die "CAN'T CLOSE FILEHANDLER TO \"$ARGV[1]\"\n";

my %seqs;

open IN, $ARGV[2] or die "CAN'T OPEN INPUT FILE WITH SEQUENCES \"$ARGV[2]\"\n";

while (<IN>) {
	
	chomp;
	my ($ensp, $seq) = split ("\t", $_);
	
	unless ($ensp && $seq) {
		die "CAN'T RECOGNIZE THIS PROTEIN SEQUENCE FORMAT \"$_\"\n";
	}
	
	$seqs{$ensp} = $seq;
}

close IN or die "CAN'T CLOSE FILEHANDLER TO \"$ARGV[2]\"\n";

my %muts;

open IN, $ARGV[0] or die "CAN'T OPEN INPUT FILE WITH MUTATION PROFILE \"$ARGV[0]\"\n";

while (<IN>) {
	
	chomp;
	my ($ensp, $prot_pos, $tum_id, $tissue) = split ("\t", $_);
	
	unless (exists $seqs{$ensp} && exists $features{$ensp}) {
		next;
	}
	
	unless ($ensp && $prot_pos && $tum_id) {
		die "CAN'T RECOGNIZE THIS MUTATION PROFILE FORMAT \"$_\"\n";
	}
	
	if ($ARGV[4]) {
		$muts{$tissue}{$ensp}{"$tum_id\t$prot_pos"} = 1;
		$muts{'PANCAN'}{$ensp}{"$tum_id\t$prot_pos"} = 1;
	}
	
	else {
		$muts{$ensp}{"$tum_id\t$prot_pos"} = 1;
	}
}

close IN or die "CAN'T CLOSE FILEHANDLER TO \"$ARGV[0]\"\n";

open OUT, ">$ARGV[3]" or die "CAN'T OPEN OUTPUT FILE \"$ARGV[3]\"\n";

if ($ARGV[4]) {

	print OUT "Tissue\tProtein\tPFR\tpval\tExpected mutations\tObserved mutations\tTotal mutations in protein\tLength feature\tLength protein\n";
	
	while (my ($tissue, $ref_tissue) = each %muts) {
		while (my ($ensp, $ref_prot) = each %$ref_tissue) {
				
			my $total_samples_with_mutations_in_this_protein = scalar (keys %$ref_prot);
				
			unless ($total_samples_with_mutations_in_this_protein > 5) {
				next;
			}
				
			my $total_mutations_in_this_protein = $total_samples_with_mutations_in_this_protein;
				
			my $prot_seq = $seqs{$ensp};
			my $prot_length = length ($prot_seq);
					
			unless ($prot_length && $prot_length > 0) { 
				next;
			}
				
			my $total_length_protein = $prot_length * $total_samples_with_mutations_in_this_protein;

			unless (exists $features{$ensp}) {
				next;
			}
				
			my $ref_feat = $features{$ensp};
	
			while (my ($feat, $x) = each %$ref_feat) {
		
				my ($feature, $start, $end) = split ("-", $feat);
					
				if ($end > $prot_length) {
					print "THERE IS SOMETHING WRONG WITH PROTEIN \"$ensp\"\nACCORDING TO THE SEQUENCE FILE IT HAS \"$prot_length\" AA, BUT IN THE FEATURE FILE IT HAS THIS FEATURE \"$feat\" WHICH ENDS AFTER THIS LENGTH!\n";
					next;
				}
					
				my $feat_length = $end-$start;
		
				my %muts_in_region;
		
				while (my ($tumor_mut, $x) = each %$ref_prot) {
			
					my ($tumor, $pos) = split ("\t", $tumor_mut);
			
					if ($pos >= $start && $pos <= $end) {
						$muts_in_region{$tumor_mut} = 1;
					}
				}
		
				my $mutations_in_region = scalar (keys %muts_in_region);
				my $total_length_region = $feat_length * $total_samples_with_mutations_in_this_protein;
		
				my $exp = int(($feat_length/$prot_length)*$total_mutations_in_this_protein);
		
				my $pval = two_tailed ($mutations_in_region, $total_mutations_in_this_protein, $total_length_region, $total_length_protein, $ensp, $feat, $tissue);

				print OUT "$tissue\t$ensp\t$feat\t$pval\t$exp\t$mutations_in_region\t$total_mutations_in_this_protein\t$feat_length\t$prot_length\n";
			}
		}
	}
	
	close OUT or die;
}

else {
	
	print OUT "Protein\tPFR\tpval\tExpected mutations\tObserved mutations\tTotal mutations in protein\tLength feature\tLength protein\n";
	
	while (my ($ensp, $ref_prot) = each %muts) {
				
		my $total_samples_with_mutations_in_this_protein = scalar (keys %$ref_prot);
				
		unless ($total_samples_with_mutations_in_this_protein > 5) {
			next;
		}
				
		my $total_mutations_in_this_protein = $total_samples_with_mutations_in_this_protein;
				
		my $prot_seq = $seqs{$ensp};
		my $prot_length = length ($prot_seq);
					
		unless ($prot_length && $prot_length > 0) { 
			next;
		}
				
		my $total_length_protein = $prot_length * $total_samples_with_mutations_in_this_protein;

		unless (exists $features{$ensp}) {
			next;
		}
				
		my $ref_feat = $features{$ensp};
	
		while (my ($feat, $x) = each %$ref_feat) {
		
			my ($feature, $start, $end) = split ("-", $feat);
					
			if ($end > $prot_length) {
				die "THERE IS SOMETHING WRONG WITH PROTEIN \"$ensp\"\nACCORDING TO THE SEQUENCE FILE IT HAS \"$prot_length\" AA, BUT IN THE FEATURE FILE IT HAS THIS FEATURE \"$feat\" WHICH ENDS AFTER THIS LENGTH!\n";
			}
					
			my $feat_length = $end-$start;
		
			my %muts_in_region;
		
			while (my ($tumor_mut, $x) = each %$ref_prot) {
			
				my ($tumor, $pos) = split ("\t", $tumor_mut);
			
				if ($pos >= $start && $pos <= $end) {
					$muts_in_region{$tumor_mut} = 1;
				}
			}
		
			my $mutations_in_region = scalar (keys %muts_in_region);
			my $total_length_region = $feat_length * $total_samples_with_mutations_in_this_protein;
		
			my $exp = int(($feat_length/$prot_length)*$total_mutations_in_this_protein);
		
			my $pval = two_tailed ($mutations_in_region, $total_mutations_in_this_protein, $total_length_region, $total_length_protein, $ensp, $feat);

			print OUT "$ensp\t$feat\t$pval\t$exp\t$mutations_in_region\t$total_mutations_in_this_protein\t$feat_length\t$prot_length\n";
		}
	}
	
	close OUT or die;
}

exit;

################################

sub two_tailed {

	my ($common, $mut, $dom, $total, $ensp, $feat) = @_;
	my $p = Text::NSP::Measures::2D::Fisher2::twotailed::calculateStatistic (n11 => $common, n1p => $mut, np1 => $dom, npp => $total);
	
	my $errorCode;
	if( ($errorCode = getErrorCode())) {
		my $params_error = join (" ", @_);
		die "COULD NOT CALC PVAL WITH THIS PARAMS: \"$params_error\"\n$errorCode\n";
	}

	return $p;
}