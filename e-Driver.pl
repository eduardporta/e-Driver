use strict;
use warnings;

unless ($ARGV[2]) {
	die "USAGE: perl e-Driver.pl MODE[1:DOM|2:INTERFACES] FILE_WITH_MUTATIONS FILE_WITH_FEATURES FILE_WITH_PROTEIN_SEQUENCES OUTPUT_FILE_ROOT\n";
}

my %features;

if ($ARGV[0] == 2) {

	open IN, $ARGV[2] or die "CAN'T OPEN INPUT FILE WITH INTERFACES \"$ARGV[2]\"\n";

	while (<IN>) {

		chomp;
		my @array = split ("\t", $_);
	
		my $ensp = shift @array;
		my $interface = shift @array;
	
		unless (scalar @array > 5) {
			next;
		}
	
		foreach my $position (@array) {
			$features{$ensp}{$interface}{$position} = 1;
		}
	}

	close IN or die "CAN'T CLOSE FILEHANDLER TO \"$ARGV[2]\"\n";
}

elsif ($ARGV[0] == 1) {
	
	open IN, $ARGV[2] or die "CAN'T OPEN INPUT FILE WITH FEATURE \"$ARGV[2]\"\n";

	while (<IN>) {
	
		chomp;
		my @array = split ("\t", $_);
	
		my $ensp = shift @array;
		
		foreach my $feature (@array) {
			$features{$ensp}{$feature} = 1;
		}
	}

	close IN or die "CAN'T CLOSE FILEHANDLER TO \"$ARGV[2]\"\n";
}

else {
	die "The first argument needs to be either '1' for linear features or '2' for three-dimensional features\n";
}

unless (scalar (keys %features) > 0) {
	die "COULDN'T GET THE FEATURES\n";
}

my %seqs;

open IN, $ARGV[3] or die "CAN'T OPEN INPUT FILE WITH SEQUENCES \"$ARGV[3]\"\n";

while (<IN>) {
	
	chomp;
	my ($ensp, $seq) = split ("\t", $_);
	
	unless ($ensp && $seq) {
		die "CAN'T RECOGNIZE THIS PROTEIN SEQUENCE FORMAT \"$_\"\n";
	}
	
	$seqs{$ensp} = $seq;
}

close IN or die "CAN'T CLOSE FILEHANDLER TO \"$ARGV[3]\"\n";

my %muts;

open IN, $ARGV[1] or die "CAN'T OPEN INPUT FILE WITH MUTATION PROFILE \"$ARGV[1]\"\n";

while (<IN>) {
	
	chomp;
	my ($ensp, $prot_pos, $tum_id, $tissue) = split ("\t", $_);
	
	unless ($_ =~ /^ENSP/) {
		next;
	}
	
	unless (exists $seqs{$ensp} && exists $features{$ensp}) {
		next;
	}
	
	unless ($ensp && $prot_pos && $tum_id) {
		die "CAN'T RECOGNIZE THIS MUTATION PROFILE FORMAT \"$_\"\n";
	}
	
	$muts{$tissue}{$ensp}{"$tum_id\t$prot_pos"} = 1;
	$muts{'PANCAN'}{$ensp}{"$tum_id\t$prot_pos"} = 1;
}

close IN or die "CAN'T CLOSE FILEHANDLER TO \"$ARGV[1]\"\n";

unless (scalar (keys %muts) > 0) {
	die "COULDN'T GET THE MUTATIONS!\n";
}

my $name_perl_out = $ARGV[4]."_perl_output.txt";
open OUT, ">$name_perl_out" or die "CAN'T OPEN OUTPUT FILE \"$name_perl_out\"\n";
print OUT "Tissue\tProtein\tPFR\tMutationsRegion\tTotalMutations\tLengthRegion\tLengthProtein\n";

my %muts_in_structure;
	
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
	
		unless (exists $features{$ensp}) {
			next;
		}
		
		if ($ARGV[0] == 2) {
				
			my $ref_interfaces = $features{$ensp};
	
			while (my ($interface, $ref_int) = each %$ref_interfaces) {
			
				my ($start, $end, $pdb, $chain1, $chain2) = split ("-", $interface);
			
				if ($end > $prot_length) {
					next;
				}
									
				my $region_covered = $end - $start + 1;
				my $interface_length = scalar (keys %$ref_int);
		
				my %muts_in_region;
				my %muts_in_covered;
		
				while (my ($tumor_mut, $x) = each %$ref_prot) {
			
					my ($tumor, $pos) = split ("\t", $tumor_mut);
				
					if ($start <= $pos && $end >= $pos) {
						$muts_in_covered{$tumor_mut} = 1;
						$muts_in_structure{"$tumor_mut\t$ensp\t$pos"} = 1;
					}
			
					if (exists ${$ref_int}{$pos}) {
						$muts_in_region{$tumor_mut} = 1;
					}
				}
		
				my $mutations_in_region = scalar (keys %muts_in_region);
				my $mutations_in_covered = scalar (keys %muts_in_covered);
			
				unless ($mutations_in_covered > 0) {
					next;
				}

				print OUT "$tissue\t$ensp\t$interface\t$mutations_in_region\t$mutations_in_covered\t$interface_length\t$region_covered\n";
			}
        }
			
        else {
				
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

                print OUT "$tissue\t$ensp\t$feat\t$mutations_in_region\t$total_mutations_in_this_protein\t$feat_length\t$prot_length\n";
			}
		}
	}
}
	
close OUT or die;

system ("Rscript", "binomial_e-Driver.r", "$name_perl_out", "$ARGV[4]");
if ($?) {
	die "I COULD GENERATE THE PERL OUTPUT, BUT THERE WAS SOME PROBLEM RUNNING THE Rscript \"binomial_e-Driver.r\"\n";
}
