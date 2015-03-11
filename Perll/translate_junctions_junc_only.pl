#!/usr/bin/perl
# This translates across junctions found by TopHat.
# need a junction file, fix or no fix reqd, then chromosome file, fix or no fix then a set of numbers
# minimum number of codons either side of the junction, minimum size of orfs to report and minimum number of times a junction was found to consider
@files = @ARGV;
if ($files[1] eq "fix"){
	print "\nOK, fixing tophat junctions file file.\n";
	&fix_mac($files[0]);
}
if ($files[3] eq "fix"){
	print "\nOK, fixing chromosomes.\n";
	&fix_mac($files[2]);
}
open(INFILE, "$ARGV[0]"); # opens tophat junc file
open(INGENOME, "$ARGV[2]"); # opens the genome as well
print "Input - minimum length of returned peptides?\n";
$min_len = <STDIN>;
chomp $min_len;
print "Input - minimum number of hits to junction?\n";
$min_hits = <STDIN>;
chomp $min_hits;
if ($min_hits < 1){ $min_hits = 1;}
if ($min_len <1){$min_len = 1;}
open(OUT, ">$ARGV[0].min_length_$min_len.min_hit$min_hits.fasta");
print "\nOK, starting.\n";
# while ($test=<INGENOME>){$partof = substr($test, 0, 5);print "XXXXX $partof YYYYY\n";}
print "\nReading genome into a hash file.\nThis can take 6 to 7 minutes or so if its the human genme.\n";
%genome = <INGENOME>;
print "\nDone.\n";
%translation = (TTT=>"F", TTC=>"F", TCT=>"S", TCC=>"S", TAT=>"Y", TAC=>"Y", TGT=>"C", TGC=>"C", TTA=>"L", TCA=>"S", TAA=>"*", TGA=>"*", TTG=>"L", TCG=>"S", TAG=>"*", TGG=>"W", CTT=>"L", CTC=>"L", CCT=>"P", CCC=>"P", CAT=>"H", CAC=>"H", CGT=>"R", CGC=>"R", CTA=>"L", CTG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", ATT=>"I", ATC=>"I", ACT=>"T", ACC=>"T", AAT=>"N", AAC=>"N", AGT=>"S", AGC=>"S", ATA=>"I", ACA=>"T", AAA=>"K", AGA=>"R", ATG=>"M", ACG=>"T", AAG=>"K", AGG=>"R", GTT=>"V", GTC=>"V", GCT=>"A", GCC=>"A", GAT=>"D", GAC=>"D", GGT=>"G", GGC=>"G", GTA=>"V", GTG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");
%rev_translation = (GGC=>"A", ACT=>"S", TCA=>"*", ACA=>"C", TCG=>"R", GAT=>"I", GTT=>"N", GCT=>"S", GTA=>"Y", TGT=>"T", CGA=>"S", CGG=>"P", CAG=>"L", TGC=>"A", CAC=>"V", CTT=>"K", AAC=>"V", GTG=>"H", TCT=>"R", GGT=>"T", TGG=>"P", CCA=>"W", GAG=>"L", GCG=>"R", CAA=>"L", TTA=>"*", CTG=>"Q", CGT=>"T", CAT=>"M", TTT=>"K", TAC=>"V", CTA=>"*", AAG=>"L", TCC=>"G", GAC=>"V", GCA=>"C", TGA=>"S", AAT=>"I", ATA=>"Y", ATT=>"N", AGT=>"T", TTG=>"Q", GTC=>"D", ACC=>"G", GGA=>"S", AAA=>"F", CCT=>"R", ACG=>"R", CCG=>"R", ATG=>"H", TAT=>"I", GGG=>"P", CCC=>"G", TAA=>"L", CTC=>"E", TAG=>"L", ATC=>"D", AGA=>"S", GAA=>"F", CGC=>"A", GCC=>"G", AGC=>"A", TTC=>"E", AGG=>"P");
%base_pair = (G=>"A", A=>"T", T=>"A", C=>"G");
while($line = <INFILE>)
{
	chomp $line;
	# split the line into an array called cells
	@cells = split(/\t/, $line);
	# columns are  0 = chr, 1 = start plus bed overlap, 2 = end minus bed, 3=junction number, 4=number of hits, 5=strand, 9=bed block length as csv
	# now we know we have an exon, lets set some variables to start off
	if ($cells[4] < $min_hits){next;}
	$chr = ">$cells[0]\n";
	print "chr $chr\n";
	$start = $cells[1];
	$stop = $cells[2];
	$csv = $cells[10];
	if(length ($genome{$chr}) < 10){next;}
	($bed_overlap1, $bed_overlap2) = split (/\D/, $csv);
	# This bit is messy but it means the intron is removed - complicated by chr numbering and string numbering and bed file conventions...
	$start_real = $start + $bed_overlap1 - 1;
	$stop_real = $stop - $bed_overlap2 - 1;
	$gap = $stop_real - $start_real;
	$temp_min_codons = 1800; #Used to be able to input a minimum number of codons to try either side but now I just set it
	$region_needed =  substr($genome{$chr}, ($start_real - $temp_min_codons), ($gap + ($temp_min_codons * 2)) );
	$temp_min_codons++;
	if (length($region_needed) < 3598){
		print "\nregion needed too short\n";
		next;
	}
	print"region needed is $region_needed\n";
	$neg_min_codons = 2 - $temp_min_codons;
	$region_translate = substr($region_needed, 0, $temp_min_codons).substr($region_needed, $neg_min_codons, $temp_min_codons);
	$junc_itself = substr($region_translate,1790,18);
	print "region translate is $region_translate\nJunction is $junc_itself\n";
	if (length($region_translate) < 3000){
		print "\nregion to translate is too short\n";
		next;
	}
	$strand = $cells[5];
	$steps_from_junc = 0;
	$start = 1800;
	%translate_now = %rev_translation;
	$step = - 3;
	if ($strand eq "+"){%translate_now = %translation; $step = 3;}
	# now to translate the forward three frames
	$left_rk = 0;
	$right_rk = 0;
	$junc_frame1 = &translate(1800);
	if (length($junc_frame1)>= $min_len){print OUT ">Junction $junc_itself at chr $cells[0] from $start_real to $stop_real with number of hits $cells[4] $strand strand frame 1\n";
		print OUT "$junc_frame1\n";}
	$junc_frame2 = &translate(1801);
	if (length($junc_frame2)>= $min_len){print OUT ">Junction $junc_itself at chr $cells[0] position $start_real to $stop_real number of hits $cells[4] $strand strand frame 2\n";
		print OUT "$junc_frame2\n";}
	$junc_frame3 = &translate(1802);
	if (length($junc_frame3)>= $min_len){print OUT ">Junction $junc_itself at chr $cells[0] position $start_real to $stop_real number of hits $cells[4] $strand strand frame 3\n";
		print OUT "$junc_frame3\n";}

	}



sub fix_mac
{
	$filefix = $_[0];
	chomp $filefix;
	print "Correcting a file.\n";
	open (OUT, ">temp");
	open (INFILE, "$filefix");
	while (<INFILE>){s/[\r\n]+/\n/g; print OUT "$_";}
	close OUT;
	close INFILE;
	system("mv temp $filefix"); 
}



sub translate
{
	$start = $_[0];
	$left_rk = 0;
	$right_rk = 0;
	$start_aa = $translate_now{substr($region_translate, ($start), 3)};
	$new_peptide = " $start_aa ";
	$steps_from_junc = $step;
	print "just before steps from junc is $steps_from_junc\n";
	do {
		$next_codon_for = substr($region_translate, ($start + $steps_from_junc), 3);
		$next_codon_rev = substr($region_translate, ($start - $steps_from_junc), 3);
		if ($left_rk < 1) {
			$new_peptide = $translate_now{$next_codon_rev}.$new_peptide;
		}
		if ($right_rk < 1) {
			$new_peptide = $new_peptide.$translate_now{$next_codon_for};
		}

		$new_peptide =~ s/RP/rp/g;
		$new_peptide =~ s/KP/kp/g;
		$new_peptide =~ s/R P/r p/g;
		$new_peptide =~ s/K P/k p/g;
		$new_peptide =~ s/\*/XX/g;
		($left,$middle,$right) = split (/ /,$new_peptide);
		print "$new_peptide , Step is $step, Steps from is $steps_from_junc , codons $next_codon_for , $next_codon_rev .\n";
		$steps_from_junc = $steps_from_junc + $step;
		$left_rk = 0;
		$right_rk = 0;
		while ($left =~ /[RKX]/g){$left_rk++;}
		while ($right =~ /[RKX]/g){$right_rk++;}
		print "$left_rk X $right_rk Y\n $left XXXXX $start_aa ZZZZZZ $right YYYYY\n";
	} until ((($left_rk >=1) && ($right_rk >=1)) || (length($new_peptide) > 1000));
	print "left: $left , right: $right .\n";
	$new_peptide = $left.$start_aa.$right;
	$new_peptide =~ s/XX//g;
	$new_peptide = uc $new_peptide;
	print "new peptide is $new_peptide\n";
	$steps_from_junc = 0;
	print "just afterwards steps from junc is $steps_from_junc\n";
	return $new_peptide;
}
	
