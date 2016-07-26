#!/usr/bin/perl

use warnings; use strict; use Getopt::Std; use foot;
use vars qw($opt_n $opt_t $opt_l $opt_r $opt_q $opt_s $opt_g $opt_c $opt_S $opt_i $opt_e $opt_f);
getopts("n:t:l:r:q:s:g:ci:S:e:f");

my $usage = "\nUsage: $YW$0$N ${GN}[options]$N -n $CY<Output name>$N -t $CY<conversion threshold>$N -l $CY<minimum R-loop length>$N -r $CY<path to reads>$N -q $CY<minimum mapping quality for each read to be considered>$N -s $CY<path to original gene sequence>$N -g $CY<path to reference genome for indexing>$N

${GN}Extra information (required):$N
1) -n: Output name
2) Either:
${CY}-r$N: must be path to reads file itself named pacbio.fastq
${CY}-S$N: supply a sam file
3) -i: index bed file (bed4 with gene name on each)
4) -e: Exon bed file
5) -g: must be path to .fa file containing reference genome used for indexing (e.g. hg19.fa)

${GN}Optional [default]:$N
-c: <default: don't include CpG> consider Cs in CpG context
-t: [0.65] percentage (0.0-1.0) of Cs that must be converted to be considered an R-loop
-l: [100] minimum length in base pairs to be considered an R-loop
-q: [0] must be phred score (e.g. score of 40 = 99.99% base call accuracy)
-s: <fasta> use this sequence instead of the one given in index

${GN}Example:$N

If you have .fq file but no SAM file:
$CY$0$N -n ${CY}myoutput$N -t ${CY}0.65$N -l$CY 100$N -q$CY 0$N -e$CY /home/mitochi/shared/footLoop_example/hg19_appris_name.exon$N -i$CY /home/mitochi/shared/footLoop_example/pacbio12index_original.bed$N -g$CY /home/mitochi/shared/footLoop_example/hg19.fa.fa$N ${YW}-r /home/mitochi/shared/footLoop_example/pacbio12ccs.fq$N [optionl: -s$BU CALM3.fa$N]

If you have SAM file, use -S instead of -r (everything is the same as above except$YW yellow$N):
$CY$0$N -n ${CY}myoutput$N -t ${CY}0.65$N -l$CY 100$N -q$CY 0$N -e$CY /home/mitochi/shared/footLoop_example/hg19_appris_name.exon$N -i$CY /home/mitochi/shared/footLoop_example/pacbio12index_original.bed$N -g$CY /home/mitochi/shared/footLoop_example/hg19.fa.fa$N ${YW}-S /home/mitochi/shared/footLoop_example/pacbio12.sam$N [optional: -s$BU CALM3.fa$N]

";

die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -n <gene_name [CALM3]> not defined\n\n${LRD}   ########## USAGE ##########\n$N$usage" if not defined($opt_n);
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -r <read.fq> and $opt_s <read.sam> both not defined\n\n${LRD}   ########## USAGE ##########\n$N$usage" if not defined($opt_r) and not defined($opt_S);
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -i <geneindex.bed> not defined\n\n${LRD}   ########## USAGE ##########\n$N$usage" if not defined($opt_i);
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -g <ref_genome.fa [hg19.fa]> not defined\n\n${LRD}   ########## USAGE ##########\n$N$usage" if not defined($opt_g);
my $read0 = defined($opt_r) ? $opt_r : "FALSE";
my $sam0 = defined($opt_S) ? $opt_S : "FALSE";
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -r $read0 and -S $sam0 both DOES NOT EXIST (provide at least one!)\n\n${LRD}########## FATAL ERROR ##########\n$N$usage" if (defined($opt_r) and not -e ($opt_r)) or (defined($opt_S) and not -e ($opt_S));
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -i $opt_i DOES NOT EXIST\n\n${LRD}########## FATAL ERROR ##########\n$N$usage" if not -e ($opt_i);
die "\n${LRD}########## FATAL ERROR ##########\n\n!!!$N -g $opt_g DOES NOT EXIST\n\n${LRD}########## FATAL ERROR ##########\n$N$usage" if not -e ($opt_g);
die "\n$LRD########## FATAL ERROR ##########\n\n!!!$N -S $CY$opt_S$N exists but seems to be empty!\n$N$usage" if defined($opt_S) and (not -e $opt_S or (-s $opt_S < 10));

#### Input ####
$opt_i = getFullpath($opt_i);
$opt_g = getFullpath($opt_g);
$opt_r = getFullpath($opt_r) if defined($opt_r); my ($readname) = getFullpath("$opt_n") . getFilename($opt_r,'full') . "\_bismark_bt2.sam" if defined($opt_r);
my $faindex = defined($opt_s) ? getFullpath($opt_s) : getFullpath("$opt_n/$opt_n\_geneIndexes.fa");
my $mysam   = defined($opt_S) ? getFullpath($opt_S) : -e $readname ? $readname : getFullpath($opt_r) . "_bismark_bt2.sam";
$opt_t = defined($opt_t) ? $opt_t : 0.65;
$opt_l = defined($opt_l) ? $opt_l : 100;
$opt_q = defined($opt_q) ? $opt_q : 0;
my $myread = defined($opt_r) ? $opt_r : "${LRD}FALSE$N (Sam File Given)";
my $usecpg = $opt_c ? "${LGN}TRUE$N" : "${LRD}FALSE$N (Don't Use C from CpG)";
my $exonFile  = defined($opt_e) ? getFullpath($opt_e) : "${LRD}FALSE$N (exon file not given)";

# Make directory
mkdir "$opt_n" if not -d "$opt_n";
my $mydir = getFullpath("$opt_n") . "/";
mkdir $mydir if not -d $mydir;
chdir $mydir;
my $logFile = "$mydir/logFile.txt";
open(my $outLog, '>', $logFile);
my $mysam2 = "$mysam"; $mysam2 .= " $GN(will be generated)$N" if not defined($opt_S);

print $outLog "\n${YW}0. Initializing... output directory = $CY$mydir$N\n";
print $outLog "
${YW}Input Parameters$N
1. -n ${CY}Out Dir$N   :$mydir
2. -r ${CY}Read$N      :$myread
3. -S ${CY}SAM$N       :$mysam2
4. -i ${CY}Index$N     :$opt_i
5. -g ${CY}Genome$N    :$opt_g
6. -s ${CY}Seq$N       :$faindex
6. -c ${CY}UseCpG?$N   :$usecpg
7. -t ${CY}Threshold$N :$opt_t
8. -l ${CY}MinLength$N :$opt_l
9. -q ${CY}MinMapQ$N   :$opt_q
10.-e ${CY}Exon$N      :$exonFile

";
print STDERR "
${YW}Input Parameters$N
1. -n ${CY}Out Dir$N   :$mydir
2. -r ${CY}Read$N      :$myread
3. -S ${CY}SAM$N       :$mysam2
4. -i ${CY}Index$N     :$opt_i
5. -g ${CY}Genome$N    :$opt_g
6. -s ${CY}Seq$N       :$faindex
6. -c ${CY}UseCpG?$N   :$usecpg
7. -t ${CY}Threshold$N :$opt_t
8. -l ${CY}MinLength$N :$opt_l
9. -q ${CY}MinMapQ$N   :$opt_q
10.-e ${CY}Exon$N      :$exonFile

";
#runs bismark (output file will be pacbio.fastq_bismark_bt2.sam) only if it hasn't been ran previously

my $outFolder     = $mydir;
my $geneIndexes   = $opt_i; die "\n$LRD!!!\t$N$opt_i doesn't exist!\n\n" if not defined($geneIndexes);
my $geneIndexesFa = "$mydir\/geneIndexes.fa";
print STDERR "\n${YW}1. Getting fasta sequence from $geneIndexes into $CY$geneIndexesFa$N\n";
print $outLog "\n${YW}1. Getting fasta sequence from $geneIndexes into $CY$geneIndexesFa$N\n";
print $outLog "\t- Running ${YW}bedtools getfasta$N -fi $opt_g -bed $geneIndexes -fo $geneIndexesFa -name\n";
system("fastaFromBed -fi $opt_g -bed $geneIndexes -fo $geneIndexesFa -name") == 0 ? print $outLog "\t${GN}SUCCESS$N: Output: $CY$geneIndexesFa$N\n" : die "Failed to run bedtools: $!\n";
system("perl -pi -e 'tr/acgt/ACGT/' $geneIndexesFa") == 0 or die "\n$LRD!!!$N\tFailed to convert atgc into ATGC\n";

print $outLog "\n${YW}1b. Running$CY bismark$YW (output file will be$CY pacbio.fastq_bismark_bt2.sam$YW) only if it hasn't been ran previously$N\n";
	
	print $outLog "\t- Running$YW bismark_genome_preparation$N --bowtie2 $outFolder\n";
	if (not -d "$outFolder/Bisulfite_Genome") {
		system("bismark_genome_preparation --bowtie2 ./") == 0 or die "Failed to run bismark genome preparation: $!\n";
	} else {print $outLog "\t${GN}SUCCESS$N: Bisulfite_Genome already exist!\n"}
	
if (not -e "$mysam" or -s $mysam <= 10) {
	print $outLog "\t- Running$YW bismark$N --bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8 ./ $opt_r\n";
	my $result = system("bismark --bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8 ./ $opt_r");
	if ($result != 0) {
		print $outLog "\t${LRD}Bisulfte_Genome seems to be corrupted so re-running:\n\t${YW}-bismark_genome_preparation$N --bowtie2 ./\n";
		system("bismark_genome_preparation --bowtie2 ./") == 0 or die "Failed to run bismark genome preparation: $!\n";
		system("bismark --bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8 ./ $opt_r") == 0 or die "$LRD!!!$N\tFailed to run bismark: $!\n";
	}
	print $outLog "\t${GN}SUCCESS$N: Output $mysam\n";
	print STDERR "\t${GN}SUCCESS$N: Output $mysam\n";
}
else {
	print $outLog "\t${GN}SUCCESS$N: Output already exist: $CY$mysam$N\n";
	print STDERR "\t${GN}SUCCESS$N: Output already exist: $CY$mysam$N\n";
}

#takes sequence of gene and splits into array of individual bases
my $seqFile = $geneIndexesFa;
my %seq;
print STDERR "\n${YW}2. Parsing in sequence for genes from sequence file $CY$seqFile$N\n";
print $outLog "\n${YW}2. Parsing in sequence for genes from sequence file $CY$seqFile$N\n";
open(SEQ, $seqFile) or die "\n$LRD!!!$N\tFATAL ERROR: Could not open $CY$seqFile$N: $!";
while (my $line = <SEQ>) {
	chomp($line);
	if ($line =~ />/) {
		my ($gene) = $line =~ /^>(.+)$/;
		$line = <SEQ>;
		die "\n$LRD!!!$N\tFATAL ERROR: Corrupted fasta file $seqFile! $gene doesn't have sequence!\n\n" if $line =~ /^>/;
		@{$seq{$gene}{seq}} = split("", $line);
		$seq{$gene}{len} = @{$seq{$gene}{seq}};
		$seq{$gene}{total} = 0;
		$seq{$gene}{badlength} = 0;
		$seq{$gene}{lowq} = 0;
		$seq{$gene}{used} = 0;
		$seq{$gene}{pos} = 0;
		$seq{$gene}{neg} = 0;
	}
}
print STDERR "\t${GN}SUCCESS$N: Sequence has been parsed from fasta file $CY$seqFile$N\n";
print $outLog "\t${GN}SUCCESS$N: Sequence has been parsed from fasta file $CY$seqFile$N\n";

if ($opt_e) {
	#0 = white
	#1 = intron line
	#2 = exon line
	print $outLog "\t${BR}2a. Parsing in exon intron data from $CY$exonFile$N:\n";
	foreach my $gene (sort keys %seq) {
	next if -e "$outFolder/exon/$gene.exon"; # DELETE THIS
		my ($genepos) = `grep $gene $geneIndexes`; chomp($genepos);
		my $length_seq = $seq{$gene}{len};
		print $outLog "\n$LRD!!!$N\tWARNING: No sequence for gene $CY$gene$N is found in -s $CY$seqFile$N!\n\n" if $length_seq == 0;
		parse_exon($exonFile, $genepos, $gene, $outFolder, $length_seq);
		print STDERR "\t${GN}SUCCESS$N: Sequence of exon $CY$gene$N has been parsed from fasta file $CY$seqFile$N (length = $CY" . $length_seq . "$N bp)\n";
		print $outLog "\t${GN}SUCCESS$N: Sequence of exon $CY$gene$N has been parsed from fasta file $CY$seqFile$N (length = $CY" . $length_seq . "$N bp)\n";
	}
}

close SEQ;


# File description, open files, etc
my $samFile = $mysam;
my $positive = "$mydir/Positive.3";
my $negative = "$mydir/Negative.3";
my $notusedfile = "$mydir/NOTUSED.3";
print STDERR "${YW}3. Parsing sam file $CY$samFile$YW and getting only high quality reads\n";
print $outLog "${YW}3. Parsing sam file $CY$samFile$YW and getting only high quality reads\n";
open(my $positiveReads, ">", $positive) or die "Could not open $positive: $!";
open(my $negativeReads, ">", $negative) or die "Could not open $negative: $!";
open(my $notused, ">", $notusedfile) or die "Cannot open $notusedfile: $!\n";
open(my $sam, $samFile) or die "$LRD!!!$N\tFATAL ERROR: Could not open $samFile: $!";



## Some stats
my $linecount = 0; my %count; ($count{total}, $count{used}, $count{diffgene}, $count{lowq}, $count{badlength}) = (0,0,0,0,0);
my %readz;

#loops through each read and writes high quality reads into two separate files <gene>Positive.txt and <gene>Negative.txt
while(my $line = <$sam>) {
	$linecount ++;
	chomp($line);

	my @fields = split("\t", $line); #line is tab separated, first column is read name, the rest is value
	my $gene = $fields[2];

	#discounts the first 4 lines of information at the top of the .sam file
	next if(@fields < 6);

	# (EXTRA) This below is to shorten read name. Unique read name is the number between last number or behind ccs. If can't be found just use the whole read name.
	my ($read) = $fields[0] =~ /^.+\/(\d+)\/\d+/i;
	($read) = $fields[0] =~ /^.+\/(\d+)\/ccs/i if not defined($read);
	($read) = $fields[0] if not defined($read);
	$read = "SEQ_$read";

	# (EXTRA) This below is to show the user an example of parsed line and tells user if it's parsed every 20k read.
	my ($readname) = $fields[0] =~ /(.{20})$/; $readname = "..." . $readname;
	$count{total} ++;
	$seq{$gene}{total} ++;
	print $outLog "\tExample at read $count{total}: name=$CY$readname$N\tstrand=$CY$fields[1]$N\tchr/gene=$CY$fields[2]$N\tpos=$CY$fields[3]$N\tmapQ=$CY$fields[4]$N\n\n" if $count{total} == 1;
	print $outLog "\tDone $GN$count{total}$N\n" if $count{total} % 20000 == 0;

	#discounts if mapping quality is not at least phred score specified
	if($fields[4] < $opt_q)
	{
		print $notused "\t$CY$readname$N quality ($CY$fields[4]$N) is less than $CY$opt_q$N\n";
		$count{lowq} ++;
		$seq{$gene}{lowq} ++;
		next;
	}

	#discounts any read that is not the gene of interest, the proper length, or at the proper location in genome(accounting for indexing)
	elsif(length($fields[9]) < 500) # buffer is obsolete
	{
		# (EXTRA) This below is just for statistics
		$count{badlength} ++;
		$seq{$gene}{badlength} ++;
		die "READNAME undef\n" if not defined($readname);
		die "fields9 undef\n" if not defined($fields[9]);
		die "length of gene $fields[2] undef\n" if not defined($seq{$fields[2]}{len});

		print $notused "\t$CY$readname$N length of seq ($CY" . length($fields[9]). "$N) is less than 500bp (length of original sequence is ($CY" . $seq{$fields[2]}{len} . "$N)!\n";
		next
	}
	
	#counts number of CT conversions (takes into account whether or not the user wants to include Cs in CpG context)
	my $CT = ($fields[13] =~ tr/xhu/xhu/);
	if($opt_c)
	{
		$CT = ($fields[13] =~ tr/zxhu/zxhu/);
	}

	#writes positive reads into <gene>Positive.txt and negative reads into <gene>Negative.txt
	if($fields[1] == 0 || $fields[1] == 16)
	{
		my $to_be_printed = "$CT\t$line\n";
		$fields[1] == 0 ? print $positiveReads $to_be_printed : print $negativeReads $to_be_printed;
		$readz{$read} ++;
		$count{used} ++;
		$seq{$gene}{used} ++;
		$seq{$gene}{pos} ++ if $fields[1] == 0;
		$seq{$gene}{neg} ++ if $fields[1] == 16;
	}
}

my $passedFilterP = `wc -l < $positive`; chomp($passedFilterP);
my $passedFilterN = `wc -l < $negative`; chomp($passedFilterN);

print STDERR "\t${GN}SUCCESS$N: Total=$count{total}, used=$count{used}, Low Map Quality=$count{lowq}, Too short=$count{badlength}\n";
print $outLog "\n\t${GN}SUCCESS$N: Total=$count{total}, used=$count{used}, Low Map Quality=$count{lowq}, Too short=$count{badlength}\n";
print $outLog "\tOutputs are two separate files:\n\t\t- $CY$positive$N ($passedFilterP reads used)\n\t\t- $CY$negative$N ($passedFilterN reads used)\n\n";

print $outLog "
Reads that passed filters:
   Positive: $passedFilterP
   Negative: $passedFilterN
	Total   : $count{total};

Per Gene:
";

foreach my $gene (keys %seq)
{
	print $outLog "
- $gene:
	Positive    = $seq{$gene}{pos}
	Negative    = $seq{$gene}{neg}
	Used        = $seq{$gene}{used}
	Total       = $seq{$gene}{total}
	Too Short   = $seq{$gene}{badlength}
	Low Quality = $seq{$gene}{lowq}
";
}

#sorts by number of CT conversions, takes top 1000, and removes the number of CT conversions in prepartion for methylation extractor
print STDERR "${YW}4. Sorting $positive and $negative by number of CT conversions and removes the number of CT conversions in prepartion for methylation extractor\n";
print $outLog "${YW}4. Sorting $positive and $negative by number of CT conversions and removes the number of CT conversions in prepartion for methylation extractor\n";
my $finalPositive = "$mydir/PositiveFinal.txt";
system("sort -k1,1rn $positive | cut -f 2- > $finalPositive");
my ($finalPositiveLine) = `wc -l $finalPositive` =~ /^(\d+) /;

my $finalNegative = "$mydir/NegativeFinal.txt";
system("sort -k1,1rn $negative | cut -f 2- > $finalNegative");
my ($finalNegativeLine) = `wc -l $finalNegative` =~ /^(\d+) /;
print $outLog "\t${GN}SUCCESS$N: Output:\n\t\t- $CY$finalPositive$N ($finalPositiveLine reads used)\n\t\t- $CY$finalNegative$N ($finalNegativeLine reads used)\n\n";

#runs bismark methylation extractor on top 1000 reads of each strand
my $CPGpos = $mydir . "CpG_context_" . "PositiveFinal.txt";
my $CPGneg = $mydir . "CpG_context_" . "NegativeFinal.txt";
my $CHGpos = $mydir . "CHG_context_" . "PositiveFinal.txt";
my $CHGneg = $mydir . "CHG_context_" . "NegativeFinal.txt";
my $CHHpos = $mydir . "CHH_context_" . "PositiveFinal.txt";
my $CHHneg = $mydir . "CHH_context_" . "NegativeFinal.txt";
my @bismarkOutput = ($CPGpos, $CPGneg, $CHGpos, $CHGneg, $CHHpos, $CHHneg);
my $bismarkOutput = $CHGpos;

print STDERR "${YW}5. Running bismark_methylation_extractor on $CY$finalPositive$YW and $CY$finalNegative$N\n";
print $outLog "${YW}5. Running bismark_methylation_extractor on $CY$finalPositive$YW and $CY$finalNegative$N\n";
if (not $opt_f and (not -e $bismarkOutput or -s $bismarkOutput <= 10))
{
	system("bismark_methylation_extractor -s --comprehensive $finalPositive");
	system("bismark_methylation_extractor -s --comprehensive $finalNegative");
}
print STDERR "\t${GN}SUCCESS$N: Output: 4-6 files of <CpG/CHG/CHH>$CY\_context_$finalPositive$N:\n";
print $outLog "\t${GN}SUCCESS$N: Output: 4-6 files of <CpG/CHG/CHH>$CY\_context_$finalPositive$N:\n";
for (my $i = 0; $i < @bismarkOutput; $i++) {
	if (not -e $bismarkOutput[$i]) {
		print $outLog "\t\t$bismarkOutput[$i]: has$LRD 0$N reads!\n"; system("touch $bismarkOutput[$i]") == 0 or die; next;
	}
	my ($linecount) = `wc -l $bismarkOutput[$i]` =~ /^(\d+) /;
	my ($readnumber) = `unique_column.pl $bismarkOutput[$i] 1 | wc -l` =~ /^(\d+)$/;
	print STDERR "\t\t- $bismarkOutput[$i]: has$LRD 0$N reads!\n" if $linecount == 0;
	print STDERR "\t\t- $bismarkOutput[$i]: $CY$linecount$N total line and $CY$readnumber$N reads\n" if $linecount > 0;
	print $outLog "\t\t- $bismarkOutput[$i]: has$LRD 0$N reads!\n" if $linecount == 0;
	print $outLog "\t\t- $bismarkOutput[$i]: $CY$linecount$N total line and $CY$readnumber$N reads\n" if $linecount > 0;
}

#pulls the CHH, CHG, and CpG sites together into methylationPos<gene>.txt and methylationNeg<gene>.txt (takes into account -c option)
print STDERR "${YW}6. Combine CHH and CHG (and CpG if -c) sites together into$CY methylationPos<gene>.txt$YW and$CY methylationNeg<gene>.txt$N\n";
print $outLog "\n${YW}6. Combine CHH and CHG (and CpG if -c) sites together into$CY methylationPos<gene>.txt$YW and$CY methylationNeg<gene>.txt$N\n";
my $methylationPos = $mydir . "methylationPos" . ".txt";
my $methylationNeg = $mydir . "methylationNeg" . ".txt";
$methylationPos = $mydir . "methylationPos" . "CG.txt" if($opt_c);
$methylationNeg = $mydir . "methylationNeg" . "CG.txt" if($opt_c);

if($opt_c)
{
	system("cat $CPGpos $CHGpos $CHHpos | sort -n > $methylationPos") if not ($opt_f) or not -e $methylationPos or -s $methylationPos < 10;
	system("cat $CPGneg $CHGneg $CHHneg | sort -n > $methylationNeg") if not ($opt_f) or not -e $methylationNeg or -s $methylationNeg < 10;
}
else
{
	system("cat $CHGpos $CHHpos | sort -n > $methylationPos") if not $opt_f or not -e $methylationPos;
	system("cat $CHGneg $CHHneg | sort -n > $methylationNeg") if not $opt_f or not -e $methylationNeg;
}

#gets the position of each conversion (conversions=1; not converted=0)
my %read; my %info;
for (my $i = 0; $i < 2; $i++) {
	
	my $analyzeFile = $i == 0 ? $methylationPos : $methylationNeg;
	
	print $outLog "\t6.$i. Parsing $CY$analyzeFile$N\n";
	open(my $analyze, $analyzeFile) or die "$LRD!!!$N\tFATAL ERROR: Could not open $analyzeFile: $!";
	
	# Initialize array @read which are zeroes with length of @seq
	while(my $line = <$analyze>) {
		chomp($line);
		next if $line  =~ /^Bismark/;
		my @fields = split("\t", $line);
		my ($read) = $fields[0] =~ /^.+\/(\d+)\/\d+/i;
		($read) = $fields[0] =~ /^.+\/(\d+)\/ccs/i if not defined($read);
		($read) = $fields[0] if not defined($read);
		$read = "SEQ_$read";
		my ($readname) = length($fields[0]) > 20 ? $fields[0] =~ /(.{20})$/ : $fields[0]; $readname = "..." . $readname;
		my ($name, $junk, $gene, $pos, $conv) = @fields; $pos = $pos - 1;
		print $outLog "\tExample: name=$CY$readname$N, junk=$CY$fields[1]$N, gene=$CY$fields[2]$N, pos=$CY$fields[3]$N, conv=$CY$fields[4]$N\n" if not defined($read{$gene}{$i});
		$read{$gene}{$i}{$read}{$pos} = $conv =~ /^[xzhu]$/ ? 1 : $conv =~ /^[XZHU]$/ ? 0 : die "$LRD!!!$N\tFATAL ERROR: conversion isn't x/z/h/u (case ins) ($CY$conv$N)in $CY$analyzeFile$N line:\n$line\n\n";
		$info{$gene}{$read}{min} = $pos if not defined($info{$gene}{$read}{min}) or $info{$gene}{$read}{min} > $pos;
		$info{$gene}{$read}{max} = $pos if not defined($info{$gene}{$read}{max}) or $info{$gene}{$read}{max} < $pos;
		#$read{$i}{$read}{$pos - $buffer} = $conv =~ /^[xzhu]$/ ? 1 : $conv =~ /^[XZHU]$/ ? 0 : die "$LRD!!!$N\tFATAL ERROR: conversion isn't x/z/h/u (case ins) ($CY$conv$N)in $CY$analyzeFile$N line:\n$line\n\n";
		#$info{$read}{min} = $pos - $buffer if not defined($info{$read}{min}) or $info{$read}{min} > $pos - $buffer;
		#$info{$read}{max} = $pos - $buffer if not defined($info{$read}{max}) or $info{$read}{max} < $pos - $buffer;
	}
}

print STDERR "\t${GN}SUCCESS$N: Done parsing methylationPos and methylationNeg files!\n";
print $outLog "\n\t${GN}SUCCESS$N: Done parsing methylationPos and methylationNeg files!\n";

print STDERR "${YW}7. Converting each gene's sequence position and methylation data\n$N";
print $outLog "${YW}7. Converting each gene's sequence position and methylation data\n$N";
foreach my $gene (sort keys %seq) {
	print $outLog "\t- Processing gene $CY$gene$N\n";
	my $filePos = $opt_c ? "$mydir/$gene\_CG_POS" : "$mydir/$gene\_POS";
	my $fileNeg = $opt_c ? "$mydir/$gene\_CG_NEG" : "$mydir/$gene\_NEG";

	my @seq = @{$seq{$gene}{seq}};
	open(my $FILEPOS, ">", "$filePos.tsv") or die "Could not open $filePos.tsv: $!";
	open(my $FILENEG, ">", "$fileNeg.tsv") or die "Could not open $fileNeg.tsv: $!";
	open(my $FILEPOSFA, ">", "$filePos.customfa") or die "Could not open $filePos.customfa: $!";
	open(my $FILENEGFA, ">", "$fileNeg.customfa") or die "Could not open $fileNeg.customfa: $!";
	open(my $FILEPOSTRANS, ">", "$filePos.trans") or die "Could not open $filePos.trans: $!";
	open(my $FILENEGTRANS, ">", "$fileNeg.trans") or die "Could not open $fileNeg.trans: $!";

	foreach my $strand (sort keys %{$read{$gene}}) {
		my $checkprint = 0;
		foreach my $name (keys %{$read{$gene}{$strand}}) {
			my @value; my @fasta;
			for (my $i = 0; $i < @seq; $i++) {
				my $nuc1 = $strand == 0 ? "C" : "G";
				my $nuc2 = $strand == 0 ? "G" : "C";
				my $add = $strand == 0 ? 1 : -1;
				if ($seq[$i] ne $nuc1 or ($seq[$i] eq $nuc1 and $i != @seq - 1 and $seq[$i+$add] eq $nuc2)) {
					if ($info{$gene}{$name}{min} > $i or $info{$gene}{$name}{max} < $i) {
						$value[$i] = "NA";
					}
					else {
						#if ($strand == 0) {
						my $currvalue = $read{$gene}{$strand}{$name}{$i}; $currvalue = "NA" if not defined($currvalue);
						$value[$i] = 2 if not $opt_c;
						$value[$i] = 2 if $seq[$i] ne $nuc1;
						$value[$i] = $currvalue + 5 if $seq[$i] eq $nuc1 and $opt_c and $currvalue ne "NA";
						$value[$i] = "NA" if $opt_c and $currvalue eq "NA";
						#}
						#else {
						#	my $currvalue2 = $read{$gene}{$strand}{$name}{$i+1}; $currvalue2 = "NA" if not defined($currvalue2);
						#	my $currvalue = $read{$gene}{$strand}{$name}{$i}; $currvalue = "NA" if not defined($currvalue);
						#	if (not $opt_c) {$value[$i] = 2}
						#	else {$value[$i] = 2; $value[$i-1] = $currvalue2 + 5
						#}
					}
				}
				else {
					my $currvalue = $read{$gene}{$strand}{$name}{$i}; $currvalue = "NA" if not defined($currvalue);
					if ($info{$gene}{$name}{min} > $i or $info{$gene}{$name}{max} < $i) {
						$value[$i] = "NA";
					}
					else {
						$value[$i] = $currvalue;
					}
					$currvalue = 0 if $currvalue eq "NA";
					my $faindice = defined($fasta[0]) ? @fasta : 0;
					push(@fasta, $currvalue);
					print $FILEPOSTRANS "$faindice\t$i\n"if $strand == 0 and $checkprint == 0;
					print $FILENEGTRANS "$faindice\t$i\n"if $strand == 1 and $checkprint == 0;
				}
			}
			print $FILEPOS "$name\t" . join("\t", @value) . "\n" if $strand == 0;
			print $FILENEG "$name\t" . join("\t", @value) . "\n" if $strand == 1;
			print $FILEPOSFA ">$name\n" . join("", @fasta) . "\n" if $strand == 0;
			print $FILENEGFA ">$name\n" . join("", @fasta) . "\n" if $strand == 1;
			$checkprint = 1;
		}
	}
	close $FILEPOS;
	close $FILENEG;
	close $FILEPOSFA;
	close $FILENEGFA;
	close $FILEPOSTRANS;
	close $FILENEGTRANS;
	# To check: fastaFromBed -fi $geneIndexesFa -bed test.bed -fo test.fa && cat test.fa
	# ALL MUST BE C
	my ($filePosLine) = `wc -l $filePos.tsv` =~ /^(\d+) /;
	my ($fileNegLine) = `wc -l $fileNeg.tsv` =~ /^(\d+) /;
	print STDERR "\t${GN}SUCCESS$N: Done parsing $YW$gene$N Output:\n\t\t- $CY$filePos.tsv$N ($filePosLine reads)\n\t\t- $CY$fileNeg.tsv$N ($fileNegLine reads)\n";
	print $outLog "\t${GN}SUCCESS$N: Done parsing $YW$gene$N Output:\n\t\t- $CY$filePos.tsv$N ($filePosLine reads)\n\t\t- $CY$fileNeg.tsv$N ($fileNegLine reads)\n\n";
	#system("mv $filePos.tsv $fileNeg.tsv $filePos.customfa $fileNeg.customfa $filePos.trans $fileNeg.trans $mydir/");
	my ($finalmd5) = "$filePos.tsv:" .`md5sum $filePos.tsv`;
	($finalmd5) = "$fileNeg.tsv:" .`md5sum $fileNeg.tsv`;
	($finalmd5) .= "$filePos.customfa:" .`md5sum $filePos.customfa`;
	($finalmd5) .= "$fileNeg.customfa:" .`md5sum $fileNeg.customfa`;
	($finalmd5) .= "$filePos.trans:" .`md5sum $filePos.trans`;
	($finalmd5) .= "$fileNeg.trans:" .`md5sum $fileNeg.trans`;
	print $outLog "$finalmd5\n";
}

foreach my $gene (sort keys %seq) {
	print $outLog "$mydir/exon/$gene.exon:" . `md5sum $mydir/exon/$gene.exon`;
}
my @seq;

print $outLog "${YW}8. determines which regions of conversion are R-loops based on conversion threshold$N\n";
print STDERR "${YW}8. determines which regions of conversion are R-loops based on conversion threshold$N\n";
#0=not converted (grey)
#1=converted (green)
#2=non-C (white)
#3=converted CpG (blue)
#4=non-converted CpG  (black)
#5=converted CpG in R-loop (purple)
#9=converted C in R-loop (red)
my $start = 1;
my $end = $opt_l;
my $conC = 0;
my $nonConC = 0;
my $conPer = 0;

foreach my $gene (sort keys %seq) {
my @seq = @{$seq{$gene}{seq}};
my $finalPos = $mydir . "$gene\_Pos" . ($opt_t*100) . ".txt";
my $finalNeg = $mydir . "$gene\_Neg" . ($opt_t*100) . ".txt";
$finalPos = $mydir . "$gene\_Pos" . ($opt_t*100) . "CG.txt" if($opt_c);
$finalNeg = $mydir . "$gene\_Neg" . ($opt_t*100) . "CG.txt" if($opt_c);
open(FINALPOS, ">", $finalPos) or die "Could not open $finalPos: $!";
open(FINALNEG, ">", $finalNeg) or die "Could not open $finalNeg: $!";
for(my $i=0; $i<2; $i++)
{
	my $filePos = $opt_c ? "$mydir/$gene\_CG_POS" : "$mydir/$gene\_POS";
	my $fileNeg = $opt_c ? "$mydir/$gene\_CG_NEG" : "$mydir/$gene\_NEG";
	my $fileLast = $i == 0 ? "$filePos.tsv" : "$fileNeg.tsv";
	print $outLog "\t- Doing gene $CY$gene$N ($fileLast)\n";
	print $outLog "$fileLast doesn't exist!\n" and next if not -e $fileLast;
	print $outLog "$fileLast has no line!\n" and next if `wc -l < $fileLast` == 0;
	my @lineFinal = `cat $fileLast`;
	for (my $p =0 ; $p < @lineFinal; $p++) 
	{
		print $outLog "\tDone $GN$p$N\n" if $p % 100 eq 0;
		my $lineFinal = $lineFinal[$p];
		chomp($lineFinal);
		my @fields = split("\t", $lineFinal);
		for(my $z=0; $z<@fields; $z++)
		{
			$fields[$z] =~ s/ //g;
		}
		for(my $z=0; $z<@seq; $z++)
		{
			#positive reads
			if($i eq 0)
			{
				if($seq[$z] ne "C")
				{
					$fields[$z+1] = 2;
				}
				if($opt_c)
				{
					if($z != (@seq-1) && $seq[$z] eq "C" && $seq[$z+1] eq "G" && $fields[$z+1] eq 1)
					{
						$fields[$z+1] = 3;
					}
					if($z != (@seq-1) && $seq[$z] eq "C" && $seq[$z+1] eq "G" && $fields[$z+1] eq 0)
					{
						$fields[$z+1] = 4;
					}	
				}
				else
				{
					if($z != (@seq-1) && $seq[$z] eq "C" && $seq[$z+1] eq "G")
					{
						$fields[$z+1] = 2;
					}
				}	
			}
			#negative reads
			if($i eq 1)
			{
				if($seq[$z] ne "G")
				{
					$fields[$z+1] = 2;
				}
				if($opt_c)
				{
					if($z != (@seq-1) && $seq[$z] eq "G" && $seq[$z+1] eq "C" && $fields[$z+1] eq 1)
					{
						$fields[$z+1] = 3;
					}
					if($z != (@seq-1) && $seq[$z] eq "G" && $seq[$z+1] eq "C" && $fields[$z+1] eq 0)
					{
						$fields[$z+1] = 4;
					}
				}
				else
				{
					if($z != (@seq-1) && $seq[$z] eq "G" && $seq[$z+1] eq "C")
					{
						$fields[$z+1] = 2;
					}
				}
			}
		}
		while($end<@fields)
		{
			for(my $temp = $start; $temp<$end; $temp++)
			{
				if($fields[$temp] eq 1 || $fields[$temp] eq 9 || $fields[$temp] eq 3 || $fields[$temp] eq 5)
				{
					$conC++;
				}
				if($fields[$temp] eq 0 || $fields[$temp] eq 4)
				{
					$nonConC++;
				}
			}
			if($conC != 0 || $nonConC != 0)
			{
				$conPer = (($conC)/($conC+$nonConC));
			}
			else
			{
				$conPer = 0;
			}
			if($conPer >= $opt_t)
			{
				for(my $k=$start; $k<=$end; $k++)
				{
					if($fields[$k] eq 1)
					{
						$fields[$k] = 9;
					}
					if($fields[$k] eq 3)
					{
						$fields[$k] = 5;
					}
				}
			}
			$start++;
			$end++;
			$conC = 0;
			$nonConC = 0;
		}
		my $newfield = join("\t", @fields);
		if($i eq 0)
		{
			print FINALPOS "$newfield\n";
		}
		if($i eq 1)
		{
			print FINALNEG "$newfield\n";
		}
		$start = 1;
		$end = 100;
		$conC = 0;
		$nonConC = 0;
		$conPer = 0;
	} 
}

my ($finalPosLine) = `wc -l $finalPos` =~ /^(\d+) /;
my ($finalNegLine) = `wc -l $finalNeg` =~ /^(\d+) /;

print $outLog "\n\t${GN}SUCCESS$N: Gene $CY$gene$N Output:\n\t\t- $CY$finalPos$N ($finalPosLine reads)\n\t\t- $CY$finalNeg$N ($finalNegLine reads)\n\n";
print STDERR "\t${GN}SUCCESS$N: Gene $CY$gene$N Output:\n\t\t- $CY$finalPos$N ($finalPosLine reads)\n\t\t- $CY$finalNeg$N ($finalNegLine reads)\n\n";

#makes heatmaps for positive and negative strand
my $finalPosPDF = $mydir . "$gene\_Pos" . ($opt_t*100) . ".pdf";
my $finalNegPDF = $mydir . "$gene\_Neg" . ($opt_t*100) . ".pdf";
$finalPosPDF = $mydir . "$gene\_Pos" . ($opt_t*100) . "CG.pdf" if($opt_c);
$finalNegPDF = $mydir . "$gene\_Neg" . ($opt_t*100) . "CG.pdf" if($opt_c);

my $Rscript = "$mydir/$gene\_MakeHeatmap.R";
open(my $out, ">", $Rscript) or die "Can't print to $Rscript: $!\n";
if($opt_c)
{
	print $out "
.libPaths()
		library(\"GMD\")
		df = read.table(\"$finalPos\", sep=\"\t\", row.names=1)
		df2 = df
		df2[df2 != 9] = 0
		df\$sum = apply(df2,1,sum)
		df = df[order(-df\$sum),]
		dimz = 200
		if(dim(df)[1] < 200)
		{
			dimz = dim(df)[1]
		}
		df = df[1:dimz,]
		df = subset(df,select=-sum)
		df2 = df
		df2[df2 != 9] = 0
		h = heatmap.3(df2,
            dendrogram=\"row\",
            Rowv=TRUE, Colv=FALSE,
            labRow=FALSE, labCol=FALSE)
		id = h\$rowInd
		id = rev(id)
		df = df[id,] 
		pdf(\"$finalPosPDF\")
		heatmap.3(
			x=df,
			dendrogram=\"none\",
			Rowv=TRUE, Colv=FALSE,
			labRow=FALSE,labCol=FALSE,
			breaks=c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5,9.5),
			color.FUN=function(x) c(\"grey\",\"green\",\"white\",\"blue\",\"black\",\"yellow\",\"red\")
		)
		dev.off()
		df3 = read.table(\"$finalNeg\", sep=\"\t\", row.names=1)
		df4 = df3
      df4[df4 != 9] = 0
      df3\$sum = apply(df4,1,sum)
      df3 = df3[order(-df3\$sum),]
		dimz = 200
		if(dim(df3)[1] < 200)
		{
			dimz = dim(df3)[1]
		}
      df3 = df3[1:dimz,]
      df3 = subset(df3,select=-sum)
      df4 = df3
      df4[df4 != 9] = 0
      h = heatmap.3(df4,
            dendrogram=\"row\",
            Rowv=TRUE, Colv=FALSE,
            labRow=FALSE, labCol=FALSE)
      id = h\$rowInd
      id = rev(id)
      df3 = df3[id,] 
		pdf(\"$finalNegPDF\")
      heatmap.3(
         x=df3,
         dendrogram=\"none\",
         Rowv=TRUE, Colv=FALSE,
         labRow=FALSE,labCol=FALSE,
         breaks=c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5,9.5),
         color.FUN=function(x) c(\"grey\",\"green\",\"white\",\"blue\",\"black\",\"yellow\",\"red\")
      )
		dev.off()
	";
}
else
{
		print $out "
.libPaths()
      library(\"GMD\")
      df = read.table(\"$finalPos\", sep=\"\t\", row.names=1)
      df2 = df
      df2[df2 != 9] = 0
      df\$sum = apply(df2,1,sum)
      df = df[order(-df\$sum),]
		dimz = 200
		if(dim(df)[1] < 200)
		{
			dimz = dim(df)[1]
		}
      df = df[1:dimz,]
      df = subset(df,select=-sum)
      df2 = df
      df2[df2 != 9] = 0
      h = heatmap.3(df2,
            dendrogram=\"row\",
            Rowv=TRUE, Colv=FALSE,
            labRow=FALSE, labCol=FALSE)
      id = h\$rowInd
      id = rev(id)
      df = df[id,] 
      pdf(\"$finalPosPDF\")
      heatmap.3(
         x=df,
         dendrogram=\"none\",
         Rowv=TRUE, Colv=FALSE,
         labRow=FALSE,labCol=FALSE,
      	breaks=c(-0.5,0.5,1.5,2.5,9.5),
			color.FUN=function(x) c(\"grey\",\"blue\",\"white\",\"red\")
		)
		dev.off()
      df3 = read.table(\"./$finalNeg\", sep=\"\t\", row.names=1)
      df4 = df3
      df4[df4 != 9] = 0
      df3\$sum = apply(df4,1,sum)
      df3 = df3[order(-df3\$sum),]
		dimz = 200
		if(dim(df3)[1] < 200)
		{
			dimz = dim(df3)[1]
		}
      df3 = df3[1:dimz,]
      df3 = subset(df3,select=-sum)
      df4 = df3
      df4[df4 != 9] = 0
      h = heatmap.3(df4,
            dendrogram=\"row\",
            Rowv=TRUE, Colv=FALSE,
            labRow=FALSE, labCol=FALSE)
      id = h\$rowInd
      id = rev(id)
      df3 = df3[id,] 
      pdf(\"$finalNegPDF\")
      heatmap.3(
         x=df3,
         dendrogram=\"none\",
         Rowv=TRUE, Colv=FALSE,
         labRow=FALSE,labCol=FALSE,
      	breaks=c(-0.5,0.5,1.5,2.5,9.5),
			color.FUN=function(x) c(\"grey\",\"blue\",\"white\",\"red\")
		)
      dev.off()
   ";
}
close $out;
my $cmd = "R --no-save < $Rscript >> $logFile\n";
print ($cmd);
system($cmd);
system("rm Rplots.pdf") if -e "Rplots.pdf";
}
