#!/usr/bin/perl

use warnings; use strict; use colorz; use Getopt::Std;
use vars qw($opt_n $opt_t $opt_l $opt_r $opt_q $opt_s $opt_g $opt_c);
getopts("n:t:l:r:q:s:g:c");

my $usage = "Usage: $YW$0$N ${GN}[options]$N -n $CY<gene name>$N -t $CY<conversion threshold>$N -l $CY<minimum R-loop length>$N -r $CY<path to reads>$N -q $CY<minimum mapping quality for each read to be considered>$N -s $CY<path to original gene sequence>$N -g $CY<path to reference genome for indexing>$N
${GN}Options:$N
-c: consider Cs in CpG context
${GN}Extra information:$N
-n: gene name must be in all caps (e.g. CALM3)
-t: percentage (0.0-1.0) of Cs that must be converted to be considered an R-loop
-l: minimum length in base pairs to be considered an R-loop
-r: must be path to reads file itself named pacbio.fastq
-q: must be phred score (e.g. score of 40 = 99.99% base call accuracy)
-s: must be path to sequence file itself
-g: must be path to .fa file containing reference genome used for indexing
";

die $usage if not ($opt_n) or not ($opt_t) or not ($opt_l) or not ($opt_r) or not defined ($opt_q) or not ($opt_s) or not ($opt_g);

#runs bismark (output file will be pacbio.fastq_bismark_bt2.sam) only if it hasn't been ran previously
unless(-e "./pacbio.fastq_bismark_bt2.sam")
{
	if( -e "./geneIndexes.bed")
   {	
     system("bedtools getfasta -fi $opt_g -bed geneIndexes.bed -fo geneIndexes.fa -name");
   }
   else
   {
      die "$YW Error:$N geneIndexes.bed does not exist in this directory";
   }
	system("bismark_genome_preparation --bowtie2 ./ > /dev/null");

   system("bismark --bowtie2 --rdg 2,1 --rfg 2,1 --score_min L,0,-0.8 ./ $opt_r > /dev/null") == 0 or die "Failed to run bismark: $!\n";
}

#takes sequence of gene and splits into array of individual bases
open(SEQ, $opt_s) or die "Could not open $opt_s: $!";
my @seq = split("", <SEQ>);
close SEQ;

my $positive = $opt_n . "Positive.txt";
my $negative = $opt_n . "Negative.txt";
open(my $positiveReads, ">", $positive) or die "Could not open $positive: $!";
open(my $negativeReads, ">", $negative) or die "Could not open $negative: $!";

my $samFile = "pacbio.fastq_bismark_bt2.sam";
open(my $sam, $samFile) or die "Could not open $samFile: $!";

my $logFile = $opt_n . "logFile.txt";
open(my $fh, '>', $logFile);

my $countSAM = "countSAM.txt";
`grep "$opt_n" $samFile > $countSAM`;

my $countS = `wc -l < $countSAM`;
`rm countSAM.txt`;
print $fh "Number of reads in SAM file (positive and negative): $countS";

#loops through each read and writes high quality reads into two separate files <gene>Positive.txt and <gene>Negative.txt
while(my $line = <$sam>)
{
	chomp($line);
   my @fields = split("\t", $line);
	#discounts the first 4 lines of information at the top of the .sam file
   if(@fields < 6)
   {
      next
   }
	#discounts if mapping quality is not at least phred score specified
   if($fields[4] < $opt_q)
   {
      next
   }
	#discounts any read that is not the gene of interest, the proper length, or at the proper location in genome(accounting for indexing)
	if($fields[2] ne "$opt_n" || length($fields[9]) < 500 || $fields[3] < 45 || $fields[3] > (@seq+55))
   {
      next
   }
	#writes positive reads into <gene>Positive.txt and negative reads into <gene>Negative.txt
   if($fields[1] == 0 && $fields[2] eq "$opt_n")
   {
      print $positiveReads "$line\n";
   }
   if($fields[1] == 16 && $fields[2] eq "$opt_n")
   {
      print $negativeReads "$line\n";
   }
}

my $passedFilterP = `wc -l < $positive`;
my $passedFilterN = `wc -l < $negative`;

print $fh "
Reads that passed filters:
   Positive: $passedFilterP
   Negative: $passedFilterN
";

my $CPGpos = "CpG_context_" . $positive;
my $CPGneg = "CpG_context_" . $negative;
my $CHGpos = "CHG_context_" . $positive;
my $CHGneg = "CHG_context_" . $negative;
my $CHHpos = "CHH_context_" . $positive;
my $CHHneg = "CHH_context_" . $negative;

#runs bismark methylation extractor on reads of each strand
my $bismarkOutput = "./" . $CPGpos;
unless(-e $bismarkOutput and -s $bismarkOutput >= 10)
{
	system("bismark_methylation_extractor -s --comprehensive $positive");
	system("bismark_methylation_extractor -s --comprehensive $negative");
}

#pulls the CHH, CHG, and CpG sites together into methylationPos<gene>.txt and methylationNeg<gene>.txt (takes into account -c option)
my $methylationPos = "methylationPos" . $opt_n . ".txt";
my $methylationNeg = "methylationNeg" . $opt_n . ".txt";
$methylationPos = "methylationPos" . $opt_n . "CG.txt" if($opt_c);
$methylationNeg = "methylationNeg" . $opt_n . "CG.txt" if($opt_c);

if($opt_c)
{
	system("cat $CPGpos $CHGpos $CHHpos | sort -n > $methylationPos");
	system("cat $CPGneg $CHGneg $CHHneg | sort -n > $methylationNeg");
}
else
{
	system("cat $CHGpos $CHHpos | sort -n > $methylationPos");
	system("cat $CHGneg $CHHneg | sort -n > $methylationNeg");
}

my $filePos = $opt_n . "PosPositions.txt";
my $fileNeg = $opt_n . "NegPositions.txt";

open(FILEPOS, ">", $filePos) or die "Could not open $filePos: $!";
open(FILENEG, ">", $fileNeg) or die "Could not open $fileNeg: $!";

#gets the position of each conversion (conversions=1; not converted=0)
my $j = 50;
for(my $i=0; $i<2; $i++)
{
	my $analyzeFile = $methylationPos if($i==0);
	$analyzeFile = $methylationNeg if($i==1);
	open(my $analyze, $analyzeFile) or die "Could not open $analyzeFile: $!";
	my (@read) = ("0\t")x(@seq);
	while(my $line = <$analyze>)
	{
		chomp($line);
		next if $line  =~ /Bismark/;
		my @fields = split("\t", $line);
		if($read[0] ne "$fields[0]\t")
   	{
      	if(@read == (@seq) && $read[0] ne "0\t" && $read[0] ne "1\t" && $read[0] ne "Bismark methylation extractor version v0.13.1\t")
      	{
				if($i==0)
				{
         		print FILEPOS "@read\n";
				}
				if($i==1)
				{
					print FILENEG "@read\n";
      		}
			}
     		@read = ();
      	@read = ("0\t")x(@seq);
      	$read[0] = "$fields[0]\t";
   	}
		if($read[0] eq "$fields[0]\t")
   	{
      	if($fields[4] eq "x" || $fields[4] eq "h" || $fields[4] eq "z" || $fields[4] eq "u" && $fields[3]>$j)
      	{
         	my $index = ($fields[3] - $j + 1);
         	$read[$index] = "1\t";
      	}
   	}
	}
}

#determines which regions of conversion are R-loops based on conversion threshold
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

my $finalPos = $opt_n . "Pos" . ($opt_t*100) . ".txt";
my $finalNeg = $opt_n . "Neg" . ($opt_t*100) . ".txt";
$finalPos = $opt_n . "Pos" . ($opt_t*100) . "CG.txt" if($opt_c);
$finalNeg = $opt_n . "Neg" . ($opt_t*100) . "CG.txt" if($opt_c);
open(FINALPOS, ">", $finalPos) or die "Could not open $finalPos: $!";
open(FINALNEG, ">", $finalNeg) or die "Could not open $finalNeg: $!";
for(my $i=0; $i<2; $i++)
{
	my $fileLast = $opt_n . "PosPositions.txt" if($i==0);
	$fileLast = $opt_n . "NegPositions.txt" if($i==1);
	my @lineFinal = `cat $fileLast`;
  	for (my $p =0 ; $p < @lineFinal; $p++) 
	{
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
			if($i==0)
			{
				if($seq[$z] ne "C")
				{
					$fields[$z+1] = 2;
				}
				if($opt_c)
				{
					if($z != (@seq-1) && $seq[$z] eq "C" && $seq[$z+1] eq "G" && $fields[$z+1] == 1)
      			{
         			$fields[$z+1] = 3;
      			}
      			if($z != (@seq-1) && $seq[$z] eq "C" && $seq[$z+1] eq "G" && $fields[$z+1] == 0)
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
			if($i==1)
			{
				if($seq[$z] ne "G")
				{
					$fields[$z+1] = 2;
				}
				if($opt_c)
				{
					if($z != (@seq-1) && $seq[$z] eq "G" && $seq[$z+1] eq "C" && $fields[$z+1] == 1)
               {
                  $fields[$z+1] = 3;
               }
               if($z != (@seq-1) && $seq[$z] eq "G" && $seq[$z+1] eq "C" && $fields[$z+1] == 0)
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
				if($fields[$temp] == 1 || $fields[$temp] == 9 || $fields[$temp] == 3 || $fields[$temp] == 5)
				{
            	$conC++;
         	}
         	if($fields[$temp] == 0 || $fields[$temp] == 4)
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
					if($fields[$k] == 1)
					{
						$fields[$k] = 9;
					}
					if($fields[$k] == 3)
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
		if($i==0)
		{
			print FINALPOS "$newfield\n";
		}
		if($i==1)
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
my $countRloopP = "countRloopP.txt";
my $countRloopN = "countRloopN.txt";
`grep "9	" $finalPos > $countRloopP`;
`grep "9	" $finalNeg > $countRloopN`;

my $countP = `wc -l < $countRloopP`;
my $countN = `wc -l < $countRloopN`;

print $fh "
Number of reads containing footprints:
   Positive: $countP
   Negative: $countN
";

`rm $countRloopP`;
`rm $countRloopN`;
#makes heatmaps for positive and negative strand
my $finalPosPDF = $opt_n . "Pos" . ($opt_t*100) . ".pdf";
my $finalNegPDF = $opt_n . "Neg" . ($opt_t*100) . ".pdf";
$finalPosPDF = $opt_n . "Pos" . ($opt_t*100) . "CG.pdf" if($opt_c);
$finalNegPDF = $opt_n . "Neg" . ($opt_t*100) . "CG.pdf" if($opt_c);

my $Rscript = "MakeHeatmap" . $opt_n . ".R";
open(my $out, ">", $Rscript) or die "Can't print to $Rscript: $!\n";
if($opt_c)
{
	print $out "
.libPaths()
		library(\"GMD\")
		df = read.table(\"./$finalPos\", sep=\"\t\", row.names=1)
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
      df = read.table(\"./$finalPos\", sep=\"\t\", row.names=1)
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
system("R --no-save < $Rscript");
system("rm Rplots.pdf");
