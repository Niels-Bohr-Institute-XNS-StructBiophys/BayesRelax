#!/usr/bin/perl
# What does this script do: a commandline version of Steen's relaxation spectrum program. Gnuplot ver. 5 required
# UCPH PLEN: Peter Ulvskov, ulvskov@plen.ku.dk
use strict;
use warnings;
use File::Basename;
use File::Copy;

my $numberOfArgs = @ARGV;

if($numberOfArgs !=2) 
{ 
	die "This utility must be fed with a three column rheometer file, (Frequency, G' and G'') and the desired number of points in the spectrum (30-100)\n"; 
}

my $numberOfPoints = $ARGV[1];
if($numberOfPoints < 30 || $numberOfPoints > 100)
{
	print "$numberOfPoints should fall in the range 30-100\n";
	exit;
}




# _____We create inputfile.d for the optional parameters of which we only handle desired number of points in spectrum for the moment

my $inputfile = "inputfile.d";

open(my $fh0 ,">", "inputfile.d") or die "could not write to inputfile.d\n";
print $fh0 "\n\n\n\n\n$numberOfPoints\n\n";
close $fh0;
#______________________

#________open rheological data file________________________________________
open (my $fh1, $ARGV[0]) or die "could not open $ARGV[0]\n";
my $rheoTable = $ARGV[0];
my($filename, $dirs, $suffix) = fileparse($rheoTable);
# convert spaces to underscores
$filename =~ tr/ /_/;
#add back file extension
$filename.= $suffix;
#copy the rheometer data into this new file (which now has a Fortran friendly name)
open(my $fh2 ,">", $filename) or die "could not write to $filename\n";
copy($fh1, $fh2) or die "copy failed: $!";
close $fh1;
close $fh2;



#________________call the error2 executable on the copy of the rhemomter file__________

`./error2 $filename > output.d`;

#error2 produces a file with the name $filename prefixed with 'er'. Standard deviations have been calculated and data has been reformatted with G' first then G''
# we copy this file back onto the original:

move "er$filename", $filename;

#________________call the bayes executable on the copy of the rhemomter file__________

`./bayes $filename >> output.d`;

#__A flurry of files are created that we need to rename to standardised names expected later

move "fit1$filename", "fit1.d";
move "fit2$filename", "fit2.d";
move "trans_me$filename", "fit.d";
move "me$filename", "estimate.d";
move "in_me$filename", "data.d";

#________________call the conv executable on the newly created files


`./conv`;


#___Clean up a bit________


unlink("prob_me$filename");
unlink("sd$filename");
unlink("$filename");
unlink("fort.55");
unlink("fort.56");
unlink("output.d");
unlink("inputfile.d");

#_________Plot data________________

`gnuplot plot1.pl`;
`gnuplot plot2.pl`;

#____delete plot scripts_________

unlink("plot1.pl");
unlink("plot2.pl");

#_____the relevant output files should be turned into tab-delimited files.

my $directory = "Results";

    unless(-e $directory or mkdir $directory) 
    {
        die "Unable to create $directory\n";
    }    
 
 my $relaxationSpectrum = "Results/RelaxationSpectrum_$filename";   
 open(my $fh3 ,">", $relaxationSpectrum) or die "could not write to $relaxationSpectrum\n";   
 
 open(my $fh4, 'estimate.d')  or die "could not estimate.d\n";
 while(<$fh4>)
 {
	if($. < 4)
	{
		$_ =~  s/^\s+//;
		print $fh3 $_;
		next;
	}	
	chomp;
	my $thisLine = $_;
	#Fortran cannot do tab - we need to get rid of the variable numbers of spaces that it inserts.
	$thisLine  =~  s/^\s+//;
	$thisLine  =~  s/\s+$//;
	$thisLine  =~  s/\s+/\t/g;
	$thisLine  .= "\n";
	print $fh3 $thisLine;
 }	
close $fh3;
close $fh4;

 my $rheoDataFit = "Results/FitOf_$filename";   

 open(my $fh5 ,">", $rheoDataFit) or die "could not write to $rheoDataFit\n";
 print $fh5 "Frequency\tMeasured G'\tMeasured G''\tG'-fit\tG''-fit\tG'-error\tG''error\n";
  

open(my $fh6, 'data.d') or die "could not open data.d\n";
my @rawData = <$fh6>;
chomp(@rawData);
my $linesInRawdata = @rawData;
close $fh6;
#the file has a seven line header before the real data begins. Remember G' is then the first half, G'' the second half
$linesInRawdata = (@rawData-7)/2;
s/^\s+// for @rawData;
s/\s+$// for @rawData;
s/\s+/\t/g for @rawData;

open(my $fh7, 'fit1.d') or die "could not open fit1.d\n";
my @storageModuli = <$fh7>;
close $fh7;
chomp(@storageModuli);
s/^\s+// for @storageModuli;
s/\s+$// for @storageModuli;
s/\s+/\t/g for @storageModuli;

open(my $fh8, 'fit2.d') or die "could not open fit2.d\n";
my @lossModuli = <$fh8>;
close $fh8;
chomp(@lossModuli);
s/^\s+// for @lossModuli;
s/\s+$// for @lossModuli;
s/\s+/\t/g for @lossModuli;

#print headers for the fit table


#and put the table together

for(my $i = 0; $i < $linesInRawdata; $i++)
{
	my ($freq, $rawGprime, $GprimeError) = split /\t/, $rawData[7+$i];
	my (undef, $rawGdoubleprime, $GdoublePrimeError) = split /\t/, $rawData[$linesInRawdata+7+$i];
	my (undef, $fitGprime) = split /\t/, $storageModuli[$i];
	my (undef, $fitGdoubleprime) = split /\t/, $lossModuli[$i];
	print $fh5 "$freq\t$rawGprime\t$rawGdoubleprime\t$fitGprime\t$fitGdoubleprime\t$GprimeError\t$GdoublePrimeError\n";
	
}
close $fh5;
#___rename the figures and put them into the Results folder

my $fig1File = "Results/fig1$filename";
$fig1File =~ s/\.txt$/.pdf/;
my $fig2File = "Results/fig2$filename";
$fig2File =~ s/\.txt$/.pdf/;
move "fig1.pdf", $fig1File;
move "fig2.pdf", $fig2File;

#__________Clean up________

unlink("fig1.pdf");
unlink("fig2.pdf");
unlink("estimate.d");
unlink("data.d");
unlink("data2.d");
unlink("dummy.d");
unlink("fit.d");
unlink("fit1.d");
unlink("fit2.d");

	
				


