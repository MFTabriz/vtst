#!/usr/bin/env perl
#;-*- Perl -*-

# Program uses the freq.dat to extract the normal modes and outputs the modes with the amplitude, to display the spectra.

@args = @ARGV;
@args == 4 || die "usage: dymspect.pl <freq.dat> <sigma: variance> <minimum frequency> <maximum frequency>\n";
$freqfile = $args[0];
$sigma = $args[1];
$x_min = $args[2];
$x_max = $args[3];

$x_min >=1 || die "Note: minimum frequency must be >= 1.\n";
$x_max >1 || die "Note: maximum frequency must be > 1.\n";
$sigma >0 || die "Note: sigma (variance) must be > 0.\n";

# Removing all spaces to one blank
open (IN,$freqfile);
while (<IN>) {
    $_ =~ s/^\s+//g;
    $file1 .= $_;
}
close (IN);

# Placing all characters separated by a new line into array
@freqfile = split(/\n/,$file1);

$numimg = 0;
@freqValue;

# Placing each freqfile array value separated by spaces into new array
for ($i=0; $i<@freqfile; $i++) {
    @line1 = split(/\s+/,$freqfile[$i]);
    
    if($line1[3] == 1) {
        $numimg++;
	next;
    }
    
    $freqValue[$i-$numimg] = $line1[0];

}

print " Note: freq.dat has $numimg imaginary frequencies.\n";
# print join("\n",@freqValue),"\n";

 # Calculating the gaussian for each frequency data point.
@xValue;
for ($i=0; $i<@freqValue; $i++) {
	for ($j=$x_min; $j<=$x_max; $j++) {
		$x = ($j - $freqValue[$i]);
                $y = ( $x / $sigma)**2;
                $z = exp(-0.5*$y);
                $xValue[$j] += $z;
	}
}

 # print the x_max values and amplitude to STDOUT
print " Frequencies \t Amplitude \n";
for ($i=1; $i<@xValue; $i++) {
    print "\t $i \t $xValue[$i]\n";
}

# print the x_max values and amplitude to spect.dat
open (MYFILE, '>>spect.dat');
print MYFILE " Note: freq.dat has $numimg imaginary frequencies.\n";
print MYFILE " Frequencies \t Amplitude \n";
for ($i=1; $i<@xValue; $i++) {
    print MYFILE "\t $i \t $xValue[$i]\n";
} 
close (MYFILE); 










