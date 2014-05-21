#!/usr/bin/env perl
#;-*- Perl -*-

print "\n";

@ARGV>0 || die "usage: dosanalyze.pl [w=number widths at quarter-height to include] [e=emin,emax] <s,p,d,a (all)> <atom num(s)>\n";

# Determine how many Widths at Half Height to include, default 2.5
$arg1 = @ARGV[0];
$e_range = grep(/e=[-+]?\d+/, $arg1);
$w = grep(/w=\d+/, $arg1);
#print $e_range."|".$w."\n";
$Emin = $Emax = 0;

if(@ARGV>0){
    $arg1 = @ARGV[0];
    if($arg1 =~ /^\d+$/ || $arg1 =~ /^(\d+)(-)(\d+)$/) {
    } elsif($w) {
        $NumStdDevs = shift(@ARGV);
        $NumStdDevs =~ s/^(w=)(0*)(\d+)(\.?)(\d*)$/$2$3$4$5/;
    } elsif($e_range) {
        $string = shift(@ARGV);
        @dummy = split(/=/, $string);
        @amy = split(/,/, $dummy[-1]);
        $Emin = $amy[0];
        $Emax = $amy[1];
#    print $Emin."\n"."$Emax\n";
    }
    if($NumStdDevs =~ m/^\D/) { die "Error in usage syntax!\n"; }
    if($e_range) {
        print "integrate from ".$Emin."eV to ".$Emax."eV.\n";
        if($Emin>$Emax) {
            die "Emin is larger than Emax!\n";
        }
    } 
}

# If no band is selected, analyze d band
$oflag = 'd';
if(@ARGV>0) {
    $arg2 = @ARGV[0];
    if($arg2 =~ /^\d+$/ || $arg2 =~ /^(\d+)(-)(\d+)$/) {
    } else {
        $oflag = shift(@ARGV)
    }
    if($oflag =~ m/^(s|p|d|a)$/i) {
    } else {
         die "Error in usage syntax!\n";
    }
}

# Determine whether file is spin polarized or not and denote the number of columns
if(@ARGV>0 && @ARGV[0]==0) {
    $DOSfile = "DOS0";
} else {
    $DOSfile = "DOS1";
}
open (DOS,$DOSfile);
<DOS>;
$in = <DOS>;
$in =~ s/^\s+//g;
@line = split(/\s+/,$in);
$spin = @line;

# Set selected atom(s)
if(@ARGV[0] =~ /^\d+$/){
    $NumAtm = @ARGV;
    @Atm = @ARGV;
} elsif(@ARGV[0] =~ /^(\d+)(-)(\d+)$/) {
    $NumAtm = $3-$1+1;
    for($i=0; $i<$NumAtm; $i++) {
        $Atm[$i] = $i+$1;
    }
}

# Find how many atom DOS there are
opendir MAINDIR, "." or die "can't open this dir!" ;
@DOSFILE = grep /^DOS\d+$/, readdir MAINDIR;
$NumDOS = @DOSFILE;
closedir MAINDIR;

# If no atom is selected, analyze all of them
if($NumAtm == 0) {
    $NumDOS = $NumDOS-1;
} else {
    $NumDOS = $NumAtm;
}

# Inform user if usage arguments not totally filled
if($arg1 =~ /^\d+$/ || $arg1 =~ m/^(s|p|d|a)$/i) {
    print "No energy range specified, integration will cover the whole range.\n";
}
if($arg2 =~ /^\d+$/) {
   print "When no band is selected, will analyze d band\n";
}


# Determine which subroutine to execute based on band
$maxDOSj = 0; # For determination of where DOS max occurs
$maxDOS = 0;
#print "spin:".$spin."\n";
if($spin == 4 || $spin == 7) { npNA(); }
if($spin == 3 || $spin == 5) { npA0(); }

# Subroutine for non-polarized analyzing a specific band
sub npNA{

    if($spin == 4) {
    # If no orbital flag, plot the d-band
        if($oflag =~ /^s$/i) { $col = 1; }
        if($oflag =~ /^p$/i) { $col = 2; }
        if($oflag =~ /^d$/i) { $col = 3; }
        if($oflag =~ /^a$/i) { $col = 4; }
    }

    if($spin == 7) {
    # If no orbital flag, analyze the d-band
        if($oflag =~ /^s$/i) { $colup = 1; $coldown = 2; }
        if($oflag =~ /^p$/i) { $colup = 3; $coldown = 4; }
        if($oflag =~ /^d$/i) { $colup = 5; $coldown = 6; }
        if($oflag =~ /^a$/i) { $colup = 7 }; # analyze all bands
    }

# Read selected DOS files
    $first = 1;
    for ($i=1; $i<=$NumDOS; $i++) {
        if($NumAtm == 0){
            $DOSFILE = "DOS"."$i";
        } else {
            $a = $Atm[$i-1];
            $DOSFILE = "DOS"."$a";
        }
        $j = 0;
        open (DOS,$DOSFILE);
        <DOS>;
        while ($in = <DOS>) {
            $in =~ s/^\s+//g;
            @line = split(/\s+/,$in);
            $eneval = $line[0];
            if($spin == 4) {
                if($col < 4) {
                    $dosval = $line[$col];
                } else {
                    $dosval = $line[1]+$line[2]+$line[3];
                }
            }
            if($spin == 7) {
                if($colup < 7) {
                    $dosval = $line[$colup]-$line[$coldown];
                } else {
                    $dosval = $line[1]-$line[2]+$line[3]-$line[4]+$line[5]-$line[6];
                }
            }
            if($first){
                $Ene[$j] = $eneval;
                $Dos[$j] = $dosval;
                if($dosval > $maxDOS) {
                    $maxDOS = $dosval;
                    $maxDOSj = $j;
                }
            } else {
                $Ene[$j] == $eneval || die "Energy $j in file $DOSFILE does not match first file.\n";
                $Dos[$j] += $dosval;
                if($Dos[$j] > $maxDOS){
                    $maxDOS = $Dos[$j];
                    $maxDOSj = $j;
                }
            }
            $j++;
        }
        close(DOS);
        if($first) {
            $numene = $j;
            $first = 0;
        } else {
            $numene == $j || die "Number of energy values in $DOSFILE does not match first file.\n";
        }
    }
# End of subroutine np
}


# Subroutine for non-polarized analyzing all DOS for all atoms (DOS0 file)
sub npA0 {

    print "WARNING: For DOS0, ALL bands are analyzed\n";
    $NumDOS = 1;
    $j = 0;
    open (DOS,"DOS0");
    while ($in = <DOS>) {
        $in =~ s/^\s+//g;
        @line = split(/\s+/,$in);
        $eneval = $line[0];
        if($spin == 3){
            $dosval = $line[1];
        }
        if($spin == 5){
            $dosval = $line[1] - $line[2];
        }
        $Ene[$j] = $eneval;
        $Dos[$j] = $dosval;
        if($dosval > $maxDOS){
            $maxDOS = $dosval;
            $maxDOSj = $j;
        }
        $j++;
    }
    close(DOS);
    $numene = $j;
# End of subroutine npA0
}


if($w){
    # Find half width at half height and corresponding indices
    $cut = 0.25*$maxDOS;
    $DOS50up = $maxDOSj;
    $DOS50down = $maxDOSj;
    while($Dos[$DOS50up] >= $cut) {
        $DOS50up++;
    }
    while($Dos[$DOS50down] >= $cut) {
        $DOS50down--;
    }
    $WHH = abs($Ene[$DOS50up] - $Ene[$DOS50down]);
    $range = $WHH*$NumStdDevs;
    $halfrange = $range/2;
    $EneUp = $Ene[$maxDOSj] + $halfrange;
    $EneDown = $Ene[$maxDOSj] - $halfrange;
    $k = $maxDOSj;
    while($k<$j && $Ene[$k]<=$EneUp){
        $numenemax = $k;
        $k++;
    }
    $h = $maxDOSj;
    while($h >= 0 && $Ene[$h] >= $EneDown) {
        $numenemin = $h;
        $h--;
    }
}

if(!$e_range) {
    $Emin = $Ene[0];
    $Emax = $Ene[$numene-1];
}

# Calculate the average energy state
$dossum = $dstsum = $dsqsum = 0;
if($w) {
    for($i=$numenemin; $i<=$numenemax; $i++) {
        $dossum += $Dos[$i];
        $dstsum += $Ene[$i]*$Dos[$i];
        $dsqsum += $Ene[$i]*$Ene[$i]*$Dos[$i]*$Dos[$i];
        $dssqsum += $Ene[$i]*$Ene[$i]*$Dos[$i];
    }
} else {
    for($i=0; $i<$numene; $i++) {
        if($Ene[$i] >= $Emin and $Ene[$i] <= $Emax) {
            $dossum += $Dos[$i];
            $dstsum += $Ene[$i]*$Dos[$i];
            $dssqsum += $Ene[$i]*$Ene[$i]*$Dos[$i];
        }
    }
}

$dossum || die "Total DOS is zero.\n";
$eneavg = $dstsum/$dossum;
$enevar = ($dssqsum*$dossum - $dstsum*$dstsum)/($dossum*$dossum);
$enestd = sqrt($enevar);

if($w || $e_range) {
    # Verify that center is not shifted preposterously far (more than 15%) otherwise warn
    for($i=0; $i<$j; $i++) {
        $testdossum += $Dos[$i];
        $testenesum += $Ene[$i]*$Dos[$i];
    }
    $testeneavg = $testenesum/$testdossum;
    $diff = $testeneavg - $eneavg;
    $norm = $diff/$eneavg;
    $absnorm = abs($norm);
    if($absnorm > 0.1){
        print "WARNING: band center has shifted by more than 10% compared to full integrated range\n";
    }
}

# Display output
print "\n";
print "Total States: $dossum\n"; 
print "Average Energy (band center): $eneavg\n";
if($w) {
    print "\n";
    print "Width at quarter-height: $WHH\n";
    print "Lower energy cutoff: $Ene[$numenemin]\n";
    print "Upper energy cutoff: $Ene[$numenemax]\n";
}
print "Standard Deviation: $enestd\n";
print "\n";
