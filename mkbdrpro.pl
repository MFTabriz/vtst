#!/usr/bin/env perl
#;-*- Perl -*-

# This program reads in sphpro.F in VASP and bdr_changes, and then combine them to create a new file called bdrpro.F.

# variables needed by the script

@args=@ARGV;
@args<2 || die "usage: mkbdrpro.pl <bdrfile>\n";
$bdrfile="bdr_changes";
if(@args==1){
  $bdrfile=$args[0];
}
$sphfile="sphpro.F";
$outfile="bdrpro.F";

# take sphpro.F and convert it to bdrpro.F
open (IN,$bdrfile) or die "Cannot open $bdrfile\n";
@bdr=<IN>;
close (IN);

open (IN,$sphfile) or die "Cannot open $sphfile\n";
@sph=<IN>;
close (IN);

$len=@sph;
$l2=0;
open(OUT,">$outfile");
for ($i=0;$i<$len;$i++) {
  if($sph[$i] =~ /MODULE msphpro/){
    if($sph[$i+1] =~ /USE/){
      $sph[$i] =~ s/msphpro/mbdrpro/;
      print OUT $sph[$i];
      next;
    }elsif($sph[$i] =~ /END/){
      next;
    }
  }
  if($sph[$i] =~ /SUBROUTINE SPHPRO\( &/) {
    $sph[$i] =~ s/SPHPRO/BDRPRO/;
    print OUT $sph[$i];
    next;
  }
  if($sph[$i] =~ /CHARACTER \(LEN=3\), ALLOCATABLE, SAVE :: LMTABLE\(:, :\)/) {
    print OUT $sph[$i];
    while($bdr[$l2] !~ /111111/) {print OUT $bdr[$l2];$l2=$l2+1;}
    $l2=$l2+1;
    next;
  } 
  if($sph[$i] =~ /REAL\(q\) WORKR\(WDES\%NRPLWV\), WORKI\(WDES\%NRPLWV\)/){
    while($bdr[$l2] !~ /111111/) {print OUT $bdr[$l2];$l2=$l2+1;}
    $l2=$l2+1;
    next;
  }
  if($sph[$i] =~ /#endif/ and $sph[$i-1] =~ /IONODE = WDES\%COMM\%IONODE/){
    print OUT $sph[$i];
    while($bdr[$l2] !~ /111111/) {print OUT $bdr[$l2];$l2=$l2+1;}
    $l2=$l2+1;
    next;
  }
  if($sph[$i] =~ /CALL CREATE_SINGLE_KPOINT_WDES\(WDES,WDES_1K,NK\)/){
    print OUT $sph[$i];
    while($bdr[$l2] !~ /111111/) {print OUT $bdr[$l2];$l2=$l2+1;}
    $l2=$l2+1;
    next;
  }
  if($sph[$i] =~ /CGDR=EXP\(CITPI/){
    print OUT $sph[$i];
    while($bdr[$l2] !~ /111111/) {print OUT $bdr[$l2];$l2=$l2+1;}
    $l2=$l2+1;
    next;
  }
  if($sph[$i] =~ /DO ISPINOR=0,WDES\%NRSPINORS-1/){
    print OUT $sph[$i];
    while($sph[$i]!~ /ENDDO/){$i=$i+1;}
    while($bdr[$l2] !~ /111111/) {print OUT $bdr[$l2];$l2=$l2+1;}
    $l2=$l2+1;
    next;
  }
  if($sph[$i] =~ /SUMI=0/ and $sph[$i-1] =~ /SUMR=0/){
    print OUT $sph[$i];
    while($sph[$i]!~ /CSUM\(LMS/){$i=$i+1;}
    $i=$i-1;
    while($bdr[$l2] !~ /111111/) {print OUT $bdr[$l2];$l2=$l2+1;}
    $l2=$l2+1;
    next;
  }
  if($sph[$i] =~ /CALLMPI\( M_sum_z\( WDES\%COMM_INB, CSUM_PHASE/){
    print OUT $sph[$i];
    while($bdr[$l2] !~ /111111/) {print OUT $bdr[$l2];$l2=$l2+1;}
    $l2=$l2+1;
    next;
  }
  if($sph[$i] =~ /ENDDO band/ and $sph[$i+1] =~ /ENDDO ion/){
    print OUT $sph[$i];
    while($bdr[$l2] !~ /111111/) {print OUT $bdr[$l2];$l2=$l2+1;}
    $l2=$l2+1;
    next;
  }
  if($sph[$i] =~ /ENDDO ion/ and $sph[$i-1] =~ /ENDDO band/){
    print OUT $sph[$i];
    while($bdr[$l2] !~ /111111/) {print OUT $bdr[$l2];$l2=$l2+1;}
    $l2=$l2+1;
    next;
  }
  if($sph[$i] =~ /DEALLOCATE\(PHAS\)/){
    print OUT $sph[$i];
    while($bdr[$l2] !~ /111111/) {print OUT $bdr[$l2];$l2=$l2+1;}
    $l2=$l2+1;
    next;
  }
  if($sph[$i] =~ m/SPHPRO_FAST/){
    while($sph[$i]!~ /END SUBROUTINE/){$i=$i+1;}
    $i=$i+1;
  }
  print OUT $sph[$i];
}
for ($i=$l2;$i<@bdr;$i++) {
  print OUT $bdr[$i];
}
close(OUT);

