#!/usr/bin/perl

$init = 2;
$hmax = 10;

@Gs	       = (0.01, 0.02, 0.03, 0.05, 0.08, 0.1, 0.2);
$Ng 		= scalar(@Gs);

 system("g77 -w -x f77-cpp-input -o 2dcollapse 2dcollapse.f specfun.f");
 
for($gcount=0; $gcount < $Ng; $gcount++)
{

$g = $Gs[$gcount];

for($hh=1; $hh <= $hmax; $hh++)
{

 $xvname     = sprintf("outxv_%i_g%1.2f_n%i.dat", $init, $g, $hh);
 $energyname     = sprintf("outenergy_%i_g%1.2f_n%i.dat", $init, $g, $hh);

 open(INFILE,">in.dat");
 printf INFILE "0.450 $g 0.1 0.01 $init\n 207. 13 0.2\n";
 close INFILE;
 
 system("./2dcollapse");
 system("mv outxv.dat $xvname");
 system("mv outenergy.dat $energyname");
  };
};