use strict;
use warnings;

our %configuration;
our %parameters;

sub build_targets
{
	my $thisVariation = $configuration{"variation"} ;

	# cad variation has two volume:
	# target container
	# and inside cell
	if($thisVariation eq "lD2")
	{
		my $nplanes = 4;

		my @oradius  =  (    50.3,   50.3,  21.1,  21.1 );
		# NOMINAL : (  -140.0,  265.0, 280.0, 280.0 );
		my @z_plane  =  (  -240.0,  265.0, 280.0, 280.0 );

		# vacuum target container
		my %detector = init_det();
		$detector{"name"}        = "target";
		$detector{"mother"}      = "root";
		$detector{"description"} = "Target Container";
		$detector{"color"}       = "22ff22";
		$detector{"type"}        = "Polycone";
		my $dimen = "0.0*deg 360*deg $nplanes*counts";
		for(my $i = 0; $i <$nplanes; $i++) {$dimen = $dimen ." 0.0*mm";}
		for(my $i = 0; $i <$nplanes; $i++) {$dimen = $dimen ." $oradius[$i]*mm";}
		for(my $i = 0; $i <$nplanes; $i++) {$dimen = $dimen ." $z_plane[$i]*mm";}
		$detector{"dimensions"}  = $dimen;
		$detector{"material"}    = "G4_Galactic";
		$detector{"style"}       = 0;
		print_det(\%configuration, \%detector);


		$nplanes = 5;
		my @oradiusT  =  (   2.5,  10.3,  7.3, 5.0,  2.5);
		my @z_planeT  =  ( -24.2, -21.2, 22.5, 23.5, 24.5);

		# actual target
		%detector = init_det();
		$detector{"name"}        = "lD2";
		$detector{"mother"}      = "target";
		$detector{"description"} = "Target Cell";
		$detector{"color"}       = "aa0000";
		$detector{"type"}        = "Polycone";
		$dimen = "0.0*deg 360*deg $nplanes*counts";
		for(my $i = 0; $i <$nplanes; $i++) {$dimen = $dimen ." 0.0*mm";}
		for(my $i = 0; $i <$nplanes; $i++) {$dimen = $dimen ." $oradiusT[$i]*mm";}
		for(my $i = 0; $i <$nplanes; $i++) {$dimen = $dimen ." $z_planeT[$i]*mm";}
		$detector{"dimensions"}  = $dimen;
		$detector{"material"}    = "LD2"; # defined in gemc database
		$detector{"style"}       = 1;
		print_det(\%configuration, \%detector);

	}
}

1;






















