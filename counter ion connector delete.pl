#!perl
# This script can be used to delete a connector (or a bond) between any two beads in a mesostructure file (.xsd).
# Make sure to go to edit sets and chose all beads which can possibly contain the bond or connector of interest. You may chose a mesomoelcule or you can chose all the beads in the entire mesostructure. After selecting, create a new set 'S' for the selected beads.

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

my $doc=$Documents{"Mesostructure.xsd"};
my $beads=$doc->UnitCell->Sets("S")->Beads;
my $connectors = $doc->UnitCell->BeadConnectors;
foreach my $connector (@$connectors) {
	my $type1 = $connector->Bead1->BeadTypeName;
    	my $type2 = $connector->Bead2->BeadTypeName;
    	$connector->Delete if ($type1 eq 'CA' || $type2 eq 'CB');
}
