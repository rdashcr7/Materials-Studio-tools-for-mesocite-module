#!perl

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

my $doc=$Documents{"BMP-2 KP with CI.xsd"};
my $beads=$doc->UnitCell->Sets("S")->Beads;
my $connectors = $doc->UnitCell->BeadConnectors;
foreach my $connector (@$connectors) {
	my $type1 = $connector->Bead1->BeadTypeName;
    	my $type2 = $connector->Bead2->BeadTypeName;
    	$connector->Delete if ($type1 eq 'A' || $type2 eq 'A');
}