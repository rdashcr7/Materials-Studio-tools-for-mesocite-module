#!perl

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

my $doc=$Documents{"KIPKASSVDTELSDISTLHH.xsd"};
my $beads=$doc->UnitCell->Sets("S")->Beads;
my $bead1 =$doc->UnitCell->Beads("CA");
my $bead2 =$doc->UnitCell->Beads("CB");
$doc->CreateBeadConnector($bead1,$bead2);