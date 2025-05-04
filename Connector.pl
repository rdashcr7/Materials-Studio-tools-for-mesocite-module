#!perl
# This code can be used to create a connector between two beads 'CA' and 'CB'. This can be useful to create ringed or closed structure within your mesomolecules.
# Go to edit sets to chose these two specific beads and group them under a new set and name the set 'S'.
# insert the name of the mesomolecule in line 10.

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

my $doc=$Documents{"mesomolecule.xsd"}; 
my $beads=$doc->UnitCell->Sets("S")->Beads;
my $bead1 =$doc->UnitCell->Beads("CA");
my $bead2 =$doc->UnitCell->Beads("CB");
$doc->CreateBeadConnector($bead1,$bead2);
