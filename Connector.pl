#!perl
# This code can be used to create a connector between two beads 'CA' and 'CB'. This can be useful to create ringed or closed structure within your mesomolecules.
# Go to edit sets to chose these two specific beads and group them under a new set and name the set 'S'.
# insert the name of the mesomolecule in line 11.
# After creating the connector, you can rename the beads CA and CB according to your needs. Make sure to declare the forcefield type for the beads after renaming them, not doing so will result in failed jobs. 

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

my $doc=$Documents{"mesomolecule.xsd"}; 
my $beads=$doc->UnitCell->Sets("S")->Beads;
my $bead1 =$doc->UnitCell->Beads("CA");
my $bead2 =$doc->UnitCell->Beads("CB");
$doc->CreateBeadConnector($bead1,$bead2);
