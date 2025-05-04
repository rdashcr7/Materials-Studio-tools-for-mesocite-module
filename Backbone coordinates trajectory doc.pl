#!perl

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);
use warnings;  


################################################

my $doc = $Documents{"0.1%.xtd"};

my $setName = "Backbone"; # define name of set

my $object = "Beads"; # Set to Atoms or Beads

################################################

my $setItems = $doc->AsymmetricUnit->Sets("$setName")->$object;

for (my $frame = 4001; $frame <= $doc->Trajectory->NumFrames; ++$frame) {

    $doc->Trajectory->CurrentFrame = $frame;
    #print "Frame $frame\n";


    foreach my $item (@$setItems) {
    
        # Write out X, Y and Z coordinate
        printf ("%0.5f %0.5f %0.5f \n", $item->X, $item->Y, $item->Z);
    }
    
}