#!perl
# Include the name of the trajectory document in line number 17. 
# Go to edit set and chose the group of beads of interest (press alt and click on a bead), define new set 'S' for the chosen beads.
# Enter the path where you wish to save the csv file containing the coordinates of the selected beads in each frame of the trajectory in line 11.

use strict;
use warnings;
use Getopt::Long;
use MaterialsScript qw(:all);
 
my $file_path = 'Coordinates.csv';

open(my $fh, '>', $file_path) or die "Could not open file '$file_path' $!";
 
print $fh "X,Y,Z\n";

my $doc = $Documents{"0.05%.xtd"};

my $setName = "S"; 

my $object = "Beads"; 

my $setItems = $doc->AsymmetricUnit->Sets("$setName")->$object;

for (my $frame = 1; $frame <= $doc->Trajectory->NumFrames; ++$frame) {

    $doc->Trajectory->CurrentFrame = $frame;
    #print "Frame $frame\n";

    foreach my $item (@$setItems) {

        my $x = $item->X;
        my $y = $item->Y;
        my $z = $item->Z;
        
        print $fh "$x,$y,$z\n";
    }
}

close $fh;

print "Coordinates file path: $file_path\n";
