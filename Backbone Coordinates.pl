#!perl

use strict;
use warnings;
use Getopt::Long;
use MaterialsScript qw(:all);

################################################

# Specify the full path where the file should be saved
# Update the path to a directory of your choice, for example, your desktop
my $file_path = 'C:\Users\rdash\Documents\Materials Studio Projects\Martini vs my Method\Blend\Charged\All mapping\Modified knuckle epitopes_Files\Documents\KP 586\0.05%\1000 ps\Backbone Coordinates.csv';

# Open the CSV file for writing
open(my $fh, '>', $file_path) or die "Could not open file '$file_path' $!";

# Write headers
print $fh "X,Y,Z\n";

my $doc = $Documents{"0.05%.xtd"};

my $setName = "Backbone"; # define name of set

my $object = "Beads"; # Set to Atoms or Beads

################################################

my $setItems = $doc->AsymmetricUnit->Sets("$setName")->$object;

for (my $frame = 1; $frame <= $doc->Trajectory->NumFrames; ++$frame) {

    $doc->Trajectory->CurrentFrame = $frame;
    #print "Frame $frame\n";

    foreach my $item (@$setItems) {
    
        # Extract X, Y and Z coordinates
        my $x = $item->X;
        my $y = $item->Y;
        my $z = $item->Z;
        
        # Write coordinates to CSV file
        print $fh "$x,$y,$z\n";
    }
}

# Close the CSV file
close $fh;

print "Coordinates successfully written to $file_path\n";
