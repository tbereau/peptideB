# Apr. 16, 2009
#
# This script is intended to be used with peptideB.
# It calculates the radius of gyration of existing PDB files
# in the current folder.
# 
# Author : Tristan Bereau (bereau@cmu.edu)
#
#
#

# Description of the script :
# It reads in PDB files, feeds the configuration in ESPResSo,
# and output the radius of gyration.
#
#
#



# ----------------------------------- Script -------------------------------------- #

# Turn the chemical resolution 'on'. This will help computing scalar products.
set ::display_O_H "on"

# Output file
set output_file "gyration2.dat"

set files [glob *.pdb]
set files [lsort $files]

set channel [open $output_file w]

foreach pdb_file $files {

	set coords [::peptideb::input::import_pdb $pdb_file]
	::peptideb::espresso::create_topology $coords 0

	analyze set chains 0 1 [expr [llength [lindex $coords 0]]/4]
	set rad_g [analyze rg 0 1 [expr [llength [lindex $coords 0]]/4]]
	
	puts $channel [lindex $rad_g 2]
	
}





