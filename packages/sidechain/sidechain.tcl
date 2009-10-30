# Oct. 30, 2009
#
# This script is intended to be used with peptideB.
# It calculates the side chain energy of existing PDB files
# in the current folder.
# 
# Author : Tristan Bereau (bereau@cmu.edu)
#
#
#

# Description of the script :
# It reads in PDB files, feeds the configuration in ESPResSo,
# and output the side chain energy
#
#
#



# ----------------------------------- Script -------------------------------------- #

# Turn the chemical resolution 'on'. This will help computing scalar products.
set ::display_O_H "on"

# Output file
set output_file "energy_sidechain.dat"

set files [glob *.pdb]
set files [lsort $files]

set channel [open $output_file w]

foreach pdb_file $files {

    set coords [::peptideb::input::import_pdb $pdb_file]
    ::peptideb::espresso::create_topology $coords 0
    
    # Calculate total energy due to side chain interactions                                                                      
    set sc_energy 0.
    for {set i 10} {$i<30} {incr i} {
	for {set j $i} {$j<30} {incr j} {
	    set sc_energy [expr $sc_energy + [analyze energy nonbonded $i $j]]
	}
    }
    puts $channel $sc_energy
    
}





