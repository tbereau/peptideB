# Parameter file for peptideb::espresso.
#
#
# List of parameters to run a simulation
#
# Polyalanine in a fake bilayer (HBond strength changes wrt z-position).
#
#
# Author: Tristan
#
#



# --------------- Protein parameters --------------------- #

# Model poly-peptide
set amino_acids {
    {ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA}

}

# Choose dihedrals randomly from an equilibrated basin
set rand_equilib 0

# Model oxygens and N-hydrogens - "on" or "off"
set display_O_H "on"


# Directory and filenames for output
set directory "bilayer"
set PDB_file "$directory"
# Starting point for the numbering of PDB files (usually 0)
set filenumber 0



# -------------- ESPResSo Simulation ----------------------- #
# Activate ESPResSo simulation. Turn "on" to activate.
set espresso "on"

# Box length. Make sure it's sufficiently big to
# contain all the beads.
set setbox_l {100.0 100.0 100.0}

# Warmup parameters
set warm_time_step 0.01
# Integration parameters
set main_time_step 0.01
set verlet_skin 0.4
set langevin_gamma 1.0
# Temperature
set systemtemp [expr 1.*$eps]

# Charges (on/off) - requires Espresso feature ELECTROSTATICS
set charges "off"

# The number of steps to integrate with capped forces
set warm_steps 50
# The number of times to call warmup
set warm_n_times 20
# The number of steps to integrate with each call to integrate
set int_steps 1000
# The number of times to call integrate
set int_n_times 1000

# DPD simulation: 0 for off, 1 for on
set dpd 0


################
# Fake Bilayer #
################

variable HB_bilayer_z0 [expr [lindex $setbox_l 2]/2.]
variable HB_bilayer_dz 40.
variable HB_bilayer_kappa 1.


###############
# Hat potential strength
variable hat_potential 1.


