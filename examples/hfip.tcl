# Parameter file for peptideb::espresso.
#
#
# List of parameters to run a simulation
#
# Author: Tristan
#
#



# --------------- Protein parameters --------------------- #

# Model poly-peptide
set amino_acids {
    {ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA}
    {ALA ALA ALA}

}

# Choose dihedrals randomly from an equilibrated basin
set rand_equilib 0

# Model oxygens and N-hydrogens - "on" or "off"
set display_O_H "on"


# Directory and filenames for output
set directory "hfip"
set PDB_file "$directory"
# Starting point for the numbering of PDB files (usually 0)
set filenumber 0


# HFIP - v/v fraction (number between 0 and 1)
set hfip_frac 0.2


# -------------- ESPResSo Simulation ----------------------- #
# Activate ESPResSo simulation. Turn "on" to activate.
set espresso "on"

# Box length. Make sure it's sufficiently big to
# contain all the beads.
set setbox_l {40.0 40.0 40.0}

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
set int_steps 500
# The number of times to call integrate
set int_n_times 1000

# DPD simulation: 0 for off, 1 for on
set dpd 0

######################
# Parallel tempering #
######################

# Set the replica exchange temperatures
# i.e. 200 239 286 342 409 489 585 700 K
#set replica_temps {
#  0.67 0.80 1.00 1.15 1.37 1.64 1.96 2.35  
#}
set replica_temps { 1.0 1.1 }
# number of times to integrate between each replica MC step
set replica_timestep 100
# number of replica exchange rounds
set replica_rounds 10000
# frequency at which pdb files are written
set replica_pdb_freq 10
 

################################
# Hamiltonian replica exchange #
################################

# side chain couplings: 
# * 1.0 is original force field
# * 0.0 is no side chain interaction
set hremd_couplings { 1.0 0.5 }
# number of times to integrate between each replica MC step
set hremd_timestep 10000
# number of replica exchange rounds
set hremd_rounds 10000
# frequency at which pdb files are written
set hremd_pdb_freq 10



