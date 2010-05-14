# Parameter file for peptideb::espresso.
#
#
# List of parameters to run a simulation
#
# Author: Tristan
#
#



# --------------- Protein parameters --------------------- #


set amino_acids {
    {GLY ALA GLY}
}

# Model oxygens and N-hydrogens - "on" or "off"
set display_O_H "on"


# Directory and filenames for output
set directory "tripeptide"
set PDB_file "$directory"
# Starting point for the numbering of PDB files (usually 0)
set filenumber 0



# -------------- ESPResSo Simulation ----------------------- #
# Activate ESPResSo simulation. Turn "on" to activate.
set espresso_switch "on"

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

# The number of steps to integrate with capped forces
set warm_steps 100
# The number of times to call warmup
set warm_n_times 20
# The number of steps to integrate with each call to integrate
set int_steps 10
# The number of times to call integrate
set int_n_times 5000

set replica_temps { 1.0 1.2 1.4 1.6 }
# number of times to integrate between each replica MC step
set replica_timestep 100
# number of replica exchange rounds
set replica_rounds 5000
# frequency at which pdb files are written
set pdb_freq 100




