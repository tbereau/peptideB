# ::peptideb::secondary
#
# read the secondary structure of a protein
# from a file
#
#
namespace eval peptideb {
    namespace eval input {
				# ::peptideb::input::read_AA_sequence
				# Read the secondary structure of a protein from
				# a TCL file of the type:
				# 
				# set amino_acids {
				#      {ALA -60 -45}
				#      {MET -59 -45}
				#        [...]
				# }
				#
				# where there can be as many amino acids as we want.
				# The two numbers following the residue name represent
				# the \phi and \psi angles. Obviously the very first
				# and last dihedral angles are useless, but we require
				# them for the sake of uniformity in the structure of
				# the arrays. The 3rd dihedral angle is 180 degrees by
				# default (encountered awfully often).
				# * If the global variable 'peptideb::rand_equilib' is
				# turned on, choose dihedral angles according to 
				# Ramachandran map (files 'rama_XXX.dat', XXX={ala, gly}).
				# This way, dihedrals will be randomly chosen from 
				# an equilibrated basin.
				# Argument: name of the file.
				# Returns an array of coordinates of all the beads.
				#
				proc read_AA_sequence { file } {
						variable coords ""
						variable ::peptideb::rand_equilib

						::mmsg::send [namespace current] "Building peptide from file '$file'..."
						# Read the sequence in a variable $amino_acids.

						if { [llength $peptideb::amino_acids]==0} {
								::mmsg::err [namespace current] "The configuration file does not contain any amino acid!\n"
						}

						# Initialize the nested list of coordinates.
						set final_coords ""

						# Now check that all the amino acids have a proper name. This will avoid bugs in the nesting of the sequence.
						check_aa $peptideb::amino_acids

						# Fill up the $coords array of beads... Do this for each peptide.
						for { set j 0 } { $j < [llength $peptideb::amino_acids] } { incr j } {
								set 3lettercode [lindex [lindex [lindex $peptideb::amino_acids $j] 0] 0]
								# Warning if a peptide contains only one amino acid. 
								# Most probably a mistake of the user in the nesting of the amino_acid sequence.
								if { [llength [lindex $peptideb::amino_acids $j]] == 1 } {
										::mmsg::warn [namespace current] "Peptide number $j has only one amino acid !"
								}
								

								# Set the variable coords to NULL. We grow the peptide and insert it into the nested list of coords.
								set coords ""
								# We start with the 1st AA by hand. We'll generalize
								# for the next ones.
								
								# Start at the origin - it's a N atom! The argument doesn't matter...symmetry!
								add_N_atom 0
								
								# Second bead has rotational symmetry, it's a Calpha.
								add_Ca_atom
								
								# Third bead has 1d rotation freedom, it's a Cbeta. Also append next C
								# last bead of 1st amino acid.
								set iso_1 [isomer $j 0]
								# if we're equilibrating the initial dihedrals, use the proper routine
								if {$rand_equilib == 1} {
										set dihedrals_list [randequilib $3lettercode]
										set phi_1 [lindex $dihedrals_list 0]
										set psi_k [lindex $dihedrals_list 1]
								} else {
										set phi_1 [::peptideb::utils::random_1d]
								}
								add_CbC_atoms $phi_1 $iso_1 $3lettercode
								
								
								
								# From now on it's all automatic, everything looks the same.
								# So we just loop over all remaining amino acids.
								for { set k 1 } { $k < [llength [lindex $peptideb::amino_acids $j]] } { incr k  } {
										set 3lettercode [lindex [lindex [lindex $peptideb::amino_acids $j] $k] 0]
										# We need the psi angle of the previous AA to orient the next N atom.
										# Next are Cbeta and C atoms. Torsion angle is phi of the current AA.
										if {[llength [lindex [lindex $peptideb::amino_acids $j ] [expr $k-1]]] >= 3 } {
												set psi_k [::peptideb::utils::deg2rad [lindex [lindex [lindex $peptideb::amino_acids $j ] [expr $k - 1]] 2]]
												set phi_k [::peptideb::utils::deg2rad [lindex [lindex [lindex $peptideb::amino_acids $j ]       $k     ] 1]]
										} else {
												# Random angles -- either they're equilibrated or uniformaly random
												if {$rand_equilib == 1} {
														set dihedrals_list [randequilib $3lettercode]
														# We need the psi angle of the *previous* AA to orien the next N atom
														set psi_k2 [lindex $dihedrals_list 1]
														set phi_k [lindex $dihedrals_list 0]
												} else {
														# Uniformly distributed random angle (in radians)
														set psi_k [expr 2.*$peptideb::pi*rand()]
														# useless variable psi_k2
														set psi_k2 0.
														set phi_k [expr 2.*$peptideb::pi*rand()]
												}
										}
										add_N_atom $psi_k
										# The next Calpha atom does not need any torsion angle. Assume 180 degrees.
										add_Ca_atom
										
										set iso_k [isomer $j $k]
										add_CbC_atoms $phi_k $iso_k $3lettercode

										# Update variables 
										if {$rand_equilib == 1} {
												set psi_k $psi_k2
										}
								}
								lappend final_coords $coords
						}

						::mmsg::send [namespace current] "Configuration was parsed successfully."
						return $final_coords
				}

				# ::peptideb::input::add_CbC_atoms
				#
				# Computes the coordinates of the Carbon beta + Carbon atoms
				# in the chain. Torsion angle is given in argument.
				# The next C is constrained from the type of isomerism.
				# Arguments: 
				#   - torsion angle phi
				#   - type of isomerism.
				#   - name of the amino acid
				# Appends the coordinates to the array.
				# Note: The function *assumes* the last atom was a
				# Calpha atom. Test on the length of $coords.
				# Determines automatically whether it's the first
				# AA of not.
				#
				# *** IMPORTANT *** : even if the routine refers to carbon-beta,
				# we are now implementing more specificity, and the atom now
				# corresponds to the center of mass of the side chain.
				proc add_CbC_atoms {phi iso name} {
						variable coords

						if { [llength $coords] %4 != 2} { 
								::mmsg::err [namespace current] "Trying to add new Cbeta and C atoms at the wrong place. Check the structure!"
						}

						# Calculate new atom Cb in local coord system
						set proj_angle [projected_angle $peptideb::angleR_CbCaC $peptideb::angleR_NCaCb $peptideb::angleR_NCaC]
						# Isomerism
						set phiCb [phi_isomer $phi $iso $proj_angle]
						set theta [expr $peptideb::pi - $peptideb::angleR_NCaCb]
						if { $name == "GLY" } { set R $peptideb::bondCaCb_GLY } else { set R $peptideb::bondCaCb_XXX }

						set Cb_vector [sphericalcoords $R $theta $phiCb]

						# Define unit vectors for transformation matrix.
						set size [llength $coords]
						set bc [::peptideb::utils::2bead_vector $coords [expr $size - 2] [expr $size - 1]]
						# $size == 2 is a special case: there is still too much freedom
						# to treat it as the general case for the ab vector.
						if { $size == 2 } {
								set ab [::peptideb::utils::randomperpvec $bc]
						} else {
								set ab [::peptideb::utils::2bead_vector $coords [expr $size - 2] [expr $size - 3]]
						}
						set n [bc_ab_to_n $bc $ab]
						set lastatom [lindex $coords end]

						lappend coords [change_frame $Cb_vector $bc $n $lastatom]

						# Calculate new atom C in local coord system
						set theta [expr $peptideb::pi - $peptideb::angleR_NCaC]
						set R $peptideb::bond_CaC
						
						set C_vector [sphericalcoords $R $theta $phi]

						# Note that this $lastatom is really the Calpha by construction.
						# It should NOT be the Cbeta atom.
						lappend coords [change_frame $C_vector $bc $n $lastatom]
				}

				# ::peptideb::input::add_N_atom
				#
				# Computes the coordinates of the next N atom in the chain.
				# Argument: torsion angle \psi.
				# Returns nothing. It appends the new coordinates in the 
				# $coords array.
				# Note: this *assumes* the previous atom was a C atom.
				# Test on the length of the array to make sure, but can't
				# make sure that the last atom was indeed C.
				proc add_N_atom { psi } {
						variable coords

						# Makes sure $coords has a number of points that's a multiple
						# of 4 (N, Calpha, Cbeta, C) and different than 0.
						if { [llength $coords] % 4 != 0 } {
								::mmsg::err [namespace current] "Trying to add a new N atom at the wrong place. Check the structure!"
						} elseif { [llength $coords] > 0 } {
								# If we're adding a N atom that's NOT on the first AA.
								
								# Calculate new atom N in local coord system
								set phi $psi
								set theta [expr $peptideb::pi - $peptideb::angleR_CaCN]
								set R $peptideb::bond_CN
								
								set N_vector [sphericalcoords $R $theta $phi]

								# Define unit vectors for transformation matrix.
								set size [llength $coords]
								set bc [::peptideb::utils::2bead_vector $coords [expr $size - 3] [expr $size - 1]]
								set ab [::peptideb::utils::2bead_vector $coords [expr $size - 3] [expr $size - 4]]
								set n [bc_ab_to_n $bc $ab]
								
								set lastatom [lindex $coords end]

								lappend coords [change_frame $N_vector $bc $n $lastatom]
						} else {
								# If we're just starting the chain, put something random, centered in the middle.
								set x [expr [lindex $peptideb::setbox_l 0]*(0.25+0.5*rand())]
								set y [expr [lindex $peptideb::setbox_l 1]*(0.25+0.5*rand())]
								set z [expr [lindex $peptideb::setbox_l 2]*(0.25+0.5*rand())]
								lappend coords "$x $y $z"
						}
				}

				# ::peptideb::input::add_Ca_atom
				#
				# Computes the coordinates of the next Ca atom in the chain.
				# Argument: none.
				# Returns nothing. It appends the new coordinates in the 
				# $coords array.
				# Note: this *assumes* the previous atom was a N atom.
				# Test on the length of the array to make sure, but can't
				# make sure that the last atom was indeed N.
				proc add_Ca_atom { } {
						variable coords

						# Makes sure $coords has a number of points that's a multiple
						# of 4 % 1 (N, Calpha, Cbeta, C, N) and different than 1.
						if { [llength $coords] % 4 != 1 } {
								::mmsg::err [namespace current] "Trying to add a new Ca atom at the wrong place. Check the structure!"
						}
						
						if { [llength $coords] > 1 } {
								# Calculate new atom Ca in local coord system
								# torsion angle is always PI.
								set phi $peptideb::pi
								set theta [expr $peptideb::pi - $peptideb::angleR_CNCa]
								set R $peptideb::bond_NCa
								
								set Ca_vector [sphericalcoords $R $theta $phi]

								# Define unit vectors for transformation matrix.
								set size [llength $coords]
								set bc [::peptideb::utils::2bead_vector $coords [expr $size - 2] [expr $size - 1]]
								set ab [::peptideb::utils::2bead_vector $coords [expr $size - 2] [expr $size - 4]]
								set n [bc_ab_to_n $bc $ab]

								set lastatom [lindex $coords end]
								
								lappend coords [change_frame $Ca_vector $bc $n $lastatom]
						} else {
								# Deal with first AA.
								set temp [::peptideb::utils::random_2d $peptideb::bond_NCa]
								set origin [lindex $coords 0]
								set x [expr [lindex $temp 0] + [lindex $origin 0]]
								set y [expr [lindex $temp 1] + [lindex $origin 1]]
								set z [expr [lindex $temp 2] + [lindex $origin 2]]
								lappend coords "$x $y $z"
						}
				}

				# ::peptideb::input::change_frame
				#
				# Transforms the local coordinate of the new atom(s)
				# to its old frame.
				# Arguments: 
				#    - D_loc is the coordinates of the new atom(s)
				# in the local frame.
				# Must be of the form:
				#    { x1 y1 z1 } { x2 y2 z2 } ...
				# If only one atom is passed in the argument, don't forget
				# extra braces "{" and "}"!
				#    - $bc is the unit vector
				#    - $n is the unit vector
				#    - $lastatom is the position of the parent atom
				# that we use for positioning the new one.
				# returns the coordinates of the new atom in the general frame.
				#
				proc change_frame { D_loc bc n lastatom} {
						set bc [::peptideb::utils::normalize $bc]
						set nxbc [::peptideb::utils::crossproduct $n $bc]
						set nxbc [::peptideb::utils::normalize $nxbc]
						set result [applymatrix $lastatom $D_loc $nxbc $n $bc]
						return $result
				}

				# ::peptideb::input::randequilib
				#
				# Choose random dihedral according to one of the two ramachandran
				# maps ('rama_XXX.dat', XXX={ala, gly}). The maps are unnormalized
				# probability distribution function. The most probable state has a
				# probability of 1. Even though these maps were calculated at T=1,
				# the same map will be used for any temperature.
				# The algorithm works the following way:
				# - choose two random dihedral angles (phi, psi)
				# - lookup the probability P of this pair from one of the 'rama_XXX' files
				# - accept this set of dihedrals if a uniformly distributed random number 
				#   r satisfies the following condition:  r < P
				# Convert dihedrals to radians
				# Arguments: name of the amino acid
				proc randequilib { name } {
						# dummy variable to control the while loop
						set accept_dihedral_set 0
						::mmsg::send [namespace current] "Randomly choose equilibrated angle for amino acid $name..." 
						while {$accept_dihedral_set == 0} {

								
								# choose two random numbers between -180 and +180
								# Values in the file are in degrees
								set phi_trial [expr 360.*rand() - 180.]
								set psi_trial [expr 360.*rand() - 180.]
								
								set filename ""
								set path [namespace eval :: {set env(PEPTIDEB_DIR)}]
								
								if {$name == "GLY"} {
										set filename "$path/src/input/rama_gly.dat"
								} else {
										set filename "$path/src/input/rama_ala.dat"
								}
								set read_file [open $filename r]
								set data [read $read_file]
								set data [split $data "\n"]

								# values of phi,psi that are closest to values in the file
								# Start with ridiculous values
								set phi_bestfit -1000
								set psi_bestfit -1000

								# Try to find the closest values of {phi,psi}_trial
								foreach line $data {
										if {[string range [lindex $line 0] 0 0]!="#"} {
												# phi is the 0th element in the lin
												set phi_line [lindex $line 0]
												if {[expr abs($phi_line-$phi_trial)] < [expr abs($phi_bestfit-$phi_trial)]} {
														set phi_bestfit $phi_line
												}
												# psi is the 1st element in the line
												set psi_line [lindex $line 1]								
												if {[expr abs($psi_line-$psi_trial)] < [expr abs($psi_bestfit-$psi_trial)]} {
														set psi_bestfit $psi_line
												}
										}
								}

								# Now read the file a second time and look for $phi_bestfit and $psi_bestfit
								set probability 0.
								foreach line $data {
										if {[lindex $line 0] == $phi_bestfit && [lindex $line 1] == $psi_bestfit} {
												set probability [lindex $line 2]
										}
								}
								
								set random_number [expr rand()]
								if {$random_number < $probability} {
										# Accept the pair phi,psi.
										set accept_dihedral_set 1
										::mmsg::send [namespace current] "Pair of dihedrals (phi,psi)=($phi_bestfit,$psi_bestfit) was accepted."
								}
						}

						return [list [::peptideb::utils::deg2rad $phi_bestfit] [::peptideb::utils::deg2rad $psi_bestfit]]

				}
    }
}
