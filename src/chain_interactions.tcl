###############################
#                             #
#      Peptide Builder        #
#                             #
#   Author: Tristan Bereau    #
#          (2008)             #
#                             #
###############################

# Contains all interaction parameters.


package require mmsg
package require peptideb::utils
package provide peptideb 1.0.0

namespace eval peptideb {
    variable bonded_parms    ""
    variable angle_parms     ""
    variable dihedral_parms  ""
    variable nb_interactions ""

    # --------------- Interactions --------------------- #
    # Bonded interactions

    # Pairs
    # bond 0 is between N and Ca
    lappend bonded_parms [list 0 harmonic $k_bond $bond_NCa]
    # bonds 15 for GLY and 16 for non-GLY - between Ca and Cb
    lappend bonded_parms [list 15 harmonic $k_bond $bondCaCb_GLY]
    lappend bonded_parms [list 16 harmonic $k_bond $bondCaCb_XXX]
    # bond 2 is between Ca and C
    lappend bonded_parms [list 2 harmonic $k_bond $bond_CaC]
    # bond 3 is between C and N
    lappend bonded_parms [list 3 harmonic $k_bond $bond_CN]
    
    
    # Angles
    
    # Angle  N - Ca - Cb
    lappend angle_parms [list 4 angle $k_angle $angleR_NCaCb]
    # Angle  N - Ca - C
    lappend angle_parms [list 5 angle $k_angle $angleR_NCaC ]
    # Angle Cb - Ca - C
    lappend angle_parms [list 6 angle $k_angle $angleR_CbCaC]
    # Angle Ca - C - N
    lappend angle_parms [list 7 angle $k_angle $angleR_CaCN ]
    # Angle C - N - Ca
    lappend angle_parms [list 8 angle $k_angle $angleR_CNCa ]
    
    # Dihedrals
    
    # psi and phi now include the dipolar interaction, that biases
    # the potential for beta-sheets rather than alpha helices.
    # Dihedral \psi and \phi 
    #lappend dihedral_parms [list 9  dihedral 3 $k_dihedral  $pi]
    # Dihedral \omega
    lappend dihedral_parms [list 10 dihedral 1 $k_dih_omega $pi]
	# Dihedral \omega for proline
	lappend dihedral_parms [list 14 dihedral 2 $k_dih_w_pro  0.]
    # Dipolar interactions
    lappend dihedral_parms [list 11 dihedral 1 $k_dipolar    0.]
    # Improper dihedral - for chirality
    lappend dihedral_parms [list 12 dihedral 1 $k_imp_dih   $angleR_Im_Cb]


    # Subtract Lennard-Jones for bonded partners    
    #lappend bonded_parms [list 13 subt_lj 0 $max_length]
    
    
    
    # Non-bonded interactions
    
    # --- Matrix of LJ interactions --- #
    # Type     Atom
    #  0          N
    #  1         Ca
    #  2          C
    #  3  Proline N (no Hbond)
    #  4 terminal-N
    #  5 terminal-C
    #  10-29     Sc
    
    # ** Backbone beads only ** 
    # interaction N  and N
    lappend nb_interactions \
	[list 0 0 lennard-jones $lj_eps [expr 2.*$rvdw_N]  [expr $cut_factor*2.*$rvdw_N ] $lj_shift $ljoffset]
    # interaction N and Ca
    lappend nb_interactions \
	[list 0 1 lennard-jones $lj_eps $lb_NCa  [expr $cut_factor*$lb_NCa ] $lj_shift $ljoffset]
    # interaction N and C
    lappend nb_interactions \
	[list 0 2 lennard-jones $lj_eps $lb_NC   [expr $cut_factor*$lb_NC  ] $lj_shift $ljoffset]
    # interaction N  and Pro-N
    lappend nb_interactions \
	[list 0 3 lennard-jones $lj_eps [expr 2.*$rvdw_N]  [expr $cut_factor*2.*$rvdw_N ] $lj_shift $ljoffset]
    # interaction Ca and Ca
    lappend nb_interactions \
	[list 1 1 lennard-jones $lj_eps [expr 2.*$rvdw_Ca] [expr $cut_factor*2.*$rvdw_Ca] $lj_shift $ljoffset]
    # interaction Ca and C
    lappend nb_interactions \
	[list 1 2 lennard-jones $lj_eps $lb_CaC  [expr $cut_factor*$lb_CaC ] $lj_shift $ljoffset]
    # interaction C  and C
    lappend nb_interactions \
	[list 2 2 lennard-jones $lj_eps [expr 2.*$rvdw_C]  [expr $cut_factor*2.*$rvdw_C ] $lj_shift $ljoffset]
    # interaction Pro-N  and Pro-N
    lappend nb_interactions \
	[list 3 3 lennard-jones $lj_eps [expr 2.*$rvdw_N]  [expr $cut_factor*2.*$rvdw_N ] $lj_shift $ljoffset]
    # interaction Pro-N and Ca
    lappend nb_interactions \
	[list 3 1 lennard-jones $lj_eps $lb_NCa  [expr $cut_factor*$lb_NCa ] $lj_shift $ljoffset]
    # interaction Pro-N and C
    lappend nb_interactions \
	[list 3 2 lennard-jones $lj_eps $lb_NC   [expr $cut_factor*$lb_NC  ] $lj_shift $ljoffset]
    # interaction term-N and N
    lappend nb_interactions \
	[list 4 0 lennard-jones $lj_eps [expr 2.*$rvdw_N]  [expr $cut_factor*2.*$rvdw_N ] $lj_shift $ljoffset]
    # interaction term-N and Ca
    lappend nb_interactions \
        [list 4 1 lennard-jones $lj_eps $lb_NCa  [expr $cut_factor*$lb_NCa ] $lj_shift $ljoffset]
    # interaction term-N and C
    lappend nb_interactions \
        [list 4 2 lennard-jones $lj_eps $lb_NC   [expr $cut_factor*$lb_NC  ] $lj_shift $ljoffset]
    # interaction term-N and Pro-N
    lappend nb_interactions \
        [list 4 3 lennard-jones $lj_eps [expr 2.*$rvdw_N]  [expr $cut_factor*2.*$rvdw_N ] $lj_shift $ljoffset]
    # interaction term-N and term-C
    lappend nb_interactions \
        [list 4 5 lennard-jones $lj_eps $lb_NC   [expr $cut_factor*$lb_NC  ] $lj_shift $ljoffset]
    # interaction term-C and N
    lappend nb_interactions \
        [list 0 5 lennard-jones $lj_eps $lb_NC   [expr $cut_factor*$lb_NC  ] $lj_shift $ljoffset]
    # interaction term-C and Ca
    lappend nb_interactions \
        [list 1 5 lennard-jones $lj_eps $lb_CaC  [expr $cut_factor*$lb_CaC ] $lj_shift $ljoffset]
    # interaction term-C and C
    lappend nb_interactions \
        [list 2 5 lennard-jones $lj_eps [expr 2.*$rvdw_C]  [expr $cut_factor*2.*$rvdw_C ] $lj_shift $ljoffset]
    # interaction term-C and Pro-N
    lappend nb_interactions \
        [list 3 5 lennard-jones $lj_eps $lb_NC   [expr $cut_factor*$lb_NC  ] $lj_shift $ljoffset]
 


    # ** Now with sidechains **
    # Sc beads are id'ed from 10 to 29.
    set index 10
    foreach name $3letter_list {
	# No interaction with GLY
	if {$name != "GLY"} {
	    # The next line is there only in case we decide to turn on the nb 
	    # interaction with GLY...
	    if { $name == "GLY" } { set rvdw_SC $rvdw_GLY } else { set rvdw_SC $rvdw_XXX }
		
	    # interaction N and Sc
	    set sigma [utils::sum $rvdw_N $rvdw_SC]
	    lappend nb_interactions \
			[list 0 $index lennard-jones $lj_eps $sigma [expr $cut_factor *$sigma] $lj_shift $ljoffset]
	    # interaction Ca and Sc
	    set sigma [utils::sum $rvdw_Ca $rvdw_SC]
	    lappend nb_interactions \
			[list 1 $index lennard-jones $lj_eps $sigma [expr $cut_factor*$sigma] $lj_shift $ljoffset]
	    # interaction C and Sc
	    set sigma [utils::sum $rvdw_C $rvdw_SC]
	    lappend nb_interactions \
			[list 2 $index lennard-jones $lj_eps $sigma [expr $cut_factor*$sigma] $lj_shift $ljoffset]
	    # interaction Pro-N and Sc
	    set sigma [utils::sum $rvdw_N $rvdw_SC]
	    lappend nb_interactions \
			[list 3 $index lennard-jones $lj_eps $sigma [expr $cut_factor *$sigma] $lj_shift $ljoffset]
	    # interaction term-N and Sc
            set sigma [utils::sum $rvdw_N $rvdw_SC]
            lappend nb_interactions \
		[list 4 $index lennard-jones $lj_eps $sigma [expr $cut_factor *$sigma] $lj_shift $ljoffset]
	    # interaction term-C and Sc
	    set sigma [utils::sum $rvdw_C $rvdw_SC]
            lappend nb_interactions \
		[list 5 $index lennard-jones $lj_eps $sigma [expr $cut_factor*$sigma] $lj_shift $ljoffset]
	    # interaction Sc and Sc - Matrix of coefficients
	    set index2 10
	    foreach partner $3letter_list {
			if {$partner != "GLY"} {
				# The following condition will never be true because of the if statement two levels above.
				if { $name == "GLY" } { set rvdw_SC2 $rvdw_GLY } else { set rvdw_SC2 $rvdw_XXX }
				set sigma [utils::sum $rvdw_SC $rvdw_SC2]
				# epsilon is the geometric mean of the two side chain energies
				set epsilon [utils::g_mean [set hp_$name] [set hp_$partner]]
				# Multiply epsilon by the side chain coupling
				set epsilon [expr $epsilon * $hp_coupling]
				# WCA shift
				set eps_rel [expr .25*(1-$epsilon/$lj_hp)]

				# WCA - repulsive part (us lj-gen because we can't have lennard-jones twice!)
				lappend nb_interactions \
					[list $index $index2 lj-gen $lj_hp $sigma [expr $cut_factor*$sigma] $eps_rel $ljoffset 12 6 1.0 1.0]

				# LJ - attractive part
				lappend nb_interactions \
					[list $index $index2 lennard-jones $epsilon $sigma $ljhp_cut \
					     [calc_lj_shift $sigma $ljhp_cut] $ljoffset 0. [expr $cut_factor*$sigma]]
			}
			incr index2
	    }
	}
	incr index
    }
    unset index
    unset sigma
    unset epsilon


    # --------------------------------- #

    # Hydrogen bonding
    # Careful with the bonded partners !  
    if {$HB_bilayer_dz > 0.} {
	lappend nb_interactions [list 0 2 lj-angle $ljangle_eps \
				     $hbond_NC $ljangle_cut 1 -1 1 \
				     -2 0 $HB_bilayer_z0 $HB_bilayer_dz \
				     $HB_bilayer_kappa $ljangle_eps_bilayer]
    } else {
	lappend nb_interactions [list 0 2 lj-angle $ljangle_eps \
				     $hbond_NC $ljangle_cut 1 -1 1 -2]
    }
    
    # Self interaction - hat potential
    if {$hat_potential > 0. && $HB_bilayer_dz > 0.} {
	lappend bonded_parms \
	    [list 17 hat $hat_potential \
		 $HB_bilayer_z0 $HB_bilayer_dz $HB_bilayer_kappa]
    }	
    

    # Electrostatics - Debye-Hueckel potential
    # Turn on if it has been set in the parameter file
    if { [info exists charges] && $charges=="on"} {
	::mmsg::send [namespace current] "Charges are turned on."
	lappend nb_interactions [list coulomb $dh_bjerrum dh $dh_kappa $dh_rcut]
    }
    
}
