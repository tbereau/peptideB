###############################
#                             #
#      Peptide Builder        #
#                             #
#   Author: Tristan Bereau    #
#          (2010)             #
#                             #
###############################

# Contains all interaction parameters for HFIP.
# Bond IDs start at 30.

package require mmsg
package require peptideb::utils
package provide peptideb 1.0.0

namespace eval peptideb {
    variable bonded_parms   
    variable angle_parms    
    variable nb_interactions

    # --------------- Interactions --------------------- #
    # Bonded interactions    
    # bond 30 is between C and F
    lappend bonded_parms [list 30 harmonic $ifp_k_bond $ifp_bond_CF]

    # Angles
    # Bond 31 Angle  F - C - F
    lappend angle_parms [list 31 angle $ifp_k_angle $ifp_angleR_FCF]
    
    # Non-bonded interactions
    # --- Matrix of LJ interactions --- #
    # Type     Atom
    #  30          C
    #  31          F
    
    # ** HFIP-HFIP ** 
    # interaction C  and C
    lappend nb_interactions \
	[list 30 30 lennard-jones $ifp_epsC [expr 2.*$ifp_sigmaC]  [expr $cut_factor*2.*$ifp_sigmaC ] $lj_shift $ljoffset]
    # interaction F  and F
    lappend nb_interactions \
	[list 31 31 lennard-jones $ifp_epsF [expr 2.*$ifp_sigmaF]  [expr $cut_factor*2.*$ifp_sigmaF ] $lj_shift $ljoffset]
    # interaction C  and F
    lappend nb_interactions \
	[list 30 31 lennard-jones $ifp_lb_e_CF $ifp_lb_s_CF   [expr $cut_factor*$ifp_lb_s_CF  ] $lj_shift $ljoffset]

    # ** HFIP-peptide **   
    # interaction N  and ifpC
    lappend nb_interactions \
	[list 0 30 lennard-jones $lj_eps \
	     [utils::sum $rvdw_N $ifp_sigmaC] [expr $cut_factor*[utils::sum $rvdw_N $ifp_sigmaC]] $lj_shift $ljoffset]
    # interaction Ca and ifpC
    lappend nb_interactions \
	[list 1 30 lennard-jones $lj_eps \
	     [utils::sum $rvdw_Ca $ifp_sigmaC] [expr $cut_factor*[utils::sum $rvdw_Ca $ifp_sigmaC]] $lj_shift $ljoffset]
    # interaction C  and ifpC
    lappend nb_interactions \
	[list 2 30 lennard-jones $lj_eps \
	     [utils::sum $rvdw_C $ifp_sigmaC] [expr $cut_factor*[utils::sum $rvdw_C $ifp_sigmaC]] $lj_shift $ljoffset]
    # interaction Pro-N and ifpC
    lappend nb_interactions \
	[list 3 30 lennard-jones $lj_eps \
	     [utils::sum $rvdw_N $ifp_sigmaC] [expr $cut_factor*[utils::sum $rvdw_N $ifp_sigmaC]] $lj_shift $ljoffset]
    # interaction N  and ifpF
    lappend nb_interactions \
	[list 0 31 lennard-jones $lj_eps \
	     [utils::sum $rvdw_N $ifp_sigmaF] [expr $cut_factor*[utils::sum $rvdw_N $ifp_sigmaF]] $lj_shift $ljoffset]
    # interaction Ca and ifpF
    lappend nb_interactions \
	[list 1 31 lennard-jones $lj_eps \
	     [utils::sum $rvdw_Ca $ifp_sigmaF] [expr $cut_factor*[utils::sum $rvdw_Ca $ifp_sigmaF]] $lj_shift $ljoffset]
    # interaction C  and ifpF
    lappend nb_interactions \
	[list 2 31 lennard-jones $lj_eps \
	     [utils::sum $rvdw_C $ifp_sigmaF] [expr $cut_factor*[utils::sum $rvdw_C $ifp_sigmaF]] $lj_shift $ljoffset]
    # interaction Pro-N and ifpF
    lappend nb_interactions \
	[list 3 31 lennard-jones $lj_eps \
	     [utils::sum $rvdw_N $ifp_sigmaF] [expr $cut_factor*[utils::sum $rvdw_N $ifp_sigmaF]] $lj_shift $ljoffset]

    # ** Now with sidechains **
    # Sc beads are id'ed from 10 to 29.
    set index 10
    foreach name $3letter_list {
	# No interaction with GLY
	if {$name != "GLY"} {
	    # The next line is there only in case we decide to turn on the nb 
	    # interaction with GLY...
	    if { $name == "GLY" } { set rvdw_SC $rvdw_GLY } else { set rvdw_SC $rvdw_XXX }
		
	    # interaction ifpC and Sc
	    set sigma [utils::sum $ifp_sigmaC $rvdw_SC]
	    lappend nb_interactions \
			[list 30 $index lennard-jones $lj_eps $sigma [expr $cut_factor *$sigma] $lj_shift $ljoffset]

	    # interaction ifpF and Sc
	    set sigma [utils::sum $ifp_sigmaF $rvdw_SC]
	    set epsilon [utils::g_mean $ifp_epsF [set hp_$name]]
	    # WCA shift
	    set eps_rel [expr .25*(1-$epsilon/$lj_hp)]
	    # WCA - repulsive part (us lj-gen because we can't have lennard-jones twice!)
	    lappend nb_interactions \
		[list 31 $index  lj-gen $lj_hp $sigma [expr $cut_factor*$sigma] $eps_rel $ljoffset 12 6 1.0 1.0]
	    
	    # LJ - attractive part
	    lappend nb_interactions \
		[list 31 $index lennard-jones $epsilon $sigma $ljhp_cut 0.00 $ljoffset 0. [expr $cut_factor*$sigma]]
	}
	incr index
    }
    unset index
    unset sigma
    unset epsilon
}
