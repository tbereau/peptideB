###############################
#                             #
#      Peptide Builder        #
#                             #
#   Author: Tristan Bereau    #
#          (2008)             #
#                             #
###############################

# Contains everything that's relevant to the sidechains (SC).

package require mmsg
package require peptideb::utils
package provide peptideb 1.0.0

namespace eval peptideb {
    
    # ------- Hydrophobicity scale ---------- #
    # Input raw hydrophobicity scale. peptideB will
    # take care of shifting the free energies in order 
    # to mimic hydrophobicity beteen sidechains.
    # Data fitted from MJ matrix to accomodate the
    # Lorentz-Berthelot mixing rule.

    # aqueous environment
    if {![info exists hp_medium] || $hp_medium == 0} {
	# with electrostatics
	if {[info exists charges] && $charges=="on"} {
	    variable dG_GLY 2.15;#2.13;# 0.0
	    variable dG_ALA 2.63;#2.73;#-0.4
	    variable dG_PRO 2.04;#1.92;#-1.0
	    variable dG_GLU 1.55;#1.36;# 0.9
	    variable dG_GLN 1.74;#1.85;# 0.3
	    variable dG_ASP 1.35;#1.40;# 1.1
	    variable dG_ASN 1.54;#1.70;# 0.8
	    variable dG_SER 1.61;#1.75;# 0.1
	    variable dG_HIS 2.81;#2.70;#-0.2
	    variable dG_LYS 0.99;#1.00;# 1.5
	    variable dG_ARG 1.98;#1.87;# 1.5
	    variable dG_THR 2.11;#2.10;#-0.3
	    variable dG_VAL 5.43;#5.40;#-2.4
	    variable dG_ILE 6.79;#6.66;#-1.6
	    variable dG_LEU 7.63;#7.76;#-2.3
	    variable dG_MET 5.59;#5.51;#-1.6
	    variable dG_PHE 7.52;#7.58;#-2.4
	    variable dG_TYR 4.29;#4.31;#-1.3
	    variable dG_CYS 4.66;#4.68;#-2.1
	    variable dG_TRP 5.30;#5.30;#-3.0
	} else {
	    variable dG_GLY 2.13;# 0.0
	    variable dG_ALA 2.73;#-0.4
	    variable dG_PRO 1.92;#-1.0
	    variable dG_GLU 1.36;# 0.9
	    variable dG_GLN 1.85;# 0.3
	    variable dG_ASP 1.40;# 1.1
	    variable dG_ASN 1.70;# 0.8
	    variable dG_SER 1.75;# 0.1
	    variable dG_HIS 2.70;#-0.2
	    variable dG_LYS 1.00;# 1.5
	    variable dG_ARG 1.87;# 1.5
	    variable dG_THR 2.10;#-0.3
	    variable dG_VAL 5.40;#-2.4
	    variable dG_ILE 6.66;#-1.6
	    variable dG_LEU 7.76;#-2.3
	    variable dG_MET 5.51;#-1.6
	    variable dG_PHE 7.58;#-2.4
	    variable dG_TYR 4.31;#-1.3
	    variable dG_CYS 4.68;#-2.1
	    variable dG_TRP 5.30;#-3.0
	}
    } else {
	# hydrophobic environment - Fauchere & Pliska
	# Table found in Protein Physics, Finkelstein.
	# This does not take charges into account
	variable dG_GLY  0.0
	variable dG_ALA -0.4
	variable dG_PRO -1.0
	variable dG_GLU  0.9
	variable dG_GLN  0.3
	variable dG_ASP  1.1
	variable dG_ASN  0.8
	variable dG_SER  0.1
	variable dG_HIS -0.2
	variable dG_LYS  1.5
	variable dG_ARG  1.5
	variable dG_THR -0.3
	variable dG_VAL -2.4
	variable dG_ILE -1.6
	variable dG_LEU -2.3
	variable dG_MET -1.6
	variable dG_PHE -2.4
	variable dG_TYR -1.3
	variable dG_CYS -2.1
	variable dG_TRP -3.0
	if {[info exists charges] && $charges=="on"} {
	    ::mmsg::err [namespace current] "Electrostatics in hydrophobic medium not currently supported."
	}			
    }
    
    # -------- Van der Waals radii ----------- #
    variable rvdw_GLY 0.0;
    variable rvdw_XXX 2.5;

    # -------- Calpha-Cb distance ------------ #
    variable bondCaCb_GLY 1.00; # where the H lies
    variable bondCaCb_XXX 1.53; # non-GLY
    


    # ---------- End of parameter list -------------- #

    # list of all the 3-letter amino acid codes
    set 3letter_list {GLY ALA PRO GLU GLN ASP ASN \
			  SER HIS LYS ARG THR VAL \
			  ILE LEU MET PHE TYR CYS TRP}
    
    variable hp_strength ""
    foreach aa $3letter_list {
	lappend hp_strength "$aa [set dG_$aa]"
    }
    set hp_strength [lsort -real -index 1 $hp_strength]
    set largest [lindex $hp_strength end]
    set smallest [lindex $hp_strength 0]
    
    for {set k 0} {$k < [llength $hp_strength]} {incr k} {
	lset hp_strength $k 1 [expr [lindex [lindex $hp_strength $k] 1] - [lindex $smallest 1]]
	lset hp_strength $k 1 [expr [lindex [lindex $hp_strength $k] 1]\
				   / ([lindex $largest 1] - [lindex $smallest 1])]
    }
    
    unset largest
    unset smallest

    # *** Hydrophobicity scale ***    
    # The depth of the LJ potential for each 3-letter code
    # is saved in the variables 'hp_XXX' where XXX represents
    # a given 3-letter amino acid code.
    # hp_GLY
    # hp_TRP
    # ...

    # $lj_hp is automatically scaled depending on the medium 
    foreach aa $hp_strength { 
	variable    hp_[lindex $aa 0] [format %5.3f [expr [lindex $aa 1]*$lj_hp]] 
    }
}



