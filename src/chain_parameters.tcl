###############################
#                             #
#      Peptide Builder        #
#                             #
#   Author: Tristan Bereau    #
#          (2008)             #
#                             #
###############################

# All the science has been put in this file !
# Contains all geometric parameters.

package require mmsg
package require peptideb::utils
package provide peptideb 1.0.0

namespace eval peptideb {

    # --- List of parameters ------------------------ #
    
    # We always need pi. 
    variable pi 3.14159265358979323846
    
    # Van der Waals radii  (Takada AA - UA)
    variable rvdw_N    1.45;    #1.32 - 1.65
    variable rvdw_Ca   1.85;    #1.48 - 1.85
    variable rvdw_C    1.75;    #1.70 - 2.00
    # Lorentz-Berthelot mixing rules
    variable lb_NCa  [utils::sum $rvdw_N $rvdw_Ca] 
    variable lb_NC   [utils::sum $rvdw_N $rvdw_C ] 
    variable lb_CaC  [utils::sum $rvdw_Ca $rvdw_C ] 
    
    
    # Bonds (in A). Read from N (start of AA) to C (end). Always ordered.
    # covalent bonds
    variable bond_NCa  1.455
    variable bond_CaC  1.510
    variable bond_CN   1.325
    # effective bonds
    variable bond_NCb  2.442
    variable bond_NC   2.444
    variable bond_CaN  2.406
    variable bond_CaCa 3.784
    variable bond_CbC  2.494
    variable bond_CCa  2.432
    # not simulated, only for vmd
    variable bond_CO   1.235
    variable bond_NH   1.000
    
    # Angles (in Degrees)
    variable angleD_NCaCb  108
    variable angleD_CbCaC  113
    variable angleD_NCaC   111
    variable angleD_CaCN   116
    variable angleD_CNCa   122
    variable angleD_CaCO   122
    variable angleD_CaNHN  115 
    # Dihedral
    variable angleD_NCaCO  135
    variable angleD_CCaNHN 120
    # Improper dihedrals - for chirality
    variable angleD_Im_Cb -122
    
    # Hydrogen bonding equilibrium distance
    variable hbond_NC 4.11
    
    # Angles (in Radians)
    variable angleR_NCaCb  [utils::deg2rad $angleD_NCaCb ]
    variable angleR_CbCaC  [utils::deg2rad $angleD_CbCaC ]
    variable angleR_NCaC   [utils::deg2rad $angleD_NCaC  ]
    variable angleR_CaCN   [utils::deg2rad $angleD_CaCN  ]
    variable angleR_CNCa   [utils::deg2rad $angleD_CNCa  ]
    variable angleR_CaCO   [utils::deg2rad $angleD_CaCO  ]
    variable angleR_NCaCO  [utils::deg2rad $angleD_NCaCO ]
    variable angleR_CaNHN  [utils::deg2rad $angleD_CaNHN ]
    variable angleR_CCaNHN [utils::deg2rad $angleD_CCaNHN]
    variable angleR_Im_Cb  [utils::deg2rad $angleD_Im_Cb ]
    
    
    
    # --------------- Interactions --------------------- #
    # Define unit of energy
    variable eps 1.0

    # Rescale all interactions
    variable rescale 1.0;#1.429
    
    # Bonded interactions
    variable k_bond     [expr $rescale*300.0]
    
    # Angles
    variable k_angle    [expr $rescale*300.0]; # From Takada
    
    # Dihedrals
    variable k_dih_omega   [expr $rescale*67.]; # dihedral for \omega - 67 From Takada (one well)
    variable k_dih_w_pro   [expr $rescale* 3.]; # Double well with peak at 6 k_BT
    variable k_dihedral    [expr $rescale* 0.]; # dihedral for \phi and \psi

    # Dipolar interaction
    variable k_dipolar     [expr $rescale*-0.3]; # best value so far is -0.3 (new -0.55)
    
    # Improper dihedral
    variable k_imp_dih     [expr $rescale*17.] ; # From Takada
    
    # Subtract Lennard-Jones for bonded partners
    variable max_length 8.0
    
    # Non-bonded interactions
    # Useful Parameters
    variable lj_eps [expr $rescale*.02*$eps]; # Takada : 0.5 kT = 0.3 kcal/mol
    variable cut_factor 1.1225;#  For 12-6 :   1.1225; # 2^{1/6}
    variable lj_shift 0.25
    variable ljoffset 0.0
    # Special case : interaction of side chains
    variable lj_hp 0.
    # Side chain strength depends on medium
    if {[info exists hp_medium] && $hp_medium==1} {
	# side chain strength in a hydrophobic environment
	variable lj_hp [expr 4.5*$eps]
    } else {
	# side chain strength in a aqueous environment
	variable lj_hp [expr 4.5*$eps]
    }
    variable ljhp_cut 10.
    # side chain coupling. Default is 1.0; 0.0 will turn the interaction off
    variable hp_coupling 1.0
    
    # Hydrogen bonding
    # Careful with the bonded partners !
    variable ljangle_eps 0.
    # HBond strength in the bilayer
    variable ljangle_eps_bilayer [expr 14.*$eps]
    # HBond strength depends on medium
    if {[info exists hp_medium] && $hp_medium==1} {
	# Hydrophobic medium
	variable ljangle_eps $ljangle_eps_bilayer
    } else {
	# aqueous environment
	variable ljangle_eps [expr 6.*$eps]; # Takada uses 4.7 kT = 2.8 kcal/mol
    }

    variable ljangle_cut 8.
    variable HB_max_cap 200.    

    # HBond in bilayer (off by default)
    if {![info exists HB_bilayer_z0]} {
	variable HB_bilayer_z0 0. }
    if {![info exists HB_bilayer_dz]} {
	variable HB_bilayer_dz 0. }
    if {![info exists HB_bilayer_kappa]} {
	variable HB_bilayer_kappa 1. }

    # Electrostatics - Debye-Hueckel
    variable dh_bjerrum   7.; # in units of \sigma
    variable dh_kappa  0.125; # 1/\kappa is 8\sigma
    variable dh_rcut  15.; # in units of \sigma

    # DPD stuff
    variable dpd_r_cut  10.; # as large as the largest cutoff in the force field: 
    # hydrophobic interactions
}

# ------------- Interactions file ------------------- #

# Source the side chain parameters
source [file join [file dirname [info script]] chain_sidechain.tcl   ]
# Source the interactions. 
source [file join [file dirname [info script]] chain_interactions.tcl]
# Source the hfip parameter file
source [file join [file dirname [info script]] hfip_parameters.tcl]
# Source the hfip interaction file
source [file join [file dirname [info script]] hfip_interactions.tcl]
