###############################
#                             #
#      Peptide Builder        #
#                             #
#   Author: Tristan Bereau    #
#          (2010)             #
#                             #
###############################

# HFIP parameters

package require mmsg
package require peptideb::utils
package provide peptideb 1.0.0

namespace eval peptideb {
    
    # --- List of parameters ------------------------ #

    # Bead sizes
    # For homogeneity, provide radius rather than diameter
    variable ifp_sigmaC  [expr 3.0/2.]; 
    variable ifp_sigmaF  [expr 5.0/2.];

    # Interaction strengths
    variable ifp_epsC    0.01;
    variable ifp_epsF    0.55;

    # Lorentz-Berthelot mixing rules
    variable ifp_lb_s_CF  [utils::sum $ifp_sigmaC $ifp_sigmaF];
    variable ifp_lb_e_CF  [utils::g_mean $ifp_epsC $ifp_epsF];

    # covalent bonds
    variable ifp_bond_CF  1.53;
    
    # Angles (in Degrees)
    variable ifp_angleD_FCF  111.;
    
    # Angles (in Radians)
    variable ifp_angleR_FCF  [utils::deg2rad $ifp_angleD_FCF ]
        
    # --------------- Interactions --------------------- #
    # Bonded interactions
    variable ifp_k_bond     300.0;
    
    # Angles
    variable ifp_k_angle    300.0;
    
}

