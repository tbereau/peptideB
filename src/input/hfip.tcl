# peptideb::input --
#
# Add HFIP to the simulation box
# Author: Tristan
#


namespace eval peptideb {
    namespace eval input {
	# Given a volume fraction of HFIP frac, determine the number of HFIP
	# molecules to add. Return the number of HFIP molecules to be added.
	proc calculate_HFIP_vv_frac {frac} {
	    # Change that!
	    return 2
	}

	# Add hfip molecules randomly in the simulation box. Return a list of
	# coordinates. Argument: v/v fraction, in order to determine the
	# number of molecules
	proc create_HFIP_molecules {frac} {
	    set hfip_coords ""
	    set index_mol 0
	    while {$index_mol < $peptideb::hfip_num_mol} {
		set hfip_mol_coords ""
		# set F1 bead randomly
		set F1x [expr [lindex $peptideb::setbox_l 0]*rand()]
		set F1y [expr [lindex $peptideb::setbox_l 1]*rand()]
		set F1z [expr [lindex $peptideb::setbox_l 2]*rand()]
		lappend hfip_mol_coords "$F1x $F1y $F1z"

		# Set Cc bead
		set Cc_bead [::peptideb::utils::random_2d $peptideb::ifp_bond_CF]
		set Ccx [expr $F1x+[lindex $Cc_bead 0]]
		set Ccy [expr $F1y+[lindex $Cc_bead 1]]
		set Ccz [expr $F1z+[lindex $Cc_bead 2]]
		lappend hfip_mol_coords "$Ccx $Ccy $Ccz"

		# Set F2 bead
		set bc [::peptideb::utils::2bead_vector $hfip_mol_coords 0 1]
		set ab [::peptideb::utils::randomperpvec $bc]
		set theta [expr $peptideb::pi - $peptideb::ifp_angleR_FCF]
		set R $peptideb::ifp_bond_CF
		set CF_vector [sphericalcoords $R $theta $peptideb::pi]
		lappend hfip_mol_coords [change_frame $CF_vector $bc \
					     [bc_ab_to_n $bc $ab] [list $Ccx $Ccy $Ccz]]

		lappend hfip_coords $hfip_mol_coords
		incr index_mol
	    }
	    return $hfip_coords
	}
    }
}

