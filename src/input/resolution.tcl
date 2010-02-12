# ::peptideb::resolution
#
# Adds details to an existing chain.
#
#

namespace eval peptideb {
    namespace eval input {
	# Adds chemical resolution defined by the parameters
	# in the configuration file.
	proc addresolution {coords} {
	    if {[llength $coords] == 0} {
		::mmsg::err [namespace current] "Add atoms to the chain before improving the resolution!"
	    }
	    if { $peptideb::display_O_H == "on" } {
		set coords [oxygen $coords] 
		set coords [hydrogenN $coords]
	    } 
	    return $coords
	}

	# Adds the oxygens to an existing chain consisting of
	# N-Ca-Cb-C beads of amino acids.
	# Returns the reconstructed chain.
	proc oxygen {coords} {
	    # As a sanity check, all lists should be a multiple of 4.
	    for { set l 0 } { $l < [expr [llength $coords] - $peptideb::hfip_flag] } { incr l } {
		if { [llength [lindex $coords $l]] % 4 != 0 } {
		    ::mmsg::err [namespace current] "Problem in trying to add oxygens: the number of beads is not a multiple of 4 (N-Ca-Cb-C)."
		}
	    }
	    
	    # Initialize the coordinate vector
	    set final_coords ""

	    for { set l 0 } { $l < [expr [llength $coords] - $peptideb::hfip_flag] } { incr l } {
		set peptide_coords [lindex $coords $l]
		set initchainlength [llength $peptide_coords]
		# From now on we assume that the chain is well constructed, 
		# and that it has the right sequence of beads N-Ca-Cb-C. 
		# We want to insert an O atom between C_i and N_{i+1}: index 4 (counting from 0). 
		
		# Here we construct the oxygen atom geometrically. We look at the atoms around and 
		# get the position of the new atom by a simple formula. Can't used that for the last AA (see below).
		for { set k 0 } { $k < [expr $initchainlength / 4] -1 } { incr k } {
		    # Reference point is C' atom
		    # The factor of 5 takes into account the previous Oxygens.
		    set pos [expr $k*5+3]
		    # The next two vectors determine the direction. The third one is just the sum of these two.
		    set r_CC [peptideb::utils::2bead_vector $peptide_coords [expr $pos - 2] $pos]
		    set r_NC [peptideb::utils::2bead_vector $peptide_coords [expr $pos + 1] $pos]
		    set r_direction [peptideb::utils::add_2_vectors $r_CC $r_NC]
		    set r_direction [peptideb::utils::normalize $r_direction]
		    set r_direction [peptideb::utils::scale_vector $peptideb::bond_CO $r_direction]
		    
		    set peptide_coords [linsert $peptide_coords [expr $pos+1] [peptideb::utils::add_2_vectors [lindex $peptide_coords $pos] $r_direction]]
		}
		# Oxygen of the last Amino Acid can't be computed the same way since 
		# we need the position of the following N atom. This method is
		# not correct, but it doesn't really matter for the last O. It's only for VMD.
		
		set k [expr $initchainlength/4 -1]
		# Same game: we transform to local coordinates and come back.
		set phi $peptideb::angleR_NCaCO
		set theta [expr $peptideb::pi - $peptideb::angleR_CaCO]
		set R $peptideb::bond_CO
		
		set O_vector [sphericalcoords $R $theta $phi]
		
		# Define unit vectors for transformation matrix.
		# The factor of 5 takes into account the previous Oxygens.
		set pos [expr $k*5+3]
		
		set bc [peptideb::utils::2bead_vector $peptide_coords [expr $pos - 2]       $pos     ]
		set ab [peptideb::utils::2bead_vector $peptide_coords [expr $pos - 2] [expr $pos - 3]]
		set n [bc_ab_to_n $bc $ab]
		
		set lastatom [lindex $peptide_coords $pos]
		
		set peptide_coords [linsert $peptide_coords [expr $pos+1] [change_frame $O_vector $bc $n $lastatom]]
		
		# At the end, add the whole peptide to the final vector of coordinates
		lappend final_coords $peptide_coords
	    }
	    # Loop over hfip molecules                                                                                                           
	    if { $peptideb::hfip_flag } {
		lappend final_coords [lindex $coords [expr [llength $coords]-1]]
	    }
	    return $final_coords
	}



	# Adds the hydrogens to an existing chain consisting of
	# N-Ca-Cb-C=O beads of amino acids. There *must* be
	# oxygen atoms.
	proc hydrogenN {coords} {
	    # As a sanity check, the list should be a multiple of 4 or 5.
	    for { set l 0 } { $l < [expr [llength $coords] - $peptideb::hfip_flag] } { incr l } {
		if { [llength [lindex $coords $l]] % 5 != 0} {
		    ::mmsg::err [namespace current] "Problem in trying to add hydrogens: the number of beads is not a multiple of 5 (N-Ca-Cb-C=O)."
		}
	    }

	    # Initialize the coordinate vector
	    set final_coords ""
	    
	    for { set l 0 } { $l < [expr [llength $coords] - $peptideb::hfip_flag]} { incr l } {
		set peptide_coords [lindex $coords $l]
		set initchainlength [llength $peptide_coords]
		# From now on we assume that the chain is well constructed, 
		# and that it has the right sequence of beads N-Ca-Cb-C=O. 
		# We want to insert an H atom between N_i and Ca_i: index 1 (counting from 0). 
		
		# Except for the first AA which doesn't have a previous C' atom, we construct
		# H atoms geometrically from a simple formula.
		
		# Same game: we transform to local coordinates and come back.
		set phi $peptideb::angleR_CCaNHN
		set theta [expr $peptideb::pi - $peptideb::angleR_CaNHN]
		set R $peptideb::bond_NH
		
		set H_vector [sphericalcoords $R $theta $phi]
		
		# Define unit vectors for transformation matrix.
		set bc [peptideb::utils::2bead_vector $peptide_coords 1 0]
		set ab [peptideb::utils::2bead_vector $peptide_coords 1 3]
		set n [bc_ab_to_n $bc $ab]
		
		set lastatom [lindex $peptide_coords 0]
		
		set peptide_coords [linsert $peptide_coords 1 [change_frame $H_vector $bc $n $lastatom]]
		
		for { set k 1 } { $k < [expr $initchainlength / 5] } { incr k } {
		    # Reference point is N atom
		    # The factor of 6 takes into account the added Oxygen+Hydrogen.
		    set pos [expr $k*6]
		    # The next two vectors determine the direction. The third one is just the sum of these two.
		    set r_CpN [peptideb::utils::2bead_vector $peptide_coords [expr $pos - 1] $pos]
		    set r_CaN [peptideb::utils::2bead_vector $peptide_coords [expr $pos + 1] $pos]
		    set r_direction [peptideb::utils::add_2_vectors $r_CpN $r_CaN]
		    set r_direction [peptideb::utils::normalize $r_direction]
		    set r_direction [peptideb::utils::scale_vector $peptideb::bond_NH $r_direction]
		    
		    set peptide_coords [linsert $peptide_coords [expr $pos+1] [peptideb::utils::add_2_vectors [lindex $peptide_coords $pos] $r_direction]]
		}
		
		lappend final_coords $peptide_coords
	    }
	    # Loop over hfip molecules                                                                                                           
	    if { $peptideb::hfip_flag } {
		lappend final_coords [lindex $coords [expr [llength $coords]-1]]
	    }
	    return $final_coords
	}


    }
}

