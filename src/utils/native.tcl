# ::peptideb::utils - native
#
# Calculate the nativeness order parameter
#


namespace eval peptideb { 
    namespace eval utils {
	# Calculate the set of r_{ij} from a given list
	# of coordinates. This will calculate all nonlocal
	# pairs of C_alpha atoms that are at least 3 amino 
	# acids apart. 
	# Arguement : coords - list of coordinates. Assume
	# there is no O or H atoms.
	# Return a list containing all the r_{ij}.
	proc native_rij {coords} {
	    set final_list ""
	    # Supports multiple peptides !
	    for { set l 0 } { $l < [llength $coords] } { incr l } {
		set chain [lindex $coords $l]
		set list_rij ""
		# Assume coords doesn't have any O or H. It should
		# then be a multiple of 4 (N-Ca-Sc-C').
		if {[llength $chain] % 4 != 0} {
		    mmsg::err [namespace current] "\[Native\] The list of coordinates does not have the right length (=[llength $chain])."
		}
		for {set i 0} {$i < [expr [llength $chain]/4]} {incr i} {
		    for {set j [expr $i + 3]} {$j < [expr [llength $chain]/4]} {incr j} {
			lappend list_rij [distance [lindex $chain [expr $i*4+1]] [lindex $chain [expr $j*4+1]]]
		    }
		}
		lappend final_list $list_rij
	    }
	    return $final_list
	}


	# Calculate the nativeness order parameter.
	# Arguments : two lists of coordinates - one
	# is the integrated structure, the other is the
	# native one.
	proc calculate_Q {coords1 coords2} {
	    set parameter 0.
	    set index 0
	    if {[llength $coords1] != [llength $coords2]} {
		mmsg::err [namespace current] "\[Native\] The two structures don't have the same number of peptides !"
	    }
	    set number_peptides [llength $coords1]
	    for {set l 0} {$l < $number_peptides} {incr l} {
		set chain1 [lindex $coords1 $l]
		set chain2 [lindex $coords2 $l]
		if {[llength $chain1] != [llength $chain2]} {
		    mmsg::err [namespace current] "\[Native\] The two peptides don't have the same length !"
		}
		foreach r1 $chain1 r2 $chain2 {
		    set parameter [expr $parameter + exp(-($r1-$r2)*($r1-$r2)/9.)]
		    incr index
		}
	    }
	    set parameter [expr $parameter/(1.*$index)]
	    return $parameter
	}
    }
}

