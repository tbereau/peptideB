#::peptideb::misc
#
# Various routines used in the 'input' namespace.
#
#

namespace eval peptideb {
    namespace eval input {
	# Takes in argument one length and two angles
	# and return an array of array of spherical coordinates
	# Args: 
	#  - bond length R
	#  - polar angle theta
	#  - azimuthal angle phi
	# Returns "{x y z}"
	proc sphericalcoords { R theta phi } {
	    set vector [expr $R*sin($theta)*cos($phi)]
	    lappend vector [expr $R*sin($theta)*sin($phi)]
	    lappend vector [expr $R*cos($theta)]

	    return "$vector"
	}


	# Returns the isomer from the amino_acids list
	# Arguments: peptide number (position), and number on the list (index starts from 0!!)
	# Returns a letter L or D
	proc isomer { peptide position } {
	    if { [llength $peptideb::amino_acids] <= $peptide } {
		if { [llength [lindex $peptideb::amino_acids $peptide]] <= $position } {
		    ::mmsg::err [namespace current] "Trying to access an amino acid that does not exist."
		}
	    }
	    if { [llength [lindex [lindex $peptideb::amino_acids $peptide] $position]] == 5 } {
		# case {MET 60. 40. 180. D}
		return [lindex [lindex [lindex $peptideb::amino_acids $peptide] $position] 4]
	    } elseif { [llength [lindex [lindex $peptideb::amino_acids $peptide] $position]] == 4 } {
		# case {MET 60. 40. D}
		if { [lindex [lindex [lindex $peptideb::amino_acids $peptide] $position] 3] == "D" } {
		    return D
		} elseif { [lindex [lindex [lindex $peptideb::amino_acids $peptide] $position] 3] == "L" } {
		    return L
		} else {
		    return L
		}
	    } else {
		if { [lindex [lindex [lindex $peptideb::amino_acids $peptide] $position] 1] == "D" } {
		    return D
		} else {
		    return L
		}
	    }
	}

	# Determines the angle phi from the original
	# dihedral angle $phi and the isomer type $iso.
	# Returns the new phi angle taking into account 
	# the projected angle $proj_a.
	# I hope the sign is right...
	proc phi_isomer { phi iso proj_a } {
	    set result 0
	    if { $iso == "L" || $iso == "l" } {
		set result [expr $phi - $proj_a]
	    } elseif { $iso == "D" || $iso == "d" } {
		set result [expr $phi + $proj_a]
	    } else {
		::mmsg::err [namespace current] "There was an error in trying to read the isomer parameter of one of the amino acids."
	    }

	    return $result
	}

	# Check that every single amino acid of a sequence has a proper name.
	# Will make the scripting a little stronger.
	# Argument : the whole sequence of amino acids
	# Returns: nothing. It just exits the program with an error if something's wrong.
	proc check_aa { aa } {
	    for { set l 0 } { $l < [llength $aa] } { incr l } {
		for { set k 0 } { $k < [llength [lindex $aa $l]] } { incr k } {
		    set res [lindex [lindex [lindex [lindex $aa $l]] $k] 0]
		    if { $res!="ASP" && $res!="GLU" && $res!="ARG" && $res!="LYS" && $res!="HIS" \
			     && $res!="ASN" && $res!="GLN" && $res!="SER" && $res!="THR" && $res!="TYR" \
			     && $res!="ALA" && $res!="GLY" && $res!="VAL" && $res!="LEU" && $res!="ILE" \
			     && $res!="PRO" && $res!="PHE" && $res!="MET" && $res!="TRP" && $res!="CYS"} {
			mmsg::err [namespace current] "Can't read amino acid : $res. Check the sequence and give a 3-letter code name. Is the sequence properly nested ?\n"
		    }
		}
	    }
	}

	# Projects an angle between two vectors on the xy plane.
	# We assume the first vector is lying in the xz plane.
	# Arguments:
	#   - angle between the two vectors alpha
	#   - polar angle of vector 1: delta
	#   - polar angle of vector 2: beta
	# returns an angle gamma in radians
	proc projected_angle { alpha delta beta } {
	    set ca [expr cos($alpha)]
	    set cd [expr cos($delta)]
	    set cb [expr cos($beta)]
	    set sd [expr sin($delta)]
	    set sb [expr sin($beta)]
	    return [expr acos(($ca-$cd*$cb)/($sd*$sb))]
	}

	# Calculate the change of frame matrix.
	# Takes in argument $lastcoords of the last atom, $local: 
	# the vector which we're working on, the $nxbc unit vector, 
	# and the $n and $bc unit vectors.
	# Returns the $resultcoords in the general frame.
	proc applymatrix { lastcoords local nxbc n bc } {
	    set coordx [lindex $lastcoords 0]
	    set coordy [lindex $lastcoords 1]
	    set coordz [lindex $lastcoords 2]

	    set coordx [expr $coordx + [lindex $nxbc 0] * [lindex $local 0] + [lindex $n 0] * [lindex $local 1] + [lindex $bc 0] * [lindex $local 2]]
	    set coordy [expr $coordy + [lindex $nxbc 1] * [lindex $local 0] + [lindex $n 1] * [lindex $local 1] + [lindex $bc 1] * [lindex $local 2]]
	    set coordz [expr $coordz + [lindex $nxbc 2] * [lindex $local 0] + [lindex $n 2] * [lindex $local 1] + [lindex $bc 2] * [lindex $local 2]]

	    return "[expr $coordx] [expr $coordy] [expr $coordz]"
	}

	# Takes two vectors bc and ab in argument and return
	# the unit vector n
	proc bc_ab_to_n { bc ab } {
	    # normalize
	    set bc [::peptideb::utils::normalize $bc]
	    set ab [::peptideb::utils::normalize $ab]
	    # The \hat{n} vector is formed by the cross product of the \hat{bc}
	    # vector with the previous 2 beads that are in the x-z plane
	    # (see doc for more information).
	    set n [::peptideb::utils::crossproduct $bc $ab]
	    set n [::peptideb::utils::normalize $n]

	    return $n
	}

    }
}
