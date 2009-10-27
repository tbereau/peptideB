# ::peptideb::utils - math.tcl
#
# Various math subroutines.
#
#

namespace eval peptideb {
    namespace eval utils {


	# ::peptideb::utils::deg2rad
	#
	# Converts an angle in degrees to radians.
	# Argument: the angle in degrees.
	# Returns an angle in radians
	#
	proc deg2rad {angleD} {
	    return [expr $angleD / 180. * $peptideb::pi]
	}



	# ::peptideb::utils::random_2d
	# 
	# Returns a random orientation on a 2D-sphere.
	# Argument: bond length.
	# Returns an array of coordinates {x,y,z}.
	# 
	#
	proc random_2d {b_length} {
	    # Main Algorithm to compute a random orientation on
	    # a 2d sphere.
	    set ran1 0.
	    set ran2 0.
	    set ransq 2.
	    while { $ransq>=1. } {
		set ran1 [expr 1.-2.*rand()]
		set ran2 [expr 1.-2.*rand()]
		set ransq [expr $ran1*$ran1+$ran2*$ran2]
	    }
	    set ranh [expr 2.*sqrt(1.-$ransq)]

	    # Store orientation in cartesian coordinates {x,y,z}
	    set x [expr $ran1*$ranh*$b_length]
	    set y [expr $ran2*$ranh*$b_length]
	    set z [expr (1.-2.*$ransq)*$b_length]

	    set coords $x
	    lappend coords $y
	    lappend coords $z
	    return $coords
	}


	#::peptideb::utils::random_1d 
	#
	# Returns a random orientation on a ring
	# Argument: none.
	# Returns an angle in radians. 
	#
	proc random_1d { } {
	    return [expr 2.*$peptideb::pi*rand()]
	}
	
	
	
	
	# ::peptideb::utils::dotproduct
	#
	# Returns the scalar product between two vectors.
	# Arguments: two arrays of length 3.
	# Returns a scalar.
	#
	proc dotproduct {A B} {
	    if { [llength $A] != 3 || [llength $B] != 3 } {
		::mmsg::err [namespace current] "Scalar product must be between two 3D vectors."
	    }
	    set result 0.
	    for { set k 0 } { $k<3 } { incr k } {
		set result [expr $result + [lindex $A $k] * [lindex $B $k]]
	    }
	    return $result
	}
	
	
	
	# ::peptideb::utils::crossproduct
	#
	# Returns the vector product between two vectors.
	# Arguments: two arrays of length 3.
	# Returns an array of coordinates {x y z}.
	#
	proc crossproduct {A B} {
	    if { [llength $A] != 3 || [llength $B] != 3 } {
		::mmsg::err [namespace current] "Vector product must be between two 3D vectors."
	    }
	    
	    # x coordinate of the resulting vector
	    set result 0.
	    set result [expr [lindex $A 1] * [lindex $B 2] - [lindex $A 2] * [lindex $B 1]]
	    # then append y and z components
	    lappend result [expr [lindex $A 2] * [lindex $B 0] - [lindex $A 0] * [lindex $B 2]]
	    lappend result [expr [lindex $A 0] * [lindex $B 1] - [lindex $A 1] * [lindex $B 0]]
	    return $result
	}
	
	# ::peptideb::utils::angle
	#
	# Returns the angle between two vectors.
	# Arguments: two arrays of length 3.
	# Returns an angle in *radians*.
	#
	proc angle {A B} {
	    if { [llength $A] != 3 || [llength $B] != 3 } {
		::mmsg::err [namespace current] "Angle must be calculated between two 3D vectors."
	    }

	    set result [expr [dotproduct $A $B]/(1.*[norm $A]*[norm $B])]
	    # Deal with numerical errors.
	    if { $result > 1.} {
		set result 1.
	    }
	    if { $result < -1.} {
		set result -1.
	    }
	    return [expr acos($result)]
	    
	}
	
	# ::peptideb::utils::norm
	#
	# Returns the norm of a vector.
	# Argument: one array of length 3.
	# Returns a scalar.
	#
	proc norm {A} {
	    # A norm is the square root of the scalar product
	    # of a vector with itself.
	    set result [expr sqrt([dotproduct $A $A])]
	    if { $result == 0} {
		::mmsg::warn [namespace current] "Found a vector of norm 0!"
	    }
	    return $result
	}
	
	
	
	
	# ::peptideb::utils::2bead_vector
	#
	# From a list of coordinates, returns the vector
	# formed by the difference in position between
	# bead1 and bead2 (i.e. vector_12= b2 - b1)
	# Arguments:
	#   - coords is an array of coordinates
	#   - bead1 is the position of bead 1 in $coords
	#   - bead2 is the position of bead 2 in $coords
	# Returns an array of 3 coordinates: {x y z}
	#
	proc 2bead_vector {coords bead1 bead2} {
	    set vector1 [lindex $coords $bead1]
	    set vector2 [lindex $coords $bead2]
	    
	    set result [expr [lindex $vector2 0] - [lindex $vector1 0]]
	    lappend result [expr [lindex $vector2 1] - [lindex $vector1 1]]
	    lappend result [expr [lindex $vector2 2] - [lindex $vector1 2]]
	    
	    return $result
	}
	
	# ::peptideb::utils::add_2_vectors
	#
	# Adds the 2 vectors A and B.
	# Arguments: 2 vectors A and B.
	# Returns an array of 3 coordinates: {x y z}
	#
	proc add_2_vectors {A B} {
	    set result [expr [lindex $A 0] + [lindex $B 0]]
	    lappend result [expr [lindex $A 1] + [lindex $B 1]]
	    lappend result [expr [lindex $A 2] + [lindex $B 2]]    
	    
	    return $result
	}
	
	
	# ::peptideb::utils::diff_2_vectors
	#
	# Subtracts one vector from the other. 
	# Arguments: 2 vectors A and B.
	# Returns array of coordinates of vector A-B.
	proc diff_2_vectors {A B} {
	    set result [expr [lindex $A 0] - [lindex $B 0]]
	    lappend result [expr [lindex $A 1] - [lindex $B 1]]
	    lappend result [expr [lindex $A 2] - [lindex $B 2]]    
	    
	    return $result
	}
	
	# ::peptideb::utils::scale_vector
	#
	# Multiplies a vector A by a scalar k 
	# Arguments : vector A and scalar k
	# Returns the new vector.
	proc scale_vector {k A} {
	    lset A 0 [expr [lindex $A 0] * $k]
	    lset A 1 [expr [lindex $A 1] * $k]
	    lset A 2 [expr [lindex $A 2] * $k]
	    
	    return $A
	}
	
	# ::peptideb::utils::normalize
	#
	# Returns a normlized vector.
	# Argument: array of 3 coordinates.
	# Returns an array of 3 coordinates.
	#
	proc normalize {A} {
	    set norm [norm $A]
	    set result [expr [lindex $A 0]/(1.*$norm)]
	    lappend result [expr [lindex $A 1]/(1.*$norm)]
	    lappend result [expr [lindex $A 2]/(1.*$norm)]
	    
	    return $result
	}
	
	# ::peptideb::utils::randomperpvec
	#
	# Returns an orthogonal, arbitrary vector
	# from another vector. 
	# Argument: one vector
	# Returns another vector.
	#
	proc randomperpvec {A} {
	    # makes sure $A is normalized
	    set A [normalize $A]
	    # Initialize the vectors to 0, that will let us
	    # enter the loop
	    set B {0 0 0}
	    set C {0 0 0}
	    # We need to avoid getting a 0 vector, we'll test that.
	    set testzero {0 0 0}
	    while {$B == $testzero || $C == $testzero } {
		set B "[expr rand()] [expr rand()] [expr rand()]"
		set B [normalize $B]
		set C [crossproduct $A $B]
		set C [normalize $C]
	    }
	    return $C
	}
	
	# ::peptideb::utils::distance
	#
	# Returns the distance between two beads
	# Arguments: the two position vectors
	# returns a scalar
	#
	proc distance {A B} {
	    set x [expr [lindex $B 0] - [lindex $A 0]]
	    set y [expr [lindex $B 1] - [lindex $A 1]]
	    set z [expr [lindex $B 2] - [lindex $A 2]]
	    return [expr sqrt($x*$x + $y*$y + $z*$z)]
	}
	
	
	# ::peptideb::utils::average
	#
	# Returns the average of two scalars a and b
	# Arguments: 2 scalars a and b
	# Returns a scalar
	#
	proc average {a b} {
	    return [expr ($a+$b)/2.]
	}
	
	# ::peptideb::utils::sum
	#
	# Returns the sum of two numbers.
	proc sum {a b} {
	    return [expr $a + $b]
	}

	# ::peptideb::utils::g_mean
	#
	# Returns the geometric mean of two numbers.
	proc g_mean {a b} {
	    return [expr sqrt($a*$b)]
	}
 
	# ::peptideb::utils::get_dihedral
	#
	# Returns the dihedral angle defined by 4
	# atoms.
	# Arguments : 4 lists of coords (of atoms)
	# Returns a scalar
	#
	#
	#              A           D
	#           r1  \         / r3
	#                B- - - -C
	#                    r2
	#
	#
	proc get_dihedral {A B C D} {
	    if { [llength $A] != 3 || [llength $B] != 3 || \
		     [llength $C] != 3 || [llength $D] != 3 } {
		::mmsg::err [namespace current] "Calculating dihedral, one of the atoms is not well defined."
	    }
	    
	    set r1 [diff_2_vectors $B $A]
	    set r2 [diff_2_vectors $C $B]
	    set r3 [diff_2_vectors $D $C]
	    
	    set p1 [normalize [crossproduct $r1 $r2]]
	    set p2 [normalize [crossproduct $r2 $r3]]
	    
	    set r2 [normalize $r2]
	    set cosphi [dotproduct $p1 $p2]
	    set sinphi [dotproduct $r2 [crossproduct $p1 $p2]]
	    
	    return [expr atan2($sinphi,$cosphi)]
	}
	
	# Returns the center of mass of chain
	# This particular routine will average over CB atoms only
	# which are positioned at indices 4*k+2, looping over k.
	# Argument : unfolded coordinates of the chain.
	# Returns : the position of the center of mass
	proc center_of_mass {coords} {
	    set result [list 0. 0. 0.]
	    set index 0
	    for {set k 0} {$k<[expr [llength $coords]/4.]} {incr k} {
		# Loop over x, y, z
		for {set i 0} {$i<3} {incr i} {
		    lset result $i [expr [lindex $result $i]+[lindex [lindex $coords [expr 4*$k+2]] $i]]
		}
		incr index
	    }
	    for {set i 0} {$i<3} {incr i} {
		lset result $i [expr [lindex $result $i]/(1.*$index)]
	    }
	    return $result
	}	
    }
}
