# ::peptideb::import
#
# read the secondary structure of a protein
# from a pdb file.
#
#

namespace eval peptideb {
    namespace eval input {
	# Import and CG the secondary structure of a protein from
	# a pdb file
	# The subroutine reads the file, keeps only useful 
	# information and coarse-grain accordingly.
	# Argument: name of the file.
	# Optional argument : whether we recreate the side chain, or read the CB atoms.
	# if exactCB == 0 we recreate the side chain, if it's 1 we read it.
	# Returns an array of coordinates of all the beads.
	# Create the sequence of amino acid.
	proc import_pdb {file {exactCB 0}} {
	    set peptideb::frozen ""
	    set peptideb::amino_acids ""
	    set return_value ""
	    # Return variable
	    set final_coords ""
	    # Bookkeeping of number of AA
	    set index_AA 0
	    # Keep the number of the first AA in index_AA_first. Do not set now.
	    # Temporary AA
	    set temp_AA ""
	    # Bookkeeping of which peptide we're working with
	    set peptideID ""
	    # New peptide
	    set newPeptide 0


	    ::mmsg::send [namespace current] "Importing peptide from pdb file '$file'..."

	    set fp [open $file r]
	    fconfigure $fp -buffering line

	    # Trick the while loop to enter it :
	    set data 0
	    # data2 is a trick used to read the last amino acid.
	    set data2 0

	    # Keep only the first snapshot of the molecule (before first ENDMDL statement).
	    while {$data != "" && [lindex $data 0]!="ENDMDL" } {

		set data2 $data
		# Reads the next line of the file
		gets $fp data

		# We're only interested in the molecule's properties, not the headers.
		if { [lindex $data 0]=="ATOM" || [lindex $data2 0]=="ATOM"} {    
		    set data [read_line $data]
		    
		    # Check that $index_AA is the same as the AA we're reading
		    if {$index_AA != [lindex $data 5]} {
			
			    
			# Now if it's not the first AA, we can update the coordinates of
			# the last AA.
			if { [info exists index_AA_first] } {
			    # Build the side chain
			    set name [lindex $temp_AA [expr $index_AA - $index_AA_first]]
			    # Check that N, CA, and C do exist !
			    if { ![info exists N_atom] || ![info exists CA_atom] || \
				     ![info exists C_atom]} {
				::mmsg::err [namespace current] "Problem in parsing $file during side chain reconstruction."
			    }
			    # For the exact reading, we need all N, CA, C, and CB.
			    if { $exactCB == 1 } {
				if { ![info exists CB_atom] } {
				    ::mmsg::err [namespace current] "Can't parse Carbon beta during side chain reconstruction."
				}
			    } else {
				# If it's the first amino acid, use random phi, otherwise
				# calculate it.
				if { [expr $index_AA - $index_AA_first] == 0 || $newPeptide==1 } { 
				    set phi [peptideb::utils::random_1d]
				} else {
				    set phi [peptideb::utils::get_dihedral [lindex $final_coords end] $N_atom $CA_atom $C_atom]
				}
				
				# Calculate new atom Cb in local coord system
				set proj_angle [projected_angle $peptideb::angleR_CbCaC $peptideb::angleR_NCaCb $peptideb::angleR_NCaC]
				set phiCb [phi_isomer $phi L $proj_angle]
				set theta [expr $peptideb::pi - $peptideb::angleR_NCaCb]

				if { $name == "GLY" } { set R $peptideb::bondCaCb_GLY } else { set R $peptideb::bondCaCb_XXX }
				
				set Cb_vector [sphericalcoords $R $theta $phiCb]
				
				set bc [peptideb::utils::diff_2_vectors $CA_atom $N_atom]
				if { [expr $index_AA - $index_AA_first] == 0 || $newPeptide==1} {
				    set ab [::peptideb::utils::randomperpvec $bc]
				} else {
				    set ab [peptideb::utils::diff_2_vectors [lindex $final_coords end] $N_atom]
				}
				set n [bc_ab_to_n $bc $ab]
				
				set CB_atom [change_frame $Cb_vector $bc $n $CA_atom]
			    }

			    if { [info exists N_atom] && [info exists CA_atom] && [info exists C_atom]} {
				lappend final_coords $N_atom
				lappend final_coords $CA_atom
				lappend final_coords $CB_atom
				lappend final_coords $C_atom
				if { $name == "GLY" } { 
				    if { [lindex $data2 9] == "0.00" } {
					lappend peptideb::frozen 1
				    } else {
					lappend peptideb::frozen 0
				    }
				}
			    } else {
				::mmsg::err [namespace current] "Problem in parsing $file. Can't recover all atoms."
			    }
			    unset N_atom
			    unset CA_atom
			    unset CB_atom
			    unset C_atom
			}
			if { ![info exists index_AA_first]} {
			    set index_AA_first [lindex $data 5]
			}
			set index_AA [lindex $data 5]

			set newPeptide 0
			# Check whether the peptideID changed. If so, append the full set
			# of coordinates to recreate the whole peptide, and start an empty set.
			if {$peptideID!=[lindex $data 11]} {
			    set peptideID [lindex $data 11]
			    # Update coordinates only if we have a non-empty array
			    if {$final_coords!=""} {
				lappend return_value $final_coords
				lappend peptideb::amino_acids $temp_AA
				set newPeptide 1
			    }
			    set final_coords ""
			    set temp_AA ""
			}
			
			# update the list of amino acids
			if {[lindex $data 0]=="ATOM"} {
			    lappend temp_AA [lindex $data 3]
			}


		    }
		    
		    set atom_type [lindex $data 2]
		    switch $atom_type "N" {
			set N_atom "[lindex $data 6] [lindex $data 7] [lindex $data 8]"
			if { [lindex $data 9] == "0.00" } {
			    lappend peptideb::frozen 1
			} else {
			    lappend peptideb::frozen 0
			}
		    } "CA" { 
			set CA_atom "[lindex $data 6] [lindex $data 7] [lindex $data 8]"
			if { [lindex $data 9] == "0.00" } {
			    lappend peptideb::frozen 1
			} else {
			    lappend peptideb::frozen 0
			}
		    } "C" { 
			set C_atom "[lindex $data 6] [lindex $data 7] [lindex $data 8]"
			if { [lindex $data 9] == "0.00" } {
			    lappend peptideb::frozen 1
			} else {
			    lappend peptideb::frozen 0
			}
		    } "CB" {
			set CB_atom "[lindex $data 6] [lindex $data 7] [lindex $data 8]"
			if { [lindex $data 9] == "0.00" } {
			    lappend peptideb::frozen 1
			} else {
			    lappend peptideb::frozen 0
			}			
		    }
		}
		

	    }
	    close $fp

	    if {$temp_AA != ""} {
		lappend peptideb::amino_acids $temp_AA
		lappend return_value $final_coords
	    }

	    ::mmsg::send [namespace current] "PDB file was parsed successfully."

	    return $return_value
	}


	# A pdb file is pretty strange. it's organized in columns, sometimes
	# there's no spaces between two different elements. This function
	# reorgranizes columns into proper elements. Each element of a line
	# becomes an element of a list that will be returned in argument.
	# Argument : data line (it's a string - as everything else in Tcl...)
	# Returns a list
	proc read_line {data} {
	    # Element #3 (4th element) is a 3-letter code. If its length is more
	    # than that, fix it ! It's because there are no spaces between elements
	    # number 2 and 3.

	    # Problem if the length of element 2 is 7. That means elements 2 and 3 have been
	    # fused and need to be separated. See below for A's and B's.
	    if { [string length [lindex $data 2]] == 7 } {
		# now separate this element in 2, and keep it only if it has label A.
		if { [string index [lindex $data 2] 3]=="A" } {
		    set data [linsert $data 3 [string range [lindex $data 2] 4 end]]
		    lset data 2 [string replace [lindex $data 2] 3 end ""]
		} elseif { [string index [lindex $data 2] 3]=="B" } {
		    set data [linsert $data 3 [string range [lindex $data 2] 4 end]]
		    # Replace whatever's on element 2 so that it won't get read later.
		    lset data 2 "__"
		} else {
		    ::mmsg::err [namespace current] "Bad 3-letter code. Can't parse $data."
		}
	    }

	    # Looks like some amino acids have recorded everything twice using some
	    # strange scheme A--- and B--- where dashes represent the usual three-
	    # letter codes. Just delete the B ones and keep the A's.
	    if { [string length [lindex $data 3]] == 4 } {
		if { [string index [lindex $data 3] 0]=="A" } {
		    lset data 3 [string replace [lindex $data 3] 0 0 ""]
		} elseif { [string index [lindex $data 3] 0]=="B" } {
		    # Replace whatever's on element 2 so that it won't get read later.
		    lset data 2 "__"
		} else {
		    ::mmsg::err [namespace current] "Bad 3-letter code. Can't parse $data."
		}
	    }

	    return $data
	}

    }
}
