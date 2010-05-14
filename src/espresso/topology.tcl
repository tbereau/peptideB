# ::peptideb::espresso
#
# Create topology of the chain
#
#

namespace eval peptideb {
    namespace eval espresso {

	# Set all bonded interactions in ESPResSo
	proc set_bonded_interactions { bonded_parms } {
	    foreach bondtype $bonded_parms {
		if { [catch {eval [concat inter $bondtype] } errmsg ] } {
		    mmsg::err [namespace current] "couldn't set interaction: $bondtype\n$errmsg"
		} else {
		    # mmsg::send [namespace current] "set interaction: $bondtype "
		}
	    }
	    return
	}

	# Set all the non-bonded interactions.
	proc set_nb_interactions { interactionlist } {
	    foreach intertype $interactionlist {
		if { [catch { eval [concat inter  $intertype ] } errmsg ] } {
		    mmsg::err [namespace current] "could not set interaction: $intertype \n$errmsg"
		}
		# mmsg::send [namespace current] "set interaction: $intertype "
	    }
	    return
	}

	# Set all the angle potentials
	proc set_angle_interactions { interactionlist } {
	    foreach intertype $interactionlist {
		if { [catch { eval [concat inter  $intertype ] } errmsg] } {
		    mmsg::err [namespace current] "could not set interaction: $intertype \n$errmsg"
		}
		# mmsg::send [namespace current] "set interaction: $intertype "
	    }
	    return
	}


	# Set all the dihedral potentials
	proc set_dihedral_interactions { interactionlist } {
	    foreach intertype $interactionlist {
		if { [catch { eval [concat inter  $intertype ] } errmsg] } {
		    mmsg::err [namespace current] "could not set interaction: $intertype\n$errmsg"
		}
		# mmsg::send [namespace current] "set interaction: $intertype "
	    }
	    return
	}



	# Takes the chain in argument and create a topology
	proc create_topology { coords {set_interactions 1}} {
	    # Source the interactions, in case a paramater was set in the config file.
	    namespace eval :: {
		source [set env(PEPTIDEB_DIR)]/src/chain_sidechain.tcl
                source [set env(PEPTIDEB_DIR)]/src/chain_interactions.tcl
		source [set env(PEPTIDEB_DIR)]/src/hfip_interactions.tcl
            }

	    # Create an id on all the beads
	    set j 0

	    # Loop over all peptides
	    for { set l 0 } { $l < [expr [llength $coords] - $peptideb::hfip_flag]} { incr l } {    
		set chain [lindex $coords $l]

		# Verify that our chain is a multiple of 4, although it's a rather weak sanity test...
		set size [llength $chain]
		if { [expr $size % 4] != 0 } {
		    ::mmsg::err [namespace current] "The chain should be of the form (N-Ca-Cb-C)."
		}
		
		set jprime 0
		for { set k 0 } { $k < [expr $size/4] } { incr k } {
		    set a_acid [lindex [lindex [lindex $peptideb::amino_acids $l] $k] 0]
		    set a_acid [aa_code_hp $a_acid]

		    # Two types of N beads : regular ones, and proline ones which do not Hbond.
		    if {$a_acid != "pro"} {
			# N atoms have type 0
			part $j pos [lindex [lindex $chain $jprime] 0] [lindex [lindex $chain $jprime] 1]\
			    [lindex [lindex $chain $jprime] 2] type 0
			part $j molecule_id 0
		    } else {
			# Pro-N atoms have type 3
			part $j pos [lindex [lindex $chain $jprime] 0] [lindex [lindex $chain $jprime] 1]\
			    [lindex [lindex $chain $jprime] 2] type 3
			part $j molecule_id 0
		    }
		    incr j
		    incr jprime
		    # Ca atoms have type 1
		    part $j pos [lindex [lindex $chain $jprime] 0] [lindex [lindex $chain $jprime] 1]\
			[lindex [lindex $chain $jprime] 2] type 1
		    part $j molecule_id 0
		    incr j
		    incr jprime
		    # Cb atoms have type ? - side chains have different types, depending on their hydrophobicity (see side_chain_hp)
		    side_chain_hp $chain $j $jprime $a_acid 0
		    incr j
		    incr jprime
		    # C atoms have type 2
		    part $j pos [lindex [lindex $chain $jprime] 0] [lindex [lindex $chain $jprime] 1]\
			[lindex [lindex $chain $jprime] 2] type 2
		    part $j molecule_id 0
		    incr j
		    incr jprime
		}

		if {[expr $l+1==$peptideb::virtual_com] && [expr $peptideb::virtual_com==1]} {
		    for { set jflag 0 } { $jflag < $j } { incr jflag } {
			part $jflag molecule_id 0
		    }
		    part $j pos 0 0 0 type 99
		    part $j molecule_id 0
		    part $j virtual 1
		    # Right now, crash the simulation if there's more than one peptide. 
		    # FIX that later
		    if { [expr [llength $coords] - $peptideb::hfip_flag]> 1} {
			::mmsg::send [namespace current] "Virtual site only supports simulations of one peptide."
			exit 1		       
		    }
		}
	    }

	    # Do we have hfip in the system?
	    if { $peptideb::hfip_flag } {
		# Loop over hfip molecules		
		set chain [lindex $coords [llength $peptideb::amino_acids]]
		foreach hfip_mol $chain {
		    if { [llength $hfip_mol] != 3 } {
			::mmsg::err [namespace current] "An HFIP molecule should contain 3 beads!"
		    }
		    # F bead
		    set bead [lindex $hfip_mol 0]
		    part $j pos [lindex $bead 0] [lindex $bead 1] [lindex $bead 2] type 31
		    incr j
		    # Cc bead
		    set bead [lindex $hfip_mol 1]
		    part $j pos [lindex $bead 0] [lindex $bead 1] [lindex $bead 2] type 30
		    incr j
		    # F bead
		    set bead [lindex $hfip_mol 2]
		    part $j pos [lindex $bead 0] [lindex $bead 1] [lindex $bead 2] type 31
		    incr j
		}
	    }

	    # Add virtual site for peptide $peptideb::virtual_com ('1' is the first peptide)	    
	    if { [expr $peptideb::virtual_com > 0] } {
		set chain [lindex $coords 0]
		set size [llength $chain]
		analyze set chains 0 1 [expr $size+1]
		incr j
	    }
	    
	    if {$set_interactions==1} {
		# Add bonded interaction.
		set_bonded_interactions $peptideb::bonded_parms
		
		# Add angle potentials
		set_angle_interactions $peptideb::angle_parms
		
		# Add dihedral potentials
		set_dihedral_interactions $peptideb::dihedral_parms
		

		# Add these interactions now that we have all the partners defined.
		set j 0
		for { set l 0 } { $l < [expr [llength $coords] - $peptideb::hfip_flag] } { incr l } {    
		    set size [llength [lindex $coords $l]]
		    for { set k 0 } { $k < [expr $size/4] } { incr k } {
			set a_acid [lindex [lindex [lindex $peptideb::amino_acids $l] $k] 0]
			set a_acid [aa_code_hp $a_acid]
			# Determine next amino acid (for PRO omega dihedral)					
			if {$k<[expr [expr $size/4]-1]} {
			    set next_a_acid [lindex [lindex [lindex $peptideb::amino_acids $l] [expr $k+1]] 0]
			    set next_a_acid [aa_code_hp $next_a_acid]
			}
			# atom N		    
			# Covalent bond N Ca
			part $j bond 0 [expr $j + 1]
			# Subt LJ N-Ca
			#		    part $j bond 13 [expr $j + 1]
			# Subt LJ N-Cb
			#		    part $j bond 13 [expr $j + 2]
			# Subt LJ N-C
			#		    part $j bond 13 [expr $j + 3]
			# previous chain test
			if { $k > 0 } {
			    # Angle C - N - Ca
			    part $j bond  8 [expr $j - 1] [expr $j + 1]
			    # Dihedral \phi : C - N - Ca - C
			    #part $j bond  9 [expr $j - 1] [expr $j + 1] [expr $j + 3] 
			    # Dipolar interaction for \phi : C- N - Ca - C
			    part $j bond 11 [expr $j - 1] [expr $j + 1] [expr $j + 3]
			}
			incr j
			
			# atom Ca

			# Covalent bond Ca-Cb (either GLY or non-GLY)
			if { $a_acid == "gly"} { set index 15 } else { set index 16 }
			part $j bond $index  [expr $j + 1]
			unset index
			# Subt LJ Ca-Cb
			#		    part $j bond 13 [expr $j + 1]
			# Covalent bond Ca C
			part $j bond 2  [expr $j + 2]
			# Subt LJ Ca-C
			#		    part $j bond 13 [expr $j + 2]
			# Angle N - Ca - Cb
			part $j bond 4  [expr $j - 1] [expr $j + 1]
			# Angle N - Ca - C
			part $j bond 5  [expr $j - 1] [expr $j + 2]
			# Angle Cb - Ca - C
			part $j bond 6  [expr $j + 1] [expr $j + 2]
			if { [expr $k*4+4] < $size } {
			    # Dihedral \psi : N - Ca - C - N
			    #part $j bond 9  [expr $j - 1] [expr $j + 2] [expr $j + 3]
			    # Dipolar interaction for \psi : N - Ca - C - N
			    part $j bond 11 [expr $j - 1] [expr $j + 2] [expr $j + 3]
			    # Subt LJ Ca-N
			    #			part $j bond 13 [expr $j + 3]
			}
			# Improper dihedral : N - Ca - Cb - C
			# This depends on the initial chirality of the amino acid
			if { [::peptideb::input::isomer $l $k] == "L" } {
			    part $j bond 12 [expr $j - 1] [expr $j + 2] [expr $j + 1]
			} else {
			    part $j bond 12 [expr $j - 1] [expr $j + 1] [expr $j + 2]
			}
			# Hat potential
			if { $peptideb::hat_potential > 0. && $peptideb::HB_bilayer_z0 > 0.} {
			    part $j bond 17 
			}
			incr j
			
			# atom Cb
			# Subt LJ Cb-C
			#                       part $j bond 13 [expr $j + 1]
			incr j
			
			# atom C
			# end of chain test
			if { [expr $k*4+4] < $size } {
			    # Covalent bond C N
			    part $j bond 3 [expr $j + 1]
			    # Subt LJ C-N
			    #			part $j bond 13 [expr $j + 1]
			    # Angle Ca - C - N
			    part $j bond 7 [expr $j - 2] [expr $j + 1]
			    # Dihedral \omega : Ca - C - N - C
			    if { $next_a_acid == "pro"} { 
				part $j bond 14 [expr $j - 2] [expr $j + 1] [expr $j + 2]
			    } else { 
				part $j bond 10 [expr $j - 2] [expr $j + 1] [expr $j + 2]						
			    }
			    
			    # Subt LJ C-Ca
			    #			part $j bond 13 [expr $j + 2]
			}
			incr j
		    }
		}
		# HFIP in the simulation?
		if { $peptideb::hfip_flag } {
		    # Loop over hfip molecules		
		    set chain [lindex $coords [llength $peptideb::amino_acids]]
		    foreach hfip_mol $chain {
			if { [llength $hfip_mol] != 3 } {
			    ::mmsg::err [namespace current] "An HFIP molecule should contain 3 beads!"
			}
			# F-Cc bond
			part $j bond 30 [expr $j+1]
			# Cc-F bond
			part [expr $j+1] bond 30 [expr $j+2]
			# F-Cc-F angle
			part [expr $j+1] bond 31 $j [expr $j+2]
			incr j 3
		    }
		}

		

		# Add non-bonded interaction.
		set_nb_interactions $peptideb::nb_interactions
		
		# Magic formula - Excludes non-bonded interactions for bonded
		# partners linked by 1 one or 2 bonds.
		part auto_exclusions 2
	    }

	    return 
	}

	# Reads the topology from espresso and output a coordinate file.
	# Assumes that the number of particles hasn't changed since the 
	# initialization of the chain.
	proc read_topology { } {
	    set final_coords ""
	    set j 0
	    for { set l 0 } { $l < [llength $peptideb::amino_acids] } { incr l } {
		set coords ""
		for { set k 0 } { $k < [expr [llength [lindex $peptideb::amino_acids $l]] * 4] } { incr k } {
		    lappend coords [part $j print pos]
		    incr j
		}
		# Now calculate the COM of the peptide $l with unfolded coordinates $coords
		set com [::peptideb::utils::center_of_mass $coords]
		# Fold the COM
		set com [fold_coords $com [setmd box_l]]
		set coords [choose_PBC_image $coords $com [setmd box_l]]

		lappend final_coords $coords
	    }
	    # hfip
	    if { $peptideb::hfip_flag } {
		set hfip_coords ""
		for { set k 0 } { $k < $peptideb::hfip_num_mol } { incr k } {
		    set coords ""
		    for { set kk 0 } { $kk < 3 } { incr kk } {
			lappend coords [part $j print pos]
			incr j
		    }
		    # Now calculate the COM of the hfip molecule $l with unfolded coordinates $coords
		    set com [::peptideb::utils::center_of_mass $coords]
		    # Fold the COM
		    set com [fold_coords $com [setmd box_l]]
		    set coords [choose_PBC_image $coords $com [setmd box_l]]
		    #set hfip_coords [concat $hfip_coords $coords]
		    lappend hfip_coords $coords
		}	       
		lappend final_coords $hfip_coords
	    }
	    return $final_coords
	}



	# Reads a set of 3 coordinates and returns the folded coordinates
	# Arguments : [list x y z] [list Lx Ly Lz]
	# Returns the folded set of coordinates
	proc fold_coords {coordinates box_size} {
	    for {set i 0} {$i<3} {incr i} {
		while {[lindex $coordinates $i] < 0.} {
		    lset coordinates $i [expr [lindex $coordinates $i]+[lindex $box_size $i]]
		}
		if {[lindex $coordinates $i] > [lindex $box_size $i]} {
		    lset coordinates $i [expr fmod([lindex $coordinates $i],[lindex $box_size $i])]
		}
	    }
	    return $coordinates
	}

	# Choose the periodic boundary image such that it minimizes the
	# distance of a set of coordinates of ONE polypeptide with a given COM.
	# Arguments : coordinates of the peptide; coords of COM; box size in 3D
	proc choose_PBC_image {chain com box} {
	    for {set k 0} {$k<[llength $chain]} {incr k} {
		for {set i 0} {$i<3} {incr i} {
		    while {[expr abs([lindex [lindex $chain $k] $i]-[lindex $com $i]\
					 +[lindex $box $i])]<[expr abs([lindex [lindex $chain $k] $i]-[lindex $com $i])]} {
			lset chain $k $i [expr [lindex [lindex $chain $k] $i]+[lindex $box $i]]
		    }
		    while {[expr abs([lindex [lindex $chain $k] $i]-[lindex $com $i]\
					 -[lindex $box $i])]<[expr abs([lindex [lindex $chain $k] $i]-[lindex $com $i])]} {
			lset chain $k $i [expr [lindex [lindex $chain $k] $i]-[lindex $box $i]]
		    }
		}
	    }
	    
	    return $chain
	}
	
	
	
	# This routine takes care of all non-bonded interactions
	# between side chains.
	# Arguments: 
	#   - chain is the set of coordinates (see chain in create
	#     topology)
	#   - j is the index of the amino acid we're looking at
	#   - jprime is some other index that lets the magic happen.
	#     see create_topo for more information.
	#   - aa is the amino acid. It can be either a 1- or 
	#     3-letter code. 
	#   - molecule ID
	# Returns : nothing. It appends the list of interactions in
	# Espresso.
	# Beads for side chain are identified from 10 to 29.
	proc side_chain_hp { chain j jprime aa_type mol_count} {
	    set index 10
	    foreach name $peptideb::3letter_list {
		set name [string tolower $name]
		if { $aa_type == $name } {
		    part $j pos [lindex [lindex $chain $jprime] 0] [lindex [lindex $chain $jprime] 1]\
			[lindex [lindex $chain $jprime] 2] type $index q [charge $name]
		    part $j molecule_id $mol_count
		}
		incr index
	    }
	}

	# Returns the charge of an amino acid side chain.
	# Takes one argument : the side chain name (lower case).
	proc charge { name } {
	    switch $name "lys" { return 1
	    } "arg" { return 1
	    } "asp" { return -1
	    } "glu" { return -1
	    } default { return 0 }
	}

	# Parses the side chain letter codes
	# It can either be a 1- or 3-letter code.
	# Arguments : amino acid string code
	# Returns : the 3-letter code in lower case.
	proc aa_code_hp { aa } {
	    set aa [string tolower $aa]
	    switch -- $aa  {
		"d" - "asp" { return "asp" }
		"e" - "glu" { return "glu" }
		"r" - "arg" { return "arg" }
		"k" - "lys" { return "lys" }
		"h" - "his" { return "his" }
		"n" - "asn" { return "asn" }
		"q" - "gln" { return "gln" }
		"s" - "ser" { return "ser" }
		"t" - "thr" { return "thr" }
		"y" - "tyr" { return "tyr" }
		"a" - "ala" { return "ala" }
		"g" - "gly" { return "gly" }
		"v" - "val" { return "val" }
		"l" - "leu" { return "leu" }
		"i" - "ile" { return "ile" }
		"p" - "pro" { return "pro" }
		"f" - "phe" { return "phe" }
		"m" - "met" { return "met" }
		"w" - "trp" { return "trp" }
		"c" - "cys" { return "cys" }
		default { ::mmsg::err [namespace current] "$aa is not a valid amino acid code." }
	    }
	}


	# Deletes all particles and interactions.
	# Useful when starting a new simulation whithin the script.
	proc reset_sim { } {
	    part deleteall
	    set peptideb::bonded_parms ""
	    set peptideb::angle_parms ""
	    set peptideb::dihedral_parms ""
	    set peptideb::nb_interactions ""
	    return 
	}
    }
}
