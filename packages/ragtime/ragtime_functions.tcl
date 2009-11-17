# RAGTIME
#
# This file contains all the routines
# required for ragtime.tcl
#
# @author : Tristan Bereau 
#          bereau@cmu.edu
#
#
# @version: 1.4
#
# @date   : 09/22/2009
#


# Call the stride program, store the output 
# in a variable, and return it.
# Takes one argument : the PDB file gotten in
# argument from ragtime
proc call_stride pdb {
    variable verbose
    variable force
    variable failed
    variable b_PBC
    variable b_PBCsize

    set without_pbc [expr $b_PBC - 1]

    # If b_PBC is turned on (periodic boundary condition), attempt to call
    # with the corresponding option. If it fails, go to the regular call.
    if { $b_PBC == 1 } {
	if { [catch {exec stride -a$b_PBCsize -h -f[concat $pdb.str] $pdb} errmsg] } {
	    puts $errmsg
	    set without_pbc 1
	    exec rm -f [concat $pdb.str]
	}
    }

    # $without_pbc will only be ==1 if:
    # 1) there is no pbc command invoked (b_PBC == 0)
    # 2) the stride command WITH pbc failed
    if { $without_pbc } {
	if { [catch {exec stride -h -f[concat $pdb.str] $pdb} errmsg] } {
	    puts "Problem with STRIDE:\n$errmsg"
	    if {$force==0} {
		exit 1
	    }
	    set failed 1
	}
	if {$verbose} {
	    puts "STRIDE analysis : $pdb"
	}
    }
}

# Read PDB file and store the number of amino acids for each chains.
# Also, store the coordinates of CA atoms.
proc read_pdb pdb {
    variable CA_coords ""
    variable alphabet
    variable number_aa ""

    set chain_id ""
    set aa_count 0
    
    set fp [open $pdb r]
    fconfigure $fp -buffering line
    gets $fp data
    while {$data != ""} {
	if {[lindex $data 0]=="ATOM"} {
	    if {[lindex $data 4]!=$chain_id} {
		if {$chain_id != ""} { 
		    lappend number_aa $aa_count
		}
		set chain_id [lindex $data 4]
		set aa_count 0
	    } else {
		if {[lindex $data 2]=="CA"} {
		    lappend CA_coords [list [lindex $data 6]\
					   [lindex $data 7]\
					   [lindex $data 8]]
		    
		    incr aa_count
		}
	    }
	}
	gets $fp data
    }
    close $fp

    # Don't forget the last amino acid
    lappend number_aa $aa_count
}

# Store all N, H, C, and O atoms in hbonds_coords.
# Assume that all N, H, C, and O atoms are present and grouped by residue.
# N-H not taken into account for PRO residue.
# Calculate the total hydrogen bond energy. Store 
# it in variable hbond--second element.
# Return hbond energy.
proc hbond_energy pdb {
    variable hbonds

    set N_coords ""
    set H_coords ""
    set C_coords ""
    set O_coords ""

    set result 0.

    set fp [open $pdb r]
    fconfigure $fp -buffering line
    gets $fp data
    while {$data != ""} {
	if {[lindex $data 0]=="ATOM"} {
	    if {[lindex $data 2]=="N" && [lindex $data 3]!="PRO"} {
		lappend N_coords [list [lindex $data 6]\
				      [lindex $data 7]\
				      [lindex $data 8]]
	    } elseif {[lindex $data 2]=="HN" && [lindex $data 3]!="PRO"} {
		lappend H_coords [list [lindex $data 6]\
				      [lindex $data 7]\
				      [lindex $data 8]]
	    } elseif {[lindex $data 2]=="C"} {
		lappend C_coords [list [lindex $data 6]\
				      [lindex $data 7]\
				      [lindex $data 8]]
	    } elseif {[lindex $data 2]=="O"} {
		lappend O_coords [list [lindex $data 6]\
				      [lindex $data 7]\
				      [lindex $data 8]]
	    }
	    
	}
	gets $fp data
    }
    close $fp
    
    # Sanity check: atom lists N-H and C-O must have the same lengths
    if {[llength $N_coords] != [llength $H_coords] || 
	[llength $C_coords] != [llength $O_coords]} {
	puts "Error -- Not a consistent number of N-H or C-O atoms. Exiting."
	exit 1
    }


    #### Calculation of HBond energy for each pair ####
    set eps 6.
    set cut 8.
    set sigma 4.11
    
    # loop over all N-H groups
    for {set j 0} {$j<[llength $N_coords]} {incr j} {
	set N_atom [lindex $N_coords $j]
	set H_atom [lindex $H_coords $j]
	# loop over all C-O groups
	for {set k 0} {$k<[llength $C_coords]} {incr k} {
	    set C_atom [lindex $C_coords $k]
	    set O_atom [lindex $O_coords $k]

	    # Are N and C sufficiently close?
	    set dist [distance $N_atom $C_atom]
	    if {$dist<=$cut} {
		set cos_jik [cos_angle [2bead_vector $N_atom $H_atom] \
				 [2bead_vector $N_atom $C_atom]]
		set cos_ikn [cos_angle [2bead_vector $C_atom $N_atom] \
				 [2bead_vector $C_atom $O_atom]]
		if {$cos_jik>0. && $cos_ikn>0.} {
		    set frac2 [expr $sigma*$sigma/($dist*$dist)]
		    set frac10 [expr $frac2*$frac2*$frac2*$frac2*$frac2]
		    set result [expr $result + ($eps * $frac10 * (5.0 * $frac2 - 6.0) \
						    * $cos_jik * $cos_jik * $cos_ikn * $cos_ikn)]
		}
	    }
	}
    }
    return $result
}

# Returns the distance between two beads
# Arguments: the two position vectors
# returns a scalar
proc distance {A B} {
    set x [expr [lindex $B 0] - [lindex $A 0]]
    set y [expr [lindex $B 1] - [lindex $A 1]]
    set z [expr [lindex $B 2] - [lindex $A 2]]
    return [expr sqrt($x*$x + $y*$y + $z*$z)]
}


# From a list of coordinates, returns the vector
# formed by the difference in position between
# A and B (i.e. vector_12= B - A)
# Arguments:
#   - A is a vector of coordinates
#   - B is a vector of coordinates
# Returns an array of 3 coordinates: {x y z}
proc 2bead_vector {A B} {
    set result [expr [lindex $B 0] - [lindex $A 0]]
    lappend result [expr [lindex $B 1] - [lindex $A 1]]
    lappend result [expr [lindex $B 2] - [lindex $A 2]]

    return $result
}

# Returns the cos between two vectors.
# Arguments: two arrays of length 3.
# Returns an angle in *radians*.
proc cos_angle {A B} {
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
    return $result

}

# Returns the scalar product between two vectors.
# Arguments: two arrays of length 3.
# Returns a scalar.
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

# Returns the norm of a vector.
# Argument: one array of length 3.
# Returns a scalar.
proc norm {A} {
    # A norm is the square root of the scalar product
    # of a vector with itself.
    set result [expr sqrt([dotproduct $A $A])]
    if { $result == 0} {
        ::mmsg::warn [namespace current] "Found a vector of norm 0!"
    }
    return $result
}


# Load the raw data from STRIDE and throw away
# anything not related to extended conformations.
# Details on secondary structure assignments ONLY
# for extended conformations; and hydrogen bond 
# information ONLY for those ones.
# One argument : the STRIDE raw data.
# Returns : Variable. Contains all HBonds; sequence; STRIDE
proc keep_only_beta pdb {
    variable s_seq
    variable verbose
    set sequence ""
    set stride ""
    set strands ""
    set hbonds ""
    set aa_length 0

    exec sed -e s/\-//g -e s/\"//g -e s/\{//g -e s/\'//g $pdb.str > $pdb.str2
    set fp [open [concat $pdb.str2] r]
    fconfigure $fp -buffering line
    gets $fp data
    while {$data != ""} {
	switch -- [lindex $data 0] "SEQ" {
	    lappend sequence [lindex $data 2]
	    set aa_length [lindex $data 3]
	} "STR" {
	    set temp_string [string range $data 10 [expr 9+$aa_length]]
	    set temp_string [string map {" " _} $temp_string]
	    lappend stride $temp_string
	} "ACC" {
	    lappend hbonds $data 
	} "DNR" {
	    lappend hbonds $data
	}
	gets $fp data
    }
    close $fp
    
    
    if {$verbose} {
	puts "Sequence : $sequence"
    }
    
    # Determine the amino acid sequence if "-s" is turned on
    if {$s_seq==1} {
	set s_seq [lindex $sequence 0]
	for {set j 0} {$j<[string length [lindex $sequence 0]]} {incr j} {
	    # Residue 1-letter code - A - P - M - E - H
	    lappend s_seq "[string index [lindex $sequence 0] $j] 0 0 0 0 0"
	}
    }
    
    # Verify that all chains have the same sequence;
    # otherwise exit
    if {$s_seq!=""} {
	foreach aa_seq $sequence {
	    if {$aa_seq!=[lindex $s_seq 0]} {
		puts "ERROR : Cannot make residue-based statistics\
on non-identical chains."
		exit 1
	    }
	}
    }
    
    set stride [string range $stride 0 \
		    [expr [string length $sequence]-1]]
    if {$verbose} {
	puts "STRIDE   : $stride\n"
    }
    
    # Delete STRIDE output now that the information
    # has been retrieved.
    exec rm -f [concat $pdb.str]
    exec rm -f [concat $pdb.str2]
    
    return [list $hbonds $sequence $stride]
}


# Store the STRIDE information into an array
# Argument : STRIDE information + pdb index to keep track
# of the element length + old residue index
# Returns : an array that contains information relative
# to secondary structure
proc create_res_array {stride index old_array} {
    variable a_conf
    variable e_conf
    variable length_stride

    array set ret_array $old_array

    set length_stride [llength $stride]

    for {set l 0} {$l<[llength $stride]} {incr l} {
	for {set k 0} {$k<[string length [lindex $stride $l]]} {incr k} {
	    set conformation [string index [lindex $stride $l] $k]
	    if {$a_conf==1 && $e_conf==0} {
		if {$conformation !="E" && $conformation !="A" \
			&& $conformation !="P" && $conformation !="M"} {
		    set conformation "_"
		}
	    } elseif {$a_conf==1 && $e_conf==1} {
		if {$conformation !="H" && $conformation !="E" \
			&& $conformation !="A" && $conformation !="P"\
			&& $conformation !="M"} {
		    set conformation "_"
		}
	    } elseif {$a_conf==0 && $e_conf==1} {
		if {$conformation !="H"} {
		    set conformation "_"
		}
	    }
	    if {$index>0} {
		set ret_array($l,$k) [concat [set ret_array($l,$k)]$conformation]
	    } else {
		set ret_array($l,$k) $conformation
	    }
	}
    }
    return [array get ret_array]
}

# Identify the correct chain ID out of the HBond information.
# A bit subtle when there's more than 26 chains because
# there's only one letter allowed to characterize the chain.
# Argument : HBond line
# Return : both chain IDs in a list
proc read_id_from_hbond hbond_line {
    variable alphabet
    variable number_aa
    set result_id_1 ""
    set result_id_2 ""

    set chain1 [lsearch $alphabet [lindex $hbond_line 2]]
    set chain2 [lsearch $alphabet [lindex $hbond_line 7]]
    set test1 [expr [lindex $hbond_line 3] - [lindex $hbond_line 4]]
    set test2 [expr [lindex $hbond_line 8] - [lindex $hbond_line 9]]

    foreach test [list $test1 $test2] chain [list $chain1 $chain2]\
	result {result_id_1 result_id_2} {
	    if {$test==1} {
		set $result $chain
	    } elseif {[expr abs($test)%([lindex $number_aa $chain]-1)]==0} {
		set $result [expr $chain+abs($test)/\
				 ([lindex $number_aa $chain]-1)*26]
	    }
	}

    return [list $result_id_1 $result_id_2]
    
}


# Perform and store the cluster analysis on peptides.
# This algorithm uses the presence of hydrogen bonds
# between peptides to detect clusters.
proc cluster_analysis {stride} {
    variable cluster_list
    variable number_aa
    variable alphabet 
    set list_id_chains ""
    
    set cluster ""

    foreach element $stride {
	set chains_id [read_id_from_hbond $element]
	set chain1 [lindex $chains_id 0]
	set chain2 [lindex $chains_id 1]
	if {$chain1 != $chain2 &&\
		[string first "$chain1 $chain2" $list_id_chains]==-1 &&\
		[string first "$chain2 $chain1" $list_id_chains]==-1} {
	    lappend list_id_chains "$chain1 $chain2"
	}
	
    }

    for {set k 0} {$k<[llength $list_id_chains]} {incr k} {
	for {set l [expr $k+1]} {$l<[llength $list_id_chains]} {incr l} {
	    set chain1 [lindex [lindex $list_id_chains $k] 0]
	    set chain2 [lindex [lindex $list_id_chains $k] 1]
	    set chain3 [lindex [lindex $list_id_chains $l] 0]
	    set chain4 [lindex [lindex $list_id_chains $l] 1]
	    if {$chain1==$chain3 || $chain1==$chain4 ||
		$chain2==$chain3 || $chain2==$chain4 } {
		set neighbor 0
		# Try to add to an existing cluster
		for {set m 0} {$m < [llength $cluster]} {incr m} {
		    foreach number {1 2 3 4} {
			if {[lsearch [lindex $cluster $m]\
				 [set chain$number]]>-1} {
			    set neighbor 1
			    foreach number {1 2 3 4} {
				if {[lsearch [lindex $cluster $m]\
					 [set chain$number]]==-1} {
				    lset cluster $m \
					[concat [lindex $cluster $m] [set chain$number]]														 
				}
			    }	
			}
		    }
		}
		# New cluster
		if {$neighbor==0} {
		    set element [list $chain1 $chain2]
		    foreach number {3 4} {
			if {[set chain$number]!=$chain1 &&\
				[set chain$number]!=$chain2} {
			    lappend element [set chain$number]
			}
		    }
		    lappend cluster $element
		}						
	    }	
	}
    }

    # It's still possible that a few lists are actually part
    # of the same cluster.
    for {set m 0} {$m < [llength $cluster]} {incr m} {
	for {set n [expr $m+1]} {$n<[llength $cluster]} {incr n} {
	    foreach c_m [lindex $cluster $m] {
		if {[lsearch [lindex $cluster $n] $c_m]>-1} {
		    foreach c_n [lindex $cluster $n] {
			if {[lsearch [lindex $cluster $m] $c_n]==-1} {
			    lset cluster $m\
				[concat [lindex $cluster $m] $c_n]
			}
		    }
		    lset cluster $n ""
		    break
		}
	    }
	}
    }

    # sort and delete empty lists
    set cluster_temp ""
    for {set m 0} {$m<[llength $cluster]} {incr m} {
	if {[lindex $cluster $m]!=""} {
	    lappend cluster_temp [lsort -integer [lindex $cluster $m]]
	}
    }
    set cluster $cluster_temp
    unset cluster_temp

    # add clusters of two-elements
    foreach 2element $list_id_chains {
	set belongs 0
	foreach c_element $cluster {
	    if {[lsearch $c_element [lindex $2element 0]]>-1} {
		set belongs 1
	    }
	}
	if {$belongs==0} {
	    # 2element is not in any cluster
	    lappend cluster $2element
	}
    }

    # add one-elements to cluster variable
    for {set k 0} {$k<[llength $number_aa]} {incr k} {
	set belongs 0
	foreach c_element $cluster {
	    if {[lsearch $c_element $k]>-1} {
		set belongs 1
	    }
	}
	if {$belongs==0} {
	    #1-element is not in any cluster
	    lappend cluster $k
	}
    }

    # Elements of cluster are now to be stored in cluster_list
    foreach element $cluster {
	set n_size [llength $element]
	while {[llength $cluster_list]<=$n_size} {
	    set cluster_list [concat $cluster_list 0]
	}
	lset cluster_list [expr $n_size-1]\
	    [expr [lindex $cluster_list [expr $n_size-1]]+1]
    }

    # return nothing
}


# Parse information on which amino acids are HBonding
# One argument : the output from the keep_only_beta
# routine. 
# Returns : list of all HBonding amino acids that form
# extended conformations
proc analyze_hbonds data {
    variable alphabet
    
    set hbonds [lindex $data 0]
    set good_hbonds ""
    set counter 1
    
    foreach element $hbonds {
	# extract ID of both amino acids involved
	set chains_id [read_id_from_hbond $element]
	# ID 1
	set chain1 [lindex $chains_id 0]
	set residue1 [expr [lindex $element 3]-1]
	# ID 2
	set chain2 [lindex $chains_id 1]
	set residue2 [expr [lindex $element 8]-1]
	if {([string index [lindex [lindex $data 2] $chain1] $residue1]
	     =="E" && 
	     [string first E [lindex [lindex $data 2] $chain2]]>-1) ||
	    ([string index [lindex [lindex $data 2] $chain2] $residue2]
	     =="E" && 
	     [string first E [lindex [lindex $data 2] $chain1]]>-1)  ||
	    ([string index [lindex [lindex $data 2] $chain1] [expr $residue1-1]]
	     =="E" && 
	     [string first E [lindex [lindex $data 2] $chain2]]>-1) ||
	    ([string index [lindex [lindex $data 2] $chain1] [expr $residue1+1]]
	     =="E" &&
	     [string first E [lindex [lindex $data 2] $chain2]]>-1) ||
	    ([string index [lindex [lindex $data 2] $chain2] [expr $residue2-1]]
	     =="E" && 
	     [string first E [lindex [lindex $data 2] $chain1]]>-1) ||
	    ([string index [lindex [lindex $data 2] $chain2] [expr $residue2+1]]
	     =="E") &&
	    ([string first E [lindex [lindex $data 2] $chain1]]>-1)} {
	    for {set k $counter} {$k<[llength $hbonds]} {incr k} {
		if {[are_partners $element [lindex $hbonds $k]]} {
		    lappend good_hbonds $element
		    break
		}
	    }
	}
	incr counter
    }

    return $good_hbonds
}

# sheet_analysis
# Version 2 of the sheet_analysis routine
# Returns an array of assigned secondary structures
# with relative orientations for sheets.
# A : antiparallel
# P : parallel
# M : mixed (parallel + antiparallel)
proc sheet_analysis {hbonds res_array index pdb} {
    variable alphabet
    variable aggregates
    variable 2agg
    variable 3agg
    variable verbose

    array set residue $res_array

    set counter 1
    set list_E  ""
    set list_agg ""
    set id_chains ""
    set hbonds_index 0

    foreach element $hbonds {
	incr hbonds_index
	set list_agg_temp ""
	set updated_agg 0

	if {$verbose==1} {
	    puts $element
	}

	set chains_id [read_id_from_hbond $element]
	# ID 1
	set chain1 [lindex $chains_id 0]
	set residue1 [expr [lindex $element 3]-1]
	# ID 2
	set chain2 [lindex $chains_id 1]
	set residue2 [expr [lindex $element 8]-1]
	
	set list_E_temp ""

	if {[string first "$chain1 $chain2" $id_chains]==-1} {
	    lappend id_chains "$chain1 $chain2"
	    
	    foreach ss [concat [array names	residue "$chain1,*"]\
			    [array names residue "$chain2,*"]] {
		if {[string index $residue($ss) end]=="E" } {
		    lappend list_E_temp $ss
		}
	    }
	}
	if {$list_E_temp!=""} {
	    lappend list_E $list_E_temp
	    # Keep track of the number of 2-aggregates
	    lappend list_agg_temp $chain1
	    lappend list_agg_temp $chain2
	    lset aggregates 0 [expr [lindex $aggregates 0]+1]
	    # Look for larger aggregates
	    for {set k $hbonds_index} {$k < [llength $hbonds]} {incr k} {
		set chains_id_2 [read_id_from_hbond\
				     [lindex [lindex $hbonds $k]]]
		set chain3 [lindex $chains_id_2 0]
		set chain4 [lindex $chains_id_2 1]
		if {(($chain1==$chain3 && $chain2!=$chain4) || 
		     ($chain2==$chain3 && $chain1!=$chain4) ||
		     ($chain2==$chain4 && $chain1!=$chain3)) &&
		    ([lsearch $list_agg $chain1]==-1 ||
		     [lsearch $list_agg $chain2]==-1 ||
		     [lsearch $list_agg $chain3]==-1 ||
		     [lsearch $list_agg $chain4]==-1)} {
		    foreach number {1 2 3 4} {
			if {[lsearch $list_agg_temp\
				 [set chain$number]]==-1} {
			    lappend list_agg_temp [set chain$number]
			}
		    }
		}
	    }
	}
	if {$list_agg==""} {
	    lappend list_agg $list_agg_temp
	} else {
	    if {[llength $list_agg_temp]>0} {
		for {set l 0} {$l<[llength $list_agg]} {incr l} {
		    foreach partner1 $list_agg_temp {
			if {[lsearch [lindex $list_agg $l] $partner1]>-1} {
			    # New aggregates belong to an existing set
			    # Update the set with new ones
			    foreach partner2 $list_agg_temp {
				if {[lsearch [lindex $list_agg $l] $partner2]==-1} {
				    lset list_agg $l [concat [lindex $list_agg $l]\
							  $partner2]
				}
			    }
			    set updated_agg 1
			    # breaking the code will make the loop
			    # go faster
			    break
			}
		    }
		}
		if {$updated_agg==0} {
		    for {set k 0} {$k<[llength $list_agg]} {incr k} {				
			set counter 1
			foreach partner1 $list_agg_temp {
			    if {[lsearch [lindex $list_agg $k] $partner1]>-1} {
				set counter 0
			    }
			}
			if {$counter==1} {
			    lappend list_agg $list_agg_temp
			    break
			}
		    }
		}
	    }
	}
    }

    # Count the number of neighboring strands and update
    # the aggregates variable. Larger aggregates must
    # decrease the number of 2-aggregates, since they-ve
    # already been counter.
    foreach agg $list_agg {
	set n_agg [expr [llength $agg]-2]
	if {$n_agg>0} {
	    lset aggregates 0 [expr [lindex $aggregates 0]-[expr $n_agg+1]]
	    while {[llength $aggregates]<=$n_agg} {
		set aggregates [concat $aggregates 0]
	    }
	    lset aggregates $n_agg [expr [lindex $aggregates $n_agg]+1]
	}
    }
    return [compute_ss $list_E $res_array $pdb $list_agg]
}


# From pdb file, determine orientation of beta sheets
# argument list_E determines the chains that should
# be looked at.
proc compute_ss {list_E res_array pdb list_agg} {
    variable alphabet
    variable CA_coords
    variable 2agg
    variable 3agg

    array set resid_array $res_array

    foreach pair $list_E {
	set pair [lsort $pair]
	set chain1 -1
	set chain2 -1
	set residues1 ""
	set residues2 ""

	# open pdb and obtain C_alpha of all atoms involved
	# in pair
	set atoms ""
	foreach element $pair {
	    set comma   [string first "," $element]
	    set chain   [string range $element 0 [expr $comma-1]]
	    set residue [string range $element [expr $comma+1] end]
	    
	    if {$chain1 ==-1} {
		set chain1 $chain
		lappend residues1 $residue
	    } elseif {$chain1 == $chain} {
		if {[string first $residue $residues1]==-1} {
		    lappend residues1 $residue
		}
	    } else {
		set chain2 $chain
		if {[string first $residue $residues2]==-1} {
		    lappend residues2 $residue
		}
	    }
	}
	
	set residues1 [lsort $residues1]
	set residues2 [lsort $residues2]

	# Now calculate the orientation of the sheets
	# by calculating vectors which go between
	# the two C_alpha at the edges.
	set x1 [expr [lindex [lindex $CA_coords 1] 0]\
		    -[lindex [lindex $CA_coords 0] 0]]
	set y1 [expr [lindex [lindex $CA_coords 1] 1]\
		    -[lindex [lindex $CA_coords 0] 1]]
	set z1 [expr [lindex [lindex $CA_coords 1] 2]\
		    -[lindex [lindex $CA_coords 0] 2]]
	set x2 [expr [lindex [lindex $CA_coords 3] 0]\
		    -[lindex [lindex $CA_coords 2] 0]]
	set y2 [expr [lindex [lindex $CA_coords 3] 1]\
		    -[lindex [lindex $CA_coords 2] 1]]
	set z2 [expr [lindex [lindex $CA_coords 3] 2]\
		    -[lindex [lindex $CA_coords 2] 2]]
	# Calculate scalar product
	set orientation [expr $x1*$x2+$y1*$y2+$z1*$z2]
	if {$orientation>0.} {
	    set orientation P
	} else {
	    set orientation A
	}

	# Plug the result in res_array
	# Note that if a sheet already has an orientation,
	# P or A will be replaced by M (if the two orientations
	# are different)
	foreach element $pair {
	    if {[string index $resid_array($element) end]=="E"} {
		set resid_array($element)\
		    [concat [string range $resid_array($element) 0 end-1]$orientation]
	    } elseif {[string index $resid_array($element) end]=="A" ||
		      [string index $resid_array($element) end]=="P"} {
		if {[string index $resid_array($element) end]!=$orientation} {
		    set resid_array($element)\
			[concat [string range $resid_array($element) 0 end-1]M]
		}
	    }

	}
    } 

    set res_array [array get resid_array]


    # Count the ratios of P and A for 2- and 3-aggregates
    foreach agg $list_agg {
	if {[llength $agg]==2} {
	    # 2-aggregates
	    # Enough to look at one strand. This
	    # will determine the orientation
	    set index_ss [lsearch -all $res_array "[lindex $agg 0],*"]
	    foreach index $index_ss {
		set code [string index\
			      [lindex $res_array [expr $index+1]] end]
		if {$code=="P"} {
		    lset 2agg 0 [expr [lindex $2agg 0]+1]
		    break
		} elseif {$code=="A"} {
		    lset 2agg 1 [expr [lindex $2agg 1]+1]
		    break
		}
	    }
	}
	if {[llength $agg]==3} {
	    # 3-aggregates
	    # Enough to look at middle strand.
	    # It will be A, P, or M.
	    set index_ss [lsearch -all $res_array "[lindex $agg 1],*"]
	    foreach index $index_ss {
		set code [string index\
			      [lindex $res_array [expr $index+1]] end]
		if {$code=="P"} {
		    lset 3agg 0 [expr [lindex $3agg 0]+1]
		    break
		} elseif {$code=="A"} {
		    lset 3agg 1 [expr [lindex $3agg 1]+1]
		    break
		} elseif {$code=="M"} {
		    lset 3agg 2 [expr [lindex $3agg 2]+1]
		    break
		}
	    }
	}
    }
    return $res_array
}


# Determine whether two lines in the STRIDE output (HBonds) are
# partners. IDs should be reversed, and all geometrical properties
# equal.
# Arguments: the two lines
# Returns : 1 or 0 for yes or no, respectively
proc are_partners {line1 line2} {
    if {[lindex $line1 2]==[lindex $line2 2]} {
	return 0
    }
    if {[lindex $line1 4]==[lindex $line2 9] &&
	[lindex $line1 10]==[lindex $line2 10] &&
	[lindex $line1 11]==[lindex $line2 11] &&
	[lindex $line1 12]==[lindex $line2 12] &&
	[lindex $line1 13]==[lindex $line2 13] &&
	[lindex $line1 14]==[lindex $line2 14]} {
	return 1
    } else {
	return 0
    }
}

# Determine whether two lines in the STRIDE output have the
# same residue. It checks the ID of the two residues and
# takes into account the two possible combinations
# Arguments : the two lines
# Returns : 0 for no, P for parallel sheet, A for antiparallel sheet
proc sheet_orientation {line1 line2} {
    set index1 [join [list [lindex $line1 2] [lindex $line1 4]]]
    set index2 [join [list [lindex $line1 7] [lindex $line1 9]]]
    set index3 [join [list [lindex $line2 2] [lindex $line2 4]]]
    set index4 [join [list [lindex $line2 7] [lindex $line2 9]]]
    if {$index1==$index3 && $index2==$index4} {
	return A
    }
    if {[string index $index1 0]==[string index $index3 0] &&
	[string index $index2 0]==[string index $index4 0] &&
	(abs([string index $index1 end]-[string index $index3 end])==2 ||
	 abs([string index $index2 end]-[string index $index4 end])==2) &&
	!(abs([string index $index1 end]-[string index $index3 end])==0 ||
	  abs([string index $index2 end]-[string index $index4 end])==0)} {
	if {[expr ([string index $index1 end]-[string index $index3 end])/ \
		 ([string index $index2 end]-[string index $index4 end])]<0.} {
	    return A
	}
    }
    
    if {$index1==$index3 && 
	[string index $index2 0]==[string index $index4 0] &&
	(abs([string index $index2 end]-[string index $index4 end])==2 ||
	 abs([string index $index2 end]-[string index $index4 end])==3)} {
	return P
    }
    if {$index2==$index4 &&
	[string index $index1 0]==[string index $index3 0] &&              
        (abs([string index $index1 end]-[string index $index3 end])==2 ||
	 abs([string index $index1 end]-[string index $index3 end])==3)} {
        return P                                                          
    }     
    if {[string index $index1 0]==[string index $index3 0] &&
	[string index $index2 0]==[string index $index4 0] &&
	abs([string index $index1 end]-[string index $index3 end])==2 &&
	abs([string index $index2 end]-[string index $index4 end])==2 &&
	[expr ([string index $index1 end]-[string index $index3 end])/ \
	     ([string index $index2 end]-[string index $index4 end])]==1} {
	return P
    }

    # All other cases
    return 0
}

# Determine the contact order parameter: N/L.
# N: number of contacting residues (distance less than $ca_threshold)
# L: number of residues in the protein
# This function ASSUMES that only one protein is considered.
proc contact_order_parameter {} {
    variable CA_coords
    variable CA_threshold
    
    set num_residues [llength $CA_coords]
    set num_contacts 0
    for {set i 0} {$i < $num_residues} {incr i} {
	for {set j [expr $i+1]} {$j < $num_residues} {incr j} {
	    if {[distance [lindex $CA_coords $i] [lindex $CA_coords $j]] < $CA_threshold} {
		incr num_contacts
	    }
	}
    }
    return [expr $num_contacts/($num_residues*1.)]
}

# Count and store number of appearing A, P, E, H, and I
proc conf_counting residue {
    variable s_seq
    variable c_switch
    variable aggregates
    variable cluster_list
    variable 2agg
    variable 3agg
    variable hbonds_energy
    variable hbonds_energy_switch
    variable length_stride

    set count_A 0
    set count_P 0
    set count_M 0
    set count_E 0
    set count_H 0
    set count_I 0; # pi-helix
    set count   0
    set longest_helix 0
    set helix_in_a_row 0
    set number_of_helices 0

    # Loop over list. Increment by 2 steps because of 
    # the structure of the list (array index is useless here)
    for {set k 1} {$k < [llength $residue]} {incr k 2} {
	set comma   [string first "," [lindex $residue [expr $k-1]]]
	set id_aa   [string range [lindex $residue [expr $k-1]] [expr $comma+1] end]

	set chain_k [lindex $residue $k]
	for {set j 0} {$j < [string length $chain_k]} {incr j} {
	    if {[string index $chain_k $j]=="_"} { incr count
	    } elseif {[string index $chain_k $j]=="A"} {incr count_A
		incr count
		if {$s_seq!=""} {
		    lset s_seq [expr $id_aa+1] 1 [expr [lindex [lindex $s_seq [expr $id_aa+1]] 1]+1]
		}
	    } elseif {[string index $chain_k $j]=="P"} {incr count_P
		incr count
		if {$s_seq!=""} {
		    lset s_seq [expr $id_aa+1] 2 [expr [lindex [lindex $s_seq [expr $id_aa+1]] 2]+1]
		}
	    } elseif {[string index $chain_k $j]=="M"} {incr count_M
		incr count
		if {$s_seq!=""} {
		    lset s_seq [expr $id_aa+1] 3 [expr [lindex [lindex $s_seq [expr $id_aa+1]] 3]+1]
		}
	    } elseif {[string index $chain_k $j]=="E"} {incr count_E
		incr count
		if {$s_seq!=""} {
		    lset s_seq [expr $id_aa+1] 4 [expr [lindex [lindex $s_seq [expr $id_aa+1]] 4]+1]
		}
	    } elseif {[string index $chain_k $j]=="H"} {incr count_H
		incr count
		if {$s_seq!=""} {
		    lset s_seq [expr $id_aa+1] 5 [expr [lindex [lindex $s_seq [expr $id_aa+1]] 5]+1]
		}
	    } elseif {[string index $chain_k $j]=="I"} {incr count_I
		incr count
		if {$s_seq!=""} {
		    lset s_seq [expr $id_aa+1] 6 [expr [lindex [lindex $s_seq [expr $id_aa+1]] 6]+1]
		}
	    } else {
		incr count
	    }
	}
    }

    # Loop over list again to find the longest helix
    set sorted_residue_list ""

    for {set l 0} {$l < $length_stride} {incr l} {
	set residue_list ""
	for {set k 1} {$k < [llength $residue]} {incr k 2} {
	    set comma   [string first "," [lindex $residue [expr $k-1]]]
	    set prefix  [string range [lindex $residue [expr $k-1]] 0 [expr $comma-1]]
	    set id_aa   [string range [lindex $residue [expr $k-1]] [expr $comma+1] end]
	    if {$prefix == $l} {
		lappend residue_list "$id_aa [lindex $residue $k]"
	    }
	}
	lappend sorted_residue_list [lsort -integer -index 0 $residue_list]
    }

    foreach prefix $sorted_residue_list {
	foreach residue $prefix {
	    if {[lindex $residue 1]=="H"} {
		incr helix_in_a_row
	    } else {
		if {$longest_helix < $helix_in_a_row} {
		    set longest_helix $helix_in_a_row
		}
		if {$helix_in_a_row > 0} {
		    incr number_of_helices
		}
		set helix_in_a_row 0			
	    }
	}
	if {$longest_helix < $helix_in_a_row} {
	    set longest_helix $helix_in_a_row
	}
    }
    

    # Cluster algorithm results
    if {$c_switch==1} {
	puts "** Cluster analysis ** :"
	puts "Size\tCount"
	set cluster_count 0
	foreach cluster_element $cluster_list {
	    incr cluster_count
	    if {$cluster_element != 0} {
		puts "$cluster_count\t$cluster_element"
	    }
	}
	puts ""
    }


    # Hbonds energy switch
    if {$hbonds_energy_switch==1} {
	puts "** Total energy from hydrogen bonds ** :"
	puts "HBonds_energy : $hbonds_energy"
	puts ""
    }

    set total [expr $count_A+$count_P+$count_M+$count_E]
    # Use the order of magnitude of the total number of aggregates to align the final output.
    set total_mag [format %e $total]
    set total_mag [string range $total_mag [expr [string first + $total_mag]+1] end]
    incr total_mag
    

    if {$s_seq!=""} {
	puts "** Residue conformation analysis ** : "
	puts "   -> column : conformation"
	puts "   -> line   : residue"
	set counter 0
	foreach index $s_seq {
	    if {$counter==0} {
		puts "   Amino acid sequence : $index"
		puts ""
		puts "  [format %[set total_mag]s A] [format %[set total_mag]s P]\
 [format %[set total_mag]s M] [format %[set total_mag]s E] [format %[set total_mag]s H] [format %[set total_mag]s I]"
	    } else {
		set letter_count 0
		foreach letter $index {
		    if { $letter_count == 0 } {
			puts -nonewline "$letter "
			incr letter_count
		    } else {
			puts -nonewline "[format %[set total_mag]s $letter] "
			incr letter_count
		    }
		}
		puts ""
	    }
	    incr counter
	}
	puts ""
    }

    puts "** Aggregation analysis ** :"
    if {$total>0} {
	puts "A ratio : A/\{A,P,M,E\} : $count_A/$total ([format %5.3f [expr 1.*$count_A/($total*1.)]]) - A/total : [format %6.4f [expr 1.*$count_A/($count*1.)]]"
	puts "P ratio : P/\{A,P,M,E\} : $count_P/$total ([format %5.3f [expr 1.*$count_P/($total*1.)]]) - P/total : [format %6.4f [expr 1.*$count_P/($count*1.)]]"
	puts "M ratio : M/\{A,P,M,E\} : $count_M/$total ([format %5.3f [expr 1.*$count_M/($total*1.)]]) - M/total : [format %6.4f [expr 1.*$count_M/($count*1.)]]"
	puts "E ratio : E/\{A,P,M,E\} : $count_E/$total ([format %5.3f [expr 1.*$count_E/($total*1.)]]) - E/total : [format %6.4f [expr 1.*$count_E/($count*1.)]]"
	puts "Aggregation percentage              : [format %6.4f [expr 1.*($count_A+$count_P+$count_M+$count_E)/($count*1.)]]"
	puts "Number of aggregates (2-,3-,4-,...) : $aggregates"
	set total_2 [expr ([lindex $2agg 0]+[lindex $2agg 1]*1.)]
	if {$total_2>0} {
	    puts "Types of 2-aggregates (P/A)         : $2agg ([format %5.3f [expr 1.*[lindex $2agg 0]/([lindex $2agg 0]+[lindex $2agg 1]*1.)]]/[format %5.3f [expr 1.*[lindex $2agg 1]/$total_2]])"
	}
	set total_3 [expr ([lindex $3agg 0]+[lindex $3agg 1]+[lindex $3agg 2]*1.)]
	if {$total_3>0} {
	    puts "Types of 3-aggregates (P/A/M)       : $3agg ([format %5.3f [expr 1.*[lindex $3agg 0]/([lindex $3agg 0]+[lindex $3agg 1] +[lindex $3agg 2]*1.)]]/[format %5.3f [expr 1.*[lindex $3agg 1]/([lindex $3agg 0]+[lindex $3agg 1]+[lindex $3agg 2]*1.)]]/[format %5.3f [expr 1.*[lindex $3agg 2]/$total_3]])"
	}
    } else {
	puts "No aggregation event recorded !"
    }
    puts "\n** Helix analysis ** :"
    puts "Helix count                         : $count_H"
    puts "Helicity percentage                 : [format %6.4f [expr 1.*$count_H/($count*1.)]]"
    puts "Length of the longest helix         : $longest_helix"
    puts "Number of helices                   : $number_of_helices"
    puts "\n** Pi-Helix analysis ** :"
    puts "Pi-helix count                      : $count_I"
    puts "Pi-helix percentage                 : [format %6.4f [expr 1.*$count_I/($count*1.)]]"
}



