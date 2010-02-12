# ::peptideb::output::pdb
#
# output the sequence of amino acids
# to a .PDB file.
#
#

namespace eval peptideb {
    namespace eval output {
	# Writes the configuration in a PDB file.
	# Arguments:
	#   - $coords is an array of the cartesian coordinates
	# of all the CG beads
	# Returns nothing. It jsut creates a .pdb file.
	proc pdb { coords {inter_folder ""}} {
	    # Open file to be read
	    catch {exec mkdir "$peptideb::directory"}
	    if {$inter_folder != ""} {
		set folder [join [list $peptideb::directory $inter_folder] "/"]
	    } else {
		set folder "$peptideb::directory"
	    }
	    catch {exec mkdir $folder}
	    set filename [peptideb::utils::get_pdbfilename $folder]
	    set filePDB [open $filename "w"]
	    
	    
	    # automatize identification of CG beads
	    set beads "{ N  } { CA } { CB } { C  }"
	    set id "{N} {C} {C} {C}"
	    # number of beads per chain
	    set numbeads 4
	    if { $peptideb::display_O_H == "on" } {
		# Model oxygens
		lappend beads " O  "
		lappend id "O"
		incr numbeads
		# Model hydrogens
		set beads [linsert $beads 1 " HN "]
		set id [linsert $id 1 "H"]
		incr numbeads
	    }

	    # Create an array containing all letters of the alphabet. 
	    # This will determine which amino acid corresponds to which peptide chain.
	    # First chain denoted by A. Second by B. etc.
	    set alphabet {A B C D E F G H I J K L M N O P Q R S T U V W X Y Z}

	    set index 0
	    # Loop over peptides
	    for { set l 0 } { $l < [llength $peptideb::amino_acids] } { incr l } {
		# Loop over amino acids
		for { set k 0 } { $k < [llength [lindex $peptideb::amino_acids $l]] } { incr k } {
		    # Make sure the beads exist in the $coords file
		    if { $k > [llength [lindex $coords $l]]  } {
			::mmsg::err [namespace current] "Trying to write the coordinates of an amino acid that hasn't been parsed."
		    }
		    # Loop over CG beads
		    for { set j 0 } { $j < $numbeads } { incr j } {
			incr index
			set x [lindex [lindex [lindex $coords $l] [expr $k*$numbeads+$j]] 0]
			set y [lindex [lindex [lindex $coords $l] [expr $k*$numbeads+$j]] 1]
			set z [lindex [lindex [lindex $coords $l] [expr $k*$numbeads+$j]] 2]
			puts -nonewline $filePDB "[format %4s "ATOM"]"
			puts -nonewline $filePDB "[format %7s $index] "
			puts -nonewline $filePDB "[format %4s [lindex $beads $j]] "
			puts -nonewline $filePDB "[format %3s [lindex [lindex [lindex $peptideb::amino_acids $l] $k] 0]] "
			puts -nonewline $filePDB "[lindex $alphabet [expr $l%26]]"
			puts -nonewline $filePDB "[format %4s [expr $k+1]]    "
			puts -nonewline $filePDB "[format %8s [8charnumber $x]]"
			puts -nonewline $filePDB "[format %8s [8charnumber $y]]"
			puts -nonewline $filePDB "[format %8s [8charnumber $z]]"
			puts -nonewline $filePDB "  1.00  0.00"
			puts -nonewline $filePDB "      [format %3s P[expr $l+1]]   [lindex $id $j]\n"
		    }
		}
	    }
	    # Loop over hfip molecules                                                                                                           
	    if { $peptideb::hfip_flag } {
		set l [llength $peptideb::amino_acids]
		set k 0
		set beads "{ F  } { C  } { F  }"
		set id "{F} {C} {F}"
		set numbeads 3
		set chain [lindex $coords [llength $peptideb::amino_acids]]
		set hfip_id [llength $peptideb::amino_acids]
		foreach hfip_mol $chain {
		    incr index
		    set j 0
		    foreach hfip_atom $hfip_mol {
			set x [lindex $hfip_atom 0]
			set y [lindex $hfip_atom 1]
			set z [lindex $hfip_atom 2]
			puts -nonewline $filePDB "[format %4s "ATOM"]"
			puts -nonewline $filePDB "[format %7s $index] "
			puts -nonewline $filePDB "[format %4s [lindex $beads $j]] "
			puts -nonewline $filePDB "[format %3s IFP] "
			puts -nonewline $filePDB "[lindex $alphabet [expr $l%26]]"
			puts -nonewline $filePDB "[format %4s [expr $k+1]]    "
			puts -nonewline $filePDB "[format %8s [8charnumber $x]]"
			puts -nonewline $filePDB "[format %8s [8charnumber $y]]"
			puts -nonewline $filePDB "[format %8s [8charnumber $z]]"
			puts -nonewline $filePDB "  1.00  0.00"
			puts -nonewline $filePDB "      [format %3s P[expr $l+1]]   [lindex $id $j]\n"
			incr j
		    }			
		    incr k
		}
	    }
	    
	    puts $filePDB "END"
	    close $filePDB

	    ::mmsg::send [namespace current] "The peptide configuration was saved in '$filename'."
	    return
	}

	# Create a startup file for VMD. 
	# This will make sure all the files get loaded correctly
	# as an animation, and the molecule display style can be changed.
	proc pdbstartup {{inter_folder ""}} {
	    # Open file to be read
	    catch {exec mkdir "$peptideb::directory"}
	    if {$inter_folder != ""} {
		set folder [join [list $peptideb::directory $inter_folder] "/"]
	    } else {
		set folder "$peptideb::directory"
	    }
	    catch {exec mkdir $folder}
	    set filename "$folder/vmd_structure.script"
	    set filePDB [open $filename "w"]

	    puts $filePDB "loadseries $peptideb::PDB_file.vmd 1 0"
	    puts $filePDB "rotate stop"
	    puts $filePDB "mol modstyle 0 0 NewCartoon 0.3 15"
	    puts $filePDB "mol modcolor 0 0 Structure"
	    puts $filePDB "mol addrep 0"
	    puts $filePDB "mol modselect 1 0 \"resname IFP\""
	    puts $filePDB "mol modstyle 1 0 Licorice"
	    close $filePDB
	    
	    set filename "$folder/vmd_backbone.script"
	    set filePDB2 [open $filename "w"]

	    puts $filePDB2 "loadseries $peptideb::PDB_file.vmd 1 0"
	    puts $filePDB2 "rotate stop"
	    puts $filePDB2 "mol modstyle 0 0 Licorice 0.3 15 15"
	    puts $filePDB2 "mol modcolor 0 0 name"
	    puts $filePDB2 "mol addrep 0"
	    puts $filePDB2 "mol modselect 1 0 \"resname IFP\""
	    puts $filePDB2 "mol modstyle 1 0 Licorice"
	    close $filePDB2


	    set filename "$folder/vmd_lowres_struct.script"
	    set filePDB3 [open $filename "w"]

	    puts $filePDB3 "loadseries $peptideb::PDB_file.vmd 1 0"
	    puts $filePDB3 "rotate stop"
	    puts $filePDB3 "mol modstyle 0 0 NewCartoon 0.3 6"
	    puts $filePDB3 "mol modcolor 0 0 Structure"
	    puts $filePDB3 "mol addrep 0"
	    puts $filePDB3 "mol modselect 1 0 \"resname IFP\""
	    puts $filePDB3 "mol modstyle 1 0 Licorice"
	    close $filePDB3



	    return 
	}

	# Inserts a 8 character number according to PDB's 
	# rather strict structure...
	# Argument: number
	# Returns string
	proc 8charnumber { number } {
	    set int [expr int(abs($number))]
	    set dec [expr int((abs($number)-int(abs($number)))*1000)]
	    if { $dec < 10 } {
		set dec "00$dec"
	    } elseif { $dec < 100 } {
		set dec "0$dec"
	    } else {
		set dec [string range $dec 0 2]
	    }
	    if { $int < 10 } {
		if { $number < 0. } {
		    set result "  -$int.$dec"
		} elseif { $number > 0.} {
		    set result "   $int.$dec"
		} else {
		    set result "   0.000"
		}
	    } elseif { $int < 100 } {
		if {$number < 0. } {
		    set result " -$int.$dec"
		} elseif { $number > 0. } {
		    set result "  $int.$dec"
		} 
	    } elseif { $int < 1000 } {
		if {$number < 0. } {
		    set result "-$int.$dec"
		} elseif { $number > 0. } {
		    set result " $int.$dec"
		} 
	    } else {
		::mmsg::err [namespace current] "Coordinate is too large to be exported to a PDB file: $number."
	    }
	    return $result
	}
    }
}
