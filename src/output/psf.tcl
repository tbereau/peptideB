# ::peptideb::output::psf
#
# output the sequence of amino acid
# to a .PSF file.
#
#

namespace eval peptideb {
    namespace eval output {
		# Writes the topology of the chain in a PSF file.
		# Argument: $coords
		proc psf { coords {inter_folder ""}} {
			# Open file to be read
			catch {exec mkdir "$peptideb::directory"}
			if {$inter_folder != ""} {
				set folder [join [list $peptideb::directory $inter_folder] "/"]
			} else {
				set folder "$peptideb::directory"
			}
			catch {exec mkdir $folder}
			set filename "$folder/$peptideb::PDB_file.vmd.psf"
			set filePSF [open $filename "w"]
			
			# Number of beads will be useful
			set nbeads 0
			set npeptides 0
			for { set l 0 } { $l < [llength $coords] } { incr l } {
				set nbeads [expr $nbeads + [llength [lindex $coords $l]]]
				incr npeptides
			}

			# automatize identification of CG beads
			set beads "{ N  } { CA } { CB } { C  }"
			# number of beads per chain
			set numbeads 4
			if { $peptideb::display_O_H == "on" } {
				# Oxygens
				lappend beads " O  "
				incr numbeads
				# Hydrogens
				set beads [linsert $beads 1 " HN "]
				incr numbeads
			}

			# Atoms
			puts $filePSF "PSF"
			puts $filePSF "[format %8s $nbeads] !NATOM"
			set index 0
			set index_aa 0
			# Loop over peptides
			for { set l 0 } { $l < [llength $peptideb::amino_acids] } { incr l } {
				# Loop over amino acids
				for { set k 0 } { $k < [llength [lindex $peptideb::amino_acids $l]] } { incr k } {
					# Make sure the beads exist in the $coords file
					if { $k > [llength [lindex $coords $l]]  } {
						::mmsg::err [namespace current] "Trying to write the coordinates of an amino acid that hasn't been converted."
					}
					incr index_aa
					# Loop over CG beads
					for { set j 0 } { $j < $numbeads } { incr j } {
						incr index
						puts -nonewline $filePSF "[format %8s $index] "
						puts -nonewline $filePSF "[format %-4s "P[expr $l+1]"] "
						puts -nonewline $filePSF "[format %4s $index_aa] "
						puts -nonewline $filePSF "[format %3s [lindex [lindex [lindex $peptideb::amino_acids $l ] $k] 0]] "
						puts -nonewline $filePSF "[format %4s [lindex $beads $j]] "
						puts -nonewline $filePSF "[format %4s [lindex $beads $j]] "
						puts $filePSF "                                     "
					}
				}
			}
			
			# Bonds
			puts $filePSF "[format %8s [expr $nbeads - $npeptides]] !NBOND"
			set index 0
			set line_check 1
			for { set l 0 } { $l < [llength $coords] } { incr l } {
				set index_char 1
				incr index
				while { $index_char < [llength [lindex $coords $l]] } {
					if { $peptideb::display_O_H == "off" } {
						# Case where neither O nor H is modeled.
						if { $index_char % 4 == 3 } {
							# Special case : Cb does NOT bind with C', but instead Ca does.
							puts -nonewline $filePSF "[format %8s [expr $index-1]]"
							puts -nonewline $filePSF "[format %8s [expr $index+1]]"
						} else {
							puts -nonewline $filePSF "[format %8s $index]"
							puts -nonewline $filePSF "[format %8s [expr $index+1]]"
						}
					} else {
						# Case where both H and O are modeled.
						if { $index_char % 6 == 0 || $index_char % 6 == 2 || $index_char % 6 == 4 } {
							# Special case : No Cb-C', or O-N, or H-Ca.
							puts -nonewline $filePSF "[format %8s [expr $index-1]]"
							puts -nonewline $filePSF "[format %8s [expr $index+1]]"
						} else {
							puts -nonewline $filePSF "[format %8s $index]"
							puts -nonewline $filePSF "[format %8s [expr $index+1]]"
						}
					}
					if { [expr $line_check % 4] == 0 } {
						puts $filePSF ""
					}
					incr index
					incr index_char
					incr line_check
				}
			}
			close $filePSF

			return
		}
    }    
}
