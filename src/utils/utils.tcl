# peptideb::utils --
#
# Contains useful routines that may be 
# useful anywhere in the code.
# Author: Tristan
#


package provide peptideb::utils 1.0.0
package require mmsg

namespace eval peptideb {
    namespace eval utils {

		# Initialize the random number generators on each processor
		proc init_random { } {
			
			set c [expr [clock seconds] - 1068130800]
			#    set c  1068130800
			
			# The number of processors is given by Espresso
			set n_procs [setmd n_nodes]
			
			::mmsg::send [namespace current]  "Setting the random seed to clock in seconds: $c "
			set seedbase $c
			for { set i 1 } { $i < $n_procs } { incr i } {
				lappend seedbase [expr $c + 2304*$i]
			}
			eval t_random seed $seedbase
			
			flush stdout
		}

		# Append observables to a file
		# coords are the coordinates of the chain
		# filename is the name of the file !
		proc append_obs {coords filename current_energy {id -1}} {

			set f [open $filename a]
			puts $f "[setmd time] " nonewline
			if { [info exists peptideb::native_file] } {
				variable ::peptideb::native_rij
				set integrated_rij [peptideb::utils::native_rij $coords]
				set param_Q [peptideb::utils::calculate_Q $integrated_rij $peptideb::native_rij]
				puts $f "\t$param_Q" nonewline
				stop_simulation $param_Q $f
			} 
			if { [info exists peptideb::sidechain_energy] } {
				# Calculate total energy due to side chain interactions
				set sc_energy 0.
				for {set i 10} {$i<30} {incr i} {
					for {set j $i} {$j<30} {incr j} {
						set sc_energy [expr $sc_energy + [analyze energy nonbonded $i $j]]
					}
				}
				puts $f "\t$sc_energy" nonewline
			}
			# Calculate z position of CoM when simulating two different environments
			if { $peptideb::2nd_environment == 1 } {
			    # loop over all chains
			    for { set l 0 } { $l < [llength $coords]} { incr l } {
				set chain [lindex $coords $l]
				set CoM [center_of_mass $chain]
				puts $f "\t[lindex $CoM 2]" nonewline
			    }
			}
			# 
			# RAGTIME - HELICITY
			#
			if { [info exists peptideb::ragtime_A] } {
				set folder [string range $filename 0 [expr [string last "/" $filename]-1]]
				set pdbfile [get_pdbfilename $folder]
				set param1 0.
				set param2 0
				set param3 0
				# run ragtime. it's possible we get an error, e.g. there's no hbonds. catch it.
				if {![catch {set output [exec $peptideb::ragtime_path $pdbfile]} errmsg]} {
					# Helicity percentage
					set aggr_list [lsearch -all $output "percentage"]
					foreach position $aggr_list {
						if {[lindex $output [expr $position-1]]=="Helicity"} {
							set param1 [lindex $output [expr $position+2]]
						}
					}
					# Length of the longest helix
					set word_length [lsearch -all $output "Length"]
					foreach position $word_length {
						if {[lindex $output [expr $position+4]]=="helix"} {
							set param2 [lindex $output [expr $position+6]]
						}
					}
					# Number of helices
					set word_length [lsearch -all $output "Number"]
					foreach position $word_length {
						if {[lindex $output [expr $position+2]]=="helices"} {
							set param3 [lindex $output [expr $position+4]]
						}
					}
				} else {
					::mmsg::send [namespace current] "ragtime encountered an error:"
					::mmsg::send [namespace current] $errmsg
					::mmsg::send [namespace current] "The parameter will be set *Arbitrarily* to 0."
					set param1 0.
				}
				puts $f "\t$param1\t$param2\t$param3" nonewline
				stop_simulation $param1 $f
			} 
			#
			# RAGTIME - SHEET PROPENSITY
			#
			if { [info exists peptideb::ragtime_B] } {
				set folder [string range $filename 0 [expr [string last "/" $filename]-1]]
				set pdbfile [get_pdbfilename $folder]
				set param 0.
				# run ragtime. it's possible we get an error, e.g. there's no hbonds. catch it.
				if {![catch {set output [exec $peptideb::ragtime_path $pdbfile]} errmsg]} {
					# Aggregation percentage
					set aggr_list [lsearch -all $output "percentage"]
					foreach position $aggr_list {
						if {[lindex $output [expr $position-1]]=="Aggregation"} {
							set param [lindex $output [expr $position+2]]
						}
					}
				} else {
					::mmsg::send [namespace current] "ragtime encountered an error:"
					::mmsg::send [namespace current] $errmsg
					::mmsg::send [namespace current] "The parameter will be set *Arbitrarily* to 0."
					set param 0.
				}
				puts $f "\t$param" nonewline
				stop_simulation $param $f
			}
			#
			# RAGTIME - HYDROGEN-BOND ENERGY
			#
			if { [info exists peptideb::ragtime_H] } {
				set folder [string range $filename 0 [expr [string last "/" $filename]-1]]
				set pdbfile [get_pdbfilename $folder]
				set param 0.
				# run ragtime. it's possible we get an error, e.g. there's no hbonds. catch it.
				# add HBonds option
				if {![catch {set output [exec $peptideb::ragtime_path "-y" $pdbfile]} errmsg]} {
					# HBond energy
					set hbondE_list [lsearch -all $output "HBonds_energy"]
					foreach position $hbondE_list {
						set param [lindex $output [expr $position+2]]
					}
				} else {
					::mmsg::send [namespace current] "ragtime encountered an error:"
					::mmsg::send [namespace current] $errmsg
					::mmsg::send [namespace current] "The parameter will be set *Arbitrarily* to 0."
					set param 0.
				}
				puts $f "\t$param" nonewline
				stop_simulation $param $f
			} 
			if { $peptideb::contact1==1 || $peptideb::contact2==1} {
				set folder [string range $filename 0 [expr [string last "/" $filename]-1]]
				set pdbfile [get_pdbfilename $folder]
				# run 'native_contacts' script for FILE1
				set param ""
				if {![catch {set output [exec $peptideb::contact_path $pdbfile $peptideb::contact_file1]} errmsg]} {
					set param [lindex $output 0]
				} else {
					::mmsg::send [namespace current] "$script_name encountered an error:"
					::mmsg::send [namespace current] $errmsg
					exit 1
				}
				puts $f " \t$param" nonewline
				stop_simulation $param $f
				# run native contacts script for optional FILE2
				if {$peptideb::contact2==1} {
					set param ""
					if {![catch {set output [exec $peptideb::contact_path $pdbfile $peptideb::contact_file2]} errmsg]} {
						set param [lindex $output 0]
					} else {
						::mmsg::send [namespace current] "$script_name encountered an error:"
						::mmsg::err $errmsg
					}
					puts $f " \t$param" nonewline
					stop_simulation $param $f
				}
			}
			puts $f " \t$current_energy" nonewline
			if { $id != -1 } {
				puts $f " \t$id"
			} else {
				puts $f "" 
			}
			close $f
		}


		# Stop simulation once the order parameter has reached a certain value
		# Delete all PDB files inside the simulation's root directory		
		proc stop_simulation {param output_file} {
			if {[info exists peptideb::stop_thresh]} {
				if {$peptideb::stop==1 && $param>=$peptideb::stop_thresh} {
					close $output_file
					::mmsg::send [namespace current] "Order parameter has reached threshold value ($param)."
					::mmsg::send [namespace current] "Deleting all .pdb files..."
					catch { eval exec rm -f [glob [set peptideb::directory]/*.pdb] }
					catch { eval exec rm -f [glob [set peptideb::directory]/*/*.pdb] }
					exit
				}
			}
			return		   
		}

		# Read histogram file from a given temperature.
		# Takes filename in argument. Returns histogram as
		# an array
		proc read_histogram {filename} {
			array set histogram ""
			set f [open $filename r]
			set data [read $f]
			set data [split $data "\n"]
			close $f
			foreach line $data {
				if {[lindex $line 0] != "\#" && [llength $line] > 0} {
					set histogram([lindex $line 0]) [lindex $line 1]
				}
			}
			return [array get histogram]
		}

		# Return the current PDB filename
		proc get_pdbfilename {{folder ""}} {
			if {$folder == ""} {
				set folder $peptideb::directory
			}
			set fnumber 0
			if {$peptideb::filenumber<10} {
				set fnumber "000$peptideb::filenumber"
			} elseif {$peptideb::filenumber<100} {
				set fnumber "00$peptideb::filenumber"
			} elseif {$peptideb::filenumber<1000} {
				set fnumber "0$peptideb::filenumber"
			} else {
				set fnumber $peptideb::filenumber
			}
			return "[set folder]/$peptideb::PDB_file.vmd$fnumber.pdb"
		}							 

		# Write histogram of energies to file.
		# The new energy is an argument, as well as the filename.
		# Returns nothing.
		proc write_histogram {filename current_energy} {
			array set histogram [read_histogram $filename]
			set f [open $filename w]
			# Add $current_energy to histogram variable.
			# The bins will have a precision of 1.
			set value [expr floor($current_energy)]
			if { [catch {incr histogram($value)} ] } {
				set histogram($value) 1
			}
			
			puts $f "\# Histogram of energies"
			puts $f "\# Potential energy \t Number of hits"
			foreach idx [array names histogram] {
				puts $f "$idx \t[set histogram($idx)]"
			}
			close $f
		}
    }
}

source [file join [file dirname [info script]] math.tcl  ]
source [file join [file dirname [info script]] warmup.tcl]
source [file join [file dirname [info script]] native.tcl]



