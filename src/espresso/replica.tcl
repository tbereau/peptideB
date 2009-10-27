# ::peptideb::espresso
#
# Various routines for the parallel tempering algorithm
#
# Author: Tristan
#

namespace eval peptideb {
    namespace eval espresso {
	# proc replica_init
	# Initializes the system prior to replica exchange.
	# Arguments: - id of the system
	# - temperature temp
	proc replica_init {id temp} {
	    set this [namespace current]
	    variable pdb_counter 0
	    variable init_temp $temp

	    # Create file that acknowledges that all clients have connected
	    if {$id == 0 } {
		set file_allthere [open "$peptideb::directory/../allthere_$::paramsfile" w]
		puts $file_allthere "allthere"
		close $file_allthere
	    }

	    # Create output files for temperature temp
	    mmsg::send $this "Starting replica exchange instance at temperature $temp"

	    # Creating observables data file
	    set folder "$peptideb::directory/temp$temp"
	    catch {exec mkdir $folder}
	    set file_f [join [list $folder /observables.dat ] ""]
	    set f [open $file_f w]
	    puts $f "\# Time"
	    set e_index 2
	    if {[info exists ::peptideb::native_file]} {
		puts $f "\# Nativeness order parameter compared to $peptideb::native_file."
		incr e_index
	    } 
	    if {[info exists ::peptideb::sidechain_energy]} {
		puts $f "\# Hydrophobicity (energy due to side chain interactions)."
		incr e_index
	    }
	    if {$peptideb::2nd_environment == 1} {
		puts $f "\# z-position of the CoM of each peptide in the box."
		incr e_index
	    }
	    if {[info exists ::peptideb::ragtime_A]} {
		puts $f "\# ragtime: helicity analysis."
		incr e_index 3
	    }
	    if {[info exists ::peptideb::ragtime_B]} {
		puts $f "\# ragtime: aggregation/beta-conformation analysis."
		incr e_index
	    }
	    if {[info exists ::peptideb::ragtime_H]} {
		puts $f "\# ragtime: Hydrogen-bond analysis."
		incr e_index
	    }
	    if {$peptideb::contact1==1 || $peptideb::contact2==1} {
		puts $f "\# contact1: nativity calculation to reference structure $peptideb::contact_file1."
		incr e_index
		if {$peptideb::contact2==1} {
		    puts $f "\# contact2: nativity calculation to reference structures $peptideb::contact_file2."
		    incr e_index
		}
	    } 
	    puts $f "\# Potential energy"

	    close $f
	    

	    # Creating histogram data file
	    set folder "$peptideb::directory/temp$temp"
	    catch {exec mkdir $folder}
	    set file_f [join [list $folder /histogram.dat ] ""]
	    set f [open $file_f w]
	    puts $f "\# Histogram of energies at temperature $temp"
	    puts $f "\# Potential Energy \t Number of hits"
	    close $f

	    # Gnuplot script to look at all histograms at once
	    set f [open "$peptideb::directory/histogram.gnu" w]
	    puts $f "\# Gnuplot script to read histograms of all temperatures"
	    puts $f ""
	    puts $f "p " nonewline
	    foreach temperature $peptideb::replica_temps {
		puts $f "\'temp$temperature/histogram.dat\' w histeps" nonewline
		if {$temperature != [lindex $peptideb::replica_temps end]} {
		    puts $f ", " nonewline
		}
	    }
	    close $f

            # Gnuplot script to look at all Energy VS time curves at once
	    set f [open "$peptideb::directory/energy_vs_time.gnu" w]
            puts $f "\# Gnuplot script to read Energy VS time for all temperatures"
            puts $f ""
            puts $f "p " nonewline
            foreach temperature $peptideb::replica_temps {
                puts $f "\'temp$temperature/observables.dat\' u 1:$e_index w histeps" nonewline
                if {$temperature != [lindex $peptideb::replica_temps end]} {
                    puts $f ", " nonewline
                }
            }
            close $f

	    if {[info exists ::peptideb::native_file]} {
		# Gnuplot script to look at all nativeness order parameters at once
		set f [open "$peptideb::directory/nativeness.gnu" w]
		puts $f "\# Gnuplot script to read nativeness (Q) order parameter of all temperatures"
		puts $f ""
		puts $f "p \[\]\[0:1\] " nonewline
		foreach temperature $peptideb::replica_temps {
		    puts $f "\'temp$temperature/observables.dat\' w l" nonewline
		    if {$temperature != [lindex $peptideb::replica_temps end]} {
			puts $f ", " nonewline
		    }
		}
		close $f
	    }


	    # Erase all existing pdb files
	    foreach file [glob -nocomplain -directory $folder $peptideb::PDB_file.vmd*] {
		file delete $file
	    }
	    
	    # Write a vmd startup file.
	    peptideb::output::pdbstartup "temp$temp"

            if { [catch {setmd box_l [lindex $peptideb::setbox_l 0] [lindex $peptideb::setbox_l 1] [lindex $peptideb::setbox_l 2] } ] } {
                set usage "Usage: Espresso peptidebuilder.tcl paramsfile \[-espresso\]"
                mmsg::err $this "Espresso + peptideB were not started correctly. $usage"
            }
            mmsg::send $this "The box size was set to: [setmd box_l]"
	    
            setmd periodic 1 1 1

	    # Integration parameters
            setmd time_step $peptideb::main_time_step
            setmd skin      $peptideb::verlet_skin

	    # Apply the thermostat
            thermostat langevin $temp $peptideb::langevin_gamma
	    
	    if { [info exists peptideb::dpd] && $peptideb::dpd == 1} {
		thermostat off
		thermostat set dpd $temp $peptideb::langevin_gamma $peptideb::dpd_r_cut
	    }


            setmd time 0


	    set bondbroken_test 1
	    while { $bondbroken_test == 1 } {
		# If a pdb file name is given as an argument in the cmd line, import it,
		# rather than construct a sequence of amino_acids
		if { [info exists peptideb::pdb_read]} {
		    # Import the pdb file.
		    set coords [peptideb::input::import_pdb $peptideb::pdb_read]
		} else {
		    # Calculate the cartesian coordinates from the sequence.
		    set coords [peptideb::input::read_AA_sequence $::paramsfile]
		}
		
		# Topology
		create_topology $coords

		# Warmup process - Make sure the peptide doesn't blow up
		if {![catch {::peptideb::utils::warmup $peptideb::warm_steps $peptideb::warm_n_times}]} {
		    set bondbroken_test 0
		}
	    }

	    # Write the psf file after adding chemical resolution
	    if { [catch {glob -directory $peptideb::folder $peptideb::PDB_file.vmd.psf} ] } {
		set outputcoords [::peptideb::input::addresolution $coords]
		::peptideb::output::psf $outputcoords "temp$temp"
	    }
	}

	# proc replica_perform
	# performs MD integration
	# Arguments : - id of system to evaluate
	# - temperature temp
	# Returns nothing
	proc replica_perform {id temp} {
	    variable pdb_counter
	    variable init_temp
	    variable ::peptideb::analysis_file
	    variable folder "$peptideb::directory"
	    variable ::peptideb::nopdb
	    
	    set folder "$peptideb::directory/temp$temp"
	    set file_f [join [list $folder /observables.dat ] ""]
	    set file_h [join [list $folder /histogram.dat   ] ""]

	    # Apply the new thermostat
	    thermostat off
	    thermostat langevin $temp $peptideb::langevin_gamma

	    if { [info exists peptideb::dpd] && $peptideb::dpd == 1} {
		thermostat off
		thermostat set dpd $temp $peptideb::langevin_gamma $peptideb::dpd_r_cut
	    }


	    integrate $peptideb::replica_timestep
	    
	    set label 0
	    foreach t $peptideb::replica_temps {
		if {$t == $init_temp} { break }
		incr label
	    }

            set current_energy [expr [analyze energy total]-[analyze energy kinetic]]
	    set coords [read_topology]

	    # Add chemical resolution and output to PDB
	    if {$pdb_counter % $peptideb::replica_pdb_freq == 0 && $pdb_counter != 0} {
		if {$peptideb::nopdb == 0 || [info exists peptideb::ragtime_letter]
		    || $peptideb::contact1==1 || $peptideb::contact2==1} {
		    set outputcoords [::peptideb::input::addresolution $coords]
		    ::peptideb::output::pdb $outputcoords "temp$temp"
		}
		if { [info exists analysis_file] } {
		    if {$analysis_file != ""} {
			source $analysis_file
		    }
		}
		
		::peptideb::utils::append_obs $coords $file_f $current_energy $label
		::peptideb::utils::write_histogram $file_h $current_energy	       
		
		# In case the user doesn't want any PDB stored but the analysis required it
		if {$peptideb::nopdb == 1 && ([info exists peptideb::ragtime_letter]
					      || $peptideb::contact1 == 1 
					      || $peptideb::contact2 == 1)} {
		    catch {exec rm -f [peptideb::utils::get_pdbfilename $folder]}
		} 
		

		set pdb_counter 0
		incr peptideb::filenumber
	    }
	    
	    incr pdb_counter
	}
	
	# proc replica_swap
	# Calculates the swap probabilities 
	# for a system.
	# Arguments : - id of the system to evaluate
	# - temperature of system 1 : temp1
	# - temperature of system 2 : temp2
	# Returns the two probabilities.
	proc replica_swap {id temp1 temp2} {
	    set epot [expr [analyze energy total] - [analyze energy kinetic]]
	    return "[expr $epot/$temp1] [expr $epot/$temp2]"
	}
    }
}



