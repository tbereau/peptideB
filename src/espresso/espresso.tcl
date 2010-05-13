# ::peptideb::espresso
#
# Feed a peptide chain in ESPResSo
#
# Author: Tristan
#

package provide peptideb::espresso 1.0.0
package require mmsg
package require peptideb::input

namespace eval peptideb {
    namespace eval espresso {
	# Initialize an espresso simulation
	# from the peptide chain we built with the script.
	# This routine will be called by the main script as well as
	# the parameter optimizer. Different paths.
	# Argument : the list of coordinates.
	# Returns nothing.
	proc espresso_init {coords} {
	    variable ::peptideb::analysis_file
	    variable folder "$peptideb::directory"
	    # variable analysis is only used as a generic variable
	    # to save information between calls to the analysis routine.
	    variable analysis ""

	    set this [namespace current]
	    mmsg::send $this "Feeding peptide chain into ESPResSo..."
	    
	    if { [catch {setmd box_l [lindex $peptideb::setbox_l 0]\
			     [lindex $peptideb::setbox_l 1]\
			     [lindex $peptideb::setbox_l 2] } ] } {
		set usage "Usage: Espresso peptidebuilder.tcl paramsfile \[-espresso\]"
		mmsg::err $this "Espresso + peptideB were not started correctly. $usage"
	    }
	    mmsg::send $this "The box size was set to: [setmd box_l]"

	    setmd periodic 1 1 1

	    # Topology
	    create_topology $coords


	    # Integration parameters
	    setmd time_step $peptideb::main_time_step
	    setmd skin      $peptideb::verlet_skin
	    thermostat langevin $peptideb::systemtemp $peptideb::langevin_gamma
	    
	    if { [info exists peptideb::dpd] && $peptideb::dpd == 1} {
		thermostat off
		thermostat set dpd $peptideb::systemtemp $peptideb::langevin_gamma $peptideb::dpd_r_cut
		mmsg::send $this "Thermostat is [thermostat]"
	    }

	    # Regular integration 
	    set loop_index 0

	    # Resume the simulation ?
	    if { [info exists peptideb::resume] } {
		foreach file [glob -directory $peptideb::directory $peptideb::PDB_file.vmd*.pdb] {
		    lappend file_list $file
		}
		# Load the latest file
		set file_list [lsort -increasing $file_list]
		set last_file [lindex $file_list end]
		
		# Now load the configuration in the program
		::mmsg::send $this "Loading file $last_file..."
		set coords [::peptideb::input::import_pdb $last_file 1]

		set string_length [string length $last_file]
		set loop_index [string range $last_file [expr $string_length - 8] [expr $string_length - 5]]

		set peptideb::filenumber $loop_index
	    } 

	    set timingstart [clock clicks -milliseconds]		

	    # Warmup process
	    ::peptideb::utils::warmup $peptideb::warm_steps $peptideb::warm_n_times
	    

	    # Write the psf file after adding chemical resolution
	    if { [catch {glob -directory $peptideb::directory $peptideb::PDB_file.vmd.psf} ] } {
		set outputcoords [::peptideb::input::addresolution $coords]
		::peptideb::output::psf $outputcoords
	    }
	    
	    ::mmsg::send $this "Main integration..."
	    # Integration
	    for { set k $loop_index } { $k < $peptideb::int_n_times } { incr k } {
		set current_energy [expr [analyze energy total]-[analyze energy kinetic]]
		# Skip first data point if we're resuming a simulation.
		# In this case, we want to integrate before writing to the file.
		if { $k != $loop_index || $k == 0 } {
		    # Write PDB file either if the user asked for it or if
		    # some of the analysis requires a file, in which case it will
		    # be deleted right after (see end of the loop).
		    if {$peptideb::nopdb == 0 || [info exists peptideb::ragtime_letter]
			|| $peptideb::contact1==1 || $peptideb::contact2==1} {
			set pdb_flag 1
			if {[info exists peptideb::pdb_freq]} {
			    set pdb_flag 0
			    if {[expr $k%$peptideb::pdb_freq==0]} {
				set pdb_flag 1
			    }
			}
			if {$pdb_flag==1} {
			    # Add chemical resolution and output to PDB
			    set outputcoords [::peptideb::input::addresolution $coords]
			    ::peptideb::output::pdb $outputcoords
			    incr peptideb::filenumber
			}
		    }
		    
		    if {$analysis_file != ""} {
			source $analysis_file
		    }
		    
		    ::peptideb::utils::append_obs $coords "$peptideb::directory/observables.dat" $current_energy
		    
		    # In case the user doesn't want any PDB stored but the analysis required it
		    if {$peptideb::nopdb == 1 && ([info exists peptideb::ragtime_letter]
						  || $peptideb::contact1 == 1 
						  || $peptideb::contact2 == 1)} {
			catch {exec rm -f [peptideb::utils::get_pdbfilename $peptideb::directory]}
		    } 
		}
		

		# Integrate in ESPResSo
		integrate $peptideb::int_steps
		
		# Recreate the new $coords file from the topology.
		set coords [read_topology]
		
	    }
	    
	    set timingcurr [clock clicks -milliseconds]
	    ::mmsg::send $this "elapsed time: [expr $timingcurr - $timingstart] ms."

	    return
	}
    }
}

source [file join [file dirname [info script]] topology.tcl]
source [file join [file dirname [info script]] replica.tcl ]
source [file join [file dirname [info script]] hremd.tcl ]

