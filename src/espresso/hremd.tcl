# ::peptideb::espresso
#
# Various routines for the Hamiltonian replica exchange algorithm
#
# Author: Tristan
#

namespace eval peptideb {
    namespace eval espresso {
	# proc hremd_init
	# Initializes the system prior to Hamiltonian replica exchange.
	# Arguments: - id of the system
	# - side chain interaction coupling 'sc_coupling' (number between 1. and 0.)
	proc hremd_init {id sc_coupling} {
	    set this [namespace current]
	    variable pdb_counter 0
	    variable init_coupling $sc_coupling

	    # Create file that acknowledges that all clients have connected
	    if {$id == 0 } {
		set file_allthere [open "$peptideb::directory/../allthere_$::paramsfile" w]
		puts $file_allthere "allthere"
		close $file_allthere
	    }

	    # Create output files for side chain coupling 'coupling'
	    mmsg::send $this "Starting Hamiltonian replica exchange instance at coupling $sc_coupling"

	    # Creating observables data file
	    set folder "$peptideb::directory/coupling$sc_coupling"
	    catch {exec mkdir $folder}
	    set file_f [join [list $folder /observables.dat ] ""]
	    set f [open $file_f w]
	    puts $f "\# Time"
	    # Where is the side chain energy located
	    set e_index 2
	    if {[info exists ::peptideb::native_file]} {
		puts $f "\# Nativeness order parameter compared to $peptideb::native_file."
		incr e_index
	    } 
	    # Always output side chain energy for hremd
	    puts $f "\# Hydrophobicity (energy due to side chain interactions)."
	    set peptideb::sidechain_energy 1

	    if {$peptideb::2nd_environment == 1} {
		puts $f "\# z-position of the CoM of each peptide in the box."
	    }
	    if {[info exists ::peptideb::ragtime_A]} {
		puts $f "\# ragtime: helicity analysis."
	    }
	    if {[info exists ::peptideb::ragtime_B]} {
		puts $f "\# ragtime: aggregation/beta-conformation analysis."
	    }
	    if {[info exists ::peptideb::ragtime_H]} {
		puts $f "\# ragtime: Hydrogen-bond analysis."
	    }
	    if {$peptideb::contact1==1 || $peptideb::contact2==1} {
		puts $f "\# contact1: nativity calculation to reference structure $peptideb::contact_file1."
		if {$peptideb::contact2==1} {
		    puts $f "\# contact2: nativity calculation to reference structures $peptideb::contact_file2."
		}
	    } 
	    puts $f "\# Potential energy"

	    close $f
	    

	    # Creating histogram data file
	    set folder "$peptideb::directory/coupling$sc_coupling"
	    catch {exec mkdir $folder}
	    set file_f [join [list $folder /histogram.dat ] ""]
	    set f [open $file_f w]
	    puts $f "\# Histogram of energies at coupling $sc_coupling"
	    puts $f "\# Side Chain Energy \t Number of hits"
	    close $f

	    # Gnuplot script to look at all histograms at once
	    set f [open "$peptideb::directory/histogram.gnu" w]
	    puts $f "\# Gnuplot script to read histograms of all side chain couplings"
	    puts $f ""
	    puts $f "p " nonewline
	    foreach instance $peptideb::hremd_couplings {
		puts $f "\'coupling$instance/histogram.dat\' w histeps" nonewline
		if {$instance != [lindex $peptideb::hremd_couplings end]} {
		    puts $f ", " nonewline
		}
	    }
	    close $f

	    # Gnuplot script to look at all Side Chain Energy VS time curves at once
	    set f [open "$peptideb::directory/energy_vs_time.gnu" w]
	    puts $f "\# Gnuplot script to read Side Chain energy VS time for all side chain couplings"
	    puts $f ""
	    puts $f "p " nonewline
	    foreach instance $peptideb::hremd_couplings {
		puts $f "\'coupling$instance/observables.dat\' u 1:$e_index w histeps" nonewline
		if {$instance != [lindex $peptideb::hremd_couplings end]} {
		    puts $f ", " nonewline
		}
	    }
	    close $f
	    

	    if {[info exists ::peptideb::native_file]} {
		# Gnuplot script to look at all nativeness order parameters at once
		set f [open "$peptideb::directory/nativeness.gnu" w]
		puts $f "\# Gnuplot script to read nativeness (Q) order parameter for all side chain couplings"
		puts $f ""
		puts $f "p \[\]\[0:1\] " nonewline
		foreach instance $peptideb::hremd_couplings {
		    puts $f "\'coupling$instance/observables.dat\' w l" nonewline
		    if {$instance != [lindex $peptideb::hremd_couplings end]} {
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
	    peptideb::output::pdbstartup "coupling$sc_coupling"

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
            thermostat langevin $peptideb::systemtemp $peptideb::langevin_gamma
	    
	    if { [info exists peptideb::dpd] && $peptideb::dpd == 1} {
		thermostat off
		thermostat set dpd $peptideb::systemtemp $peptideb::langevin_gamma $peptideb::dpd_r_cut
	    }

	    # The correct coupling will be set here
	    set peptideb::hp_coupling $sc_coupling
    	    
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
		::peptideb::output::psf $outputcoords "coupling$sc_coupling"
	    }
	}

	# proc hremd_perform
	# performs MD integration
	# Arguments : - id of system to evaluate
	# - side chain coupling 'sc_coupling'
	# Returns nothing
	proc hremd_perform {id sc_coupling} {
	    variable pdb_counter
	    variable init_coupling
	    variable ::peptideb::analysis_file
	    variable folder "$peptideb::directory"
	    variable ::peptideb::nopdb
	    
	    set folder "$peptideb::directory/coupling$sc_coupling"
	    set file_f [join [list $folder /observables.dat ] ""]
	    set file_h [join [list $folder /histogram.dat   ] ""]

	    # Set new coupling
	    set peptideb::hp_coupling $sc_coupling
	    # Now apply new nonbonded interactions
	    set peptideb::nb_interactions ""
	    namespace eval :: {
		source [set env(PEPTIDEB_DIR)]/src/chain_interactions.tcl
            }
	    set_nb_interactions $peptideb::nb_interactions
	    
	    integrate $peptideb::hremd_timestep
	    
	    set label 0
	    foreach c $peptideb::hremd_couplings {
		if {$c == $init_coupling} { break }
		incr label
	    }

            set current_energy [expr [analyze energy total]-[analyze energy kinetic]]
	    set coords [read_topology]

	    set sc_energy 0.
	    for {set i 10} {$i<30} {incr i} {
		for {set j $i} {$j<30} {incr j} {
		    set sc_energy [expr $sc_energy + [analyze energy nonbonded $i $j]]
		}
	    }


	    # Add chemical resolution and output to PDB
	    if {$pdb_counter % $peptideb::hremd_pdb_freq == 0 && $pdb_counter != 0} {
		if {$peptideb::nopdb == 0 || [info exists peptideb::ragtime_letter]
		    || $peptideb::contact1==1 || $peptideb::contact2==1} {
		    set outputcoords [::peptideb::input::addresolution $coords]
		    ::peptideb::output::pdb $outputcoords "coupling$sc_coupling"
		}
		if { [info exists analysis_file] } {
		    if {$analysis_file != ""} {
			source $analysis_file
		    }
		}
		
		::peptideb::utils::append_obs $coords $file_f $current_energy $label
		::peptideb::utils::write_histogram $file_h $sc_energy	       
		
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
	
	# proc hremd_swap
	# Calculates the swap probabilities 
	# for a system.
	# Arguments : - id of the system to evaluate
	# - side chain coupling of system 1 : coupling1
	# - side chain coupling of system 2 : coupling2
	# Returns the two probabilities.
	proc hremd_swap {id coupling1 coupling2} {
	    # At this point, calculate the total side chain energy
	    # it's the total nonbonded energy between side chain particles
	    # NOTE: side chain-backbone interactions are not included in the sum.
	    # These interactions will stay as strong for any coupling since they are
	    # required to obtain physical Ramachandran plots.
	    # indices i and j only sum on side chain particles.
	    set sc_energy 0.
	    for {set i 10} {$i<30} {incr i} {
		for {set j $i} {$j<30} {incr j} {
		    set sc_energy [expr $sc_energy + [analyze energy nonbonded $i $j]]
		}
	    }
	    return "[expr $sc_energy*$coupling1] [expr $sc_energy*$coupling2]"
	}
    }
}


 
