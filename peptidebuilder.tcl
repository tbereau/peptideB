#!/usr/bin/env tclsh
#
# Might have to be changed depending on the location
# of tclsh. This can be figured out by typing 
# 'which tclsh'

#
# Main script for peptideB.
#


# Load all the libraries in Tcl.
if { [catch { set lib_path "$env(PEPTIDEB_DIR)" } ] } {
    puts "Error. Can't find the PEPTIDEB_DIR environment variable."
    exit
} else {
    lappend auto_path "$lib_path/src/"
}


# ---- Make sure the different modules are present -------- #
set this [namespace current]
package require mmsg
::mmsg::send $this ""
::mmsg::send $this "               ***********************  "
::mmsg::send $this "              *                       * "
::mmsg::send $this "              *        peptideB       * "
::mmsg::send $this "              *     Tristan Bereau    * "
::mmsg::send $this "              *         (2008)        * "
::mmsg::send $this "              *                       * "
::mmsg::send $this "               ***********************  "
::mmsg::send $this ""
::mmsg::send $this "loaded version [package require peptideb] of peptideb."
::mmsg::send $this "loaded version [package require peptideb::input] of peptideb::input."
::mmsg::send $this "loaded version [package require peptideb::espresso] of peptideb::espresso."
::mmsg::send $this "loaded version [package require peptideb::output] of peptideb::output."
::mmsg::send $this "loaded version [package require peptideb::utils] of peptideb::utils."
::mmsg::send $this "loaded version [package require mmsg] of mmsg."

::mmsg::send $this ""
# ---------------------------------------------------------- #

# ---- Process Command Line Args --------------------------- #
set options {
}

set usage "Usage: 
\tpeptidebuilder.tcl <CONFIG_FILE>
     or 
\tEspresso peptidebuilder.tcl <CONFIG_FILE> 
     * Possible options :
\t-pdb PDB_FILE \t\t load a PBD file
\t-espresso \t\t turns on espresso integration
\t-norandom \t\t turns off the random number generator
\t-script FILE \t\t loads a script instead of the regular program
\t-analysis FILE \t\t runs an analysis script during integration
\t-replica \[-connect\] \t starts a parallel tempering simulation
\t\t\t\t-connect starts a client simulation
\t-hremd \[-connect\] \t starts a Hamilitonian replica exchange simulation
\t\t\t\t-connect starts a client simulation
\t-native FILE \t\t compare the simulated structure to FILE using an order parameter
\t-sidechain \t\t calculate energy contribution coming from sidechains
\t-ragtime CONF_LETTER \t calculate alpha (letter A) or beta (B) content and store in observables.dat
\t\t\t\tletter H will provide hydrogen bond calculation (any combination of A, B, and H can be invoked)
\t-stop PARAMETER_VALUE \t stops the simulation once the order parameter value has reached PARAMETER_VALUE
\t-contact1 FILE1 \t calculate nativity of a peptide with reference to a topology file FILE1
\t-contact2 FILE1 FILE2    Calculates two nativity parameters corresponding to each topology file
\t-resume \t\t resume the simulation where the last checkpoint is (depending on CONFIG_FILE information)
\t-nopdb \t\t\t the simulation performs measurements and store them, but doesn't save any PDB
\t-hpmedium \t\t simulation in a hydrophobic environment
\t-randequilib \t\t choose dihedrals randomly from an equilibrated basin (no high free energy configurations)
\t-virtual_com N \t\t add virtual site: center of mass of peptide number N (starting from 1)
\t-freeze \t\t freeze atoms that have 0.00 occupancy in input PDB
"

# ----------------------------------------------------------- #

if { $argc<1} {
    ::mmsg::err $this "Parameter file missing. $usage"
}

set paramsfile [lindex $argv 0]

namespace eval peptideb {

    set espresso on
    set norandom 0
    set script 0
    variable analysis_file ""
    set native 0
    set ragtime 0
    set replica 0
    set hremd 0
    set contact1 0
    set contact2 0
    set replica_connect 0
    set hremd_connect 0
    set stop 0
    set nopdb 0
    set hp_medium 0
    set rand_equilib 0
    set 2nd_environment 0
    set hfip_frac 0.
    set virtual_com -1
    set freeze 0
    set frozen ""

    # Read the amino acids sequence, plus other parameters. Sequence can be accessed by $amino_acids.
    source $paramsfile

    # Cmdline arguments
    if { $argc >= 2 } {
	for { set k 1 } { $k < $argc } { incr k } {
	    switch -- [lindex $argv $k] "-pdb" {
		set pdb_read [lindex $argv [expr $k+1]]
		::mmsg::send $this "Importing pdb file : $pdb_read."
		incr k
	    } "-espresso" {
		set espresso on
	    } "-norandom" {
		set norandom 1
		::mmsg::send $this "Random seed was not reinitialized."
	    } "-script" {
		set script 1
		set script_name [lindex $argv [expr $k+1]]
		::mmsg::send $this "Reading script $script_name."
		incr k
	    } "-analysis" {
		# Can't have both script and analysis on
		if {$script==1} {
		    mmsg::err $this "Can't run both a custom script and custom analysis."
		}
		set analysis_file [lindex $argv [expr $k+1]]
		::mmsg::send $this "Loading analysis script $analysis_file."
		incr k
	    } "-replica" { 
		set replica 1
		::mmsg::send $this "Replica exchange algorithm has been turned on."
		if {[lindex $argv [expr $k+1]]== "-connect"} {
		    # Read file for hostname
		    while { [catch {set file_hn [open "hostfile" r] } ] } {
			after 1000
		    } 
		    set replica_connect [lindex [read $file_hn] 0]
		    close $file_hn
		    incr k
		} else {
		    if { ![info exists tcp_port] } {
			set tcp_port 12000
		    }
		    # Create file for hostname with its port
		    set file_hn [open "hostfile" w]
		    puts $file_hn "[exec hostname] $tcp_port"
		    close $file_hn
		    # delete any allthere file
		    catch { exec rm -f allthere_$::paramsfile }
		}
	    } "-hremd" {
		set hremd 1
		::mmsg::send $this "Replica exchange algorithm has been turned on."
		if {[lindex $argv [expr $k+1]]== "-connect"} {
		    # Read file for hostname
		    while { [catch {set file_hn [open "hostfile" r] } ] } {
			after 1000
		    } 
		    set hremd_connect [lindex [read $file_hn] 0]
		    close $file_hn
		    incr k
		} else {
		    if { ![info exists tcp_port] } {
			set tcp_port 12000
		    }
		    # Create file for hostname
		    set file_hn [open "hostfile" w]
		    puts $file_hn "[exec hostname] $tcp_port"
		    close $file_hn
		    # delete any allthere file
		    catch { exec rm -f allthere_$::paramsfile }
		}
	    } "-native" {
		set native 1
		set native_file [lindex $argv [expr $k+1]]
		mmsg::send $this "Nativeness order parameter will be calculated."		
		incr k
	    } "-sidechain" {
		set sidechain_energy 1
		mmsg::send $this "Side chain energies will be calculated."
	    } "-ragtime" {
		# Check that the ragtime package exists - should be in PEPTIDEB_DIR/packages
		set ragtime_path "$env(PEPTIDEB_DIR)/packages/ragtime/ragtime.tcl"
		if { ![file exists $ragtime_path ] } {
		    ::mmsg::err $this "Ragtime package wasn't found - should be\n$ragtime_path"
		}
		# Also check that STRIDE exists - should be in the PATH
		if { [ catch { set stride_exists [exec which stride] } ] || $stride_exists=="" } {
		    ::mmsg::err $this "'stride' program wasn't found in your \$PATH variable. Please install & setup. See README for more details."
		}

		set ragtime 1
		set ragtime_letter [string tolower [lindex $argv [expr $k+1]]]
		if {$ragtime_letter != "a" && $ragtime_letter != "b" && $ragtime_letter != "h"} {
		    mmsg::err $this "1-letter code for residue analysis cannot be recognized (use A for alpha-helix, B for beta-sheet, and H for hydrogen bond energy calculation)."
		}
		if {$ragtime_letter == "a"} {
		    set ragtime_A 1
		} elseif {$ragtime_letter == "b"} {
		    set ragtime_B 1
		} elseif {$ragtime_letter == "h"} {
		    set ragtime_H 1
		}
		mmsg::send $this "Residue analysis will be calculated for conformation type [string toupper $ragtime_letter]."
		incr k
	    } "-contact1" {
		# Check that the 'native_contacts' package exists - should be in peptideB/packages/
		set contact_path "$env(PEPTIDEB_DIR)/packages/native_contacts/native_contacts.sh"
		if { ![file exists $contact_path ] } {
		    ::mmsg::err $this "'native_contacts.sh' package wasn't found - should be\n$contact_path"
		}

		set contact1 1
		set contact_file1 [lindex $argv [expr $k+1]]
		incr k
	    } "-contact2" {
		# Check that the 'native_contacts' package exists - should be in peptideB/packages/
		set contact_path "$env(PEPTIDEB_DIR)/packages/native_contacts/native_contacts.sh"
		if { ![file exists $contact_path ] } {
		    ::mmsg::err $this "'native_contacts.sh' package wasn't found - should be\n$contact_path"
		}

		set contact2 1
		# Can't have contact2 on with contact1
		if {$contact1==1} {
		    mmsg::err $this "Can't turn on -contact1 and -contact2."
		}
		set contact_file1 [lindex $argv [expr $k+1]]
		set contact_file2 [lindex $argv [expr $k+2]]
		incr k 2
	    } "-stop" {
		set stop 1
		set stop_thresh [lindex $argv [expr $k+1]]
		::mmsg::send $this "The simulation will be stopped automatically once an order parameter reaches $stop_thresh."
		incr k
	    } "-resume" {
		set resume 1
		::mmsg::send $this "The simulation will be resumed at the last checkpoint."
		incr k
	    } "-nopdb" {
		set nopdb 1
		::mmsg::send $this "No PDB will be written during the simulation."
	    } "-hpmedium" {
		set hp_medium 1
		::mmsg::send $this "Environment set to hydrophobic."
	    } "-randequilib" {
		set rand_equilib 1
		::mmsg::send $this "Dihedrals will be chosen randomly from an equilibrated basin (no high free energy configs.)."
	    } "-virtual_com" {
		set virtual_com [lindex $argv [expr $k+1]]
		::mmsg::send $this "Virtual site will be added at the center of mass of the peptide"
		incr k
	    } "-freeze" {
		set freeze 1
		::mmsg::send $this "Atoms with 0.00 occupancy will be fixed"
	    } default {
		mmsg::err $this "Wrong argument [lindex $argv $k] -- exiting."
	    }
	}
    }
    
    
    # At this point, check whether the user wants a hydrophobic medium.
    # If so, resource the parameter files, and then the configuration file
    if {$hp_medium == 1} {
	namespace eval :: {
	    source [set env(PEPTIDEB_DIR)]/src/chain_parameters.tcl
        }
	source $paramsfile

    }

    # Are we simulating two environments? (HBond strength)
    if { [info exists peptideb::HB_bilayer_dz] } {
	if { $peptideb::HB_bilayer_dz > 0. } {
	    ::mmsg::send $this "z-position analysis"
	    set 2nd_environment 1
	}
    }

    
    if { $norandom==0 && $espresso=="on"} {
	if { [catch {utils::init_random} ] } {
	    set usage "Usage: Espresso peptidebuilder.tcl paramsfile \[-espresso\]"
	    mmsg::err $this "Espresso + peptideB were not started correctly. $usage"
	}
    } 


    if { $espresso=="on" } {
	::mmsg::send $this "Espresso integration is turned on."
    } else {
	::mmsg::send $this "Espresso integration is turned off ($espresso)."
    }

    # Thermostat
    if {[info exists dpd]} {
	if {$dpd==1} {
	    mmsg::send $this "DPD thermostat has been turned on."
	} else {
	    mmsg::send $this "Langevin thermostat has been turned on."
	}
    } else {
	mmsg::send $this "Langevin thermostat has been turned on."
    }

    # HFIP stuff
    set hfip_flag 0
    set hfip_num_mol 0
    if { $hfip_frac > 0.} {
	set hfip_flag 1
	set hfip_num_mol [input::calculate_HFIP_vv_frac $hfip_frac]
    }

    
    # ---------------------------------------------------------- #
    # Allow children namespaces
    if { [catch {::mmsg::setnamespaces ":: [namespace children ::peptideb] [namespace children ::parallel_tempering]"} ] } {
	# No parallel tempering
	::mmsg::setnamespaces ":: [namespace children ::peptideb]"
    }


    # Enable debug messages
    ::mmsg::enable debug


    # ----------------- Start the script ----------------------- #

    if {[info exists native_file]} {
	set native_coords [input::import_pdb $native_file]
	variable native_rij [utils::native_rij $native_coords]
    }

    if { [info exists script_name] } {
	# Load the script $script_name
	source $script_name
    } elseif { $replica == 1 } {
	# Replica exchange (parallel tempering) is turned on
	catch {exec rm -rf $directory}
	catch {exec mkdir $directory}

	# Copy parameter file in folder
	catch {exec cp "$paramsfile" "$peptideb::directory/$paramsfile"}
	
	::mmsg::send $this "Starting parallel tempering at temperatures : $replica_temps"
	
	if {$replica_connect == 0} {
	    # attempt to run PT 3 times if it catches an error
	    set n_tries 0
	    set max_n_tries 3
	    while { $n_tries <= $max_n_tries } {
		::mmsg::send $this "host [exec hostname], port $tcp_port"
		if { [ catch { parallel_tempering::main -values $replica_temps -rounds $replica_rounds \
				   -init peptideb::espresso::replica_init -swap peptideb::espresso::replica_swap \
				   -perform peptideb::espresso::replica_perform -port $tcp_port -info comm } errmsg ] } {
		    ::mmsg::send $this $errmsg
		} else {
		    break
		}
		incr tcp_port
		set file_hn [open "hostfile" w]
		puts $file_hn "[exec hostname] $tcp_port"
		close $file_hn
		incr n_tries
	    }
	} else {
	    # attempt to run PT 3 times if it catches an error
	    set n_tries 0
	    set max_n_tries 2
	    while { $n_tries <= $max_n_tries } {
		::mmsg::send $this "Connection attempt [expr $n_tries+1]/[expr $max_n_tries]..."
		# Attempt to connect 10 seconds after the hostname file has been written
		set long_enough 0
		while { $long_enough == 0 } {
		    if { [file exists hostfile] } {
			after 10000
			set file_hn [open "hostfile" r]
			set data [read $file_hn]
			set data [split $data " "]
			close $file_hn
			set tcp_port [lindex $data 1]
			set long_enough 1
		    }
		}
		::mmsg::send $this "host [exec hostname], port $tcp_port"
		if { [catch { parallel_tempering::main -connect $replica_connect -init peptideb::espresso::replica_init \
				  -swap peptideb::espresso::replica_swap -perform peptideb::espresso::replica_perform \
				  -port $tcp_port -info comm } errmsg ] } {
		    mmsg::send $this $errmsg
		} else {
		    break
		}
		incr n_tries
	    }
	}
    } elseif { $hremd == 1} {
	# Hamiltonian replica exchange is turned on
	catch {exec rm -rf $directory}
	catch {exec mkdir $directory}
	
	# Copy parameter file in folder
	catch {exec cp "$paramsfile" "$peptideb::directory/$paramsfile"}
	
	if { ![info exists tcp_port] } {
	    set tcp_port 12000
	}
	
	::mmsg::send $this "Starting Hamiltonian replica exchange with side chain coupling : $hremd_couplings"
	
	if {$hremd_connect == 0} {
	    # attempt to run HREMD 3 times if it catches an error
	    set n_tries 0
	    set max_n_tries 3
	    while { $n_tries <= $max_n_tries } {
		::mmsg::send $this "host [exec hostname], port $tcp_port"
		if { [catch { parallel_tempering::main -values $hremd_couplings -rounds $hremd_rounds \
				  -init peptideb::espresso::hremd_init -swap peptideb::espresso::hremd_swap \
				  -perform peptideb::espresso::hremd_perform -port $tcp_port -info comm } errmsg ] } {
		    ::mmsg::send $this $errmsg
		} else {
		    break
		}
		incr tcp_port
		set file_hn [open "hostfile" w]
		puts $file_hn "[exec hostname] $tcp_port"
		close $file_hn
		incr n_tries
	    }
	} else {
	    # attempt to run PT 2 times if it catches an error
	    set n_tries 0
	    set max_n_tries 2
	    while { $n_tries <= $max_n_tries } {
		::mmsg::send $this "Connection attempt [expr $n_tries+1]/[expr $max_n_tries]..."
		# Attempt to connect only 2 seconds after the hostname file has been written
		set long_enough 0
		while { $long_enough == 0 } {
		    if { [file exists hostfile] } {
			after 2000
			set file_hn [open "hostfile" r]
			set data [read $file_hn]
			set data [split $data " "]
			close $file_hn
			set tcp_port [lindex $data 1]
			set long_enough 1
		    }
		}
		::mmsg::send $this "host [exec hostname], port $tcp_port"
		if { [catch { parallel_tempering::main -connect $hremd_connect -init peptideb::espresso::hremd_init \
				  -swap peptideb::espresso::hremd_swap -perform peptideb::espresso::hremd_perform \
				  -port $tcp_port -info comm } errmsg ] } {
		    mmsg::send $this $errmsg
		} else {
		    break
		}
		incr n_tries
	    }
	}
    } else {
	# If a pdb file name is given as an argument in the cmd line, import it,
	# rather than construct a sequence of amino_acids
	if { [info exists pdb_read]} {
	    # Import the pdb file.
	    set coords [input::import_pdb $pdb_read]
	    # Add HFIP if turned on
	    if { $hfip_flag } {
		lappend coords [input::create_HFIP_molecules $hfip_frac]
	    }
	} else {
	    # Calculate the cartesian coordinates from the sequence.
	    set coords [input::read_AA_sequence $paramsfile]    
	    # Add HFIP if turned on
	    if { $hfip_flag } {
		lappend coords [input::create_HFIP_molecules $hfip_frac]
	    }
	}


	# Data file for nativeness + energy
	# Open file to be read
	catch {exec mkdir "$peptideb::directory"}
	set filename "$peptideb::directory/observables.dat"
	set file_nat [open $filename w]
	puts $file_nat "\# Time"
	if {[info exists native_file]} {
	    puts $file_nat "\# Nativeness order parameter compared to $native_file."
	} 
	if {[info exists sidechain_energy]} {
	    puts $file_nat "\# Hydrophobicity (energy due to side chain interactions)."
	}
	if {$2nd_environment == 1} {
	    puts $file_nat "\# z-position of the CoM of each peptide in the box."
	}
	if {[info exists ragtime_A]} {
	    puts $file_nat "\# ragtime: Helicity fraction."
	    puts $file_nat "\# ragtime: Length of the longest helix."
	    puts $file_nat "\# ragtime: Number of helices."
	}
	if {[info exists ragtime_B]} {
	    puts $file_nat "\# ragtime: Aggregation/beta-conformation fraction."
	}
	if {[info exists ragtime_H]} {
	    puts $file_nat "\# ragtime: Hydrogen-bond energy."
	}
	if {$contact1==1 || $contact2==1} {
	    puts $file_nat "\# contact1: nativity calculation to reference structure $contact_file1."
	    if {$contact2==1} {
		puts $file_nat "\# contact2: nativity calculation to reference structures $contact_file2."
	    }
	} 
	puts $file_nat "\# Potential energy"
	
	close $file_nat
	
	# If nativeness must be calculated, store information regarding
	# the native structure itself - unless we're resuming a simulation
	if {[info exists native_file]} {
	    if { [info exists resume] } {
		if { [catch {glob -directory $peptideb::directory observables.dat}]} {
		    ::mmsg::err $this "Can't find the 'observables.dat' file"
		} else {
		    ::mmsg::send $this "The 'observables.dat' file was found."
		}
	    }
	}
	
	# Determine if resolution should be added.
	if { ![info exists display_O_H]} {
	    set display_O_H "off"
	} elseif { $display_O_H=="on" || $display_O_H=="off"} {
	    # Do nothing special.
	} else {
	    set display_O_H "off"
	}
	::mmsg::send $this "Chemical resolution was turned $display_O_H."
	
	# Erase all existing pdb files if we're not resuming
	if { ![info exists resume] } {
	    foreach file [glob -nocomplain -directory $directory $PDB_file.vmd*] {
		file delete $file
	    }
	}
	

	# Copy parameter file in folder
	catch {exec cp "$paramsfile" "$peptideb::directory/$paramsfile"}


	# Write a vmd startup file.
	output::pdbstartup
	
	
	if { $espresso == "on"} {
	    # Output the initial configuration to ESPResSo without
	    # adding chemical resolution.
	    espresso::espresso_init $coords
	} else {
	    # Add chemical resolution
	    set coords [input::addresolution $coords]
	    # Write the psf file after adding chemical resolution
	    output::psf $coords
	    # Output the coordinates in a VMD file.
	    output::pdb $coords
	}
    }


    # Delete hostfile for replica exchange
    catch { exec rm -f hostfile }
    # Delete 'allthere' file for replica exchange
    catch { exec rm -f allthere_$::paramsfile }

    # Quit
    ::mmsg::send $this "\n"

    
    exit
}





