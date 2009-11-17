#!/usr/bin/env tclsh
#
# *** RAGTIME ***
#
#
# Protein aggregation analysis script.
# Relies on STRIDE to determine the secondary
# structure of proteins, parse the result, and 
# provides further analysis.
# Uses scalar product to determine relative 
# orientation of beta-sheets.
# 
# Script is robust against high (>26) number of peptides.
# New -c option to do cluster analysis.
#
# New: analysis of hbonds.
#
# @author: Tristan BEREAU (bereau@cmu.edu)
# @version 1.4
# @date: 09/22/2009


# Source the routines file
source [file join [file dirname [info script]] ragtime_functions.tcl  ]


set usage "Usage : ragtime \[OPTIONS\] pdb_file
           -a           Only display aggregation related conformations
           -bSize       use periodic boundary condition with box size: Size
           -c           Cluster analysis
           -e           Only display helix conformations
           -f           stride error messages will be ignored
           -h           Display this usage message
           -o           Contact order parameter
           -y           Hydrogen bonds
           -s           Statistics on residues (only for identical chains)
           -v           Verbose mode
           Any number of pdb files can be used as input.
"

if {$argc<1} {
    puts "PDB file is missing.\n$usage"
    exit 1
}

if {[lindex $argv 0]=="-h" && [lindex $argv 0]=="--help"} {
    puts $usage
    exit
}

variable a_conf 0
variable b_PBC 0
variable b_PBCsize 0.
variable c_switch 0
variable o_switch 0
variable o_values ""
variable e_conf 0
variable force 0
variable failed 0
variable s_seq ""
# hbonds: - first element is activation switch
#         - second element is list of hbond energies (one for each PDB)
variable hbonds_energy_switch 0
variable hbonds_energy ""
variable verbose 0
variable alphabet {A B C D E F G H I J K L M \
		       N O P Q R S T U V W X Y Z}
variable number_aa ""
variable CA_coords ""
variable CA_threshold 9.0;# contact distance between residues [in Angstrom]
variable aggregates [list 0 0 0]
variable cluster_list ""
# 2agg counts P and A 2-aggregates
variable 2agg [list 0 0]
# 3agg counts P, A, and Mixed 3-aggregates
variable 3agg [list 0 0 0]

if {$argc>=1} {
    for { set k 0 } { $k < $argc } { incr k } {
	if {[string index [lindex $argv $k] 0]=="-"} {
	    for { set j 1 } { $j < [string length [lindex $argv $k]] } { incr j } {
		switch -- [string index [lindex $argv $k] $j] "a" {
		    set a_conf 1
		} "b" {
		    set b_PBC 1
		    set b_PBCsize [string range [lindex $argv $k] 2 end]
		} "c" {
		    set c_switch 1
		} "e" {
		    set e_conf 1
		} "f" {
		    set force 1
		} "h" {
		    puts $usage
		    exit
		} "o" {
		    set o_switch 1
		} "y" {
		    set hbonds_energy_switch 1
		} "s" { 
		    set s_seq 1
		} "v" {
		    set verbose 1
		}
	    }
	}
    }
}

set index 0
set res_array ""
#########################################
foreach pdb $argv {
    if {[string index $pdb 0]!="-"} {
	# (re)initialize variable
	set failed 0
	# Call STRIDE and save output to file
	call_stride $pdb
	if {$failed==0} {
	    # Read PDB line by line to determine the number of
	    # amino acids per chain. Also, store CA positions.
	    read_pdb $pdb
	    # Optional hydrogen bond energy analysis.
	    if {$hbonds_energy_switch==1} {
		lappend hbonds_energy [hbond_energy $pdb]
	    }
	    # Contact order parameter
	    if {$o_switch == 1} {
		lappend o_values [contact_order_parameter]
	    }
	    # Read file into memory and delete file
	    set data [keep_only_beta $pdb]
	    # Create (or append) an array containing the STRIDE data
	    set res_array [create_res_array [lindex $data 2] $index $res_array]
	    
	    if {$c_switch==1} {
		cluster_analysis [lindex $data 0]
	    }

	    # Analyze hydrogen bonds from STRIDE output and throw irrelevant ones
	    set hbonds [analyze_hbonds $data]
	    # Analyzes extended conformations and assign
	    # Antiparallel or Parallel sheet
	    set res_array [sheet_analysis $hbonds $res_array $index $pdb]
	    
	    # Store data in array
	    array set residue $res_array
	    
	    incr index
	}
    }
}

if {$index>0} {
    if {$verbose} {
	# Display array
	parray residue
    }

    if {$o_switch == 1} {
	set o_average 0.
	foreach value $o_values {
	    set o_average [expr $o_average + $value]
	}
	set o_average [expr $o_average / [llength $o_values]]
	puts "Contact order parameter:\t $o_average"
	puts ""
    }
    # Store number of appearing A, P, E, and H
    set confs [conf_counting [array get residue]]
} else {
    puts "Input PDB is missing!"
}




exit



