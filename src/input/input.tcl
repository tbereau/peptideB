# peptideb::input --
#
# This package reads a protein configuration 
# either from a list of amino acids + angles
# or from a PDB file (later...).
# Author: Tristan
#


package provide peptideb::input 1.0.0
package require mmsg
package require peptideb::utils

namespace eval peptideb {
    namespace eval input {

    }
}

source [file join [file dirname [info script]] misc.tcl      ]
source [file join [file dirname [info script]] resolution.tcl]
source [file join [file dirname [info script]] secondary.tcl ]
source [file join [file dirname [info script]] import.tcl    ]