# peptideb::output --
#
# This package outputs the sequence of amino
# acids to different formats. Examples:
# - .PDB file 
# - ESPResSo input
# Author: Tristan
#


package provide peptideb::output 1.0.0
package require peptideb::input
package require mmsg

namespace eval peptideb {
    namespace eval output {
    
    }
}

source [file join [file dirname [info script]] pdb.tcl     ]
source [file join [file dirname [info script]] psf.tcl     ]