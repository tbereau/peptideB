# Messages
#
#  Utility for printing messages of various types
#

package provide mmsg 1.0.0
namespace eval mmsg {
    # The following variables determine what mmsg will print.  To
    # disable printing from a particular module delete it from the
    # list of allowable namespaces or vice versa if you want to
    # enable.  To enable printing of error warnings etc then set the
    # value of the appropriate variable to 1 otherwise set it to zero.
    variable allowablenamespaces { :: ::peptideb::input ::peptideb::output ::peptideb::utils }
    variable enablesend 1
    variable enableerr 1
    variable enablewarn 1
    variable enabledebug 1

    # For keeping track of whether we are printing again on the
    # sameline or not.
    variable sameline 0

    # Proceedure for specifying which namespaces to allow output from.
    # This might not work on all systems (eg macintosh) if multiple copies
    # of mmsg get loaded for some reason.
    proc setnamespaces {namespacelist} {

	variable allowablenamespaces 
	set allowablenamespaces $namespacelist

    }

    # Procedure for enabling messages from of particular types
    proc enable {type} {
	variable enablesend 
	variable enableerr 
	variable enablewarn 
	variable enabledebug 

	switch $type {
	    "send" {
		set enablesend 1
	    }
	    "err" {
		set enableerr 1
	    }
	    "warn" {
		set enablewarn 1
	    }
	    "debug" {
		set enabledebug 1
	    }
	    "default" {
		puts [namespace current ] "no message type called $type"
	    }
	}
	return
    }


    # Procedure for disabling messages from of particular types
    proc disable {type} {
	variable enablesend 
	variable enableerr 
	variable enablewarn 
	variable enabledebug 

	switch $type {
	    "send" {
		set enablesend 0
	    }
	    "err" {
		set enableerr 0
	    }
	    "warn" {
		set enablewarn 0
	    }
	    "debug" {
		set enabledebug 0
	    }
	    "default" {
		puts [namespace current ] "no message type called $type"
	    }
	}
	return
    }

    # Wrapper for the ::mmsg::print command which prints messages with no prefix
    proc send {namespace string {newline "yes"} } {
	variable enablesend
	if { $enablesend } {
	    print $namespace $string "" $newline
	}
    }

    # Wrapper for the ::mmsg::print command which prints messages with the
    # prefix Error and which calls exit
    proc err {namespace string {newline "yes"} } {
	variable enableerr
	if { $enableerr } {
	    print $namespace $string "Error: " $newline -force
	}
	exit
    }

    # ::mmsg::warn --
    #
    # Wrapper for the ::mmsg::print command which prints messages with the
    # prefix Warning 
    #
    proc warn {namespace string {newline "yes"} } {
	variable enablewarn
	if { $enablewarn } {
	    print $namespace $string "Warning: " $newline
	}
    }

    # Wrapper for the ::mmsg::print command which prints debug messages
    proc debug {namespace string {newline "yes"} } {
	variable enabledebug
	if { $enabledebug } {
	    print $namespace $string "Debug: " $newline
	}
    }

    # Check the namespace provided against the list of allowable
    # namespaces return 1 if there is a match or 0 if not
    proc checknamespace {namespace allowable} {
	foreach name $allowable {
	    if { $name == $namespace } {
		return 1
	    }
	}
	return 0
    }

    # Print a message provided it is from an allowable namespace
    proc print {namespace string prefix newline args} {
		variable allowablenamespaces
		variable sameline
		set options {
			{force "override namespace restrictions"}
		}
		set usage "Usage: print \[force] "
		
	    set namespace_OK [checknamespace $namespace $allowablenamespaces]

		if { $namespace_OK } {
			if { $newline != "yes" } {
				if { $sameline } {
					puts -nonewline  "$string"
				} else {
					puts -nonewline  "$namespace > $prefix $string"
					set sameline 1
				}
			} else {
				if { $sameline } {
					puts  "$string"
					set sameline 0
				} else {
					puts  "$namespace > $prefix $string"
				}
			}
		}
    }
	
}
