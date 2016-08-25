# ::peptideb::utils - warmup
#
# Warmup process
#


namespace eval peptideb {
    namespace eval utils {
	# Warmup routine
	# takes in argument the number of steps and times to call warmup
	proc warmup { steps times args } {
	    ::mmsg::send [namespace current] "Warming up $times times $steps timesteps."
		set startcap 0.01
		set capgoal  1000
	    set usage "Usage : warmup steps times"

	    if { $steps==0 || $times==0 } {
		::mmsg::warn [namespace current] "Warmup steps are zero."
		return
	    }

	    # Cap increment is determined here
	    set capincr [expr ($capgoal - $startcap)/($times*1.)]

	    # Set the initial force cap
	    set cap $startcap

	    for { set i 0 } { $i < $times } { incr i } {
		::mmsg::debug [namespace current] "Warmup step $i of $times at time=[setmd time] (cap=$cap)."
		inter forcecap $cap
		integrate $steps
		set cap [expr $cap + $capincr]

		flush stdout
	    }

	    # Uncapping forces
	    ::mmsg::send [namespace current] "Uncapping forces."
	    inter forcecap individual

      setmd time 0


	    return
	}
    }
}
