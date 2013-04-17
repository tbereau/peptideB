peptideB
========
Tristan Bereau 

Model
-----
This software simulates a peptide model described in:

> T. Bereau, M. Deserno, J Chem Phys 130, 235106 (2009) [Journal][jpc]

Program
-------

peptideB is a peptide builder. From a given sequence of amino acids, it builds a three-dimensional structure of (almost) all beads. The program can build as many peptides as necessary.

The sequence of amino acids must be set in a configuration file. The dihedral angles are optional, if not set random values will be given.

There are two main outputs to the program :

* ESPResSo - one can feed the peptide in the molecular dynamics (MD) software  package. peptideB takes care of defining bonds, interactions, positions,  etc. All the relevant ESPResSo parameters may be set in the peptideB parameter file which will be sourced in the MD package. Then, ESPResSo runs a simulation on the peptide(s).
* VMD - Visual Molecular Dynamics is a GUI that lets one visualize molecules in 2D. peptideB is able to write .pdb & .psf files that are read by VMD.

peptideB is able to combine these two characteristics by applying MD simulations on a peptide and output it in VMD as a function of simulation time.


[jpc]: http://link.aip.org/link/doi/10.1063/1.3152842
[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/1adc9e622abb43b177d23674afdc5238 "githalytics.com")](http://githalytics.com/tbereau/peptideB)
