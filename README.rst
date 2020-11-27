Copyright 2020, Camila Machadod de Araújo

LNLS - Brazilian Synchrotron Light Laboratory

DCMC - Scientific division of codensed matter

EMA beamline

Project objectif: Automate the use of Quantum Espresso suite (QE) for 
eletronic-structure calculations of diferent crystals, specially the 
calculus of critical superconductivity temperature


Automated Crystal-based Electronic Structure Calculator
=======================================================

ACE is a set of tools and Python modules for calculating materials 
propeties with minimun input information, only a cristalographic
information file is needed. It is based on Quantum ESPRESSO suite
for electronic structures calculations and ase library for atomistic
manipulations.

Requirements
------------

* Python 3.6 or later
* NumPy (base N-dimensional array package)
* SciPy (library for scientific computing)
* Ase (library for atomistic simulations)
* Quantum-ESPRESSO 6.5 or later


Installation
------------

Instructions in how to use ACE
->Unpack .tar.gz
->Open terminal in unpacked folder 
->Add Quantum Espresso and ace/bin to Linux path(*)
->Enable ace/bin/ace.py as executable(**)

call in terminal:

$ ace.py <routine> <input_file> <op_mode> <input_var>

- Routine: (optional)
    'tc' - calculate superconductivity critical temperature
    'ph' - calculate phonons in Gamma
    'eos' - calculate equation of state
    default - 'tc'

- Input file: (mandatory)
    formats - .cif 

- Inputvar: (optional)  
    - TC.in 
    - PH.in  
    - EOS.in
    default - internal pre-set variables respectively to the routine chosen

- Operation modes: (optional)
    'n' - (new) create a new directory for output files
    'w' - (overwrite) in case of name conflict, overwrite the data of a 
          previous calculation
    'c' - (continue) continue a calculation interrupted after some steps  
    default - 'n'
    
#Enjoy your TC's!!!

(*) You can run the commands in shell Linux environment:
$ export PATH = $PATH:<quantum_espresso_bin_dir>
$ export PATH = $PATH:<ace_bin_dir>
$ echo $PATH (to check result)
obs: temporary solution, valid only in current shell

or

add the command lines to the .bashrc file for permanet solution, and run:
$ souce .bashrc (to execut alterations and make them valid)

(**) Run the commands:
$ chmod +x bin/ace.py

Testing
-------
(in construction)

Contact
-------

* Mailing: camila.araujo@lnls.br

Please send us bug-reports, patches, code, ideas and questions.

Example
-------
(in construction)

.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: http://matplotlib.org/
.. _ase-users: https://listserv.fysik.dtu.dk/mailman/listinfo/ase-users
.. _Quantum-ESPRESSO: https://www.quantum-espresso.org/
