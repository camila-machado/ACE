Copyright 2020, Camila Machadod de Araújo

LNLS - Brazilian Synchrotron Light Laboratory
DCMC - Scientific division of codensed matter
EMA beamline

Objectif
--------

Automate the use of Quantum Espresso suite (QE) for eletronic-structure
calculations of diferent crystals, specially the calculus of critical 
superconductivity temperature


Automated Crystal-based ELectronic Structure Calculator
=======================================================

ACE is a set of tools and Python modules for calculating materials 
propeties with minimun input information, only a cristalographic
information file is needed. It is based on Quantum ESPRESSO suite
for electronic structures calculations and ase library for atomistic
manipulations.

Requirements
------------

* Python_ 3.6 or later
* NumPy_ (base N-dimensional array package)
* SciPy_ (library for scientific computing)
* Ase_ (library for atomistic simulations)
* Quantum-ESPRESSO 6.5 or later


Installation
------------

Instructions in how to use TC_calculus

->Unpack .tar.gz
->Open terminal in unpacked folder 
->Add Quantum Espresso and program.py to Linux path(*)
->Enable program.py as executable(**)

call in terminal:

$ program.py <directory/input_file> <command> <directory/input_var>

- Program:
    'tc.py' - calculate superconductivity critical temperature
    'phonon.py' - calculate phonons in Gamma
    'eos.py' - calculate equation of state
- Input file format: '.cif'
- Inputvar format: '.TC.in', 'PH.in',  'EOS.in'
- Commands:
    'n' - (new) create a new directory for output files
    'w' - (overwrite) in case of name conflict, overwrite the data of a 
          previous calculation
    'c' - (continue) continue a calculation interrupted after some steps  

obs1: Command is optional, the default operation mode is 'n' (new)
obs2: Inputvar file is optional, the default variables are defined in
      the program

#Enjoy your TC's!!!

(*) You can run the commands in shell Linux environment:
$ export PATH = $PATH:<quantum_espresso_bin_dir>
$ export PATH = $PATH:<Program.py_dir>
$ export PATH = $PATH:<Write.py_dir>
$ echo $PATH (to check result)
obs: temporary solution, valid only in current shell

or

add the command lines to the .bashrc file for permanet solution, and run:
$ souce .bashrc (to execut alterations and make them valid)

(**) Run the commands:
$ chmod +x calculations_qe.py
$ chmod +x controllers.py
$ chmod +x models_program.py
$ chmod +x models_qe.py
$ chmod +x tc.py
$ chmod +x tools.py


Testing
-------


Contact
-------

* Mailing: camila.araujo@lnls.br

Please send us bug-reports, patches, code, ideas and questions.


Example
-------



.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: http://matplotlib.org/
.. _ase-users: https://listserv.fysik.dtu.dk/mailman/listinfo/ase-users
.. _Quantum-ESPRESSO: https://www.quantum-espresso.org/