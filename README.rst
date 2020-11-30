Automated Crystal-based Electronic Structure Calculator
=======================================================

ACE is a set of tools and Python modules for calculating materials 
propeties with minimun input information, only a cristalographic
information file is needed. It is based on Quantum ESPRESSO suite
for electronic structures calculations and ase library for atomistic
manipulations. This is a project of the Brazilian Synchrotron Light Laboratory (LNLS).

Requirements
------------

* Python_ 3.6 or later
* NumPy_ (base N-dimensional array package)
* SciPy_ (library for scientific computing)
* Ase_ (library for atomistic simulations)
* Quantum-ESPRESSO_ 6.5 or later


Installation
------------

Add ~/ace to your $PYTHONPATH environment variable, add
~/ace/bin and ~/QE/bin to $PATH (assuming ~/ase is where your ACE folder is
and ~/QE is where Quantum ESPRESSO folder is).


* 'add to path temporary' tutorial: 
   - $ export PYTHONPATH =  ~/ace
   - $ export PATH = $PATH: ~QE/bin
   - $ export PATH = $PATH: ~ace/bin
   - $ echo $PATH (to check result)

* 'add to path permanent' tutorial:     
   - add the same command lines to the .bashrc file 
   - run in terminal: $ souce .bashrc (to execut alterations and make them valid)

* allow ace.py as executable tutorial:
   - $ chmod +x ~/ace/bin/ace.py


Execution
------------

$ ace.py <routine> <input_file> <op_mode> <input_var>

- Routine: (optional)
    - 'tc' - calculate superconductivity critical temperature
    - 'ph' - calculate phonons in Gamma
    - 'eos' - calculate equation of state
    - default - 'tc'

- Input file: (mandatory)
    - formats - .cif 

- Inputvar: (optional)  
    - formats: TC.in, PH.in, EOS.in
    - default - internal pre-set variables respectively to the routine chosen

- Operation modes: (optional)
    - 'n' (new) create a new directory for output files
    - 'w' (overwrite) in case of name conflict, overwrite the data of a previous calculation
    - 'c' (continue) continue a calculation interrupted after some steps  
    - default - 'n'


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
.. _Ase: https://listserv.fysik.dtu.dk/mailman/listinfo/ase-users
.. _Quantum-ESPRESSO: https://www.quantum-espresso.org/