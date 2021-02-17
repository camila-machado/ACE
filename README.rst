Automated Crystal-based Electronic Structure Calculator
=======================================================

ACE is a set of tools and Python modules for calculating materials 
propeties with minimun input information, only a cristalographic
information file is needed. It is an interface to the Quantum ESPRESSO suite
for electronic structures calculations and uses the ase library for atomistic
manipulations. This is a project of the Brazilian Synchrotron Light Laboratory (LNLS).

Requirements
------------

* Quantum-ESPRESSO_ 6.5 or later
* Python_ 3.6 or later

Python libraries:

* NumPy_ (base N-dimensional array package)
* SciPy_ (library for scientific computing)
* Ase_ (library for atomistic simulations)
* pickleDB_ (library for key-value store on JASON files)

Installation
------------

For performing the calculations, ACE and Quantum ESPRESSO bin directories must
be added to your PATH environment variable so that the contents of both can be
called by the system at any working directory. Also, the ACE folder must be
added to your $PYTHONPATH environment variable and the bin/ace.py file must be
set as an executable.

Assuming that ~/ase is your ACE folder and ~/QE your Quantum ESPRESSO folder, run on your teminal:

.. code-block:: shell

    export PYTHONPATH=$PYTHONPATH:~/ACE/ace
    export PATH=$PATH:~/qe-6.7/bin
    export PATH=$PATH:~/ACE/ace/bin
    # Note that you should replace ~/QE and ~/ace by the
    # actual Quantum-ESPRESSO and ACE paths on your system.
    
This will update PATH and PYTHONPATH for your current session. For making the update permanent, add lines above to your ~/.bashrc file and run:

.. code-block:: shell

    source ~/.bashrc

For checking the results:

.. code-block:: shell

    echo $PATH
    echo $PYTHONPATH

Finally, allow ace.py to be an executable:

.. code-block:: shell

    chmod +x ~/ace/bin/ace.py

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
