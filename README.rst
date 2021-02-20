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
    export PATH=$PATH:~/QE/bin
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

After including the ACE directories to your system PATH and setting ace.py as an executable, the program may be called from anywhere on the system by:

.. code-block::

    ace.py <routine> <input_file> <op_mode> <input_var>
    
With the up to three parameters are:

- <routine> (optional)
    The calculation you want to perform. Currently available values are:
        - tc - calculate superconductivity critical temperature
        - ph - calculate phonons in Gamma
        - eos - calculate equation of state
        Default - tc
    
- <input_file>: (mandatory)
    File in .cif format containing the desired crystal structure
        Default - None

- <input_var>: (optional)  
    Input file containing calculation parameters. Every calculation has its own input file type (TC.in, PH.in, EOS.in), but all lf than have the same. format
        Default - If no <input_var> file is used, ACE will mount its own input file based on internally defined default parameters. Such parameters are usually enough for quick caltulations. The <input_var> file thus generated is stored in the output folder and may be modified for performing another calculation.

- Operation modes: (optional)
    Defines behaviour in terms of previously performed calculation with the same .cif file. Currently available options are:
        - n - (new) create a new directory for output files
        - w - (overwrite) in case of name conflict, overwrite the data of a previous calculation
        - c - (continue) continue a calculation interrupted after some steps  
        Default - n
    
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
.. _Ase: https://wiki.fysik.dtu.dk/ase/
.. _pickleDB: https://pythonhosted.org/pickleDB/
.. _Quantum-ESPRESSO: https://www.quantum-espresso.org/
