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
    
Up to three inputs (described bellow) are supplied after calling ace.py. Their order does not matter.

- <routine> (optional)
    The calculation you want to perform. Currently available values are:
        - tc - calculate superconductivity critical temperature
        - ph - calculate phonons in Gamma
        - eos - calculate equation of state
        Default - tc
    
- <op_mode> (optional)
    Defines behaviour in terms of previously performed calculation with the same .cif file. Currently available options are:
        - n - (new) create a new directory for output files
        - w - (overwrite) in case of name conflict, overwrite the data of a previous calculation
        - c - (continue) continue a calculation interrupted after some steps  
        Default - n
        
- <input_file>: (mandatory)
    File in .cif format containing the desired crystal structure
        Default - None

- <input_var>: (optional)  
    Input file containing calculation parameters. Every calculation has its own input file type (.TC.in, .PH.in or .EOS.in), but all lf than have the same. format
        Default - If no <input_var> file is used, ACE will mount its own input file based on internally defined default parameters. Such parameters are usually enough for quick caltulations. The <input_var> file thus generated is stored in the output folder and may be modified for performing another calculation.
        
When performing the calculation, a new directory named after the .cif file used and the type of calculation performed will be created on the same directory as you run ACE. After the calculation completion, such directory will contain:
    - .in file containing all the calculation parameters used. THis is the file that might be used as <input_var> input
    - .db file, which is used internally by ACE not relevant for the user since it contains the same information present on the .in file
    - /calc directory containing internal calculation files of Quantum ESPRESSO
    - /input directory containing input files for the respective Quantum ESPRESSO calculation. These files might also be used for performing the same calculaiton, or a similar one by slighttly modifying the files, directly by running Quantum ESPRESSO, without ACE
    - /output directory containing the output Quantum ESPRESSO files and a 'result'file containing the usually most relevant result of the calculation
    
*Exemple*

A user may run:

.. code-block::

    ace.py ph Si.cif
    
For performing a Gamma point phonon calculation (ph) of silicon using the structure file Si.cif as input. In such case, the standard behaviour 'n' will be used and a new calculation will be performed, regardless of any previous calculation.

This will create a Si_PH directory (or Si_PH_new, if a previous calculation was performed, since we are not overwriting 'w' nor continuing 'c' a prevous calculation). Such directory will contain a Si.PH.in file, a Si_config.db file and the directories /calc, /input and /output. The /output directory will contain a result file with the phonon frequencies and base vectors of the normal mode.

If a user wants to repeat this calculation with, for example, a denser K-mesh, the Si.PH.in file option "kpoints_div' must be modified from [9,9,9] to a denser mesh, such as [12,12,12]. The modified file must be placed on the same directory as the original .cif file and the calculation is repeated by:

.. code-block::

    ace.py ph Si.cif Si.PH.in

Pseudopotentials
----------------

This package accompanies PAW and Ultrasoft non-relativistic pseudopotentials with exchange-correlation functionals of type GGA (ace/pseudo directory). By default, ultrasoft pseudopotentials are used.

All the included pseudopotentials are from PSlibrary (DOI: 10.1016/j.commatsci.2014.07.043, WEB: http://www.quantum-espresso.org/pseudopotentials, LICENSE: GNU General Public License (version 2 or later))

Please cite the pseudopotentials used and give proper credit to their authors. More infomration on the PSlibrary library, including citation information, may be found at: https://dalcorso.github.io/pslibrary/

Testing
-------
(in construction)

Contact
-------

* Mailing: camila.araujo@ee.ufcg.edu.br ; lucas.francisco@lnls.br

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
