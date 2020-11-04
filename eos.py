#!/usr/bin/env python3
"""
##------------------------------------------------------------------------------
# CNPEM - National Brazilian Center for Research in Energy and Materials
# Sirius - Research Group EMA
#
# Code project: TC_caculus
# 
# Objectif:
# Automate the use of Quantum Espresso (QE) for the calculation of 
# Superconductivity Critical Temperature (Tc) of diferent molecules and 
# cell structures
# 
##------------------------------------------------------------------------------

"""
import time
from controllers import ControllerProgram

start_time = time.time()

EOS = ControllerProgram(clmode = 'EOS')
EOS.calculate()

print('\nTotal execution time:', time.time()-start_time)