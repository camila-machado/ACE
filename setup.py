#!/usr/bin/env python3

"set environment variables"

import os
import sys

from os.path import join
from os.path import isdir

project = 'ace'
packages = ['bin', 'gui', 'util', 'core']

proj_repository = os.path.join(os.getcwd(),project)

assert os.path.basename(proj_repository) == project, 'This is not the setup directory'

for pack in packages:
    assert isdir(join(proj_repository, pack)), 'Missing {pack} package'.format(pack = pack)

import ace 