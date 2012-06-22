#Copyright 2012 Manor Askenazi
#This is a "light" version of the multiplierz platform.
#
# Copyright 2008 Dana-Farber Cancer Institute
# multiplierz is distributed under the terms of the GNU Lesser General Public License
#
# This file is part of multiplierz.
#
# multiplierz is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# multiplierz is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with multiplierz.  If not, see <http://www.gnu.org/licenses/>.



"""Core multiplierz functions and classes

The mzTools package contains functions, classes and modules that enable scripting
against multiplierz reports. In general, mzTools can be used for the following:
1. Interacting with reports
2. Interacting with a user for inputs, such as text input and file browsing
3. Creating reports by automatically downloading/formatting Mascot search results / Genbank data
4. Appending precursor information such as peak area to search results (using mzAPI)
5. Generating plots and images, and inserting them in reports (XLS and MZD only)

To use mzTools classes and functions, simply include the following line at the beginning of your script:
import mzTools
>>> import multiplierz.mzTools as mzTools

mzTools functions can be called in the following manner:
mzTools.function(parameters)

mzTools classes can be instantiated in the following manner:
instance = mzTools.class(parameters)

To call the methods of the instantiated class:
instance.method(parameters)

To use a sub-module of the package:
import multiplierz.mzTools.module as module
module.function(parameters)

Examples provided within each class
Some examples require file input and outputs.
These files can be found in the install directory in the 'Sample Files' folder.
The variable 'samples_dir' refers to the 'Sample Files' folder.

"""

__author__ = 'James Webber, Jignesh Parikh'

__all__ = ['mz_image', 'multifile', 'precursor_peaks']


import os

from tempfile import mkstemp
from collections import defaultdict

# import these so they'll always be present in mzTools
import mz_image
import multifile
import precursor_peaks

import mz.report as Report

from mz import myHome, myData, myTemp, logger_message
from mz.genbank import genbank_report

