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


from distutils.core import setup

import os

from mz import __version__, _get_home

packages = ['mz',
            'mz.API',
            'mz.report',
            'mz.report.formats',
            'mz.tools']

mz_data = ['COPYING','COPYING.LESSER','unimod.xml']

data_files = ['mz/unimod.xml']

setup(
    name='mz',
    version=__version__,
    author='James Webber, Jignesh Parikh, Manor Askenazi',
    author_email='manor.askenazi@gmail.com',
    packages=packages,
    package_data={ '': ['LICENSE.txt'],
                   'mz': mz_data },

    data_files = [(os.path.join(_get_home(), '.mz'), data_files)],
    url='https://github.com/manor/mz',
    description='"Lightweight" Mass Spectrometry/Proteomics Toolkit'
)
