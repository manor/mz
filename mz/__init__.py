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

__author__ = 'Jignesh Parikh, James Webber, Manor Askenazi'
__version__ = '1.2.0'

__all__ = ['API', 'Functions', 'Tools', 'Report',
           'myHome', 'myData', 'logger_message']

import logging
import os


# # Based on matplotlib's _get_home function (in matplotlib/__init__.py)
# def _get_home():
#     """
#     Find user's home directory if possible. Otherwise raise an error.
#     """
#
#     path = ''
#
#     try:
#         path=os.path.expanduser("~")
#     except:
#         pass
#
#     if not os.path.isdir(path):
#         for evar in ('HOME', 'USERPROFILE', 'TMP'):
#             try:
#                 path = os.environ[evar]
#                 if os.path.isdir(path):
#                     break
#             except:
#                 pass
#
#     if path:
#         return path
#     else:
#         raise RuntimeError('No home environment variable found')


class Logger(object):
    '''A class for the mz logger.'''

    def __init__(self, name, fmt, level=30):
        self._logger = logging.getLogger(name)
        logger_handler = logging.StreamHandler()
        logger_handler.setFormatter(logging.Formatter(fmt))
        self._logger.addHandler(logger_handler)
        self._logger.setLevel(level)

    def __call__(self, level=30, message=''):
        """Sets a new logger message at a specified level

        Logger messages are displayed in a console if they have a higher level than the current logger level.

        Example:
        If current logger level = 30, any message greater than 30 will be displayed
        >>> settings.logger_level = 30
        >>> logger_message(40, 'This message will be displayed because the message level of 40 > 30')
        """
        self._logger.log(level, str(message))

    def set_level(self, level):
        self._logger.setLevel(level)


# can be called as a function, e.g. mz.logger_message(30, 'Message at level 30')
logger_message = Logger('mz', '%(message)s')


# # these directories should be cross-platform:
#
# # user's home directory
# myHome = _get_home()
#
# # multiplierz data folder
# myData = os.path.join(myHome, ".mz")
# myTemp = os.path.join(myData, "TEMP")
#
#
# # if the .multiplierz folder does not exist we will create it
# if not os.path.exists(myData):
#     logger_message(50, 'Could not find .mz folder, creating one...')
#     os.mkdir(myData)
#
# if not os.path.exists(myTemp):
#     logger_message(50, 'Could not find .mz%sTEMP folder, creating one...' % os.sep)
#     os.mkdir(myTemp)
