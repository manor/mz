from distutils.core import setup
import py2exe

#
# This import is here because of:
# http://stackoverflow.com/questions/15478230/pack-a-software-in-python-using-py2exe-with-libiomp5md-dll-not-found
#
import numpy

setup(console=['liberate.py'],options={"py2exe":{"includes":["pythoncom"]}})