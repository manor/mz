__author__ = 'manor'

#Copyright 2012 Manor Askenazi
#This is a "light" version of the multiplierz platform.
#
# This module is used to generate source code for mod.py which is a hardcoded version of unimod.py
# The idea is to eliminate the need for lxml except for true unimod.xml or mzML access...
#
f = open("mod.py",'w')

from mz import logger_message
from unimod import UnimodDatabase
import os

myData = os.path.dirname(os.path.abspath(__file__))

unimod = None
try:
    # load the unimod database
    unimod = UnimodDatabase(os.path.join(myData, 'unimod.xml'))
except IOError:
    # this shouldn't prevent multiplierz in general from functioning,
    # although of course anything that relies on unimod will be broken
    logger_message(50, ('Warning: %s not found.\n'
                        'Please download a copy from http://www.unimod.org\n'
                        'or use the copy on your local Mascot server')
    % (os.path.join(myData, 'unimod.xml')))
    import sys
    sys.exit()

AW = unimod.elements # atomic masses
amino_acids = unimod.amino_acids # amino acid masses
name = unimod.mods()
delta = {}
for n in name:
    delta[n] = unimod.get_mod_delta(n)
sites = {}
for n in name:
    #The list,set,list move is just to reduce duplicates if they occur (as I do not quite understand the specificity "groups")...
    sites[n] = list(set([s for specs in unimod.get_mod_specificities(n).values() for (s,p) in specs]))
neutral_losses = {}
for n in name:
    for spec in sites[n]:
        neutral_losses[(n,spec)] = unimod.get_mod_neutral_loss(n,spec)

def export(name,stream):
    print >> stream, "%s = %s ;" % (name,repr(globals()[name]))

export("AW",f)
export("amino_acids",f)
export("name",f)
export("delta",f)
export("sites",f)
export("neutral_losses",f)

f.close()

