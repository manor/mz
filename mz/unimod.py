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

# this module is for accessing the unimod database
# (which we'll include with Multiplierz)

# customizations are a bit tricky--they require modifying the unimod.xml file.
# the easiest way to do that is to use a Mascot server's configuration editor
# to add the mods to its own unimod file, and then copy that file into the
# multiplierz directory.

import lxml.etree as ET
from collections import defaultdict

# namespace dictionary for the unimod schema. 'f' is just to distinguish
# the 'lower' function, using the Blais URL is not really necessary
NS = {'umod': 'http://www.unimod.org/xmlns/schema/unimod_2',
      'f': 'http://blais.dfci.harvard.edu'}

def lower(context, a):
    return a[0].lower()

class LookupError(Exception):
    '''Exception thrown when a modification does not exist. Also
    thrown if a mod does not have a given specificity.'''

    def __init__(self, mod, site=None):
        self.mod = mod
        self.site = site

    def __str__(self):
        if self.site:
            return 'Modification %s at site %s was not found' % (self.mod, self.site)
        else:
            return 'Modification %s was not found' % self.mod

class UnimodDatabase(object):
    '''Class to represent the Unimod XML file, with some methods
    for looking up modification deltas, specificities, etc, as well
    as things like element masses and amino acid masses. All masses
    are monoisotopic by default.'''

    def __init__(self, file_name):
        self.file_name = file_name

        # this is a special class for evaluating XPath expressions quickly.
        # we have to give it the extra function 'lower' because lxml doesn't
        # have xpath 2.0 functions (yet)
        self.utree = ET.XPathEvaluator(ET.parse(self.file_name),
                                       namespaces=NS,
                                       extensions={(NS['f'], 'lower'): lower})

    @property
    def elements(self):
        '''Gets element monoisotopic masses'''

        elem_list = self.utree('.//umod:elements/umod:elem[@title][@mono_mass]')

        return dict((elem.get('title'), float(elem.get('mono_mass'))) for elem in elem_list)

    @property
    def amino_acids(self, mass='monoisotopic'):
        '''Gets amino acid monoisotopic or average masses. Also includes the N and C terminals.

        'mass' should be 'monoisotopic' or 'average'. Defaults to monoisotopic mass.'''

        if mass == 'monoisotopic':
            aa_list = self.utree('.//umod:amino_acids/umod:aa[@title][@mono_mass]')

            return dict((aa.get('title'), float(aa.get('mono_mass'))) for aa in aa_list)
        elif mass == 'average':
            aa_list = self.utree('.//umod:amino_acids/umod:aa[@title][@avge_mass]')

            return dict((aa.get('title'), float(aa.get('avge_mass'))) for aa in aa_list)
        else:
            raise ValueError("%s is not a valid argument for mass type, "
                             "should be 'monoisotopic' or 'average'." % mass)

    def get_amino_acid_composition(self, aa):
        '''Gets the amino acid composition. 'aa' should be the single-letter abbreviation
        for an amino acid.'''

        aa_list = self.utree('.//umod:amino_acids/umod:aa[@title="%s"]' % aa)

        if aa_list:
            return dict((e.get('symbol'), int(e.get('number')))
                        for e in aa_list[0].xpath('.//umod:element[@symbol][@number]',
                                                  namespaces=NS))
        else:
            return None

    def mod_exists(self, mod_name, site=None):
        '''Searches for a mod by name--which name applies depends on the mod
        (see unimod.org help for details). Returns True is the mod is found.

        'site' is an optional argument--if given, the function only returns true
        if the modification exists at the given site (AA code or 'N|C-term').
        '''

        mod_list = self.utree('.//umod:mod[f:lower(@title) = "%s"]' % mod_name.lower())

        if mod_list:
            if site:
                spec_list = [mod.xpath('./umod:specificity[@site="%s"]' % site, namespaces=NS)
                             for mod in mod_list]
                return any(spec_list)
            else:
                return mod_list[0].get('title')
        else:
            return False

    def mods(self):
            """Returns all possible mods.
            """
            mod_list = self.utree('.//umod:mod')
            return map(lambda x: x.get('title'), mod_list)

    def get_mod_delta(self, mod_name, mass='monoisotopic'):
        '''Searches for a mod by name--which name applies depends on the mod
        (see unimod.org help for details). Returns monoisotopic or average mass.

        'mass' should be 'monoisotopic' or 'average'. Defaults to monoisotopic mass.'''

        if mass == 'monoisotopic':
            mod_list = self.utree(('.//umod:mod[f:lower(@title) = "%s"]/'
                                   'umod:delta[@mono_mass]') % mod_name.lower())

            if mod_list:
                return float(mod_list[0].get('mono_mass'))
            else:
                raise LookupError(mod_name)
        elif mass == 'average':
            mod_list = self.utree(('.//umod:mod[f:lower(@title) = "%s"]/'
                                   'umod:delta[@avge_mass]') % mod_name.lower())

            if mod_list:
                return float(mod_list[0].get('avge_mass'))
            else:
                raise LookupError(mod_name)
        else:
            raise ValueError("%s is not a valid argument for mass type, "
                             "should be 'monoisotopic' or 'average'." % mass)

    def get_mod_specificities(self, mod_name):
        '''Searches for a mod by name and returns a dictionary of specificities, organized by
        group number. Keys are integers and values are lists of AAs'''

        spec_list = self.utree(('.//umod:modifications/'
                                'umod:mod[f:lower(@title) = "%s"]/'
                                'umod:specificity[@site][@position]') % mod_name.lower())

        if spec_list:
            spec_dict = defaultdict(list)
            for s in spec_list:
                spec_dict[int(s.get('spec_group'))].append((s.get('site'), s.get('position')))
            return spec_dict
        else:
            raise LookupError(mod_name)

    def get_mod_neutral_loss(self, mod_name, site, mass='monoisotopic'):
        '''For a given modification and a site specificity, returns the neutral
        loss list if possible. This returns None due to the mod not existing,
        or the site not existing, and it returns an empty tuple if the site has
        no neutral losses.

        'mass' should be 'monoisotopic' or 'average'. Defaults to monoisotopic mass.'''

        specificity_list = self.utree(('.//umod:modifications/'
                                       'umod:mod[f:lower(@title) = "%s"]/'
                                       'umod:specificity[@site="%s"]') % (mod_name.lower(), site))

        if not specificity_list and site in ('N-term', 'C-term'):
            specificity_list = self.utree(('.//umod:modifications/'
                                           'umod:mod[f:lower(@title) = "%s"]/'
                                           'umod:specificity[fn:ends-with(f:lower(@position), "%s")]') % (mod_name.lower(), site))

        if specificity_list:
            if mass == 'monoisotopic':
                return tuple(float(nl.get('mono_mass')) for s in specificity_list
                             for nl in s.xpath('.//umod:NeutralLoss[@mono_mass]',
                                               namespaces=NS))
            elif mass == 'average':
                return tuple(float(nl.get('avge_mass')) for s in specificity_list
                             for nl in s.xpath('.//umod:NeutralLoss[@avge_mass]',
                                               namespaces=NS))
            else:
                raise ValueError("%s is not a valid argument for mass type, "
                                 "should be 'monoisotopic' or 'average'." % mass)
        else:
            raise LookupError(mod_name, site)

    def get_mod_composition(self, mod_name):
        '''Searches for a mod, and returns the atomic composition of the delta
        as a dictionary of element symbol -> number pairs'''

        mod_list = self.utree(('.//umod:mod[f:lower(@title)="%s"]/'
                               'umod:delta[@mono_mass]') % mod_name.lower())

        if mod_list:
            return dict((e.get('symbol'), int(e.get('number')))
                        for e in mod_list[0].xpath('.//umod:element[@symbol][@number]',
                                                   namespaces=NS))
        else:
            raise LookupError(mod_name)
