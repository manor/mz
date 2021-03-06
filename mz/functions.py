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

import os
import re

from random import randint, sample
from collections import defaultdict

from mz import logger_message
#from mz import myData, logger_message
#from unimod import UnimodDatabase, LookupError
#
#try:
#    # load the unimod database
#    unimod = UnimodDatabase(os.path.join(myData, 'unimod.xml'))
#
#    AW = unimod.elements # atomic masses
#    amino_acids = unimod.amino_acids # amino acid masses
#except IOError:
#    # this shouldn't prevent multiplierz in general from functioning,
#    # although of course anything that relies on unimod will be broken
#    logger_message(50, ('Warning: %s not found.\n'
#                        'Please download a copy from http://www.unimod.org\n'
#                        'or use the copy on your local Mascot server')
#                   % (os.path.join(myData, 'unimod.xml')))
#
#    unimod = None
#    AW = None
#    amino_acids = None

import mod
AW = mod.AW
amino_acids = mod.amino_acids

# dictionary of modification shorthand
mod_dictionary = {'p': 'Phospho',
                  'i': 'iTRAQ4plex',
                  'o': 'Oxidation',
                  'c': 'Carbamidomethyl',
                  'd': 'Deamidated'}


def digest(protein, enzyme="Trypsin", missed_cleavages=0):
    '''Produces digested peptide sequences for a given protein sequence and an enzyme.

    Takes a protein sequence and digests it with the chosen enzyme.
    Optionally select # of missed cleavages
    Returns a list of tuples: [(peptide, (start,end), num_missed_cleavages)]

    The following enzymes are available:
        Arg-C
        Asp-N
        Bromelain
        CNBr_HSer
        CNBr_HSerLac
        Cathepsin B
        Cathepsin D
        Cathepsin G
        Chymotrypsin
        Clostripain
        Elastase
        Glu-C_Bic
        Glu-C_Phos
        Hydroxylamine
        Lys-C
        Lys-N
        Papain
        Pepsin
        Proteinase K
        Subtilisin
        Thermolysin
        Trypsin

    If the enzyme name starts with [ , the enzyme is used as the regular expression

    Example:
        >>> protein = 'MWKASAGHAVSIAQDDAGADDWETDPDFVNDVSEKEQRWGAKTVQGSGHQEHINIHKLRENVFQEHQTLKEKELETGPKASHGYGGKF'
        >>> enzyme = 'Trypsin'
        >>> digests = mzFunctions.digest(protein, enzyme)
        >>> for digest in digests:
        ... 	print digest[0],digest[1]
        ...
        MWK [1, 3]
        ASAGHAVSIAQDDAGADDWETDPDFVNDVSEK [4, 35]
        EQR [36, 38]
        WGAK [39, 42]
        TVQGSGHQEHINIHK [43, 57]
        LR [58, 59]
        ENVFQEHQTLK [60, 70]
        EK [71, 72]
        ELETGPK [73, 79]
        ASHGYGGK [80, 87]
        F [88, 88]
    '''

    protein = re.sub("\s+","",protein)
    protein = re.sub("\n+","",protein)

    #Enzyme Specificities
    # list of enzymes
    enzSpec = {'Arg-C':          '[R][A-Z]',
               'Asp-N':          '[A-Z][D]',
               'Bromelain':      '[KAY][A-Z]',
               'CNBr_HSer':      '[M][A-Z]',
               'CNBr_HSerLac':   '[M][A-Z]',
               'Cathepsin B':    '[R][A-Z]',
               'Cathepsin D':    '[LF][^VAG]',
               'Cathepsin G':    '[YWF][A-Z]',
               'Chymotrypsin':   '[YWFL][A-Z]',
               'Clostripain':    '[R][P]',
               'Elastase':       '[AVLIGS][A-Z]',
               'Glu-C_Bic':      '[E][A-Z]',
               'Glu-C_Phos':     '[ED][A-Z]',
               'Hydroxylamine':  '[N][G]',
               'Lys-C':          '[K][A-Z]',
               'Lys-N':          '[A-Z][K]',
               'Papain':         '[RK][A-Z]',
               'Pepsin':         '[LF][^VAG]',
               'Proteinase K':   '[YWF][A-Z]',
               'Subtilisin':     '[^RHK][A-Z]',
               'Thermolysin':    '[LFIVMA][^P]',
               'Trypsin':        '[KR][^P]'
               }

    digest_array = []

    if enzyme[0] == "[":
        reg_exp = enzyme
    else:
        reg_exp = enzSpec[enzyme]

    expr = re.compile(reg_exp)

    start = 0
    while True:
        m = expr.search(protein, start)
        if m is None:
            break
        else:
            digest_array.append((protein[start:(m.end()-1)], (start, m.end()-1), 0))
            start = 1 + m.start()

    if start != len(protein):
        digest_array.append((protein[start:], (start, len(protein)), 0))

    #Add missed cleavages
    if missed_cleavages:
        dig_all = digest_array[:]
        for i,p in enumerate(digest_array[:-1]):
            pep = p[0]
            start = p[1][0]
            for j in range(min((len(digest_array) - i - 1), missed_cleavages)):
                pep += digest_array[i+j+1][0]
                end = digest_array[i+j+1][1][1]
                dig_all.append((pep, (start,end),j+1))
        return dig_all

    return digest_array


def fragment(peptide, ions=('b', 'y'), labels=False):
    '''Calculates the ion masses for a given peptide, ion types, and residue modifications.

    Takes:
    - a peptide sequence (with possible modifications)
    - a list of ion types (defaults to b and y)
    - 'labels' flag: if True, function returns two lists: ion masses and corresponding labels

    If 'labels' is not set, the return is a list of lists of ion masses, in the same order as
    provided to the function. The length of the lists will vary--most ions types will not return
    the final ion (e.g. the b7 ion of 'PEPTIDE') because it is almost never observed.

    The function will calculate neutral losses

    The following ions are available:
        imm (immonium)
        a, a0 (H20 loss), a* (NH3 loss)
        b, b0 (H20 loss), b* (NH3 loss)
        c
        x
        y, y0 (H20 loss), y* (NH3 loss)
        z, z+1, z+2
        intya (internal ya), intyb (internal yb)

        Add a '++' after any of the ions to make it doubly charged

    Also, passing 'MH+', 'MH++', etc will calculate the peptide m/z values

    The following modifications are available:
        p: phospho = H2PO3 -H
        i: iTRAQ = H12C4C*3NN*O
        o: Oxidation = O
        c: Carbamidomethyl = CH2CONH2 -H
        d: Deamidated = O -HN

    Alternatively, modifications can be specified using brackets with either a mass
    or Unimod name inside.

    Example:
        >>> peptide = '[iTRAQ4plex]-PEPpTIDE'
        >>> ions = ['b','y','b++','y++','b0','y0']
        >>> frags = mzFunctions.fragment(peptide, ions)
        >>> b_ions = frags[0]
        >>> print b_ions
        [242.16210299999997, 371.20469600000001, 468.25746000000004, 649.27147000000002, 762.35553400000003, 877.38247699999999]
    '''

    calc_masses = []
    label_dict = {}

    ion_list = ions
    include = set(ion_list)

    sequence, mods, vm_masses, neutral_loss = mz_pep_decode(peptide)

    neutral_loss_list = []
    neutral_loss_C_term = 0
    neutral_loss_N_term = 0

    net_NL_list = []
    left_label_list = []
    right_label_list = []
    var_mod_string_list = []

    # polarity is + for now
    polarity = '+'
    double = polarity * 2

    # amino acid masses
    masses = amino_acids.copy()

    min_internal_mass = 0.0
    max_internal_mass = 700.0

    running_sum = []

    # first calculate running sum of residues, including variable mods
    # ie calculate the approximate mass at each fragmentation point

    running_sum.append(masses[sequence[0]])
    if mods[1]:
        running_sum[0] += sum(vm_masses[mod] for mod in mods[1])
        neutral_loss_list.extend(neutral_loss[mod][0] for mod in mods[1])
    else:
        neutral_loss_list.append(0.0)

    for i,aa in enumerate(sequence[1:]):
        running_sum.append(running_sum[-1] + masses[aa])
        if mods[i+2]:
            running_sum[-1] += sum(vm_masses[mod] for mod in mods[i+2])
            for mod in mods[i+2]:
                neutral_loss_list.append(neutral_loss_list[-1] + neutral_loss[mod][0])
        else:
            neutral_loss_list.append(neutral_loss_list[-1])

    # n-term and c-term modifications:
    if mods[0]:
        masses['N-term'] += sum(vm_masses[mod] for mod in mods[0])
        neutral_loss_N_term = sum(neutral_loss[mod][0] for mod in mods[0])
    else:
        neutral_loss_N_term = 0.0

    if mods[-1]:
        masses['C-term'] += sum(vm_masses[mod] for mod in mods[-1])
        neutral_loss_C_term = sum(neutral_loss[mod][0] for mod in mods[-1])
    else:
        neutral_loss_C_term = 0.0

    calc_MR = running_sum[-1] + masses['N-term'] + masses['C-term']

    hydrogen = AW['H']
    carbon = AW['C']
    nitrogen = AW['N']
    oxygen = AW['O']
    electron = AW['e']

    CO = carbon + oxygen
    NH3 = nitrogen + 3 * hydrogen
    H2O = 2 * hydrogen + oxygen
    if (polarity == "+"):
        charge_mass = hydrogen - electron
    else:
        charge_mass = - hydrogen + electron

    # nl_label = the net neutral as a text string for a peak label
    # frag_seq = the sequence of a fragment
    # mod_string = substring of mods corresponding to the fragment
    #   used for figuring out where neutral losses need to be permuted
    #   for high-energy fragment, string excludes new terminal residue
    # net_nl = total neutral loss for a fragment

    ions = defaultdict(list)

    RKNQ_re = re.compile(r'[RKNQ]')
    STED_re = re.compile(r'[STED]')

    len_seq = len(sequence)

    # n-term fragments
    for i in range(len_seq - 1):
        frag_seq = sequence[:i+1]

        mod_string = mods[:i+2]

        net_nl = neutral_loss_list[i] + neutral_loss_N_term

        # label is blank if the loss is ~zero
        nl_label = '%s' % (int(round(-net_nl)) or '')

        ions['a'].append(running_sum[i] + masses['N-term'] - CO
                         - net_nl - hydrogen + charge_mass)

        ions['b'].append(running_sum[i] + masses['N-term']
                         - net_nl - hydrogen + charge_mass)

        ions['c'].append(ions['b'][-1] + NH3)

        for k in ('a', 'b', 'c'):
            if k in include:
                calc_masses.append(ions[k][-1])
                net_NL_list.append(net_nl)
                left_label_list.append('%s(%d)' % (k, i+1))
                right_label_list.append('')
                var_mod_string_list.append(mod_string)
                label_dict[ions[k][-1]] = '%s(%d)%s' % (k, i+1, nl_label)

        # NH3-loss if fragment includes [RKNQ]
        if RKNQ_re.search(frag_seq):
            ions['a*'].append(ions['a'][-1] - NH3)
            ions['b*'].append(ions['b'][-1] - NH3)

            for k in ('a*', 'b*'):
                if k in include:
                    calc_masses.append(ions[k][-1])
                    net_NL_list.append(net_nl)
                    left_label_list.append('%s(%d)' % (k, i+1))
                    right_label_list.append('')
                    var_mod_string_list.append(mod_string)
                    label_dict[ions[k][-1]] = '%s(%d)%s' % (k, i+1, nl_label)

                    if ('%s++' % k) in include:
                        ions['%s++' % k].append((ions[k][-1] + charge_mass) / 2)

                        calc_masses.append(ions['%s++' % k][-1])
                        net_NL_list.append(net_nl)
                        left_label_list.append('%s(%d)' % (k, i+1))
                        right_label_list.append(double)
                        var_mod_string_list.append(mod_string)
                        label_dict[ions['%s++' % k][-1]] = '%s(%d)%s%s' % (k, i+1, nl_label, double)

        # water-loss if fragment includes [STED]
        if STED_re.search(frag_seq):
            ions['a0'].append(ions['a'][-1] - H2O)
            ions['b0'].append(ions['b'][-1] - H2O)

            for k in ('a0', 'b0'):
                if k in include:
                    calc_masses.append(ions[k][-1])
                    net_NL_list.append(net_nl)
                    left_label_list.append('%s(%d)' % (k, i+1))
                    right_label_list.append('')
                    var_mod_string_list.append(mod_string)
                    label_dict[ions[k][-1]] = '%s(%d)%s' % (k, i+1, nl_label)

                    if ('%s++' % k) in include:
                        ions['%s++' % k].append((ions[k][-1] + charge_mass) / 2)

                        calc_masses.append(ions['%s++' % k][-1])
                        net_NL_list.append(net_nl)
                        left_label_list.append('%s(%d)' % (k, i+1))
                        right_label_list.append(double)
                        var_mod_string_list.append(mod_string)
                        label_dict[ions['%s++' % k][-1]] = '%s(%d)%s%s' % (k, i+1, nl_label, double)

        # doubly-charged ions
        for k in ('a', 'b', 'c'):
            if ('%s++' % k) in include:
                ions['%s++' % k].append((ions[k][-1] + charge_mass) / 2)

                calc_masses.append(ions['%s++' % k][-1])
                net_NL_list.append(net_nl)
                left_label_list.append('%s(%d)' % (k, i+1))
                right_label_list.append(double)
                var_mod_string_list.append(mod_string)
                label_dict[ions['%s++' % k][-1]] = '%s(%d)%s%s' % (k, i+1, nl_label, double)

    # c-term fragments
    for i in range(len_seq - 1):
        frag_seq = sequence[-i-1:]

        mod_string = mods[-i-2:]

        net_nl = (neutral_loss_list[len_seq - 1]
                  - neutral_loss_list[len_seq - 2 - i]
                  + neutral_loss_C_term)

        nl_label = '%s' % (int(round(-net_nl)) or '')

        ions['y'].append(running_sum[len_seq - 1]
                         - running_sum[len_seq - 2 - i]
                         + masses['C-term'] + hydrogen
                         - net_nl + charge_mass)

        ions['x'].append(ions['y'][-1] - 2 * hydrogen + CO)

        ions['z'].append(ions['y'][-1] - NH3)

        ions['z+1'].append(ions['z'][-1] + hydrogen)

        ions['z+2'].append(ions['z+1'][-1] + hydrogen)

        for k in ('y', 'x', 'z', 'z+1', 'z+2'):
            if k in include:
                calc_masses.append(ions[k][-1])
                net_NL_list.append(net_nl)
                left_label_list.append('%s(%d)' % (k, i+1))
                right_label_list.append('')
                var_mod_string_list.append(mod_string)
                label_dict[ions[k][-1]] = '%s(%d)%s' % (k, i+1, nl_label)

        if RKNQ_re.search(frag_seq):
            ions['y*'].append(ions['y'][-1] - NH3)

            if 'y*' in include:
                calc_masses.append(ions['y*'][-1])
                net_NL_list.append(net_nl)
                left_label_list.append('y*(%d)' % (i+1))
                right_label_list.append('')
                var_mod_string_list.append(mod_string)
                label_dict[ions['y*'][-1]] = 'y*(%d)%s' % (i+1, nl_label)

                if 'y*++' in include:
                    ions['y*++'].append((ions['y*'][-1] + charge_mass) / 2)

                    calc_masses.append(ions['y*++'][-1])
                    net_NL_list.append(net_nl)
                    left_label_list.append('y*(%d)' % (i+1))
                    right_label_list.append(double)
                    var_mod_string_list.append(mod_string)
                    label_dict[ions['y*++'][-1]] = 'y*(%d)%s%s' % (i+1, nl_label, double)

        if STED_re.search(frag_seq):
            ions['y0'].append(ions['y'][-1] - H2O)

            if 'y0' in include:
                calc_masses.append(ions['y0'][-1])
                net_NL_list.append(net_nl)
                left_label_list.append('y0(%d)' % (i+1))
                right_label_list.append('')
                var_mod_string_list.append(mod_string)
                label_dict[ions['y0'][-1]] = 'y0(%d)%s' % (i+1, nl_label)

                if 'y0++' in include:
                    ions['y0++'].append((ions['y0'][-1] + charge_mass) / 2)

                    calc_masses.append(ions['y0++'][-1])
                    net_NL_list.append(net_nl)
                    left_label_list.append('y0(%d)' % (i+1))
                    right_label_list.append(double)
                    var_mod_string_list.append(mod_string)
                    label_dict[ions['y0++'][-1]] = 'y0(%d)%s%s' % (i+1, nl_label, double)

        # doubly-charged ions
        for k in ('y', 'x', 'z', 'z+1', 'z+2'):
            if ('%s++' % k) in include:
                ions['%s++' % k].append((ions[k][-1] + charge_mass) / 2)

                calc_masses.append(ions['%s++' % k][-1])
                net_NL_list.append(net_nl)
                left_label_list.append('%s(%d)' % (k, i+1))
                right_label_list.append(double)
                var_mod_string_list.append(mod_string)
                label_dict[ions['%s++' % k][-1]] = '%s(%d)%s%s' % (k, i+1, nl_label, double)

    # immonium
    if 'imm' in include:
        for i,aa in enumerate(sequence):
            ions['imm'].append(masses[aa] - CO + charge_mass)
            net_nl = 0
            if mods[i+1]:
                ions['imm'][-1] += vm_masses[mods[i+1]]
                net_nl = neutral_loss[mods[i+1]][0]
            #elif neutral_loss[aa]:
                #net_nl = neutral_loss[aa][0]
            ions['imm'][-1] -= net_nl

            nl_label = '%s' % (int(round(-net_nl)) or '')

            calc_masses.append(ions['imm'][-1])
            net_NL_list.append(net_nl)
            left_label_list.append(aa)
            right_label_list.append('')
            var_mod_string_list.append(mods[i+1])
            label_dict[ions['imm'][-1]] = '%s%s' % (aa, nl_label)

    # internals
    if 'intya' in include or 'intyb' in include:
        for i in range(len_seq - 3):
            for j in range(i+2, len_seq - 1):
                frag_seq = sequence[i+1:j+1]

                mod_string = mods[i+2:j+2]

                net_nl = neutral_loss_list[j] - neutral_loss_list[i]

                nl_label = '%s' % (int(round(-net_nl)) or '')

                ions['intyb'].append(running_sum[j] - running_sum[i]
                                     - net_nl + charge_mass)

                ions['intya'].append(ions['intyb'][-1] - CO)

                for k,lbl in (('intyb',''), ('intya','-CO')):
                    if k in include:
                        if min_internal_mass <= ions[k][-1] <= max_internal_mass:
                            calc_masses.append(ions[k][-1])
                            net_NL_list.append(net_nl)
                            left_label_list.append('%s%s' % (frag_seq, lbl))
                            right_label_list.append('')
                            var_mod_string_list.append(mod_string)
                            label_dict[ions[k][-1]] = '%s%s%s' % (frag_seq, lbl, nl_label)

    # going to ignore d/d', v, w/w' ions
    #if 'd' in include:
        #if sequence[0] == 'R':
            #ions['d'].append(masses['N-term'] - neutral_loss_N_term + charge_mass
                             #+ carbon * 2 + hydrogen * 4 + nitrogen)
            #calc_masses.append(ions['d'][-1])
            #net_NL_list.append(neutral_loss_N_term)
            #left_label_list.append('d(1)')
            #right_label_list.append('')
            #var_mod_string_list.append(mods[0]

    # if 'v' in include:
    # ...
    # if 'w' in include:
    # ...

    # If any variable mods have multiple neutral losses, we now need to
    # permute out additional calculated values
    for mod in range(1, len(vm_masses) + 1):
        if len(neutral_loss[mod]) > 1:
            for i in range(len(calc_masses)):
                if mod in reduce(set.union, var_mod_string_list[i], set()):
                    count = [m for ms in var_mod_string_list[i] for m in ms].count(mod)

                    for nl_mod in neutral_loss[mod][1:]:
                        delta = count * (nl_mod - neutral_loss[mod][0])
                        if delta:
                            net_NL_list.append(net_NL_list[i] + delta)
                            left_label_list.append(left_label_list[i])
                            right_label_list.append(right_label_list[i])
                            var_mod_string_list.append(var_mod_string_list[i])

                            nl_label = '%s' % (int(round(-net_NL_list[-1])) or '')

                            if right_label_list[i]:
                                # allow for charge
                                calc_masses.append(calc_masses[i] - delta / 2)
                            else:
                                calc_masses.append(calc_masses[i] - delta)

                            label_dict[calc_masses[-1]] = (left_label_list[i]
                                                           + nl_label
                                                           + right_label_list[i])


    # this code never worked and never ran...
    #for pep_nl_list in pep_neutral_loss.values():
        #for nl in pep_nl_list:
            #calc_masses.append((calc_MR - nl + getCharge() * charge_mass)
                               #/ pep.getCharge())
            #label_dict[calc_masses[-1]] = "M%.0f%s" % (-nl, polarity * pep.getCharge())

    #for pep_nl_list in req_pep_neutral_loss.values():
        #for nl in pep_nl_list:
            #calc_masses.append((calc_MR - nl + pep.getCharge() * charge_mass)
                               #/ pep.getCharge())
            #label_dict[calc_masses[-1]] = "M%.0f%s" % (-nl, polarity * pep.getCharge())

    mh_list = [mh for mh in include if mh.startswith('MH')]
    for mh in mh_list:
        charge = mh.count(polarity)
        if charge > 0:
            calc_masses.append((calc_MR + charge_mass * charge) / charge)
            label_dict[calc_masses[-1]] = mh

    if not labels:
        return [ions[i] for i in ion_list]
    else:
        return calc_masses, label_dict


def generate_labels(scan, peptide, ions, charge=None, tolerance=0.6, **settings):
    '''Takes an MS2 scan, and a peptide (in modification format as described
    in the fragment function above) and generates a list of mass-label pairs
    by matching theoretical ion masses to experimental values within a
    tolerance.

    - scan should be a list (or other sequence) of mz,intensity pairs to be matched.
      This should probably NOT be an entire scan--that will result in a lot of
      likely-erroneous matches (to noise). One option is to use the top 50 most
      intense ions
    - peptide should be in fragment-ready format (with modification masses)
    - ions should be a list or set of ions to look for
    - tolerance is in Daltons, 0.6 is a common value
    - settings is a dictionary of optional arguments (default in brackets):
      - show_theor_mz = [ True ] | False - Display theoretical ion mass
      - ms2_mz_figs = integer >= 0, default 2 - Significant figures for mass error
      - show_mass_error = True | [ False ] - Display mass error of ions
      - mass_error_figs = integer >= 0, default 2 - Significant figures for mass error
      - mass_error_units = ['ppm'] | 'Da' - Units to display mass error

    Returns a tuple of mass-label pairs--the experimental peaks that matched
    the theoretical ion masses within tolerance.
    '''

    # initialize defaults and update with optional arguments
    _settings = dict(show_theor_mz=True, ms2_mz_figs=2, show_mass_error=False,
                     mass_error_figs=2, mass_error_units='ppm')
    _settings.update(settings)

    if _settings['mass_error_units'] == 'ppm':
        calc_error = lambda exp_mz, theor_mz: (abs(theor_mz - exp_mz) / theor_mz) * 1E6
    else:
        calc_error = lambda exp_mz, theor_mz: abs(theor_mz - exp_mz)

    if _settings['show_theor_mz'] and _settings['show_mass_error']:
        label_text = '%%s [%%.%df - %%.%df %s]' % (_settings['ms2_mz_figs'],
                                                   _settings['mass_error_figs'],
                                                   _settings['mass_error_units'])
    elif _settings['show_theor_mz']:
        label_text = '%%s [%%.%df]' % _settings['ms2_mz_figs']
    elif _settings['show_mass_error']:
        label_text = '%%s [%%.%df %s]' % (_settings['mass_error_figs'],
                                          _settings['mass_error_units'])

    if charge and 'MH' in ions:
        ions = list(ions)
        ions.remove('MH')

        for z in range(int(charge), 0, -1):
            ions.append('MH%s' % ('+'*z))

    calc_masses, label_dict = fragment(peptide, ions, labels=True)

    scan = sorted(scan)

    match_count = defaultdict(lambda: -1)

    matched_calc = []
    matched_exp = []
    matched_int = []

    # for each mass/int peak in the scan
    for j,(mass,inte) in enumerate(scan):
        # go through each calculated mass and try to match
        for i,cmass in enumerate(calc_masses):
            if abs(mass - cmass) <= tolerance:
                if match_count[i] > -1:
                    if matched_int[match_count[i]] < inte:
                        matched_exp[match_count[i]] = mass
                        matched_int[match_count[i]] = inte
                    continue

                matched_calc.append(cmass)
                matched_exp.append(mass)
                matched_int.append(inte)
                match_count[i] = len(matched_int) - 1

    matched_exp.sort()
    matched_calc.sort()

    if _settings['show_theor_mz'] and _settings['show_mass_error']:
        return tuple((e, (label_text % (label_dict[c], c, calc_error(e,c))))
                     for e,c in zip(matched_exp,matched_calc))
    elif _settings['show_theor_mz']:
        return tuple((e, (label_text % (label_dict[c], c)))
                     for e,c in zip(matched_exp,matched_calc))
    elif _settings['show_mass_error']:
        return tuple((e, (label_text % (label_dict[c], calc_error(e,c))))
                     for e,c in zip(matched_exp,matched_calc))
    else:
        return tuple((e, label_dict[c])
                     for e,c in zip(matched_exp,matched_calc))


def mz(peptide, charge=1):
    """Returns the mz for an amino acid sequence.

    Takes a peptide sequence (including possible mods) and returns the mz

    Example:
        >>> peptide = 'PEPTIDE'
        >>> totalMass = mzFunctions.mz(peptide, charge=1)
        >>> print totalMass
        799.359964319

    """

    sequence, mods, vm_masses, neutral_loss = mz_pep_decode(peptide)

    masses = amino_acids.copy()

    # randomly choose AAs for 'Z' and 'B'
    # Is this a good idea? Makes the mz function randomized
    masses['Z'] = masses[['E','Q'][randint(0,1)]]
    masses['B'] = masses[['N','D'][randint(0,1)]]

    total_mass = sum(masses[AA] for AA in sequence) + masses['N-term'] + masses['C-term']
    total_mass += sum(vm_masses[m] for mod in mods for m in mod)

    if charge:
        total_mass += (AW["H"] - AW["e"]) * charge
        total_mass /= charge

    return total_mass


def mw(peptide):
    """Returns the monoisotopic mass for an amino acid sequence.

    Calls mz(), with charge = 0

    """

    return mz(peptide, charge=0)


def parse_fasta(file_name):
    '''
    Parses a FASTA file and provides a generator of (header,sequence) tuples.
    Can be used efficiently even on very large FASTA files.

    Usage:
    for header,sequence in parse_fasta(file_name):
        print header
        process(sequence)
    '''

    fasta_file = open(file_name)
    header = ''
    sequence = []

    for i,line in enumerate(fasta_file):
        if line[0] == '>':
            if i > 0:
                yield (header, ''.join(sequence))
            sequence = []
            header = line[1:].strip()
        else:
            sequence.append(line.strip())

    yield (header, ''.join(sequence))

    fasta_file.close()


def write_fasta(fasta_dict, save_file, write_mode='w'):
    '''
    Creates a fasta file given a dictionary with prot: sequence.
    This method adds ">" at the beginning of each header.
    The proteins will not be written in any particular order.
    '''

    fasta_file = open(save_file, write_mode)

    for prot in fasta_dict:
        seq = fasta_dict[prot]
        fasta_file.write(">%s\n" % prot)
        fasta_file.write("%s\n\n" % seq)

    fasta_file.close()


def randAA(length):
    '''
    Returns a random amino acid sequence of given length. Note:
    all amino acids will appear with equal probability.

    Example:
        >>> randomSequence = mzFunctions.randAA(10)
        >>> print randomSequence
        QGSQHQYRGD

    '''

    AA_list = ['A','C','D','E','F',
               'G','H','I','K','L',
               'M','N','P','Q','R',
               'S','T','V','W','Y']

    return ''.join(sample(AA_list, length))


def mz_pep_format(peptide, modifications='', iTRAQ=False, carb=False):
    '''Converts a given sequence and a modifications string to multiplierz format

    The modifications string is in the following format:
    "Ri: Modification; " where R = residue, i = position, Modification = modification name
    Ex: "M3: Oxidation; "

    A boolean option for iTRAQ is provided, where True indicates that all Lysines (K)
    and the N-Terminus is labeled with an 'i' modification.

    A boolean option for carbamidomethyl is provided, where True indicates that all Cysteines (C)
    are labeled with a 'c' modification.

    >>> mz_pep_format('ATLPRTLK', 'T2: Phospho; ', iTRAQ = True)
    'i-ApTLPRTLiK'

    '''

    sorted_var_mods_list = []
    n_term_mod = None
    c_term_mod = None

    if peptide.startswith('i-'):
        peptide = peptide[2:]
    peptide = re.sub('iK', 'K', peptide)
    peptide = re.sub('cC', 'C', peptide)

    if not modifications:
        modifications = ''

    vm_re = re.compile(r'\s*([NC]-term|\w(\d+))\s*:\s*(.+)')
    for item in modifications.split(';'):
        m = vm_re.match(item)
        #m = re.search('[A-Z](\d+)\s*\:\s*(.+)', item)
        if m:
            if m.group(2):
                pos = int(m.group(2)) - 1
                mod_type = m.group(3)[0].lower()
                if mod_type in mod_dictionary:
                    sorted_var_mods_list.append((pos, mod_type))
            elif m.group(1).startswith('N'):
                mod_type = m.group(3)[0].lower()
                if mod_type in mod_dictionary:
                    n_term_mod = mod_type
            elif m.group(1).startswith('C'):
                mod_type = m.group(3)[0].lower()
                if mod_type in mod_dictionary:
                    c_term_mod = mod_type

    sorted_var_mods_list.sort()
    for (index,(pos, mod_type)) in enumerate(sorted_var_mods_list):
        peptide = peptide[:pos + index] + mod_type + peptide[pos + index:]

    if n_term_mod:
        peptide = '%s-%s' % (n_term_mod, peptide)
    if c_term_mod:
        peptide = '%s-%s' % (peptide, c_term_mod)

    if iTRAQ:
        peptide = 'i-' + 'iK'.join(peptide.split('K'))

    if carb:
        peptide = 'cC'.join(peptide.split('C'))

    return peptide


def mz_pep_decode(peptide):
    '''Close to the opposite of mz_pep_format, this takes a multiplierz-format
    peptide and extracts out the basic sequence and information about modifications.

    The primary use of this function is for the mz and fragment functions above. It
    is *not* the inverse of mz_pep_format, as it does not return modification names.
    This is because it is possible for the name to be unknown, e.g. if there was a
    bare delta value in the peptide.

    The output is a 4-tuple:
    - The peptide sequence, with all mods removed
    - An array of mods, one for each amino acid as well as the N and C terminii
      (the first and last elements, respectively). Each element of the mods array
      is a set of integers, each integer being a unique mod of this peptide.
    - A dictionary of mod deltas, with keys being the integers used above.
    - A dictionary of neutral losses, keys as above and values being tuples
      of neutral loss values, (0.0,) being the minimum possibility.

    '''

    # empty mod brackets are not allowed
    peptide = re.sub(r'\[\]', '', peptide)

    pep_cleaner = re.compile(r'-?(\[.+?\]|[a-z])-?')

    # peptide without any modifications
    sequence = pep_cleaner.sub('', peptide)

    # variable modifications we can handle at the moment
    vm_masses = {}
    neutral_loss = {}

    mods = [set() for i in range(len(sequence) + 2)]

    # first regex finds mods per residue (either a bracketed group or an abbreviation)
    mod_finder_A = re.compile(r'\[(.+?)\]|([cdiop])')
    # second regex parses a bracketed group (a '; '-delimited list of names/masses)
    mod_finder_B = re.compile(r'([-+]?\d+\.?(?:\d+)?|\.\d+)|(.+)')

    mod_set = dict()

    for mA in mod_finder_A.finditer(peptide):
        if mA.end() < len(peptide):
            AA = peptide[mA.end()]
            if AA == '-':
                AA = 'N-term'
        else:
            AA = 'C-term'

        if mA.group(1): # bracketed group
            m = []
            for g in mA.group(1).split('; '):
                mB = mod_finder_B.match(g) # this regex can't fail
                if mB.group(1): # number, no neutral losses known
                    m.append((mB.group(1),
                              float(mB.group(1)),
                              (0.0,)))
                else: # mod name (this will raise an exception if invalid)
                    m.append((mB.group(2),
                              mod.delta[mB.group(2)],
                              (0.0,) + mod.neutral_losses[(mB.group(2), AA)]))
                    #m.append((mB.group(2),
                    #          unimod.get_mod_delta(mB.group(2)),
                    #          (0.0,) + unimod.get_mod_neutral_loss(mB.group(2), AA)))
        else: # abbreviation
            m = [(mA.group(2),
                  mod.delta[mod_dictionary[mA.group(2)]],
                  (0.0,) + mod.neutral_losses[mod_dictionary[mA.group(2)], AA])]
            #m = [(mA.group(2),
            #      unimod.get_mod_delta(mod_dictionary[mA.group(2)]),
            #      (0.0,) + unimod.get_mod_neutral_loss(mod_dictionary[mA.group(2)], AA))]

            # note: the first neutral loss is always 0.0, so that the mod delta is reflected in
            # the default set of ions. additional neutral losses (with corresponding labels) are
            # available if you set labels=True
            for a_mod,mod_delta,mod_neutral_losses in m:
                if a_mod not in mod_set:
                    vm_masses[len(vm_masses) + 1] = mod_delta
                    neutral_loss[len(vm_masses)] = mod_neutral_losses
                    mod_set[a_mod] = len(vm_masses)

                if mA.end() < len(peptide):
                    mods[len(pep_cleaner.sub('', peptide[:mA.end() + 1]))].add(mod_set[a_mod])
                else:
                    mods[len(sequence) + 1].add(mod_set[a_mod])

    mods = [frozenset(m) for m in mods]

    return sequence, mods, vm_masses, neutral_loss
