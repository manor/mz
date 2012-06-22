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
from scipy import log2, array

import functions
import report

#mz spreadsheet IO
def get_entrez_XLS(f):
    #Store Entrez gene names
    with report.reader(f, sheet_name='GenBank_Info') as genbank:
        entrez = dict((row.acc, (row.gene).upper()) for row in genbank)

    return entrez


def get_iTRAQ_prot_raw(iTRAQ_pep_raw):
    iTRAQ_prot_raw = {}

    for pep in iTRAQ_pep_raw:
        gene = ''.join(pep.split("_")[:-1])
        #Raw Values
        raw_counts = array(iTRAQ_pep_raw[pep])

        #Sum intensities per peptide per protein
        if gene in iTRAQ_prot_raw:
            iTRAQ_prot_raw[gene] += raw_counts
        else:
            iTRAQ_prot_raw[gene] = raw_counts


    return iTRAQ_prot_raw


def get_iTRAQ_pep_raw_XLS(f, iTRAQ_reporter_ions, entrez):
    iTRAQ_data_raw = {}

    with report.reader(f) as rdr:
        for row in rdr:
            seq, mods, acc = (row.seq, row.var_mods, row.acc)
            acc = acc.replace(" ", "")

            gene = entrez[acc].upper()
            pep = functions.mz_pep_format(seq, mods)

            #Raw Values
            raw_counts = []
            for ion in iTRAQ_reporter_ions:
                raw_counts.append(float(row["Rep%d" % ion]))

            raw_counts = array(raw_counts)

            key = "%s_%s" % (gene, pep)
            #Sum intensities per peptide per protein
            if key in iTRAQ_data_raw:
                iTRAQ_data_raw[key] += raw_counts
            else:
                iTRAQ_data_raw[key] = raw_counts

    return iTRAQ_data_raw


def iTRAQ_log2_ratios(iTRAQ_data_raw, iTRAQ_reporter_ions, den):
    #Calculate iTRAQ ratios
    iTRAQ_ratios = {}
    ratio_defins = []
    gene_counter = 0

    for gene in iTRAQ_data_raw:
        gene_counter += 1
        iTRAQ_ratios[gene] = []
        den_index = iTRAQ_reporter_ions.index(den)
        den_value = iTRAQ_data_raw[gene][den_index]
        for j in range(len(iTRAQ_data_raw[gene])):
            if j == den_index:
                continue
            ratio = log2(iTRAQ_data_raw[gene][j] / den_value)
            if gene_counter == 1:
                ratio_defins.append("%d:%d" % (iTRAQ_reporter_ions[j], den))
            iTRAQ_ratios[gene].append(ratio)

    return (iTRAQ_ratios, ratio_defins)


def write_annots_XLS(annotations, p_values, ftype="KEGG_Path", f=None):
    """ Write annotations to xls file f. If no f then f = mz_annotations.xls
    annotations is a dictionary {gene:annotations}
    p_values is a dictionary {annotation: (m, M, n, N, p_plus, p_minus)}

    m = # genes of type A in list
    M = # all genes in list
    n = # genes of type A in entire population
    N = # all genes in population
    p_plus = P_value for enrichment of annotation
    p_minus = P_value for depletion of annotation

    ftype is the sheet name, which will have _Annots and _P_Values appended to end for two seperate sheets

    """

    if not f:
        f = os.path.join(os.getcwd(), "mz_annotations.xls")

    wtr = report.writer(f, sheet_name='%s_Annots' % ftype,
                          columns=['Gene', ftype])

    for k,v in annotations.items():
        wtr.write([k,v])

    wtr.close()

    pcols = [ftype,
             "# of Genes in Subset with Annotation",
             "# of Genes in Subset",
             "# of Genes in Entire Population with Annotation",
             "# of Genes in Entire Population",
             "P_+_Value",
             "P_-_Value"]

    wtr = report.writer(f, sheet_name='%s_P_Values' % ftype,
                          columns=pcols)

    for (k,(m, M, n, N, p_plus, p_minus)) in p_values.items():
        if not p_plus and not p_minus:
            continue
        wtr.write((k, m, M, n, N, p_plus or 99999, p_minus or 99999))

    wtr.close()
