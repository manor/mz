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

from decimal import Decimal
import scipy.stats as stats
from scipy.misc import comb
from collections import defaultdict
import re

def factorial(n):

    a = 1
    for i in range(1,n+1):
        a *= i
    return a


def combination(n, k):

    assert n >= k
    num = factorial(n)
    den = factorial(k) * factorial(n-k)
    return num / den


def hypergeom_pmf(m,M,n,N):
    """Returns probability mass function value for hypergeometric distribution with:
    m = # genes of type A in list
    M = # all genes in list
    n = # genes of type A in entire population
    N = # all genes in population

    """

##    num = Decimal(combination(n, m) * combination(N-n, M-m))
##    num = Decimal(comb(n, m, exact=1) * comb(N-n, M-m, exact=1))
##    den = Decimal(combination(N, M))
##    den = Decimal(comb(N, M, exact=1))
##    return float(num / den)
    return float(Decimal(comb(n, m, exact=1) * comb(N-n, M-m, exact=1)) / Decimal(comb(N, M, exact=1)))


def hypergeom_p_plus(m, M, n, N, alpha=1, midP=False):
    """Returns p value for hypergeometric distribution with:
    m = # genes of type A in list
    M = # all genes in list
    n = # genes of type A in entire population
    N = # all genes in population

    """


    p_plus = 0
    start = m
    stop = min(M,n)
    assert start <= stop
    if midP and start + 1 <= stop:
        start += 1

    for a in range(start, stop+1):
        val = hypergeom_pmf(a,M,n,N)
        p_plus += val
        if p_plus > alpha:
            return None

    if midP:
        p_plus += 0.5*hypergeom_pmf(m,M,n,N)
        if p_plus > alpha:
            return None

    return p_plus


def hypergeom_p_minus(m, M, n, N, alpha=1, midP=False):
    """Returns p value for hypergeometric distribution with:
    m = # genes of type A in list
    M = # all genes in list
    n = # genes of type A in entire population
    N = # all genes in population

    """

    p_minus = 0
    start = 0
    stop = m
    assert stop <= min(M,n)
    if midP and stop - 1 >= start:
        stop -= 1
    for a in range(stop, start-1, -1):
        val = hypergeom_pmf(a,M,n,N)
        p_minus += val
        if p_minus > alpha:
            return None

    if midP:
        p_minus += 0.5*hypergeom_pmf(m,M,n,N)
        if p_minus > alpha:
            return None

    return p_minus


def parse_KEGG_path(map_title_tab_file, gene_map_tab_file, synonym_file):
    """Returns dictionary with {gene:pathway}

    """

    path_title = {}
    fh = open(map_title_tab_file)
    for (i,line) in enumerate(fh):
        line = re.sub("\n+", "", line)

        a_l = line.split()
        id = a_l[0]
        defin = " ".join(a_l[1:])

        path_title[id] = defin

    fh.close()

    #Mouse gene to pathway annotations with kegg gene identifier as key
    annots = {}
    fh = open(gene_map_tab_file)
    for (i,line) in enumerate(fh):
        line = re.sub("\n+", "", line)

        a_l = line.split()
        id = a_l[0]
        annotation = a_l[1:]
        annots[id] = annotation

    #Definitions
    id_convert_kegg = defaultdict(list)     #kegg_id:[genes] (2 names can have same kegg_id, same name can have 2 kegg_ids)

    fh = open(synonym_file)
    for (i,line) in enumerate(fh):
        line = re.sub("\n+", "", line)

        a_l = line.split()
        name = a_l[0].upper()
        id = a_l[1]

        id_convert_kegg[id].append(name)

    annotations = None
    annotations = {}
    for (i,kegg_id) in enumerate(annots):
        pathway_ids = annots[kegg_id]
        pathway_ids = list(set(pathway_ids))
        pathway_defs = []
        for path_id in pathway_ids:
            pathway_defs.append(path_title[path_id])

        pathway_defs = list(set(pathway_defs))

        if not pathway_defs:
            continue
        genes = id_convert_kegg[kegg_id]
        for gene in genes:
            try:
                a = eval(gene)
            except:
                if gene in annotations:
                    annotations[gene].extend(pathway_defs)
                    annotations[gene] = list(set(annotations[gene]))
                else:
                    annotations[gene] = list(set(pathway_defs))

    return annotations


def parse_GO_MGI(gene_association_file):
    "Returns gene_association_dict with {gene:[go_id]} given MGI gene_association file"

    #Association
    gene_association_dict = defaultdict(list)
    fh = open(gene_association_file)
    for (i,line) in enumerate(fh):
        line = re.sub("\n+", "", line)
        if line[0] == "!":
            continue

        a_l = line.split()
        gene = a_l[2].upper()
        annot_id = a_l[3]
        gene_association_dict[gene].append(annot_id)

    return gene_association_dict


def parse_GO(ontology_file, gene_association_dict):
    """Returns dictionary with {gene:go_annotation}
    gene_association_dict is a dictionary with {gene:[go_id]}
    For MGI gene_association_file, the dict can be created using parse_GO_MGI

    """


    go_parents = {} #child:[parents]
    fh = open(ontology_file)
    parents = ["GO:0000000"]
    annot_defs = {"GO:0000000":"Fake Annotation"}
    for (i,line) in enumerate(fh):
        line = re.sub("\n+", "", line)
        if line[0] in ["!", "$"]:
            continue
        if line[1] == "%":
            continue
        a_l = line.split(";")
        annot_id = a_l[1].replace(" ","").split(",")[0]
        m = re.match("(GO:\d+)", annot_id)
        annot_id = m.groups()[0]
        space_num = a_l[0].find("%")
        if space_num < 1:       #DO not include "< - IS A PART OF" relationships
            continue
        annotation = a_l[0][space_num+1:-1]
        if len(parents) <= space_num - 1:
            parents.append(annot_id)
        else:
            parents[space_num - 1] = annot_id


        go_parents[annot_id] = parents[1:space_num - 1]
        annot_defs[annot_id] = annotation

    fh.close()


    #Association
    annots = {}
    for (i,gene) in enumerate(gene_association_dict):
        annot_ids = gene_association_dict[gene]
        for annot_id in annot_ids:
            if annot_id in annot_defs:
                if gene in annots:
                    annots[gene].append(annot_defs[annot_id])
                else:
                    annots[gene] = [annot_defs[annot_id]]
                parents = go_parents[annot_id]
                for parent in parents:
                    annots[gene].append(annot_defs[parent])

                annots[gene] = list(set(annots[gene]))


    return annots


def generic_annotate(genes, annotations_dict, superset_genes=None, exclude_unknown=True, **kwds):
    """Annotate genes with annotation and determine p values
    Takes a dictionary with {gene: [annotations]}

    returns two dictionaries
    1. {gene: annotation}
    2. {annotation: (m, M, n, N, p_plus, p_minus)}

    m = # genes of type A in list
    M = # all genes in list
    n = # genes of type A in entire population
    N = # all genes in population
    p_plus = P_value for enrichment of annotation
    p_minus = P_value for depletion of annotation

    """



    #Make list of annotation: gene count
    if superset_genes:
        total_genes = 0
    else:
        total_genes = len(annotations_dict)
    total_genes_per_annot = defaultdict(int)
    for gene in annotations_dict:
        if superset_genes:
            if gene not in superset_genes:
                continue
            total_genes += 1
        for annotation in annotations_dict[gene]:
            total_genes_per_annot[annotation] += 1

    if not exclude_unknown:
        total_genes = len(superset_genes)


    #for calculating p values
    gene_annotation = {}
    query_annotation_counts = defaultdict(int)

    unknown_count = 0
    for (i,gene) in enumerate(genes):
        if superset_genes:
            assert gene in superset_genes
        if gene not in annotations_dict:
            annots = []
        else:
            annots = annotations_dict[gene]

        #for calculating p values
        for annot in annots:
            query_annotation_counts[annot] += 1

        if not annots:
            annotation = "Unknown"
            unknown_count += 1
        else:
            annotation = "; ".join(annots)

        gene_annotation[gene] = annotation


    annotation_p_values = {}
    for (i,annotation) in enumerate(query_annotation_counts):
        m = query_annotation_counts[annotation]        # genes of type A in list
        M = len(genes)                          # all genes in list
        if exclude_unknown:
            M -= unknown_count
        n = total_genes_per_annot[annotation]               # genes of type A in entire population
        N = total_genes                             # all genes in population

        p_plus = hypergeom_p_plus(m, M, n, N, **kwds)
        p_minus = hypergeom_p_minus(m, M, n, N, **kwds)

        annotation_p_values[annotation] = (m, M, n, N, p_plus, p_minus)

    return (gene_annotation, annotation_p_values)
