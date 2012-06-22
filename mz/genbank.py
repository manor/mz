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

import lxml.etree as ET

import pycurl
import cStringIO

import tools
import report

from multiplierz import myData

class GenBankProtein(object):
    '''A class to represent the INSDSeq XML tree of a protein from the NCBI database.
    This class is primarily here for scripts that may want additional acesss to the
    XML tree, or for one-off instances that don't need to be cached'''

    def __init__(self, xml_string):
        '''Constructor takes the XML string, not an accession number. The XML should be
        downloaded separately.'''

        # construct an XPathEvaluator of the XML
        self.etree = ET.XPathEvaluator(ET.XML(xml_string))

    @property
    def protein_description(self):
        '''The protein description'''

        prot_desc = self.etree('.//INSDSeq_definition/text()')
        if len(prot_desc) > 1:
            print "More than one protein definition, that's weird"

        if prot_desc:
            return prot_desc[0]
        else:
            return ''

    @property
    def region_note(self):
        '''Region note(s)'''

        # complex XPath expression because this XML is not intuitively designed.
        # I first find the 'Region' feature (based on a child tag) and then find any
        # note qualifiers (again based on a child tag), then return the value text
        region_notes = self.etree(".//INSDFeature[child::INSDFeature_key/text() = 'Region']"
                                  "/INSDFeature_quals"
                                  "/INSDQualifier[child::INSDQualifier_name/text() = 'note']"
                                  "/INSDQualifier_value/text()")

        return '; '.join(region_notes)

    @property
    def go_annotations(self):
        '''GO annotation information. Possibly deprecated in NCBI, but here for completeness'''
        go_annotes = {'GO_component': [], 'GO_function': [], 'GO_process': []}

        cds_note = self.etree(".//INSDFeature[child::INSDFeature_key/text() = 'CDS']"
                              "/INSDFeature_quals"
                              "/INSDQualifier[child::INSDQualifier_name/text() = 'note']"
                              "/INSDQualifier_value/text()")

        if cds_note:
            for note in cds_note[0].strip().split('; '):
                s = note.split(':', 1)
                if s[0] in go_annotes:
                    go_annotes[s[0]].append(s[1])

        return go_annotes

    @property
    def gene(self):
        '''Gene ID'''

        gene = self.etree(".//INSDFeature[child::INSDFeature_key/text() = 'CDS']"
                          "/INSDFeature_quals"
                          "/INSDQualifier[child::INSDQualifier_name/text()='gene']"
                          "/INSDQualifier_value/text()")

        if len(gene) > 1:
            print "More than one gene?"

        if gene:
            return gene[0]
        else:
            return ''

    @property
    def xrefs(self):
        '''Gene cross-references, returned as a dictionary of name,value pairs'''

        xrefs = self.etree(".//INSDFeature[child::INSDFeature_key/text() = 'CDS']"
                           "/INSDFeature_quals"
                           "/INSDQualifier[child::INSDQualifier_name/text()='db_xref']"
                           "/INSDQualifier_value/text()")

        xref_dict = dict((k,'') for k in ('GeneID','HGNC','HPRD','MIM','CCDS'))
        for xref in xrefs:
            x = xref.split(':')
            if x[0] in xref_dict:
                xref_dict[x[0]] = x[-1]

        return xref_dict


class CacheProtein(object):
    '''A class to represent a single cached protein--a very simple XML tree'''

    def __init__(self, accession, root=None, genbank_protein=None):
        '''Takes either the root of cache protein tree, or takes a
        GenBankProtein instance and extracts the data we cache'''

        if root is not None:
            self.etree = ET.XPathEvaluator(root)
        elif genbank_protein:
            protein = ET.Element('protein', attrib={'gi': accession})

            ET.SubElement(protein, 'description').text = genbank_protein.protein_description
            ET.SubElement(protein, 'region_note').text = genbank_protein.region_note

            go = ET.SubElement(protein, 'GO_annotations')
            go_annotes = genbank_protein.go_annotations
            for k in go_annotes:
                for g in go_annotes[k]:
                    ET.SubElement(go, k).text = g

            ET.SubElement(protein, 'gene').text = genbank_protein.gene

            xref = ET.SubElement(protein, 'cross_references')
            xrefs = genbank_protein.xrefs
            for k in xrefs:
                ET.SubElement(xref, k).text = xrefs[k]

            self.etree = ET.XPathEvaluator(protein)
        else:
            raise ValueError("Must provide either a root node or a GenBankProtein instance")

    @property
    def root(self):
        return self.etree('/*')[0]

    @property
    def protein_description(self):
        '''Protein description'''

        prot_desc = self.etree('.//description/text()')

        if prot_desc:
            return prot_desc[0]
        else:
            return ''

    @property
    def region_note(self):
        '''Region note'''

        region_note = self.etree(".//region_note/text()")

        if region_note:
            return region_note[0]
        else:
            return ''

    @property
    def go_annotations(self):
        '''GO annotation information. Possibly deprecated in NCBI, but here for completeness'''
        go_annotes = {'GO_component': [], 'GO_function': [], 'GO_process': []}

        for k in go_annotes:
            go_annotes[k].extend(self.etree("./GO_annotations/%s/text()" % k))

        return go_annotes

    @property
    def gene(self):
        '''Gene ID'''

        gene = self.etree('.//gene/text()')

        if gene:
            return gene[0]
        else:
            return ''

    @property
    def xrefs(self):
        '''Gene cross-references, returned as a dictionary of name,value pairs'''

        xrefs = self.etree('./cross_references/*')

        xref_dict = {}
        for xref in xrefs:
            xref_dict[xref.tag] = xref.text

        return xref_dict


class GenBankCache(object):
    '''A class to hold the GenBank cache, which is information about proteins that a
    user has downloaded in the past. If they need information for the same proteins
    again, the cache allows them to skip the download step'''

    def __init__(self, file_name=None):
        if file_name:
            self.gtree = ET.XPathEvaluator(ET.parse(file_name))
        else:
            self.gtree = ET.XPathEvaluator(ET.XML('<genbank_cache/>'))

    def has_protein(self, accession):
        '''Returns True if a protein is in cache, False otherwise'''
        p = self.gtree('./protein[@gi="%s"]' % accession)

        if p:
            return True
        else:
            return False

    def get_protein(self, accession):
        '''Return a protein from the cache, or None if it is not found'''
        p = self.gtree('./protein[@gi="%s"]' % accession)

        if p:
            return CacheProtein(accession, root=p[0])
        else:
            return None

    def add_protein(self, cache_protein):
        '''Add a protein to the cache. The input should be an instance
        of the CacheProtein class'''
        self.gtree('/*')[0].append(cache_protein.root)

    def write(self, file_name):
        ET.ElementTree(self.gtree('/*')[0]).write(file_name, pretty_print=True)


def download_genbank(acc):
    acc = str(acc)

    crl = pycurl.Curl()

    crl.setopt(pycurl.HTTPGET, True)
    crl.setopt(pycurl.URL,
               ('http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?'
                'db=protein&val=%s&dopt=gpc_xml&sendto=on' % acc))

    response = cStringIO.StringIO()
    crl.setopt(pycurl.WRITEFUNCTION, response.write)

    crl.perform()
    crl.close()

    return GenBankProtein(response.getvalue())


def genbank_report(report_file):
    tools.logger_message(30, 'Fetching GenBank Information...')

    genbank_cache_file = os.path.join(myData, 'genbank_cache.xml')
    if os.path.exists(genbank_cache_file):
        genbank_cache = GenBankCache(genbank_cache_file)
    else:
        genbank_cache = GenBankCache()

    is_XLS = os.path.splitext(report_file.lower())[1] in ('.xls', '.xlsx')

    rdr = report.reader(report_file)

    acc_set = set()

    subre = re.compile('["\s]')
    accre = re.compile('gi\|(\d+)')

    for row in rdr:
        acc = str(row['Accession Number'])
        desc =  str(row['Protein Description'])
        rank = str(row['Protein Rank'])

        if not acc:
            acc = ''
        if not desc:
            desc = ''
        if not rank:
            rank = None
        else:
            rank = int(rank)

        acc = subre.sub('', acc)

        acc_m = accre.match(acc)
        if acc_m:
            acc_set.add((rank, acc, desc, acc_m.group(1)))
        else:
            acc_set.add((rank, acc, desc, None))

    rdr.close()

    cols = ['Protein Rank', 'Accession Number', 'Protein Description',
            'REGION NOTE', 'GO FUNCTION', 'GO PROCESS',
            'GO COMPONENT', 'GENE', 'ENTREZ GENE LINK',
            'HUGO GENE LINK', 'HPRD LINK', 'OMIM LINK',
            'CONSENSUS CDS LINK']

    if report_file.lower().endswith('.mzd'):
        wtr = report.writer(report_file, columns=cols, table_name='GenBank_Info')
    elif is_XLS:
        wtr = report.writer(report_file, columns=cols, sheet_name='GenBank_Info')
    else:
        return

    gi_url = '=HYPERLINK("http://www.ncbi.nlm.nih.gov/protein/%s", "%s")'

    for i,(rank,acc,desc,gi) in enumerate(sorted(acc_set)):
        if gi is None:
            row = (rank, acc, re.sub('"', '', desc), 'Accession is not in correct format.') + (None,) * 9
            wtr.write(row)
            continue

        prot = genbank_cache.get_protein(gi)

        if prot is None:
            tools.logger_message(20, 'Downloading %s from GenBank' % acc)
            prot = CacheProtein(gi, genbank_protein=download_genbank(gi))
            genbank_cache.add_protein(prot)

        if desc == '':
            desc = prot.protein_description

        row = [rank,
               None if is_XLS else acc,
               re.sub('"', '', desc) or None,
               prot.region_note or None,
               '; '.join(prot.go_annotations['GO_function']) or None,
               '; '.join(prot.go_annotations['GO_process']) or None,
               '; '.join(prot.go_annotations['GO_component']) or None,
               prot.gene or None]

        if is_XLS:
            metadata = [('Accession Number', 'formula', (gi_url % (gi, acc)))]
        else:
            metadata = []

        xrefs = prot.xrefs

        if is_XLS:
            row.extend([None] * 5)
            for c,x,url in zip(cols[8:],
                               ('GeneID','HGNC','HPRD','MIM','CCDS'),
                               ('http://www.ncbi.nlm.nih.gov/gene/%s',
                                'http://www.genenames.org/data/hgnc_data.php?hgnc_id=%s',
                                'http://www.hprd.org/protein/%s',
                                'http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=%s',
                                'http://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=%s')):
                xref = xrefs.get(x, None)
                if xref:
                    metadata.append((c, 'formula', '=HYPERLINK("%s","%s")' % (url % xref, xref)))
        else:
            for x in ('GeneID','HGNC','HPRD','MIM','CCDS'):
                row.append(xrefs.get(x, None))

        wtr.write(row, metadata)

    wtr.close()
    genbank_cache.write(genbank_cache_file)

    if is_XLS:
        report.Spreadsheet.genbank_sheet_format(report_file, i+1)


def download_entrez(entrez_id):
    """
    Returns dictionary {official_full_name: name,
                        AKA: [aka_names]}
    """

    try:
        entrez_id = int(entrez_id)
    except:
        pass

    entrez_id = str(entrez_id)

    crl = pycurl.Curl()

    crl.setopt(pycurl.URL,
               'http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=full_report&list_uids=%s&retmode=xml' % entrez_id)

    response = cStringIO.StringIO()
    crl.setopt(pycurl.WRITEFUNCTION, response.write)

    crl.perform()

    data = response.getvalue()

    crl.close()

    entrez_tree = ET.XML(data)

    x = entrez_tree.xpath("//*[normalize-space(text()) = 'Official Full Name']")
    official_name = x[0].getnext().text if x else ''

    x = entrez_tree.xpath("//*[normalize-space(text()) = 'Also known as']")
    aka_names = x[0].getnext().text if x else ''

    aka_names = [a.strip() for a in aka_names.split(';')]

    return {'official_full_name': official_name,
            'AKA': aka_names}


def eutils_entrez(acc):
    acc = str(acc)

    crl = pycurl.Curl()

    crl.setopt(pycurl.URL,
               r'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=gene&id=%s&retmode=xml' % acc)

    response = cStringIO.StringIO()
    crl.setopt(pycurl.WRITEFUNCTION, response.write)

    crl.perform()

    gi_xml = ET.XML(response.getvalue())

    g = gi_xml.xpath('.//Link/Id/text()')
    if not g:
        return (0,'')
    else:
        geneID = g[0]

    crl.setopt(pycurl.URL,
               r'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s&retmode=xml' % geneID)

    response = cStringIO.StringIO()
    crl.setopt(pycurl.WRITEFUNCTION, response.write)

    crl.perform()

    gene_xml = ET.XML(response.getvalue())

    g = gene_xml.xpath('.//Item[@Name="Name"]/text()')
    if g:
        geneName = g[0]
    else:
        geneName = ''

    crl.close()

    return (geneID, geneName)


if __name__ == '__main__':
    acc = 171455
    print 'Testing MyGenBankHash, acc#: %d' % acc
    MyGenBankHash = download_genbank(acc)
    print MyGenBankHash['FEATURES'].keys()

    gid = 12345
    print 'Testing download_entrez, id#: %d' % gid
    d = download_entrez(gid)
    print d

    print 'Testing eutils_entrez, id#: %d' % gid
    e = eutils_entrez(gid)
    print e
