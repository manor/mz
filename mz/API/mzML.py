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

import lxml.etree as etree
import os
import re

import base64
import cPickle
import cStringIO
import gzip
import struct
import zlib

# namespace for mzML. Shows up as prefix before every tag in the tree
NSd = {'mz': "http://psi.hupo.org/ms/mzml"}
NS = ("{http://psi.hupo.org/ms/mzml}",)

from multiplierz.mzAPI import mzScan, mzInfoFile, mzFile as mzAPImzFile

class FileRange:
    '''Wraps a file-like object and present a subset of the file--
    the idea is to allow accessing a slice of an mzML without reading
    the whole thing into memory as a StringIO instance.

    Some trepidation about not implementating the whole file-object
    interface here. It seems to work in its intended application,
    though. etree.iterparse only calls the read() method.
    '''
    def __init__(self, file_ob, start, stop):
        self.fo = file_ob
        self.start = start
        self.stop = stop
        self.fo.seek(start)
        self.cur = start - 1

    def read(self, size=None):
        '''Reads data within the given range, truncating if it
        gets to the end, and wrapping the whole thing with a
        'spectrumList' element because the parser wants that.
        '''
        if self.start <= self.cur < self.stop:
            r = self.fo.read(size)
            if len(r) + self.cur > self.stop:
                r = r[:(self.stop - self.cur)]
                self.cur = self.stop
            else:
                self.cur = self.fo.tell()
            return r
        elif self.cur < self.start:
            self.cur += 1
            return '<spectrumList>\n'
        elif self.cur == self.stop:
            self.cur += 1
            return '</spectrumList>\n'
        else:
            return ''


class InfoBuilderTarget:
    '''Builds a list of spectrum offsets as it reads the file.

    Takes a file-like object handle on creation, which it uses
    to get the current location whenever it encounters the end
    of a spectrum element. Meant to be fed a file line by line
    using the parser.feed(data) interface.
    '''
    offsets = []

    def __init__(self, fh):
        self.fh = fh

    def end(self, tag):
        if tag == 'spectrum':
            self.offsets.append(self.fh.tell())

    def close(self):
        offsets, self.offsets = self.offsets[:], []
        return offsets


def make_info_file(data_file):
    '''Makes an info file for an mzML file, so that we don't have
    to iterate through the whole thing every time.
    '''

    # building our own table is about as fast as reading through
    # the file and extracting theirs, so to simplify code we'll
    # just do it our way (because we need that function anyway)
    if os.path.exists(data_file + '.mzi'):
        os.remove(data_file + '.mzi')

    m = mzFile(data_file)
    m._build_info_scans()

    fh = open(data_file + '.mzi', 'wb')
    cPickle.dump(m._info_scans, fh)
    fh.close()


class mzFile(mzAPImzFile):
    """Class for access to mzML data files. This class doesn't load the
    file into memory--instead, it uses lxml's 'iterparse' class to
    create a specialized file-reader, which only loads one scan at a
    time. This uses very little memory (an insignificant amount when
    compared to multiplierz in general) and is not too slow.

    For a large number of accesses (extracting XICs, for instance), it's
    usually worth building a cheat-sheet of scan locations--for a one-time
    use, there's the _build_info_scans method, but in general it makes
    more sense to create an .mzi file.
    """

    _xp_time = etree.XPath(('./mz:scanList/mz:scan/'
                            'mz:cvParam[@name="scan start time"]/'
                            'attribute::value'),
                           namespaces=NSd, smart_strings=False)
    _xp_ms2 = etree.XPath('./mz:cvParam[@name="ms level"][@value="2"]',
                          namespaces=NSd, smart_strings=False)
    _xp_pre = etree.XPath(('.//mz:selectedIonList/mz:selectedIon/'
                           'mz:cvParam[@name="selected ion m/z"]/'
                           'attribute::value'),
                          namespaces=NSd, smart_strings=False)
    _xp_prof = etree.XPath('./mz:cvParam[@name="profile spectrum"]',
                           namespaces=NSd, smart_strings=False)

    def __init__(self, data_file, **kwargs):
        self.file_type = 'mzml'
        self.data_file = data_file

        if data_file.lower().endswith('.mzml.gz'):
            self.fileobj = gzip.GzipFile(data_file, mode='rb')
        else:
            self.fileobj = open(data_file, mode='rb')

        if os.path.exists(data_file + '.mzi'):
            self._info_file = data_file + '.mzi'
            info_fh = open(self._info_file)
            self._info_scans = cPickle.load(info_fh)
            info_fh.close()
        else:
            self._info_file = None

    def close(self):
        self.fileobj.close()

    def scan_list(self, start_time=None, stop_time=None, start_mz=0, stop_mz=99999):
        if start_time is None:
            start_time = -1.0
        # rather than guess the stop_time, if it's None we'll just ignore it

        if self._info_file:
            return [ (i['time'], i['mz']) for i in self._info_scans
                     if (start_time <= i['time'] and ((not stop_time) or i['time'] <= stop_time))
                     and (i['scan_mode'] == 'MS1' or start_mz <= i['mz'] <= stop_mz) ]

        scan_list = []

        self.fileobj.seek(0)
        context = etree.iterparse(self.fileobj, events=('end',),
                                  tag='%sspectrum' % NS)
        for event, elem in context:
            xt = self._xp_time(elem)
            if xt:
                time = float(xt[0])
                if ((start_time <= time)
                    and ((not stop_time)
                         or time <= stop_time)):
                    if self._xp_ms2(elem):
                        p = self._xp_pre(elem)
                        if p:
                            mz = float(p[0])
                            if start_mz <= mz <= stop_mz:
                                scan_list.append((time,mz))
                        else:
                            print "this ms2 didn't have precursor mz...", elem.get("id")
                    else:
                        scan_list.append((time,0.0))
            else:
                print "this spectrum didn't have a scan time...", elem.get("id")

            elem.clear()

        return scan_list

    def scan_info(self, start_time, stop_time=0, start_mz=0, stop_mz=99999):
        # scan_name is the spectrum id attribute which looks like:
        #    controllerType=0 controllerNumber=1 scan=1
        # (for RAW-file-derived mzML files)

        if stop_time == 0:
            stop_time = start_time

        if self._info_file:
            keys = ('time', 'mz', 'scan_name', 'scan_type', 'scan_mode')
            return [ tuple(i[k] for k in keys) for i in self._info_scans
                     if (start_time <= i['time'] <= stop_time)
                     and (i['scan_mode'] == 'MS1' or start_mz <= i['mz'] <= stop_mz) ]

        scan_info = []

        self.fileobj.seek(0)
        context = etree.iterparse(self.fileobj, events=('end',),
                                  tag='%sspectrum' % NS)
        for event, elem in context:
            xt = self._xp_time(elem)
            if xt:
                time = float(xt[0])
                if start_time <= time <= stop_time:
                    if self._xp_ms2(elem):
                        p = self._xp_pre(elem)
                        if p:
                            mz = float(p[0])
                            if start_mz <= mz <= stop_mz:
                                if self._xp_prof(elem):
                                    scan_mode = 'p'
                                else:
                                    scan_mode = 'c'

                                scan_name = elem.get('id')
                                scan_info.append((time, mz, scan_name, 'MS2', scan_mode))
                        else:
                            print "this ms2 didn't have precursor mz...", elem.get("id")
                    else:
                        if self._xp_prof(elem):
                            scan_mode = 'p'
                        else:
                            scan_mode = 'c'

                        scan_name = elem.get('id')
                        scan_info.append((time, 0.0, scan_name, 'MS1', scan_mode))
            else:
                print "this spectrum didn't have a scan time...", elem.get("id")
            elem.clear()

        return scan_info

    def scan_time_from_scan_name(self, scan_name):
        # scan_name is the spectrum id attribute which looks like:
        #    controllerType=0 controllerNumber=1 scan=1
        # (for RAW-file-derived mzML files)

        if self._info_file:
            return self._info_scans.filter(list_key='scan_name', value_list=[scan_name])[0]['time']

        self.fileobj.seek(0)
        context = etree.iterparse(self.fileobj, events=('end',),
                                  tag='%sspectrum' % NS)
        for event, elem in context:
            if elem.get('id') == scan_name:
                x = self._xp_time(elem)
                if x:
                    return float(x[0])
                else:
                    print "this spectrum didn't have a scan time"
            elem.clear()
        else:
            print "scan not found:", scan_name

    def scan(self, time):
        # this implementation takes a few seconds and uses very little memory
        # (by keeping only the current closest scan)

        if self._info_file:
            closest_item = self._info_scans.closest(key='time', value=time)
            spec_start, spec_size = closest_item['offset'], closest_item['size']
            self.fileobj.seek(spec_start)
            spec = etree.XML(self.fileobj.read(spec_size))
            return self._scan_from_spec_node(spec, closest_item['time'], prefix=False)

        self.fileobj.seek(0)
        context = etree.iterparse(self.fileobj, events=('end',),
                                  tag='%sspectrum' % NS)
        event,elem = context.next()
        min_dtime = abs(time - float((self._xp_time(elem) or (0.0,))[0]))
        closest_elem = elem

        for event, elem in context:
            x = self._xp_time(elem)
            if x:
                x = float(x[0])
                if x == time:
                    return self._scan_from_spec_node(elem, x)
                elif abs(time - x) < min_dtime:
                    min_dtime = abs(time - x)
                    closest_elem.clear()
                    closest_elem = elem
                else:
                    elem.clear()
            else:
                elem.clear()
        else:
            return self._scan_from_spec_node(closest_elem, float((self._xp_time(closest_elem) or (0.0,))[0]))

    def xic(self, start_time, stop_time, start_mz, stop_mz, filter=None):
        # ignoring filter parameter, for now at least.

        if self._info_file:
            items = [ item for item in self._info_scans if start_time <= item['time'] <= stop_time ]
            if not items:
                return [] # no scans in time range
            items.sort(key=lambda i: i['offset'])
            start_offset = items[0]['offset']
            stop_offset = items[-1]['offset'] + items[-1]['size']

            fileobj = FileRange(self.fileobj, start_offset, stop_offset)
            p = ('',)
        else:
            self.fileobj.seek(0)
            fileobj = self.fileobj
            p = NS

        xic_data = []

        context = etree.iterparse(fileobj, events=('end',),
                                  tag='%sspectrum' % p)
        for event, elem in context:
            if elem.find('./%scvParam[@name="ms level"][@value="1"]' % p) is not None:
                t = elem.find('./%sscanList/%sscan/%scvParam[@name="scan start time"]' % (p*3))
                if t is not None:
                    time = float(t.get('value'))
                    if start_time <= time <= stop_time:
                        scan = self._scan_from_spec_node(elem, time, prefix=(p==NS))
                        xic_data.append((time, sum(i for mz,i in scan if start_mz <= mz <= stop_mz)))
                else:
                    print "this spectrum didn't have a scan time...", s.get("id")
            elem.clear()

        xic_data.sort()
        return xic_data

    def time_range(self):
        if self._info_file:
            times = [ i['time'] for i in self._info_scans ]
        else:
            self.fileobj.seek(0)
            context = etree.iterparse(self.fileobj, events=('end',),
                                      tag='%sspectrum' % NS)
            parse_e = (lambda e: (self._xp_time(e), e.clear()))

            times = [ float(x[0]) for x in (parse_e(elem)[0]
                                            for event,elem in context)
                      if x ]

        return (min(times), max(times)) if times else (0.0, 0.0)

    def _build_info_scans(self):
        '''Build the info_scans object in memory--this allows the speed-up
        of having an .mzi file. For a file that's going to be accesssed many
        times, it makes more sense to just build the .mzi, but this might be
        useful if a lot of data is being pulled from a file, but only once.
        '''
        if self._info_file:
            return # self._info_scans is already present

        (start, end) = self.time_range()
        scan_info = self.scan_info(start, end)

        self.fileobj.seek(0)

        line = ''
        specre = re.compile(r'<\s*spectrumList')
        while not specre.search(line):
            line = self.fileobj.readline()

        offsets = [ self.fileobj.tell() ] # first offset

        # parser takes lines and accumulates offset list
        parser = etree.XMLParser(target = InfoBuilderTarget(self.fileobj))
        parser.feed(line)

        # stop when we reach the end of the spectrum list
        specre = re.compile(r'</\s*spectrumList\s*>')
        while not specre.search(line):
            line = self.fileobj.readline()
            parser.feed(line)

        #  get results from the parser
        offsets.extend(parser.close()) # has one extra offset
        sizes = [ (n - offsets[i]) for i,n in enumerate(offsets[1:]) ]
        # get rid of last offset (which was the end of the spectrumList)
        offsets.pop()

        keys = ('time', 'mz', 'scan_name', 'scan_type',
                'scan_mode', 'offset', 'size')

        self._info_file = ":memory:"
        self._info_scans = mzInfoFile(dict(zip(keys, si + (offsets[i], sizes[i])))
                                      for i,si in enumerate(scan_info))

    def _scan_from_spec_node(self, spec, scan_time, prefix=True):
        '''Gets the mzScan data from a 'spectrum' mzML node. Does some error-checking
        due to the ambiguous nature of mzML's spec--no guarantees about what's in the
        tree.
        '''
        if prefix:
            p = NS
        else:
            p = ('',)

        if spec.find('./%scvParam[@name="profile spectrum"]' % p) is not None:
            scan_mode = 'p'
        else:
            scan_mode = 'c'

        array_length = int(spec.get('defaultArrayLength'))

        mz_array = None
        int_array = None

        for bin in spec.iterfind('%sbinaryDataArrayList/%sbinaryDataArray' % (p*2)):
            if bin.find('%scvParam[@name="64-bit float"]' % p) is not None:
                fmt = '%dd' % array_length
            else:
                fmt = '%df' % array_length

            if bin.find('%scvParam[@name="no compression"]' % p) is not None:
                compression = False
            else:
                compression = True

            if bin.find('%scvParam[@name="m/z array"]' % p) is not None:
                if mz_array:
                    print "Overwriting m/z array!?"
                if compression:
                    mz_array = struct.unpack(fmt,
                                             zlib.decompress(base64.standard_b64decode(bin.find('%sbinary' % p).text)))
                else:
                    mz_array = struct.unpack(fmt,
                                             base64.standard_b64decode(bin.find('%sbinary' % p).text))
            elif bin.find('%scvParam[@name="intensity array"]' % p) is not None:
                if int_array:
                    print "Overwriting intensity array!?"
                if compression:
                    int_array = struct.unpack(fmt,
                                              zlib.decompress(base64.standard_b64decode(bin.find('%sbinary' % p).text)))
                else:
                    int_array = struct.unpack(fmt,
                                              base64.standard_b64decode(bin.find('%sbinary' % p).text))
            else:
                print ("Found some other kind of binary array in here",
                       [b.get('name') for b in bin.iterfind('%scvParam[@name]' % p)])

        if (mz_array and int_array):
            return mzScan(zip(mz_array, int_array), scan_time, mode=scan_mode)
        else:
            if mz_array:
                print "No intensity values"
            elif int_array:
                print "No m/z values"
            else:
                print "No m/z or intensity values"
            return None


class mzFileInMemory:
    """Class for access to mzML data files. Loads the entire tree into memory,
    which takes a while and uses a lot of space. Access is then faster. For
    small files and/or machines with huge amounts of RAM, this class might be
    a better choice.
    """

    def __init__(self, data_file, **kwargs):
        self.file_type = 'mzml'
        self.data_file = data_file

        if data_file.lower().endswith('.mzml.gz'):
            self.fileobj = etree.parse(gzip.GzipFile(data_file, mode='rb'))
        else:
            self.fileobj = etree.parse(data_file)

    def close(self):
        self.fileobj = None
        # now there's nothing to do, because there's no handle on the mzML file
        #pass

    def scan_list(self, start_timeNone, stop_time=None, start_mz=0, stop_mz=99999):
        if start_time is None or stop_time is None:
            (file_start_time, file_stop_time) = self.time_range()
        if start_time is None:
            start_time = file_start_time
        if stop_time is None:
            stop_time = file_stop_time

        scan_list = []

        for s in self.fileobj.iterfind('.//%srun/%sspectrumList/%sspectrum' % (NS*3)):
            st = s.find('%sscanList/%sscan/%scvParam[@name="scan start time"]' % (NS*3))
            if st is not None:
                time = float(st.get("value"))
                if start_time <= time <= stop_time:
                    m = s.find('%scvParam[@name="ms level"][@value="2"]' % NS)
                    if m is not None:
                        p = s.find('.//%sselectedIonList/%sselectedIon/%scvParam[@name="selected ion m/z"]' % (NS*3))
                        if p is not None:
                            mz = float(p.get("value"))
                            if start_mz <= mz <= stop_mz:
                                scan_list.append((time, mz))
                        else:
                            print "this ms2 didn't have precursor mz...", s.get("id")
                    else:
                        scan_list.append((time, 0.0))
            else:
                print "this spectrum didn't have a scan time...", s.get("id")

        return scan_list

    def scan_info(self, start_time, stop_time=0, start_mz=0, stop_mz=99999):
        # scan_name is the spectrum id attribute which looks like:
        #    controllerType=0 controllerNumber=1 scan=1
        # (for RAW-file-derived mzML files)

        if stop_time == 0:
            stop_time = start_time

        scan_info = []

        for s in self.fileobj.iterfind('.//%srun/%sspectrumList/%sspectrum' % (NS*3)):
            st = s.find('%sscanList/%sscan/%scvParam[@name="scan start time"]' % (NS*3))
            if st is not None:
                time = float(st.get("value"))
                if start_time <= time <= stop_time:
                    m = s.find('%scvParam[@name="ms level"][@value="2"]' % NS)
                    if m is not None:
                        p = s.find('.//%sselectedIonList/%sselectedIon/%scvParam[@name="selected ion m/z"]' % (NS*3))
                        if p is not None:
                            mz = float(p.get("value"))
                            if start_mz <= mz <= stop_mz:
                                if s.find('%scvParam[@name="profile spectrum"]' % NS) is not None:
                                    scan_mode = 'p'
                                else:
                                    scan_mode = 'c'

                                scan_name = s.get('id')
                                scan_info.append((time, mz, scan_name, 'MS2', scan_mode))
                        else:
                            print "this ms2 didn't have precursor mz...", s.get("id")
                    else:
                        if s.find('%scvParam[@name="profile spectrum"]' % NS) is not None:
                            scan_mode = 'p'
                        else:
                            scan_mode = 'c'

                        scan_name = s.get('id')
                        scan_info.append((time, 0.0, scan_name, 'MS1', scan_mode))
            else:
                print "this spectrum didn't have a scan time...", s.get("id")

        return scan_info

    def scan_time_from_scan_name(self, scan_name):
        spec = self.fileobj.find('.//%srun/%sspectrumList/%sspectrum[@id="%s"]' % (NS*3+(scan_name,)))
        if spec is not None:
            st = s.find('%sscanList/%sscan/%scvParam[@name="scan start time"]' % (NS*3))
            if st is not None:
                return float(st.get("value"))
            else:
                print "this spectrum didn't have a scan time"
        else:
            print "scan not found:", scan_name

    def scan(self, time):
        search_term = ('.//%srun/%sspectrumList/%sspectrum/%sscanList/%sscan/'
                       '%scvParam[@name="scan start time"]') % (NS*6)

        min_cv = min((cv for cv in self.fileobj.iterfind(search_term)),
                     key = lambda c: abs(time - float(c.get('value'))))

        spec = min_cv.getparent().getparent().getparent()

        return self._scan_from_spec_node(spec, float(min_cv.get('value')))

    def xic(self, start_time, stop_time, start_mz, stop_mz, filter=None):
        # something is wrong with this, it returns values that are consistently different from
        # raw.xic and the proteowizard TICs (which agree with each other).

        # Answer: raw.xic is using the centroid value (which it can calculate), but
        # mzML just records the raw data and we sum it (as does ProteoWizard).

        # ignoring filter parameter, for now at least.

        xic_data = []

        search_term = './/%srun/%sspectrumList/%sspectrum/%scvParam[@name="ms level"][@value="1"]' % (NS*4)

        for cv in self.fileobj.iterfind(search_term):
            s = cv.getparent()
            st = s.find('%sscanList/%sscan/%scvParam[@name="scan start time"]' % (NS*3))
            if st is not None:
                time = float(st.get("value"))
                if start_time <= time <= stop_time:
                    scan = self._scan_from_spec_node(s, time)
                    xic_data.append((time, sum(i for mz,i in scan if start_mz <= mz <= stop_mz)))
            else:
                print "this spectrum didn't have a scan time...", s.get("id")

        return xic_data

    def time_range(self):
        search_term = './/%sspectrum/%sscanList/%sscan/%scvParam[@name="scan start time"]' % (NS*4)
        times = [ float(t.get("value")) for t in self.fileobj.iterfind(search_term) ]

        return (min(times), max(times))

    def _scan_from_spec_node(self, spec, scan_time):
        '''Gets the mzScan data from a 'spectrum' mzML node. Right now, does some
        error-checking due to the ambiguous nature of mzML's spec--no guarantees
        about what's in the tree.
        '''
        if spec.find('%scvParam[@name="profile spectrum"]' % NS) is not None:
            scan_mode = 'p'
        else:
            scan_mode = 'c'

        array_length = int(spec.get('defaultArrayLength'))

        mz_array = None
        int_array = None

        for bin in spec.iterfind('%sbinaryDataArrayList/%sbinaryDataArray' % (NS*2)):
            if bin.find('%scvParam[@name="64-bit float"]' % NS) is not None:
                fmt = 'd' * array_length
            else:
                fmt = 'f' * array_length

            if bin.find('%scvParam[@name="no compression"]' % NS) is not None:
                compression = False
            else:
                compression = True

            if bin.find('%scvParam[@name="m/z array"]' % NS) is not None:
                if mz_array:
                    print "Overwriting m/z array!?"
                if compression:
                    mz_array = struct.unpack(fmt,
                                             zlib.decompress(base64.standard_b64decode(bin.find('%sbinary' % NS).text)))
                else:
                    mz_array = struct.unpack(fmt,
                                             base64.standard_b64decode(bin.find('%sbinary' % NS).text))
            elif bin.find('%scvParam[@name="intensity array"]' % NS) is not None:
                if int_array:
                    print "Overwriting intensity array!?"
                if compression:
                    int_array = struct.unpack(fmt,
                                              zlib.decompress(base64.standard_b64decode(bin.find('%sbinary' % NS).text)))
                else:
                    int_array = struct.unpack(fmt,
                                              base64.standard_b64decode(bin.find('%sbinary' % NS).text))
            else:
                print ("Found some other kind of binary array in here",
                       [b.get('name') for b in bin.iterfind('%scvParam[@name]' % NS)])

        if (mz_array and int_array):
            return mzScan(zip(mz_array, int_array), scan_time, mode=scan_mode)
        else:
            if mz_array:
                print "No intensity values"
            elif int_array:
                print "No m/z values"
            else:
                print "No m/z or intensity values"
            return None
