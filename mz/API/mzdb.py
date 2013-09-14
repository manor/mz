#Copyright 2012 Manor Askenazi

from mz.API import Scan, File as mzAPImzFile

import sqlite3
import struct

def _to_float(x):
    try :
        out = float(x)
    except ValueError :
        out = str(x)
    return out


class mzScan(list):
    """Subclass of list object for custom scan methods
    The mode can be 'p' or 'c' for profile or centroid respectively
    """

    def __init__(self, s, time, mode='p', mz=0.0, z=0):
        list.__init__(self, s)
        self.time = time
        self.mode = mode
        self.mz = mz
        self.z = z

    def peak(self, mz, tolerance):
        return max([i for m,i in self if abs(m-mz) <= tolerance] or [0])


class File(mzAPImzFile):
    """mzdb-backed mzFile"""

    def filters(self):
        conn = sqlite3.connect(self.data_file)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        if not self._filters :
            answer = []
            c.execute("""select time,filter from scans;""")
            for row in c:
                answer.append( (float(row[0]),row[1]) )
            self._filters = answer
        return self._filters

    def headers(self):
        '''Doesn't actually store the full headers. Generates the full scan_info
        list by looking at filter and header values. The results are cached.
        '''
        conn = sqlite3.connect(self.data_file)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        if not self._headers:
            self._headers = []
            c.execute("""select ID,time,precursor,filter,header from scans order by ID""") ;
            for (scan_num,rt,precursor,filter,header) in c:
                scan_mode = 'MS1'
                if filter.find("ms2") > -1:
                    scan_mode = 'MS2'

                self._headers.append((rt, # scan time
                                      precursor, # scan m/z from header, or 0 if MS1
                                      scan_num, # scan name == scan number
                                      scan_mode, # MS1 or MS2, referred to as 'scan type' in our API)
                                      'c')) # data mode, 'p' or 'c'
        return self._headers

    def __init__(self, data_file, **kwargs):
        """Initializes mzdb-based mzFile

        >>> dataFile = 'C:\\Documents and Settings\\User\\Desktop\\file.mzdb'
        >>> myFile = mzdb.mzFile(dataFile)

        """

        # At some point this should be set to 'mzdb' rather than pretending to be a raw file.
        # This is currently being done to protect from any multiplierz-y code that checks raw/wiff/mzML...
        self.file_type = 'raw'
        self.data_file = data_file

        self._filters = None
        self._headers = None

    def close(self):
        """Closes the open MS data file

        Example:
        >>> myPeakFile.close()

        """
        pass

    def scan_list(self, start_time=None, stop_time=None, start_mz=0, stop_mz=99999):
        """Gets a list of [(time,mz)] in the time and mz range provided

        All full MS scans that fall within the time range are included.
        Only MS/MS scans that fall within the mz range (optional) are included

        Example:
        >>> scan_list = my_peakfile.scan_list(30.0, 35.0, 435.82, 436.00)

        """
        if not start_time or not stop_time:
            (file_start_time, file_stop_time) = self.time_range()
        if not start_time:
            start_time = file_start_time
        if not stop_time:
            stop_time = file_stop_time

        return [(t,mz) for t,mz,sn,st,sm in self.headers()
                if start_time <= t <= stop_time and (st == 'MS1' or start_mz <= mz <= stop_mz)]

    def scan_info(self, start_time, stop_time=0, start_mz=0, stop_mz=99999):
        """Gets a list of [(time, mz, scan_name, scan_type, scan_mode)] in the time and mz range provided

        scan_name = number for RAW files, (cycle, experiment) for WIFF files.

        All full MS scans that fall within the time range are included.
        Only MS/MS scans that fall within the mz range (optional) are included

        Example:
        >>> scan_info = my_peakfile.scan_info(30.0, 35.0, 435.82, 436.00)
        """
        if stop_time == 0:
            stop_time = start_time

        return [(t,mz,sn,st,sm) for t,mz,sn,st,sm in self.headers()
                if start_time <= t <= stop_time and (st == 'MS1' or start_mz <= mz <= stop_mz)]

    def scan_time_from_scan_name(self, scan_name):
        """Essentially, gets the time for a raw scan number

        Example:
        >>> #raw
        >>> scan_time = myPeakFile.scan_time_from_scan_name(2165)

        """
        #At some point this becomes an anachronism... For mzdb files it is by fiat!!!
        return None

    def scan(self, time):
        """Gets scan based on the specified scan time

        The scan is a list of (mz, intensity) pairs.

        Example:
        >>> scan = myPeakFile.scan(20.035)

        """
        conn = sqlite3.connect(self.data_file)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        #Sneaky ability to access scans directly by scan_num... Shhhhh....
        if isinstance(time, float):
            scan_num = self.scanForTime(time)
        else:
            scan_num = time

        (scan_time,mz,sn,st,scan_mode) = self.headers()[scan_num - 1]
        c.execute("""select peak_num,mzs,ys from scans where time==?""",(scan_time,))
        r = c.fetchone()
        peak_num = r["peak_num"]
        mzs = struct.unpack("!%df"%peak_num,r["mzs"])
        ys = struct.unpack("!%df"%peak_num,r["ys"])
        sdata = zip(mzs,ys)
        z = self.extra_info(scan_num)["Charge State"]
        return Scan(sdata, scan_time, scan_mode, mz, z)

    def centroid(self, time, peakwidth=None):
        """Gets centroided scan based on the specified scan time

        The scan is a list of (mz, intensity) pairs.

        Example:
        >>> scan = myPeakFile.centroid(20.035)

        """
        return self.scan(time)

    def xic(self, start_time, stop_time, start_mz, stop_mz, filter=None):
        """Generates eXtracted Ion Chromatogram (XIC) for given time and mz range

        The function integrates the precursor intensities for given time and mz range.
        The xic is a list of (mz,intensity) pairs.

        Example:
        >>> xic = myPeakFile.xic(31.4, 32.4, 435.82, 436.00)

        """
        conn = sqlite3.connect(self.data_file)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        #s_start =self.scanForTime(start_time)
        #s_stop =self.scanForTime(stop_time)
        c.execute("""select time,peak_num,mzs,ys from scans where time>=? and time<=? and precursor = 0.0 ORDER BY time""",(start_time,stop_time)) ;
        xic = []
        for r in c:
            peak_num = r["peak_num"]
            s_time = r["time"]
            mzs = struct.unpack("!%df"%peak_num,r["mzs"])
            ys = struct.unpack("!%df"%peak_num,r["ys"])
            sdata = zip(mzs,ys)
            xic.append( (s_time,sum([y for (m,y) in sdata if m >= start_mz and m <= stop_mz])) )
        return xic

    def bpc(self):
        """Generates a Base Peak Chromatogram (BPC)

        Example:
        >>> xic = myPeakFile.bpc()

        """
        conn = sqlite3.connect(self.data_file)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        #s_start =self.scanForTime(start_time)
        #s_stop =self.scanForTime(stop_time)
        c.execute("""select time,bpc from scans where precursor = 0.0 ORDER BY time""") ;
        bpc = []
        for r in c:
            bpc.append( (r["time"],r["bpc"]) )
        return bpc

    def time_range(self):
        """Returns a pair of times corresponding to the first and last scan time

        Example:
        >>> time_range = myPeakFile.time_range()

        """
        start = self.headers()[0][0]
        stop = self.headers()[-1][0]
        return (start,stop)

    def scanForTime(self,the_time):
        bestDelta = 1000000.0
        scan_num = 0
        for (rt,precursor,scan_num,scan_mode,scan_type) in self.headers():
            if abs(rt-the_time) < bestDelta :
                bestDelta = abs(rt-the_time)
                bestScan = scan_num
            else:
                break
        return bestScan

    def timeForScan(self,the_scan):
        if the_scan < len(self.headers())+1 :
            return self.headers()[the_scan-1][0]
        else:
            return None

    def scan_range(self):
        start = self.headers()[0][2]
        stop = self.headers()[-1][2]
        return (start,stop)

    def lscan(self,scanNum):
        conn = sqlite3.connect(self.data_file)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        s_time = self.timeForScan(scanNum)
        sdata = []
        # no signal to noise...
        c.execute("""select peak_num,mzs,ys,zs from scans where time=?""" , (s_time,))
        r = c.fetchone()
        peak_num = r["peak_num"]
        mzs = struct.unpack("!%df"%peak_num,r["mzs"])
        ys = struct.unpack("!%df"%peak_num,r["ys"])
        zs = struct.unpack("!%dB"%peak_num,r["zs"])
        sdata = zip(mzs,ys,peak_num*[-1],zs)
        return sdata

    def cscan(self,scanNum):
        return self.lscan(scanNum)

    def extra_info(self,scanNum):
        conn = sqlite3.connect(self.data_file)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()
        str = c.execute("select header from scans where ID = ?",(scanNum,)).fetchone()[0]
        vals = str.split("\t")
        labels = vals[0::2]
        values = map(_to_float,vals[1::2])
        return dict(zip( labels, values) )

    def scanInjectionTime(self,scanNum):
        return self.extra_info(scanNum)["Ion Injection Time (ms)"]

    def scanPrecursor(self,scanNum):
        if isinstance(scanNum, float):
            the_scan = self.scanForTime(scanNum)
        else:
            the_scan = scanNum
        extra = self.extra_info(the_scan)
        print the_scan
        print extra["Monoisotopic M/Z"]
        return (extra["Monoisotopic M/Z"],int(extra["Charge State"]))
