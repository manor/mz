#Copyright 2012 Manor Askenazi
#This is a "light" version of the multiplierz platform.
#
import comtypes
from comtypes.client import CreateObject
from ctypes import *
import re

debug = True

def _to_float(x):
    try :
        out = float(x)
    except ValueError :
        out = str(x)
    return out

from mz.API import Scan, File as mzAPImzFile

# In multiplierz, we'll use the mzAPI.mzScan class.
# For stand-alone use of raw.py, uncomment this code

#class mzScan(list):
    #"""Subclass of list object for custom scan methods

    #The mode can be 'p' or 'c' for profile or centroid respectively

    #"""

    #def __init__(self, s, time, mode='p', mz=0.0, z=0):
        #list.__init__(self, s)
        #self.time = time
        #self.mode = mode
        #self.mz = mz
        #self.z = z

    #def peak(self, mz, tolerance):
        #return max([i for m,i in self if abs(m-mz) <= tolerance] or [0])

class File(mzAPImzFile):
    """msfilereader-based implementation of mzAPI's mzFile class"""

    def filters(self):
        if not self._filters :
            answer = []
            first = c_long()
            last  = c_long()
            scan  = c_long()

            retval = self.source.GetFirstSpectrumNumber(byref(first))
            if retval :
                if debug :
                    print retval
                return None

            retval = self.source.GetLastSpectrumNumber(byref(last))
            if retval :
                if debug :
                    print retval
                return None

            scan_time = c_double()
            for scan in range(first.value,last.value+1):
                filter = comtypes.automation.BSTR(None)
                retval = self.source.GetFilterForScanNum(scan,byref(filter))
                if retval :
                    if debug :
                        print retval
                    return None
                retval = self.source.RTFromScanNum(scan,byref(scan_time))
                if retval :
                    if debug :
                        print retval
                    return None
                answer.append( (scan_time.value, filter.value) )

            self._filters = answer
        return self._filters

    def headers(self):
        '''Doesn't actually store the full headers. Generates the full scan_info
        list by looking at filter and header values. The results are cached.
        '''

        if not self._headers:
            self._headers = []

            first = c_long()
            last  = c_long()
            scan_time = c_double()

            mode_re = re.compile(r'\s+(ms(?:\d+)?)\s+')
            data_re = re.compile(r'\s+([cp])\s+')
            mz_re = re.compile(r'(\d+\.\d+)@')

            retval = self.source.GetFirstSpectrumNumber(byref(first))
            if retval:
                if debug:
                    print retval
                return None

            retval = self.source.GetLastSpectrumNumber(byref(last))
            if retval:
                if debug:
                    print retval
                return None

            for scan in xrange(first.value, last.value + 1):
                filter_str = comtypes.automation.BSTR(None)

                retval = self.source.GetFilterForScanNum(scan, byref(filter_str))
                if retval:
                    if debug:
                        print retval
                    return None

                retval = self.source.RTFromScanNum(scan, byref(scan_time))
                if retval:
                    if debug:
                        print retval
                    return None

                data_m = data_re.search(filter_str.value)
                mode_m = mode_re.search(filter_str.value)

                scan_mode = str(mode_m.group(1)).upper()
                if scan_mode == 'MS':
                    scan_mode = 'MS1'

                if scan_mode == 'MS1':
                    mz = 0.0
                else:
                    header_fields = comtypes.automation.VARIANT()
                    num_fields = c_long()
                    retval = self.source.GetTrailerExtraLabelsForScanNum(scan, header_fields, comtypes.byref(num_fields))
                    if retval:
                        if debug:
                            print retval
                        return None

                    if 'Monoisotopic M/Z:' in header_fields.value:
                        mz_value = comtypes.automation.VARIANT()
                        retval = self.source.GetTrailerExtraValueForScanNum(scan, u'Monoisotopic M/Z:', mz_value)
                        if retval:
                            if debug:
                                print retval
                            return None
                        mz = mz_value.value
                    else:
                        mz = float(mz_re.search(filter_str.value).group(1))
                self._headers.append((scan_time.value, # scan time
                                      mz, # scan m/z from header, or 0 if MS1
                                      scan, # scan name == scan number
                                      scan_mode, # MS1 or MS2, referred to as 'scan type' in our API)
                                      str(data_m.group(1)).lower() if data_m else 'p')) # data mode, 'p' or 'c'
        return self._headers

    def __init__(self, data_file, **kwargs):
        """Initializes mzAPI and opens a new file

        >>> dataFile = 'C:\\Documents and Settings\\User\\Desktop\\rawFile.RAW'
        >>> myPeakFile = mzAPI.mzFile(dataFile)

        """

        self.file_type = 'raw'
        self.data_file = data_file
        self.source = None

        try:
            obj = CreateObject("MSFileReader.XRawfile")
        except WindowsError:
            obj = CreateObject("XRawfile.XRawfile")

        self.source = obj

        retval = obj.Open(data_file)
        if retval:
            if debug :
                print retval

        retval = obj.SetCurrentController(c_long(0),c_long(1))
        if retval:
            obj.Close()
            if debug :
                print retval

        self._filters = None
        self._headers = None

    def close(self):
        """Closes the open MS data file

        Example:
        >>> myPeakFile.close()

        """
        self.source.Close()

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

        scan = c_long(scan_name)
        scantime = c_double()
        retval = self.source.RTFromScanNum(scan,byref(scantime))
        if retval :
            if debug :
                print retval
            return None

        return scantime.value

    def scan(self, time):
        """Gets scan based on the specified scan time

        The scan is a list of (mz, intensity) pairs.

        Example:
        >>> scan = myPeakFile.scan(20.035)

        """

        #Sneaky ability to access scans directly by scan_num... Shhhhh....
        if isinstance(time, float):
            scan_num = self.scanForTime(time)
        else:
            scan_num = time

        (scan_time,mz,sn,st,scan_mode) = self.headers()[scan_num - 1]
        the_scan = c_long(scan_num)

        cpw = c_double()
        peaknum = c_long()
        ms = comtypes.automation.VARIANT()
        flags = comtypes.automation.VARIANT()
        retval = self.source.GetMassListFromScanNum(the_scan,None,0,0,0,False,cpw,ms,flags,peaknum)
        if retval :
            if debug :
                print retval
            return None

        z_value = comtypes.automation.VARIANT()
        retval = self.source.GetTrailerExtraValueForScanNum(the_scan, u'Charge State:', z_value)
        if retval:
            if debug:
                print retval
            return None
        z = z_value.value

        return Scan(zip(ms.value[0],ms.value[1]), scan_time, scan_mode, mz, z)

    def centroid(self, time, peakwidth=None):
        """Gets centroided scan based on the specified scan time

        The scan is a list of (mz, intensity) pairs.

        Example:
        >>> scan = myPeakFile.centroid(20.035)

        """

        #Sneaky ability to access scans directly by scan_num... Shhhhh....
        if isinstance(time, float):
            scan_num = self.scanForTime(time)
        else:
            scan_num = time

        (scan_time,mz,sn,st,scan_mode) = self.headers()[scan_num - 1]
        the_scan = c_long(scan_num)

        if peakwidth:
            cpw = c_double(peakwidth)
        else:
            cpw = c_double()

        peaknum = c_long()
        ms = comtypes.automation.VARIANT()
        flags = comtypes.automation.VARIANT()
        retval = self.source.GetMassListFromScanNum(the_scan,None,0,0,0,True,cpw,ms,flags,peaknum)
        if retval :
            if debug :
                print retval
            return None

        z_value = comtypes.automation.VARIANT()
        retval = self.source.GetTrailerExtraValueForScanNum(the_scan, u'Charge State:', z_value)
        if retval:
            if debug:
                print retval
            return None
        z = z_value.value

        return Scan(zip(ms.value[0],ms.value[1]), scan_time, scan_mode, mz, z)

    def xic(self, start_time, stop_time, start_mz, stop_mz, filter=None):
        """Generates eXtracted Ion Chromatogram (XIC) for given time and mz range

        The function integrates the precursor intensities for given time and mz range.
        The xic is a list of (mz,intensity) pairs.

        Example:
        >>> xic = myPeakFile.xic(31.4, 32.4, 435.82, 436.00)

        """

        if filter is None:
            filter =  "Full ms "
        massRange = "%s-%s" % (start_mz, stop_mz)

        ftLB = c_double(start_time)
        ftUB = c_double(stop_time)
        cdata = comtypes.automation.VARIANT()
        flags = comtypes.automation.VARIANT()
        val_num = c_long()

        retval = self.source.GetChroData(0,0,0,filter,massRange,"",0.0,byref(ftLB),byref(ftUB),0,0,cdata,flags,byref(val_num))
        if retval :
            if debug :
                print retval
            return None
        return tuple(zip(cdata.value[0],cdata.value[1]))

    def time_range(self):
        """Returns a pair of times corresponding to the first and last scan time

        Example:
        >>> time_range = myPeakFile.time_range()

        """

        start = c_double()
        stop = c_double()

        retval = self.source.GetStartTime(byref(start))
        if retval :
            if debug :
                print retval
            return None

        retval = self.source.GetEndTime(byref(stop))
        if retval :
            if debug :
                print retval
            return None
        return (start.value,stop.value)

    def scanForTime(self,the_time):
        the_scan = c_long()
        retval = self.source.ScanNumFromRT(c_double(the_time),byref(the_scan))
        if retval :
            if debug :
                print retval
            return None
        return the_scan.value

    def timeForScan(self,the_scan):
        the_time = c_double()
        retval = self.source.RTFromScanNum(c_long(the_scan),byref(the_time))
        if retval :
            if debug :
                print retval
            return None
        return the_time.value


    def scan_range(self):
        start = c_long()
        stop = c_long()

        retval = self.source.GetFirstSpectrumNumber(byref(start))
        if retval :
            if debug :
                print retval
            return None

        retval = self.source.GetLastSpectrumNumber(byref(stop))
        if retval :
            if debug :
                print retval
            return None
        return (start.value,stop.value)

    def lscan(self,scanNum):
        the_scan = c_long(scanNum)
        ms = comtypes.automation.VARIANT()
        flags = comtypes.automation.VARIANT()
        retval = self.source.GetLabelData(ms,flags,the_scan)
        if retval :
            if debug :
                print retval
            return None
        return zip(ms.value[0],ms.value[1],ms.value[4],ms.value[5])

    def cscan(self,scanNum):
        return self.lscan(scanNum)

    def extra_info(self,scanNum):
        the_scan = c_long(scanNum)
        labels = comtypes.automation.VARIANT()
        values = comtypes.automation.VARIANT()
        val_num = c_long()

        retval = self.source.GetTrailerExtraForScanNum(the_scan,labels,values,val_num)
        if retval :
            if debug :
                print retval
            return None
        return dict(zip( map(lambda x: str(x[:-1]), labels.value) , map(_to_float, values.value) ))

    def scanInjectionTime(self,scanNum):
        the_scan = c_long(scanNum)
        labels = comtypes.automation.VARIANT()
        values = comtypes.automation.VARIANT()
        val_num = c_long()

        retval = self.source.GetTrailerExtraForScanNum(the_scan,labels,values,val_num)
        if retval :
            if debug :
                print retval
            return None
        vals = dict(zip( map(lambda x: str(x[:-1]), labels.value) , map(_to_float, values.value) ))
        return vals["Ion Injection Time (ms)"]

    def scanPrecursor(self,scanNum):
        if isinstance(scanNum, float):
            the_scan = c_long(self.scanForTime(scanNum))
        else:
            the_scan = c_long(scanNum)

        labels = comtypes.automation.VARIANT()
        values = comtypes.automation.VARIANT()
        val_num = c_long()

        retval = self.source.GetTrailerExtraForScanNum(the_scan,labels,values,val_num)
        if retval :
            if debug :
                print retval
            return None
        vals = dict(zip( map(lambda x: str(x[:-1]), labels.value) , map(_to_float, values.value) ))
        return (vals["Monoisotopic M/Z"],int(vals["Charge State"]))
