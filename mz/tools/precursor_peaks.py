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

from tempfile import mkstemp

from scipy import mean, median, pi, sqrt

import mz.fit as fit
import mz.API as API

from mz import logger_message
#from mz import myData, myTemp, logger_message
from mz.tools import mz_image
from mz.functions import mz_pep_format


def extract_peaks(peak_data_path, time_window=(0.5, 0.5), mz_window=(0.1, 0.1),
                  plot_ms1=False, plot_xic=False, peak_area=False,
                  reporter_ions=False, peakfilter=None, plot_ms2=False,
                  ion_list=('b', 'y'), instrument='ESI-TRAP', isMZD=False, im_size=(8.0,6.0)):

    mz_file = API.mzFile(peak_data_path)
    peak_file_type = mz_file.file_type

    if instrument == 'ETD-TRAP':
        ion_list = ('c', 'z')

    logger_message(10, 'MS File opened...')

    (first_time,last_time) = mz_file.time_range()

    logger_message(30, 'Accessing Precursor Peaks ...')

    # all the MS2 scans
    ms2_scans = [(t,mz) for t,mz,sn,st,sm in mz_file.scan_info(first_time, last_time)
                 if st == 'MS2']

    raw_reg = re.compile(r'.+\.(\d+)\.(\d+)\.(.+)\.dta.*', flags=re.I)
    wiff_reg = re.compile(r'.*\(sample number (\d+)\), Elution: (.+) min, '
                          'Period: (\d+), Cycle\(s\): (\d+).*? \(Experiment (\d+)\).*?', flags=re.I)
    wiff_reg2 = re.compile(r'(?:Locus:)?(\d+)\.(\d+)\.(\d+)\.(\d+)\.(\d+)', flags=re.I)

    # doubles as a 'generic' scan regex {file_name}/scans/{scan_time}
    url_reg = re.compile(r'((?:http://.+/files/)?(?:.+)/scans/(\d+(?:\.\d+)?))(?:/.*)?', flags=re.I)

    row = (yield None) # generator initialization

    try:
        while row:
            image_tuples = []

            mz = row['Experimental mz']
            scan_time = row.get('ms2 time', None)

            if scan_time:
                scan = scan_time
            else:
                spec_desc = row['Spectrum Description']

                raw_m = raw_reg.match(spec_desc)
                wiff_m = wiff_reg.match(spec_desc)
                wiff_m2 = wiff_reg2.match(spec_desc)
                url_m = url_reg.match(spec_desc)
                if raw_m:
                    scan = int(raw_m.group(1))
                    if scan == 0:
                        row = yield (row, image_tuples)
                    else:
                        scan_time = mz_file.scan_time_from_scan_name(scan)
                elif wiff_m:
                    scan = (int(wiff_m.group(4)), int(wiff_m.group(5)))
                    scan_time = mz_file.scan_time_from_scan_name(scan)
                elif wiff_m2:
                    scan = (int(wiff_m2.group(4)), int(wiff_m.group(5)))
                    scan_time = mz_file.scan_time_from_scan_name(scan)
                elif url_m:
                    scan = url_m.group(1)
                    scan_time = float(url_m.group(2))
                else:
                    row = yield (row, image_tuples)

            logger_message(10,'Scan = %s, Scan Time = %s, mz = %s' % (str(scan), str(scan_time), str(mz)))

            (time_first_half, time_second_half) = time_window
            (mz_first_half, mz_second_half) = mz_window

            floatmz = float(mz)

            start_mz = floatmz - mz_first_half
            end_mz = floatmz + mz_second_half

            start_time = max(scan_time - time_first_half, first_time + 0.00001)
            end_time = min(scan_time + time_second_half, last_time - 0.00001)

            logger_message(20, 'Making XIC..')
            logger_message(10, '%s %s' % (str(start_time), str(end_time)))

            xic = mz_file.xic(start_time, end_time, start_mz, end_mz, peakfilter)
            (max_time, max_int) = max(xic, key=lambda x: x[1])
            max_xic_index = xic.index((max_time,max_int))

            # plot XIC
            if plot_xic:
                scan_dot = (scan_time, xic[-1][1])
                for (i,r) in enumerate(xic[1:]):
                    if r[0] > scan_time:
                        # indexes are shifted because of slice: xic[i] is the previous xic
                        Slope = float(xic[i][1] - r[1]) / float(xic[i][0] - r[0])
                        scan_dot = (scan_time,
                                    Slope * (scan_time - r[0]) + r[1])
                        break

                xic_times = [i[0] for i in xic] # X axis
                xic_ints = [i[1] for i in xic] # Y axis


                bin_times = [a[0] for a in ms2_scans
                             if (a[0] != scan_time
                                 and start_time < a[0] <= end_time
                                 and start_mz <= a[1] <= end_mz)]
                bin_ints = []

                for bt in bin_times:
                    for (i,r) in enumerate(xic[1:]):
                        if r[0] > bt:
                            # indexes are shifted because of slice: xic[i] is the previous xic
                            Slope = float(xic[i][1] - r[1]) / float(xic[i][0] - r[0])
                            bin_ints.append(Slope * (bt - r[0]) + r[1])
                            break
                        if i == len(xic) - 2:
                            bin_ints.append(r[1])

                if isMZD:
                    image_tuples.append(('Peak Width (sec)', 'xic',
                                         (mz, xic_times, xic_ints, scan_dot, bin_times, bin_ints)))
                else:
                    (h, pngpath) = mkstemp(suffix='.png', prefix='xic', dir=myTemp)
                    os.close(h)

                    logger_message(20, 'Drawing XIC Plot...')

                    mz_image.make_xic_im(pngpath, mz, xic_times, xic_ints, scan_dot, bin_times, bin_ints, im_size=im_size)

                    logger_message(20, 'Inserting XIC Plot into Spreadsheet...')

                    image_tuples.append(('Peak Width (sec)', 'image', pngpath))

            #precursor mass graph info
            if plot_ms1:
                precursorScanTime = scan_time
                for (i,r) in enumerate(xic[1:]):
                    if r[0] > scan_time:
                        precursorScanTime = xic[i][0]
                        break

                data = mz_file.scan(precursorScanTime)

                scan_mode = data.mode

                xy = [ (x,y) for (x,y) in data if abs(x - floatmz) <= 2.0 ]

                pm_scanDot = (floatmz, max(i[1] for i in xy) if xy else 0.0)

                # don't want to use last entry, and can't use first entry.
                for (i,(x,y)) in enumerate(xy[1:-1]):
                    if x > floatmz:
                        # indexes are shifted because of slice: xy[i] is the previous xy
                        Slope = float(xy[i][1] - y) / float(xy[i][0] - x)
                        pm_scanDot = (floatmz, Slope * (floatmz - x) + y)
                        break

                if isMZD:
                    image_tuples.append(('Experimental mz', 'ms1',
                                         (mz, xy, scan_mode, pm_scanDot)))
                else:
                    (h, pngpath) = mkstemp(suffix='.png', prefix='pm', dir=myTemp)
                    os.close(h)

                    logger_message(20, 'Drawing Precursor Mass Plot...')

                    mz_image.make_ms1_im(pngpath, mz, xy, scan_mode, pm_scanDot, im_size=im_size)

                    logger_message(20, 'Inserting Precursor Mass Plot into Spreadsheet...')

                    image_tuples.append(('Experimental mz', 'image', pngpath))

            if plot_ms2:
                ms_ms_scan = mz_file.scan(scan_time)
                scan_mode = ms_ms_scan.mode

                peptide = mz_pep_format(row['Peptide Sequence'], row['Variable Modifications'] or '')
                charge = row['Charge']
                score = row['Peptide Score']

                if isMZD:
                    image_tuples.append(('Peptide Sequence', 'ms2',
                                         (ms_ms_scan, scan_mode, peptide,
                                          None, ion_list, charge, score)))
                else:
                    (h, pngpath) = mkstemp(suffix='.png', prefix='ms2', dir=myTemp)
                    os.close(h)

                    logger_message(20, 'Drawing MS MS Mass Plot...')

                    mz_image.make_ms2_im(pngpath, ms_ms_scan, scan_mode, peptide,
                                                 None, ion_list, charge, score, im_size=im_size)

                    logger_message(20, 'Inserting MS MS Plot into Spreadsheet...')

                    image_tuples.append(('Peptide Sequence', 'image', pngpath))

            #Calculate Peak Width at half max (FWHM)
            halfIntensity = max_int / 2.0

            reverseScans = xic[::-1]

            beforeHalfTime = start_time
            afterHalfTime = end_time
            peak_comment = 'Good'

            ind_off = len(reverseScans) - max_xic_index - 1
            if max_xic_index > 0:
                for i,r in enumerate(reverseScans[-max_xic_index:]):
                    if r[1] < halfIntensity:
                        beforeSlope = (float(reverseScans[i+ind_off][1] - r[1])
                                       / float(reverseScans[i+ind_off][0] - r[0]))
                        beforeHalfTime = (halfIntensity - r[1]) / beforeSlope + r[0]
                        break
                else:
                    # reached end (beginning) without falling below half intensity
                    beforeHalfTime = r[0]
                    peak_comment = 'Half Intensity Not Reached Before Peak'
            else:
                # if the peak is at the beginning, then that'll be the 'before time'
                beforeHalfTime = xic[0][0]

            ind_off = max_xic_index - 1
            for i,r in enumerate(xic[max_xic_index:]):
                if r[1] < halfIntensity:
                    afterSlope = (float(r[1] - xic[i+ind_off][1])
                                  / float(r[0] - xic[i+ind_off][0]))
                    afterHalfTime = (halfIntensity - r[1]) / afterSlope + r[0]
                    break
            else:
                # reached end without falling below half intensity,
                # including the case where max is at end of the XIC
                afterHalfTime = r[0]
                if peak_comment != 'Half Intensity Not Reached Before Peak':
                    peak_comment = 'Half Intensity Not Reached After Peak'
                else:
                    peak_comment = 'Half Intensity Not Reached at Either End'

            if max_xic_index == 0:
                peak_comment = 'Peak at Beginning of XIC'
            elif max_xic_index == len(xic)-1:
                peak_comment = 'Peak at End of XIC'

            peak_width = (afterHalfTime - beforeHalfTime) * 60.0

            row['MS2 Time'] = scan_time
            row['Peak Time'] = max_time
            row['Peak Intensity'] = max_int
            row['Peak Width (sec)'] = peak_width
            row['Peak Comment'] = peak_comment

            if peak_area:
                logger_message(20,'Calculating Peak Area...')

                row['Peak Area'] = calc_peak_area(xic)

            if reporter_ions:
                if not plot_ms2:
                    ms_ms_scan = mz_file.scan(scan_time)

                halfwindow = 0.02

                row['Rep114'] = max([m[1] for m in ms_ms_scan if abs(m[0] - 114.11) <= halfwindow] or [0])
                row['Rep115'] = max([m[1] for m in ms_ms_scan if abs(m[0] - 115.11) <= halfwindow] or [0])
                row['Rep116'] = max([m[1] for m in ms_ms_scan if abs(m[0] - 116.11) <= halfwindow] or [0])
                row['Rep117'] = max([m[1] for m in ms_ms_scan if abs(m[0] - 117.11) <= halfwindow] or [0])


            row = yield (row, image_tuples)
    finally:
        mz_file.close()


def calc_peak_area(data, multiplier=60.0):
    # Get Gaussian Fit Peak Area
    params = (mean([y for (x,y) in data]),
              median([x for (x,y) in data]),
              0.3, 0)
    (f,p,R2) = fit.fit_data(data=data, parameters=params, function=fit.gauss)

    # Use integral of gaussian
    a = float(p[0])
    c = float(p[2])
    return abs(a * c * multiplier * sqrt(2*pi))


def parse_generic_spec_desc(spec_desc):
    m = re.match("^\((.+),(.+)\).*", spec_desc)
    g = m.groups()
    return (g[0], float(g[1]))


def parse_raw_spec_desc(spec_desc):
    m = re.match('.+\.(\d+)\.(\d+)\.(.+)\.dta',spec_desc)
    if m:
        return (int(m.group(1)),int(m.group(2)))
    else:
        return (0,0)


def parse_wiff_spec_desc(spec_desc):
    m = re.match((r'.*\(sample number (\d+)\), Elution: (.+) min,'
                  ' Period: (\d+), Cycle\(s\): (\d+).*? \(Experiment (\d+)\).*?'),
                 spec_desc, flags=re.I)
    m2 = re.match(r'(?:Locus:)?(\d+)\.(\d+)\.(\d+)\.(\d+)\.(\d+)',
                  spec_desc, flags=re.I) # alternate format

    if m:
        sample = int(m.group(1))
        elution = m.group(2)
        em = re.match('(.+)\sto\s(.+)',elution)
        if em:
            elution = (float(em.group(1)) + float(em.group(2))) / 2
        elution = float(elution)
        period = int(m.group(3))
        cycle = int(m.group(4))
        experiment = int(m.group(5))
    elif m2:
        sample = int(m.group(2))
        elution = 0.0 # elution information not present in this version
        period = int(m.group(3))
        cycle = int(m.group(4))
        experiment = int(m.group(5))
    else:
        sample = 0
        elution = 0.0
        period = 0
        cycle = 0
        experiment = 0

    return (sample,elution,period,cycle,experiment)


def add_centroid_scan_points(scan):
    scan_data = []
    for sc in scan:
        # add data point with zeros before and after
        scan_data.extend(((sc[0],0), sc, (sc[0],0)))

    return scan_data


def add_analyst_scan_points(scan):
    # Add fake zeros if next data point is > 0.03
    scan.sort()
    scan_data = [(scan[0][0]-0.00001,0)] # First Point
    for i,sc in enumerate(scan[:-1]):
        if (scan[i+1][0] - sc[0]) > 0.03:
            # Add data point, zero right after it, zero right before next point
            scan_data.extend((sc,(sc[0]+0.00001,0), (scan[i+1][0]-0.00001,0)))
        else:
            scan_data.append(sc)

    # Add last data point
    scan_data.extend((scan[-1],(scan[-1][0]+0.00001,0)))

    return scan_data


def trim_nearest_points(scan, image_width, dpi=100):
    '''
    This function removes redundant peaks from a scan if they will not be
    displayed in an image of the specified width and dpi.
    '''

    # sort decreasing by intensity
    scan.sort(key=lambda i: i[1], reverse=True)

    # granularity of the ms2 plot: ~ max-min mz / image-width, doubled because it seems we can?
    grain = max(( 2.0 * (max(scan, key=lambda i: i[0])[0] - min(scan, key=lambda i: i[0])[0])
                  / (image_width * dpi) ), 1.0)

    scan_bins = set([])
    new_scan = []
    for (x,y) in scan:
        xbin = int(x/grain)
        if not xbin in scan_bins:
            scan_bins.add(xbin)
            new_scan.append((x,y))

    return new_scan
