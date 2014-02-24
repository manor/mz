import mz.API
import mz.functions
import time
import sys

print sys.argv

fname = ""
target_scan_list = set()
if len(sys.argv) > 1 :
    fname = sys.argv[1]
if len(sys.argv) == 3 :
    f = open(sys.argv[2])
    for l in f:
        target_scan_list.add(int(l.strip()))
    f.close()

if fname == "" :
    print "Please specify a raw-file!"
    sys.exit(-1)

start_time = time.time()
f = mz.API.File(fname)
out = open(fname[:-4]+"_precursors.txt",'w')


(first_scan_time, last_scan_time) = f.time_range()
prev_ms1_time = 0.0
prev_ms1_snum = 0
scans = {}
for (stime,precursor,snum,stype,smode) in f.scan_info(first_scan_time,last_scan_time):
    charge = 0
    half_width = 0.0
    if precursor > 0.0 :
        charge = f.extra_info(snum)["Charge State"]
        half_width = f.extra_info(snum)["MS2 Isolation Width"] / 2.0
    scans[snum-1] = (stime,precursor,charge,half_width)

proton = mz.functions.AW["H"]
def calc_envelope(precursor,charge):
    delta = proton / charge
    masses = set()
    masses.add(precursor - 2.0 * delta)
    masses.add(precursor - delta)
    masses.add(precursor)
    masses.add(precursor + delta)
    masses.add(precursor + 2.0 * delta)
    masses.add(precursor + 3.0 * delta)
    masses.add(precursor + 4.0 * delta)
    return masses

peak_tolerance = 0.01
def calc_s2i_and_p2t(file,snum,precursor,charge,half_width):
    total_signal = 0.0
    precursor_signal = 0.0
    p2t = 1.0
    masses = calc_envelope(precursor,charge)
    for (p_mz,p_intensity,p_threshold,p_charge) in file.lscan(snum):
        if abs(p_mz-precursor) > half_width :
            continue
        else:
            total_signal += p_intensity
            if min(map(lambda x: abs(p_mz - x),masses)) < peak_tolerance:
                precursor_signal += p_intensity
            if abs(precursor-p_mz) < peak_tolerance :
                p2t = max(p2t,p_intensity/p_threshold)
    if total_signal == 0.0:
        return (0.0,p2t)
    else:
        return (precursor_signal/total_signal,p2t)
            
    
def next_ms1(offset,scans):
    max_scan = len(scans)
    next_offset = offset + 1
    while next_offset < max_scan :
        if scans[next_offset][1] > 0.0 :
            next_offset += 1
        else :
            return next_offset
    return -1

def prev_ms1(offset,scans):
    prev_offset = offset - 1
    while prev_offset > -1 :
        if scans[prev_offset][1] > 0.0 :
            prev_offset -= 1
        else :
            return prev_offset
    return -1

#Make XICs range a half minute in either direction...
xic_time_half_window = 1.0
#Make XICs mz-range a "centiDalton" in either direction...
xic_mz_half_window = 0.01
print >> out, "SCAN\tTIME_min\tPRECURSOR\tCHARGE\tS2I\tP2T\tT2APEX_sec\tFWHM_sec\tCOMMENT"
counter = 0
for offset in range(len(scans)):
    snum = offset + 1
    if target_scan_list :
        if not snum in target_scan_list :
            continue
    precursor = scans[offset][1]
    if precursor > 0.0 :
        charge = scans[offset][2]
        half_width = scans[offset][3]
        E_offset = prev_ms1(offset,scans)
        #For now, ignore MS2 without a preceding MS1...
        if E_offset == -1:
            continue
        E_time = scans[E_offset][0]
        L_offset = next_ms1(offset,scans)
        #For now, ignore MS2 without a followup MS1...
        if L_offset == -1:
            continue
        L_time = scans[L_offset][0]
        M_time = scans[offset][0]
        xic_start = max(M_time - xic_time_half_window, first_scan_time + 0.00001)
        xic_end = min(M_time + xic_time_half_window, last_scan_time - 0.00001)
        xic = f.xic(M_time-xic_time_half_window, M_time+xic_time_half_window, precursor-xic_mz_half_window, precursor+xic_mz_half_window)
        (apex_time, apex_int) = max(xic, key=lambda x: x[1])
        T2APEX = (M_time - apex_time) * 60.0

        max_xic_index = xic.index((apex_time, apex_int))
        halfIntensity = apex_int / 2.0
        reverseScans = xic[::-1]
        beforeHalfTime = xic_start
        afterHalfTime = xic_end
        peak_comment = 'FWHM_OK'
        ind_off = len(reverseScans) - max_xic_index - 1
        if max_xic_index > 0:
            for i,r in enumerate(reverseScans[-max_xic_index:]):
                if r[1] < halfIntensity:
                    beforeSlope = (float(reverseScans[i+ind_off][1] - r[1])/float(reverseScans[i+ind_off][0] - r[0]))
                    beforeHalfTime = (halfIntensity - r[1]) / beforeSlope + r[0]
                    break
            else:
                beforeHalfTime = r[0]
                peak_comment = 'Half Intensity Not Reached Before Peak'
        else:
            beforeHalfTime = xic[0][0]

        ind_off = max_xic_index - 1
        for i,r in enumerate(xic[max_xic_index:]):
            if r[1] < halfIntensity:
                afterSlope = (float(r[1] - xic[i+ind_off][1])/ float(r[0] - xic[i+ind_off][0]))
                afterHalfTime = (halfIntensity - r[1]) / afterSlope + r[0]
                break
        else:
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
	if charge == 0 :
	    peak_comment = 'Charge state unknown at acquisition time!'
	    counter += 1
	    print >> out, "%d\t%.3f\t%.2f\t\t\t\t%.2f\t%.2f\t%s" % (snum,M_time,precursor,T2APEX,peak_width,peak_comment)
	else:
            (S2I_E,P2T_E) = calc_s2i_and_p2t(f,E_offset+1,precursor,charge,half_width)
            (S2I_L,P2T_L) = calc_s2i_and_p2t(f,L_offset+1,precursor,charge,half_width)
            S2I = S2I_E + (M_time - E_time)*(S2I_L-S2I_E)/(L_time-E_time)
            P2T = P2T_E + (M_time - E_time)*(P2T_L-P2T_E)/(L_time-E_time)
            counter += 1
            print >> out, "%d\t%.3f\t%.2f\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s" % (snum,M_time,precursor,charge,S2I,P2T,T2APEX,peak_width,peak_comment)
out.close()
f.close()

stop_time = time.time()

print "Characterized %d scans from %s in %.2f seconds." % (counter,fname.split("\\")[-1],stop_time-start_time)
