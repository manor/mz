import sqlite3
import struct
import sys
import os
import time
import mz
from mz import API

def _to_float(x):
    try :
        out = float(x)
    except ValueError :
        out = str(x)
    return out

def extract_file_into_db(fname,dbname=None):
    if not dbname:
        dbname = fname[:-4]+".mzdb"
    if os.path.exists(dbname):
        os.remove(dbname)
    conn = sqlite3.connect(dbname)
    c = conn.cursor()
    c.execute('''create table scans(
    ID INTEGER PRIMARY KEY,
    time REAL,
    precursor REAL,
    filter text,
    header text,
    peak_num INTEGER,
    tic REAL,
    bpc REAL,
    mzs BLOB,
    ys BLOB,
    zs BLOB)''')
    c.execute("""create index timeindex ON scans (time)""") ;
    c.execute("""create index precursorindex ON scans (precursor)""") ;
    f = API.File(fname)
    filters = f.filters()
    scan_list = f.scan_list(0,100000,0,100000)
    for (offset,(scan_time,scan_precursor)) in enumerate(scan_list):
        scan_num = offset+1
        (scan_time_again,scan_filter) = filters[offset]
        detector = scan_filter[:4]
        #scan_header = cPickle.dumps(f.extra_info(scan_num), cPickle.HIGHEST_PROTOCOL)
        extra_info = f.extra_info(scan_num)
        scan_header = "\t".join(["\t".join([str(a),str(_to_float(v))]) for (a,v) in extra_info.items()])
        if not detector in ["ITMS","FTMS"]:
            print detector
            import sys
            sys.exit()
        if detector == "FTMS" :
            scan_data  = [(d1,d2,d4) for (d1,d2,d3,d4) in f.lscan(scan_num)]
            peak_num = len(scan_data)
            scan_tic = sum([d2 for (d1,d2,d4) in scan_data])
            scan_bpc = max([d2 for (d1,d2,d4) in scan_data])
            #There is apparently no loss of accuracy encoding as float rather than double
            #I don't know whether this is because raw.py drops the ball and loses info or
            #because there never was more than a float in there to begin with...
            scan_mzs = struct.pack("!%df"%peak_num,*[d1 for (d1,d2,d3) in scan_data])
            scan_ys = struct.pack("!%df"%peak_num,*[d2 for (d1,d2,d3) in scan_data])
            scan_zs = struct.pack("!%dB"%peak_num,*[int(d3) for (d1,d2,d3) in scan_data])
            #scan = cPickle.dumps([(d1,d2,d4) for (d1,d2,d3,d4) in f.lscan(scan_num)], cPickle.HIGHEST_PROTOCOL)
        else:
            #***NB*** scan_lists already gives the correct precursor...
            #alt_precursor = extra_info['Monoisotopic M/Z']
            scan_data  = [(d1,d2) for (d1,d2) in f.scan(scan_num)]
            peak_num = len(scan_data)
            scan_tic = sum([d2 for (d1,d2) in scan_data] + [0.0])
            scan_bpc = max([d2 for (d1,d2) in scan_data] + [0.0])
            scan_mzs = struct.pack("!%df"%peak_num,*[d1 for (d1,d2) in scan_data])
            scan_ys = struct.pack("!%df"%peak_num,*[d2 for (d1,d2) in scan_data])
            scan_zs = ""
            #scan = cPickle.dumps(f.scan(scan_num), cPickle.HIGHEST_PROTOCOL)
        #***NB*** here is the test that proves it...
        #if detector == "ITMS" and scan_precursor > 0.0:
        #    if abs(scan_precursor-alt_precursor) > 0.001:
        #        print scan_precursor
        #        print alt_precursor
        c.execute("insert into scans values (?,?,?,?,?,?,?,?,?,?,?)", (scan_num,scan_time,scan_precursor,scan_filter,scan_header,peak_num,scan_tic,scan_bpc,buffer(scan_mzs),buffer(scan_ys),buffer(scan_zs)) )

    c.close()
    conn.commit()
    conn.close()
 
    f.close()
    stop_time = time.time()

    print "Extracted %s into .mzdb in %.2f seconds." % (fname,stop_time - start_time)


start_time = time.time()

if len(sys.argv) < 2 :
    import os
    fnames = [f for f in os.listdir(".") if f.endswith(".RAW") or f.endswith(".raw")]
    for fname in fnames:
    os_start_time = time.time()
        os.system("python liberate.py %s" % (fname))
    os_stop_time = time.time()
    print "%s extracted in %.0f seconds" % (fname,os_stop_time - os_start_time)
else:
    fname = sys.argv[1]
    extract_file_into_db(fname)

