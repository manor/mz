import os
import mz
from mz.API import mzdb

tolerance = 0.01

def get_precursor(pMZ, ms):
    return max([(mz,s,n,z) for (mz,s,n,z) in ms if abs(mz-pMZ) <= tolerance] or [(pMZ,0,0.0,0.0)],
               key=lambda x: (x[1],x[0]))

# raw_files is a list of file names.
# mgf_name is the mgf to create (or None to use the RAW name)
# multiple_files is True if each RAW file should be made into its own MGF
def extract_mgf(raw_files, mgf_name=None, multiple_files=False):
    if not multiple_files:
        if not mgf_name:
            mgf_name = raw_files[0][:-5]+".mgf"
        MGF = open(mgf_name, 'w')
        print >> MGF, "MASS=Monoisotopic"
        print >> MGF, "SEARCH=MIS"

    for raw_file in raw_files:
        if multiple_files:
            MGF = open(raw_file[:-5]+".mgf", 'w')
            print >> MGF, "MASS=Monoisotopic"
            print >> MGF, "SEARCH=MIS"

        base_name = os.path.basename(raw_file)[:-5]
        mz_file = mzdb.File(raw_file)
        titles = dict((x+1,y) for (x,y) in enumerate(mz_file.filters()))

        (start,stop) = mz_file.scan_range()
        lastMS1 = None

        for sid in range(start,stop):
            (time,title) = titles[sid]
            if title.find("Full ms ") > -1:
                lastMS1 = mz_file.lscan(sid)
            else:
                a = title.find("Full ms2 ")
                b = title.find("@cid")
                if a > -1 and b > -1:
                    putativeMZ = float(title[(a+9):b])
                    putativeZ = 0
                    (accmz,acccharge) = mz_file.scanPrecursor(sid)
                    if accmz > 0.0 :
                        z  = acccharge
                        mz = accmz
                    else:
                        (mz,s,n,z) = get_precursor(putativeMZ,lastMS1)
                        if z > 0 :
                            (c13mz,c13s,c13n,c13z) = get_precursor(mz - 1.0/float(z) , lastMS1)
                            if c13z == z and c13s > .7 * s :
                                mz = c13mz
                                s = c13s
                                n = c13n
                                z = c13z
                    # Mascot doesn't like charge state 9, skip those scans
                    if z < 9 and len(mz_file.scan(sid)):
                        print >> MGF, "BEGIN IONS"
                        print >> MGF, "TITLE=%s.%d.%d.%d.dta" % (base_name,sid,sid,z)
                        print >> MGF, "PEPMASS=%.5f" % mz
                        if z:
                            print >> MGF, "CHARGE=%d+" % z
                        for (ms2_mz,ms2_s) in mz_file.scan(sid):
                            print >> MGF, "%.5f\t%.2f" % (ms2_mz,ms2_s)
                        print >> MGF, "END IONS\n"

        mz_file.close()
        if multiple_files:
            MGF.close()

    if not multiple_files:
        MGF.close()

    return mgf_name

if __name__ == "__main__":
    mzdb_files = [fname for fname in os.listdir(".") if fname.endswith(".mzdb")]
    extract_mgf(mzdb_files, None, True)
