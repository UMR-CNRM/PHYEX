#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import os
import shutil
import epygram
import numpy
import matplotlib.pyplot as plt
epygram.init_env()

def comp_DDH(filename1, filename2, output_fig, tol_ad=3E-7, tol_rd=1.E-6, verbose=False):
    r1 = epygram.formats.resource(filename1, 'r')
    r2 = epygram.formats.resource(filename2, 'r')

    l1 = set(r1.listfields())
    l2 = set(r2.listfields())
    pb_var = len(l1.symmetric_difference(l2)) != 0

    def comp(fid, v1, v2):
        t = numpy.all(v1 == v2)
        toplt = False
        if not t:
            if verbose: print(fid, ':')
            if numpy.array(v1).ndim == 0:
                v1 = numpy.array([v1])
                v2 = numpy.array([v2])
            for i in range(len(v1)):
                if v1[i] - v2[i] != 0.:
                    ad = v1[i] - v2[i]
                    rd = 200 * (v1[i] - v2[i]) / (v1[i] + v2[i])
                    if verbose: print("   v1={v1}, v2={v2}, diff={ad}, rdiff={rd}".format(v1=v1[i], v2=v2[i], ad=ad, rd=rd))
                    if abs(ad) > tol_ad and abs(rd) > tol_rd:
                      if verbose: print("  ==> plot")
                      toplt = True
        return fid if toplt else None
    toplt = []
    for fid in [fid for fid in l1.intersection(l2) if fid != 'DOCFICHIER']:
        v1 = r1.readfield(fid)
        v2 = r2.readfield(fid)
        if isinstance(v1, epygram.base.FieldSet):
            for i in range(len(v1)): #fieldset
                toplt.append(comp(fid, v1[i].getdata(), v2[i].getdata()))
        else:
            toplt.append(comp(fid, v1.getdata(), v2.getdata()))
    toplt = [fid for fid in toplt if fid is not None]
    pb_val = len(toplt) > 0
    if pb_val:
        figure, ax = plt.subplots(ncols=len(toplt), figsize=(5 * len(toplt), 10), squeeze=False)
        ax = ax[0, :]
        figure.suptitle(filename1 + ' ' + filename2)
        for ifid, fid in enumerate(toplt):
            v1 = r1.readfield(fid)
            v2 = r2.readfield(fid)
            assert(len(v1) == len(v2))
            for i in range(len(v1)): #fieldset
                ad = v1[i].getdata() - v2[i].getdata()
                ax[ifid].plot(v1[i].getdata(), v1[i].geometry.vcoordinate.levels, label='v1')
                ax[ifid].plot(v2[i].getdata(), v2[i].geometry.vcoordinate.levels, label='v2')
                ax[ifid].legend()
                ax[ifid].twiny().plot(ad, v1[i].geometry.vcoordinate.levels, label='diff', color='black', ls=':')
                ad = numpy.abs(ad)
                rd = (200 * numpy.abs(v1[i].getdata() - v2[i].getdata()) / numpy.abs(v1[i].getdata() + v2[i].getdata()))
                rd = rd[ad != 0.].max()
                ad = ad.max()
                ax[ifid].set_title("{fid}:\nmax_ad={ad}\nmax_rd={rd}%".format(fid=fid, ad=ad, rd=rd))
        figure.savefig(output_fig[0])
        for filename in output_fig[1:]:
             shutil.copyfile(output_fig[0], filename)
    if pb_var and pb_val:    
        message = "Variables are different and values of common variables are also different"
    elif pb_var:
        message = "Variables are different but values of common variables are equal"
    elif pb_val:
        message = "Values are different"
    else:
        message = ""
    if pb_val:
        message += ", plot is available in the folowing file(s): " + ', '.join(output_fig)

    print(message)
    return 1 if pb_var or pb_val else 0

if __name__ == '__main__':
    import sys
    sys.exit(comp_DDH(sys.argv[1], sys.argv[2], sys.argv[3:]))
