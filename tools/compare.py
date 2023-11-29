#!/usr/bin/env python3

import os
os.environ['NUMEXPR_MAX_THREADS'] = '1'
import shutil
import numpy
import sys
import subprocess
import difflib
import re

#List of budgtes groups to compare
avail_groups=['Stations/sta1',
              'LES_budgets/Miscellaneous/Cartesian/Not_time_averaged/Not_normalized/cart/',
              'LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart/',
              'LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart/',
              'LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart/',
              'LES_budgets/Surface/Cartesian/Not_time_averaged/Not_normalized/cart/',
              'LES_budgets/BU_KE/Cartesian/Not_time_averaged/Not_normalized/cart/',
              'LES_budgets/BU_THL2/Cartesian/Not_time_averaged/Not_normalized/cart/',
              'LES_budgets/BU_WTHL/Cartesian/Not_time_averaged/Not_normalized/cart/',
              'LES_budgets/BU_RT2/Cartesian/Not_time_averaged/Not_normalized/cart/',
              'LES_budgets/BU_WRT/Cartesian/Not_time_averaged/Not_normalized/cart/',
              'LES_budgets/BU_THLR/Cartesian/Not_time_averaged/Not_normalized/cart/',
              'Budgets/TH','Budgets/UU','Budgets/WW',
              'Budgets/RV','Budgets/RI','Budgets/RC',
              'Budgets/RG','Budgets/RS','Budgets/RH','Budgets/TK']

def compareBACKUPFiles(file_user, file_ref):
  import xarray as xr
  status = 0
  da = xr.open_dataset(file_user)
  da2 = xr.open_dataset(file_ref)
  JPHEXT=1
  JPVEXT=1
  ni=len(da['ni'])
  nj=len(da['nj'])
  nk=len(da['level'])
  variables = list(da.keys())
  for var in [var for var in variables if da[var].dtype.char != 'S']:
    try:
      if da[var].ndim == 4: #Variables time, level, nj, ni
        ecart_min=float(da2[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].min())-float(da[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].min())
        ecart_moy=float(da2[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].mean())-float(da[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].mean())
        ecart_max=float(da2[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].max())-float(da[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].max())
      elif  da[var].ndim == 3 and da['L2D'] == 0: #Variables time, nj, ni
        ecart_min=float(da2[var][0,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].min())-float(da[var][0,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].min())
        ecart_moy=float(da2[var][0,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].mean())-float(da[var][0,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].mean())
        ecart_max=float(da2[var][0,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].max())-float(da[var][0,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].max())
      elif  da[var].ndim == 3 and da['L2D'] == 1: #Variables time, level, nj or ni (2D simulation)
        if len(da['ni']) > len(da['nj']):
          nij=len(da['ni'])
        else:
          nij=len(da['nj'])
        ecart_min=float(da2[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nij-1-JPHEXT].min())-float(da[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nij-1-JPHEXT].min())
        ecart_moy=float(da2[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nij-1-JPHEXT].mean())-float(da[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nij-1-JPHEXT].mean())
        ecart_max=float(da2[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nij-1-JPHEXT].max())-float(da[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nij-1-JPHEXT].max())
      else:
        ecart_min=float(da2[var].min())-float(da[var].min())
        ecart_moy=float(da2[var].mean())-float(da[var].mean())
        ecart_max=float(da2[var].max())-float(da[var].max())
      if (ecart_min !=0 or ecart_moy !=0 or ecart_max !=0):
        status += 1
        print(var, ecart_min, ecart_moy, ecart_max)
    except:
      #raise
      pass
  return status

def compareTSERIESFiles(file_user, file_ref, tol_ad=1E-12):
  import xarray as xr
  status = 0
  da = xr.open_dataset(file_user)
  da2 = xr.open_dataset(file_ref)
  variables = list(da.keys())
  JPVEXT=1
  try: 
    nk=len(da['level_les'])
  except:
    pass
  for var in variables:
    try:
      ecart_min = float(da2[var].min())-float(da[var].min())
      ecart_moy = float(da2[var].mean())-float(da[var].mean())
      ecart_max = float(da2[var].max())-float(da[var].max())
      if (ecart_min !=0 or ecart_moy !=0 or ecart_max !=0):
        status += 1
        print(var, ecart_min, ecart_moy, ecart_max)
    except:
      pass
  # Groups comparison
  for grp in avail_groups:
    try: # LES or Stations variables in 1D/2D
      da = xr.open_dataset(file_user, group=grp)
      da2 = xr.open_dataset(file_ref, group=grp)
      variables = list(da.keys())
      for var in variables:
        try: #LES variables in 2D
          ecart_min = float(da2[var][:,:nk-JPVEXT].min())-float(da[var][:,:nk-JPVEXT].min())
          ecart_moy = float(da2[var][:,:nk-JPVEXT].mean())-float(da[var][:,:nk-JPVEXT].mean())
          ecart_max = float(da2[var][:,:nk-JPVEXT].max())-float(da[var][:,:nk-JPVEXT].max())
          if (ecart_min !=0 or ecart_moy !=0 or ecart_max !=0):
            status += 1
            print(var, ecart_min, ecart_moy, ecart_max)
        except:
          try: # Sations or Budgets variables (Budget box without HALO points)
            ecart_min = float(da2[var][:].min())-float(da[var][:].min())
            ecart_moy = float(da2[var][:].mean())-float(da[var][:].mean())
            ecart_max = float(da2[var][:].max())-float(da[var][:].max())
            if (abs(ecart_min) >=tol_ad or abs(ecart_moy) >=tol_ad or abs(ecart_max) >=tol_ad):
              status += 1
              print(grp, var, ecart_min, ecart_moy, ecart_max)
          except:
            pass

    except:
      pass
  return status

def comp_DDH(filename1, filename2, output_fig, tol_ad=3E-7, tol_rd=1.E-6, verbose=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import epygram
    epygram.init_env()

    if filename1 != filename2:
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
        if pb_val and output_fig is not None and len(output_fig) > 0:
            ncols, nrows = min(10, len(toplt)), 1+(len(toplt)-1)//10
            figure, ax = plt.subplots(ncols=ncols, nrows=nrows,
                                      figsize=(5 * ncols, 10 * nrows), squeeze=False)
            ax = ax.flatten()
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
        if pb_val and output_fig is not None and len(output_fig) > 0:
            message += ", plot is available in the folowing file(s): " + ', '.join(output_fig)
    else: #same file
        pb_var = False
        pb_val = False
        message = ""

    print(message)
    return 1 if pb_var or pb_val else 0

default_spnorms='VORTICITY,DIVERGENCE,TEMPERATURE,KINETIC ENERGY'
default_gpnorms='HUMI.SPECIFI,CLOUD_WATER,ICE_CRYSTAL,CLOUD_FRACTI'
def comp_NODE(f1, f2, spnorms=None, gpnorms=None, norm_max_diff=0.050000):
    """
    Adapted from the Ryad's diffNODE.001_01
    :param f1, f2: filenames to compare
    :param spnorms: coma separated list of spectral variables. '' to not
                    compare spectral norms. None to use the default
                    (""" + default_spnorms + """)
    :param gpnorms: coma separated list of gridpoint variables. '' to not
                    compare gridpoint norms. None to use the default
                    (""" + default_gpnorms + """)
    :param norm_max_diff: maximum difference allowed
    """
    
    from collections import defaultdict

    if spnorms is None: spnorms = default_spnorms
    if gpnorms is None: gpnorms = default_gpnorms    
    opts = {'spnorms': spnorms.split(','),
            'gpnorms': gpnorms.split(','),
            'norm_max_diff': norm_max_diff}

    def xave(f):
        with open(f, 'r') as fh:
            gpregs = []
            if opts['gpnorms']:
                if len(opts['gpnorms']) == 1 and opts['gpnorms'][0] == '*':
                    gpregs = [re.compile(r'^\s*GPNORM\s+\b(\S+)\b')]
                else:
                    gpregs = [re.compile(rf'^\s*GPNORM\s+\b({re.escape(gp)})\b') for gp in opts['gpnorms']]
            x = []
            for line in fh:
                line = line.strip()
                if re.match(r'^\s*GPNORM\s+', line):
                    for gpreg in gpregs:
                        match = gpreg.match(line)
                        if match:
                            F = match.group(1)
                            line = next(fh)
                            if not re.match(r'\s*AVE\s+', line):
                                continue
                            line = line.strip().lstrip('AVE')
                            values = line.split()
                            x.extend([[F, value] for value in values])
                            break
                elif re.match(r'^\s*SPECTRAL\s+NORMS\s+-\s+', line):
                    spnormk = []
                    spnormv = []
                    while True:
                        line = next(fh)
                        if not re.match(r'\s+LEV\s+', line):
                            break
                        index = {spnorm: line.find(spnorm) for spnorm in opts['spnorms']}
                        index = {k:v for (k, v) in index.items() if v > 0}
                        spnormk = sorted([sp for sp in opts['spnorms'] if sp in index.keys()],
                                         key=lambda spnorm: index[spnorm] if index[spnorm] >= 0 else float('inf'))
                        line = next(fh)
                        if not re.match(r'\s+AVE\s+', line):
                            break
                        spnormv = line.split()[1:] #'1:' to remove 'AVE'
                        for spnormk, spnormv in zip(spnormk, spnormv):
                            x.append([spnormk, spnormv])
            return x

    def center(s, n):
        i = 0
        while len(s) < n:
            s = f" {s}" if i % 2 else f"{s} "
            i += 1
        return s

    def main(f1, f2):
    
        fx1 = xave(f1)
        fx2 = xave(f2)

        title = "************* NORMS DIFFERENCES *************"
        print(center(title, 118))
        print(center('=' * len(title), 118))
        print()
    
        x = [[]]
        diff = {}
        zero = 0
        numb = 0
    
        tag1 = "NORMDIFF"
        tag2 = "NORMSTAT"
    
        nout = 0
        while fx1 and fx2:
            (field1, x1) = fx1.pop(0)
            (field2, x2) = fx2.pop(0)
    
            if field1 != field2:
                return 1
                #sys.exit("Field mismatch {} != {}".format(field1, field2))
    
            x1, x2 = x1.strip(), x2.strip()
            if x1 and x2:
                dx = float(x1) - float(x2)
                dr = (2 * dx) / (float(x1) + float(x2)) if float(x1) + float(x2) > 0 else 0.0
    
                sdx = '{:17.9e}'.format(dx)
                sdr = '{:17.9e}'.format(dr)
    
                dx = float(sdx)
                dr = float(sdr)
    
                x[-1].append(" {} | {:20} | {:17.9e}  |  {:17.9e}  |  {:17}  |  {:17} {}\n".format(
                    tag1, center(field1, 20), float(x1), float(x2), sdx, sdr, '*' if dr > opts['norm_max_diff'] else ''))
    
                nout += 1 if abs(dr) > opts['norm_max_diff'] else 0
    
                if abs(dr) > 0:
                    n = int(numpy.floor((numpy.log(abs(dr)) / numpy.log(10))))
                    diff[n] = diff.get(n, 0) + 1
                else:
                    zero += 1
    
                numb += 1
            else:
                x.append([])
    
        print(" {} |                      |{:19} | {:19} | {:19} | {:19}".format(
            tag1, center("NORM(REF)", 19), center("NORM(EXP)", 19),
            center("NORM(REF)-NORM(EXP)", 19), center("(NORM(REF)-NORM(EXP))", 19)))
    
        print(" {} |                      |{:19} | {:19} | {:19} | {:19}".format(
            tag1, '', '', '', center("/NORM(REF)", 19)))
    
        for i in range(len(x)):
            if not x[i]:
                break
            print("".join(x[i]))
    
        print("\n")
    
        diff_cumul = 0
        perc_cumul = 0
        for n1 in sorted(diff.keys()):
            n2 = n1 + 1
            diff_val = diff[n1]
            perc = 100 * diff_val / numb
            diff_cumul += diff_val
            perc_cumul += perc
            print(" {} | {:3d} .. {:3} | {:3d} / {:3d} | {:3d} / {:3d} | {:6.2f} %, {:6.2f} %\n".format(
                  tag2, n1, n2, diff_val, numb, diff_cumul, numb, perc, perc_cumul))
    
        if nout:
            print("\n")
            text = "WARNING : SOME NORMS DIFFERENCES ARE OUTSIDE ALLOWED LIMIT OF {:6.2f} %\n".format(
                100 * opts['norm_max_diff'])
            print(text * 5)
        return nout
    
    return main(f1, f2)

def comp_binary(f1, f2, offset):
    #pyhton filcmp does not allow to specify an offset
    offset = [str(offset), str(offset)] if offset != 0 else []
    p = subprocess.run(['cmp', f1, f2] + offset, capture_output=True, encoding='UTF8') 
    if p.returncode != 0:
        print(p.stdout)
    return p.returncode

def comp_ncdump(f1, f2, nbytes):
    ncdumps = [subprocess.run(['ncdump', f], capture_output=True, encoding='UTF8').stdout[:int(nbytes)]
               for f in (f1, f2)]
    diff = ''.join(difflib.unified_diff(*[ncdumps[i].splitlines(keepends=True) for i in (0, 1)]))
    if diff != '':
        print(diff)
    return 0 if ncdumps[0] == ncdumps[1] else 1

def comp_testprogs(f1, f2):
    def read(f):
        with open(f, 'r') as fd:
            s = fd.read()
            for p in ('\.\.', '~=', '!='): s = re.sub(p, '', s)
            s = re.sub(r'\-0.00000E\+00([|\- ])', r' 0.00000E+00\1', s)  
            s = re.sub(r'\n\sTotal time:.*\n', '\n', s)
            s = re.sub(r'IBL =[ ]*', 'IBL = ', s)
        return s[s.index('IBL'):]
    r = read(f1) == read(f2)
    if not r:
        print("{f1} {f2} differ".format(f1=f1, f2=f2))
    return 0 if r else 1

if __name__ == "__main__":
   import argparse
   import sys
   parser = argparse.ArgumentParser(description='Compare all the variables in the different files')
   parser.add_argument('--backup', metavar=('BACKUP_USER', 'BACKUP_REF'), nargs=2,
                       type=str, help="Backup files (user and reference)")
   parser.add_argument('--diac', metavar=('DIAC_USER', 'DIAC_REF'), nargs=2,
                       type=str, help="Diachronic .000 files (user and reference)")
   parser.add_argument('--ddh', nargs=2, metavar=('DDH_USER', 'DDH_REF'),
                       type=str, help="DDH files (user and reference)")
   parser.add_argument('--ddhplots', metavar='DDHPLOT', nargs='+',
                       help="Plot filenames for DDH differences")
   parser.add_argument('--node', nargs=2, metavar=('DDH_USER', 'DDH_REF'),
                       type=str, help="NODE files (user and reference) to compare norms")
   parser.add_argument('--binary', nargs=3, metavar=('BIN_USER', 'BIN_REF', 'OFFSET'),
                       type=str, help="Binary files (user and reference) and offset")
   parser.add_argument('--ncdump', nargs=3, metavar=('NC_USER', 'NC_REF', 'BYTES'),
                       type=str, help="Netcdf files (user and reference) whose ncdump output" + \
                                      "first bytes must be compared, with number of bytes")
   parser.add_argument('--testprogs', nargs=2, metavar=('LST_USER', 'LST_REF'),
                       type=str, help="Testprogs listing (user and ref)")
   args = parser.parse_args()

   totalstatus = 0
   if args.backup:
       status = compareBACKUPFiles(*args.backup)
       totalstatus += status
       print('status for backup files = ' + str(status))
   if args.diac:
       status = compareTSERIESFiles(*args.diac)
       totalstatus += status
       print('status for diachronic files = ' + str(status))
   if args.ddh:
       status = comp_DDH(*args.ddh, args.ddhplots)
       totalstatus += status
       #print('status for ddh files = ' + str(status))
   if args.node:
       status = comp_NODE(*args.node, norm_max_diff=0.)
       totalstatus += status
       #print('status for NODE files = ' + str(status))
   if args.binary:
       status = comp_binary(*args.binary)
       totalstatus += status
       #print('status for binary files = ' + str(status))
   if args.ncdump:
       status = comp_ncdump(*args.ncdump)
       totalstatus += status
       print('status for ncdump of files = ' + str(status))
   if args.testprogs:
       status = comp_testprogs(*args.testprogs)
       totalstatus += status
       #print('status for testprogs listings = ' + str(status))
   sys.exit(totalstatus)
