#!/usr/bin/env python3

import xarray as xr

def compareFiles(path_user, path_ref):
  status = 0
  
  filen = '16JAN.1.12B18.001.nc'
  da = xr.open_dataset(path_user + '/' + filen)
  da2 = xr.open_dataset(path_ref + '/' + filen)
  JPHEXT=1
  JPVEXT=1
  ni=len(da['ni'])
  nj=len(da['nj'])
  nk=len(da['level'])
  variables = list(da.keys())
  for var in variables:
    try:
      if da[var].ndim == 4: #Variables time, level, nj, ni
        ecart_min=float(da2[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].min())-float(da[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].min())
        ecart_moy=float(da2[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].mean())-float(da[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].mean())
        ecart_max=float(da2[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].max())-float(da[var][0,JPVEXT:nk-1-JPVEXT,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].max())
      elif  da[var].ndim == 3: #Variables time, nj, ni
        ecart_min=float(da2[var][0,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].min())-float(da[var][0,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].min())
        ecart_moy=float(da2[var][0,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].mean())-float(da[var][0,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].mean())
        ecart_max=float(da2[var][0,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].max())-float(da[var][0,JPHEXT:nj-1-JPHEXT,JPHEXT:ni-1-JPHEXT].max())
      else:
        ecart_min=float(da2[var].min())-float(da[var].min())
        ecart_moy=float(da2[var].mean())-float(da[var].mean())
        ecart_max=float(da2[var].max())-float(da[var].max())
      if (ecart_min !=0 or ecart_moy !=0 or ecart_max !=0):
        print(var,ecart_min,ecart_moy,ecart_max)
    except:
      pass
  
  filen = '16JAN.1.12B18.000.nc'
  da = xr.open_dataset(path_user + '/' + filen)
  da2 = xr.open_dataset(path_ref + '/' + filen)
  variables = list(da.keys())
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

  return status

if __name__ == "__main__":
   import argparse
   import sys
   parser = argparse.ArgumentParser(description='Compare toutes les variables si trouvées dans les deux fichiers')
   value = argparse.ArgumentParser()
   parser.add_argument('file1', metavar='file1', type=str, help="file1 user ")
   parser.add_argument('file2', metavar='file2', type=str, help="file2 reference")
   args = parser.parse_args()
   sys.exit(compareFiles(args.file1, args.file2))

