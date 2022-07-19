#!/usr/bin/env python3

import xarray as xr

avail_groups=['Stations/sta1',
'LES_budgets/Miscellaneous/Cartesian/Not_time_averaged/Not_normalized/cart/',
'LES_budgets/Miscellaneous/Cartesian/Time_averaged/Not_normalized/cart/',
'LES_budgets/Mean/Cartesian/Not_time_averaged/Not_normalized/cart/',
'LES_budgets/Mean/Cartesian/Time_averaged/Not_normalized/cart/',
'LES_budgets/Resolved/Cartesian/Not_time_averaged/Not_normalized/cart/',
'LES_budgets/Resolved/Cartesian/Time_averaged/Not_normalized/cart/',
'LES_budgets/Subgrid/Cartesian/Not_time_averaged/Not_normalized/cart/',
'LES_budgets/Subgrid/Cartesian/Time_averaged/Not_normalized/cart/',
'LES_budgets/Surface/Cartesian/Not_time_averaged/Not_normalized/cart/',
'LES_budgets/Surface/Cartesian/Time_averaged/Not_normalized/cart/',
'LES_budgets/BU_KE/Cartesian/Not_time_averaged/Not_normalized/cart/',
'LES_budgets/BU_KE/Cartesian/Time_averaged/Not_normalized/cart/',
'LES_budgets/BU_THL2/Cartesian/Not_time_averaged/Not_normalized/cart/',
'LES_budgets/BU_THL2/Cartesian/Time_averaged/Not_normalized/cart/',
'LES_budgets/BU_WTHL/Cartesian/Not_time_averaged/Not_normalized/cart/',
'LES_budgets/BU_WTHL/Cartesian/Time_averaged/Not_normalized/cart/',
'LES_budgets/BU_RT2/Cartesian/Not_time_averaged/Not_normalized/cart/',
'LES_budgets/BU_RT2/Cartesian/Time_averaged/Not_normalized/cart/',
'LES_budgets/BU_WRT/Cartesian/Not_time_averaged/Not_normalized/cart/',
'LES_budgets/BU_WRT/Cartesian/Time_averaged/Not_normalized/cart/',
'LES_budgets/BU_THLR/Cartesian/Not_time_averaged/Not_normalized/cart/',
'LES_budgets/BU_THLR/Cartesian/Time_averaged/Not_normalized/cart/']



def compareBACKUPFiles(file_user, file_ref):
  status = 0
  da = xr.open_dataset(file_user)
  da2 = xr.open_dataset(file_ref)
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
      pass
  return status

def compareTSERIESFiles(file_user, file_ref):
  status = 0
  da = xr.open_dataset(file_user)
  da2 = xr.open_dataset(file_ref)
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
  # Groups comparison
  for grp in avail_groups:
    try:
      nk=len(da['level_les'])
      da = xr.open_dataset(file_user, group=grp)
      da2 = xr.open_dataset(file_ref, group=grp)
      variables = list(da.keys())
      print(grp)
      for var in variables:
        try:
          ecart_min = float(da2[var][:,:nk-JPVEXT].min())-float(da[var][:,:nk-JPVEXT].min())
          ecart_moy = float(da2[var][:,:nk-JPVEXT].mean())-float(da[var][:,:nk-JPVEXT].mean())
          ecart_max = float(da2[var][:,:nk-JPVEXT].max())-float(da[var][:,:nk-JPVEXT].max())
          if (ecart_min !=0 or ecart_moy !=0 or ecart_max !=0):
            status += 1
            print(var, ecart_min, ecart_moy, ecart_max)
        except:
          pass
    except:
      pass
  return status

if __name__ == "__main__":
   import argparse
   import sys
   parser = argparse.ArgumentParser(description='Compare toutes les variables si trouvées dans les fichiers backup et time series')
   value = argparse.ArgumentParser()
   parser.add_argument('--f1', metavar='file1', type=str, help="Backup file1 user ")
   parser.add_argument('--f2', metavar='file2', type=str, help="Backup file2 reference")
   parser.add_argument('--f3', metavar='file3', type=str, help=".000 file1 user ")
   parser.add_argument('--f4', metavar='file4', type=str, help=".000 file2 reference")
   args = parser.parse_args()
   status1=compareBACKUPFiles(args.f1, args.f2)
   print('status1 = ' + str(status1))
   if args.f3:
     status2=compareTSERIESFiles(args.f3, args.f4)
     print('status2 = ' + str(status2))

