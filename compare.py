import xarray as xr
import os

#REFDIR a renseigner
#REFDIR="/home/rodierq/"
REFDIR = os.environ['REFDIR']

def compareFiles(file1,file2):
 path_user=REFDIR+'MNH-V5-5-0/MY_RUN/KTEST/007_16janvier/008_run2_'+file1
 if file2 == "ref":
  path_ref=REFDIR+'MNH-V5-5-0/MY_RUN/KTEST/007_16janvier/008_run2'
 else:
  path_ref=REFDIR+'MNH-V5-5-0/MY_RUN/KTEST/007_16janvier/008_run2_'+file2
 
 filen='16JAN.1.12B18.001.nc'
 da = xr.open_dataset(path_user+'/'+filen)
 da2 = xr.open_dataset(path_ref+'/'+filen)
 variables=list(da.keys())
 for var in variables:
  try:
   ecart_min=float(da2[var].min())-float(da[var].min())
   ecart_moy=float(da2[var].mean())-float(da[var].mean())
   ecart_max=float(da2[var].max())-float(da[var].max())
   if (ecart_min !=0 or ecart_moy !=0 or ecart_max !=0):
    print(var,ecart_min,ecart_moy,ecart_max)
  except:
   pass
 
 filen='16JAN.1.12B18.000.nc'
 da = xr.open_dataset(path_user+'/'+filen)
 da2 = xr.open_dataset(path_ref+'/'+filen)
 variables=list(da.keys())
 for var in variables:
  try:
   ecart_min=float(da2[var].min())-float(da[var].min())
   ecart_moy=float(da2[var].mean())-float(da[var].mean())
   ecart_max=float(da2[var].max())-float(da[var].max())
   if (ecart_min !=0 or ecart_moy !=0 or ecart_max !=0):
    print(var,ecart_min,ecart_moy,ecart_max)
  except:
   pass

if __name__ == "__main__":
 import argparse
 parser = argparse.ArgumentParser(description='Compare toutes les variables si trouvées dans les deux fichiers')
 value = argparse.ArgumentParser()
 parser.add_argument('file1', metavar='file1', type=str, help="file1 user ")
 parser.add_argument('file2', metavar='file2', type=str, help="file2 reference; ref for MNH-V5-5-0 reference")
 args = parser.parse_args()
 compareFiles(args.file1,args.file2)

