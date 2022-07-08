# INSTALLATION NEEDED THE prep\_code.sh script

## ABOUT THIS DOCUMENT 

This document is intended for persons who want to use the prep\_code.sh script (directly or through a check\_commit\_\*.sh script).

This document is written using the markdown language. With pandoc, it can be converted to HTML (pandoc -s \<filename\>.md -o \<filename\>.html) or PDF (pandoc -s \<filename\>.md -o \<filename\>.pdf).

Two packages must be installed:

  - [filepp](#filepp)
  - [MNH\_Expand\_Array](#mnh\_expand\_array)

## filepp
In the \<git repository\>/tools/mnh\_expand directory:

```
wget https://www-users.york.ac.uk/~dm26/filepp/filepp-1.8.0.tar.gz
tar xvf filepp-1.8.0.tar.gz
cd filepp-1.8.0
./configure --prefix=$PWD
make install
cd ..
ln -s filepp-1.8.0 filepp
```

## MNH\_Expand\_Array
In the \<git repository\>/tools/mnh\_expand directory, clone the MNH\_Expand\_Array repository:

```
git clone https://github.com/JuanEscobarMunoz/MNH_Expand_Array.git
```


