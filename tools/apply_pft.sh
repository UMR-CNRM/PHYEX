#!/bin/bash

#This script apply pft.py transformations to arrange the code for GPU testing
PATH_SRC=$(pwd)/../src/common/


#List of files to apply transformations
turb_hor=$(find $PATH_SRC/turb -name '*turb_hor*')
lima=$(find $PATH_SRC/micro -name '*lima*')
rain_ice_old=$(find $PATH_SRC/micro -name '*rain_ice_old*')
other_micro=$(find $PATH_SRC/micro -name '*xker*')
other_turb="$PATH_SRC/turb/mode_ibm_mixinglength.F90 $PATH_SRC/turb/mode_tridiag_w.F90 $PATH_SRC/turb/mode_sbl.F90 $PATH_SRC/turb/mode_coefj.F90 $PATH_SRC/turb/mode_rotate_wind.F90"
micro_rrcol="$PATH_SRC/micro/mode_nrcolss.F90 $PATH_SRC/micro/mode_rrcolss.F90 $PATH_SRC/micro/mode_rscolrg.F90 $PATH_SRC/micro/mode_nscolrg.F90"
black_list="$turb_hor $lima $rain_ice_old $other_micro $other_turb $micro_rrcol $PATH_SRC/micro/mode_ice4_rsrimcg_old.F90"
turb=$(find $PATH_SRC/turb -name \*90)
micro=$(find $PATH_SRC/micro -name \*90)
all="$turb $micro $PATH_SRC/aux/shuman_phy.F90 $PATH_SRC/aux/mode_gradient_m_phy.F90 $PATH_SRC/aux/mode_gradient_u_phy.F90 $PATH_SRC/aux/mode_gradient_v_phy.F90 $PATH_SRC/aux/mode_gradient_w_phy.F90"
for file in $black_list
do
	all=${all//$file/}
done

echo $all
# Link all .h files localy for --addIncludes transformation
includes=$(find $PATH_SRC -name \*.h)
for file in $includes
do
	ln -sf $file .
done
# (Temporary) Use of temp file to save where STACK must be written (--addStack writes then --checkStackArginCall reads it)
rm -f subroutines_wth_stack.txt
touch subroutines_wth_stack.txt

# Apply transformations
for file in $all
do
	echo $file
	pft.py $file $file --addIncludes
	pft.py $file $file --applyCPP
	pft.py $file $file --deleteNonColumnCalls --simplify
	pft.py $file $file --inlineContainedSubroutines
	pft.py $file $file --changeIfStatementsInIfConstructs 
	pft.py $file $file --expandAllArrays
	pft.py $file $file --attachArraySpecToEntity
	pft.py $file $file --removeIJLoops
	pft.py $file $file --reDimKlonArrayToScalar
	#pft.py $file $file --addStack "{kind}, DIMENSION({doubledotshape}), ALLOCATABLE :: {name}#ALLOCATE({name}({shape}))" "MESONH"
	pft.py $file $file --addStack "temp ({kind}, {name}, ({shape}))#alloc ({name})" "AROME"
done


# Add stack (for AROME)
# Recursive checks (maximum found at PHYEX 0.5.0 is 2 mode_prandtl>psi3 then mode_turb_ver_>psi3
for i in {1..3}
do
	for file in $all
	do
		echo $file
		pft.py $file $file --checkStackArginCall
	done

done
