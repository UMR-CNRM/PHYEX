#!/bin/bash

#This script apply pft.py transformations to arrange the code for GPU testing
#To use this script : 
# link it in PHYEX/src/common
# link all .h files in src/common

#List of files to apply transformations
turb_hor=$(find turb -name '*turb_hor*')
lima=$(find micro -name '*lima*')
rain_ice_old=$(find micro -name '*rain_ice_old*')
other_micro=$(find micro -name '*xker*')
other_turb="turb/mode_ibm_mixinglength.F90 turb/mode_tridiag_w.F90 turb/mode_sbl.F90 turb/mode_coefj.F90 turb/mode_rotate_wind.F90"
other_aux="aux/mode_thermo.F90 aux/shuman.F90 aux/gradient_u.F90 aux/gradient_v.F90 aux/gradient_w.F90 aux/gradient_m.F90 aux/mode_msg.F90"
micro_rrcol="micro/mode_nrcolss.F90 micro/mode_rrcolss.F90 micro/mode_rscolrg.F90 micro/mode_nscolrg.F90"
black_list="$turb_hor $lima $rain_ice_old $other_micro $other_turb $other_aux $micro_rrcol micro/mode_ice4_rsrimcg_old.F90"
turb=$(find turb -name \*90)
micro=$(find micro -name \*90)
all="$turb $micro aux/shuman_phy.F90 aux/mode_gradient_m_phy.F90 aux/mode_gradient_u_phy.F90 aux/mode_gradient_v_phy.F90 aux/mode_gradient_w_phy.F90"
for file in $black_list
do
	all=${all//$file/}
done

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
