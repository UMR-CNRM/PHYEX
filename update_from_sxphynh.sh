#!/bin/bash

touch /scratch/work/riette/202005_externalisation_physique/update_from_sxphynh.sh

for file in check_commit_ial.sh comp_DDH.py diffNODE.001_01 Tools prep_code.sh conf_tests update_from_sxphynh.sh env.sh; do
  rsync -rltp --delete --timeout=30 \
      sxphynh.cnrm.meteo.fr:/cnrm/phynh/data1/riette/DATA/202005_externalisation_physique/$file \
      /scratch/work/riette/202005_externalisation_physique/
done
