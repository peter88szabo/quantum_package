#!/bin/bash

if [[ -z ${QP_ROOT} ]]
then
  echo "Please source quantum_package.rc file"
  exit 1
fi

cat << EOF > FeO4.xyz
5
Td 
Fe      -0.00000     -0.00000      0.00000
 O       0.90149      0.90149     -0.90149
 O       0.90149     -0.90149      0.90149
 O      -0.90149      0.90149      0.90149
 O      -0.90149     -0.90149     -0.90149
EOF

qp_create_ezfio_from_xyz FeO4.xyz -b def2-svpd -o FeO4
qp_run SCF FeO4
qp_edit -c FeO4
echo 0.95  > FeO4/determinants/threshold_generators
echo 1000000 > FeO4/determinants/n_det_max 
qp_run fci_zmq FeO4 
echo T > FeO4/determinants/read_wf



