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
qp_set_frozen_core.py FeO4
echo F > FeO4/perturbation/do_pt2_end
echo 32 > FeO4/davidson/n_states_diag
echo 20 > FeO4/davidson/davidson_sze_max
echo 1.e-3 > FeO4/davidson/threshold_davidson
echo F > FeO4/determinants/read_wf
echo 0.95  > FeO4/determinants/threshold_generators
echo 65408 > FeO4/determinants/n_det_max 
qp_run fci_zmq FeO4 
qp_run save_natorb FeO4
echo 300000 > FeO4/determinants/n_det_max 
echo 0.92  > FeO4/determinants/threshold_generators
qp_run fci_zmq FeO4 
echo T > FeO4/determinants/read_wf
echo 600000 > FeO4/determinants/n_det_max 
echo 0.90 > FeO4/determinants/threshold_selectors
echo 0.90  > FeO4/determinants/threshold_generators

#qp_run fci_zmq_nosave FeO4 > FeO4.out




