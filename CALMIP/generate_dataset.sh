#!/bin/bash
#SBATCH -N 1 -n 1 -c 64 --qos=mesca --mem=1000000 --time=10:00:00

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
echo 0.91  > FeO4/determinants/threshold_generators
echo 1000000 > FeO4/determinants/n_det_max 
qp_run mp2_wf FeO4 
qp_run save_natorb FeO4 
qp_set_mo_class -core "[1-9]" -act "[10-50]" -del "[51-126]" FeO4
echo F > FeO4/perturbation/do_pt2_end
qp_run fci_zmq FeO4 
qp_set_frozen_core.py FeO4
echo T > FeO4/determinants/read_wf
echo 1.e-3 > FeO4/davidson/threshold_davidson
echo 20 > FeO4/davidson/davidson_sze_max
echo 20 > FeO4/davidson/n_states_diag
echo 2000000 > FeO4/determinants/n_det_max 
echo 0.91 > FeO4/determinants/threshold_selectors

QP_PREFIX=time qp_run fci_zmq_nosave FeO4 > FeO4.out




