#!/bin/bash
# closed fluid model: 2D Euler + kinetic Landau damping in each vortex's local frame + momentum deposited in the jets (3 Oct 2026, night)
cd /rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING/scripts/rayleigh
PY=/home/ir-kenn3/rds/rds-ukaea-ap001/ir-kenn3/MINICONDA_PYROKINETICS/INSTALL/envs/pyrokinetics/bin/python3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run() { $PY euler.py "$@" quiet=1 > logs/$1.log 2>&1; }
run kin_hxoff run=hxoff frame=686 hypx=0 hypy=0.0024 kin=1 T=3500 tsamp=4
run kin_vort05 run=hxoff frame=2098 hypx=0 hypy=0.0024 kin=1 eps=0.5 T=2500 tsamp=4
run kin_all05 run=hxoff frame=2098 hypx=0 hypy=0.0024 kin=1 eps=0.5 zeps=0.5 T=2500 tsamp=4
run kin_ref run=leg4 frame=0 hypx=0.0024 hypy=0.0024 kin=1 T=5000 tsamp=4
