#!/bin/bash
# fluid model with the Landau damping of the vortices (2 Oct 2026): the drain
cd /rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING/scripts/rayleigh
PY=/home/ir-kenn3/rds/rds-ukaea-ap001/ir-kenn3/MINICONDA_PYROKINETICS/INSTALL/envs/pyrokinetics/bin/python3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run() { $PY euler.py "$@" quiet=1 > logs/$1.log 2>&1; }
mkdir -p logs
run dr_L18 run=hxoff frame=686 hypx=0 hypy=0.0024 landau=0.0018 T=3600 tsamp=4
run dr_L20 run=hxoff frame=686 hypx=0 hypy=0.0024 landau=0.0020 T=3600 tsamp=4
run dr_L0 run=hxoff frame=686 hypx=0 hypy=0.0024 landau=0 T=3600 tsamp=4
run dr_ref run=leg4 frame=0 hypx=0.0024 hypy=0.0024 landau=0.0018 T=5600 tsamp=4
run dr_vort05 run=hxoff frame=2098 hypx=0 hypy=0.0024 landau=0.0018 eps=0.5 T=3000 tsamp=4
run dr_all05 run=hxoff frame=2098 hypx=0 hypy=0.0024 landau=0.0018 eps=0.5 zeps=0.5 T=3000 tsamp=4
wait
