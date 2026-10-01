#!/bin/bash
# fluid-model runs, set 1 (1 Oct 2026): long base run, sink scan, resolution scan, long restarts
cd /rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING/scripts/rayleigh
PY=/home/ir-kenn3/rds/rds-ukaea-ap001/ir-kenn3/MINICONDA_PYROKINETICS/INSTALL/envs/pyrokinetics/bin/python3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run() { $PY euler.py "$@" quiet=1 > logs/$1.log 2>&1; }
mkdir -p logs
( run base T=5600 tsamp=2 ; run hx4 hyp=0.0096 T=3000 tsamp=4; run hd4 hyp=0.0006 T=3000 tsamp=4 ) &
( run a03_long eps=0.3 T=5600 tsamp=4 ; run a003_long eps=0.03 T=5600 tsamp=4 ) &
( run nx128p nx=128 T=3000 tsamp=4 ; run nx128g nx=128 dxfac=0.5 T=3000 tsamp=4 ; run vort sink=vort T=3000 tsamp=4 ) &
( run nx256p nx=256 T=3000 tsamp=4 ; run nx256g nx=256 dxfac=0.25 T=3000 tsamp=4 ) &
wait
