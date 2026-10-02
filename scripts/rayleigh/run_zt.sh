cd /rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING/scripts/rayleigh
for zt in 0.55 0.45 0.65; do env -u LD_LIBRARY_PATH python3 euler.py kin_zt$zt run=hxoff frame=686 hypx=0 hypy=0.0024 kin=1 ztot=$zt T=3400 tsamp=4 quiet=1 > logs/kin_zt$zt.log 2>&1; done
