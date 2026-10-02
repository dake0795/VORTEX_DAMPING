cd /rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING/scripts/rayleigh
for zc in 0.5 1.0; do env -u LD_LIBRARY_PATH python3 euler.py kin_zc$zc run=hxoff frame=686 hypx=0 hypy=0.0024 kin=1 zcut=$zc T=3400 tsamp=4 quiet=1 > logs/kin_zc$zc.log 2>&1; done
