cd /rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING/scripts/rayleigh
for C in 4.85 7.0; do env -u LD_LIBRARY_PATH python3 euler.py kin_lst$C run=hxoff frame=686 hypx=0 hypy=0.0024 kin=1 lst=$C T=3400 tsamp=4 quiet=1 > logs/kin_lst$C.log 2>&1; done
