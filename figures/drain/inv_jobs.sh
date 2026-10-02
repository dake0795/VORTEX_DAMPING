cd /rds/project/rds-aSo1XX0UOlw/ir-kenn3/VORTEX_DAMPING/figures/drain
P=/home/ir-kenn3/rds/rds-ukaea-ap001/ir-kenn3/MINICONDA_PYROKINETICS/INSTALL/envs/pyrokinetics/bin/python3
$P invariants.py ref 64 0.0024 0.0024 4 /rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy/leg_0001/out /rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy/leg_0003/out /rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy/leg_0004/out /rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy/leg_0005/out > log_ref.txt 2>&1
$P invariants.py hxoff 64 0 0.0024 4 /rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxoff/leg_0001/out > log_hxoff.txt 2>&1
$P invariants.py nx128 128 0.0024 0.0024 4 /rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy_sinktest/nx128/out > log_nx128.txt 2>&1
$P invariants.py hxy_d4 64 0.0006 0.0006 4 /rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy_sinktest/hxy_d4/out /rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy_sinktest/hxy_d4_leg2/out > log_hxyd4.txt 2>&1
$P invariants.py ctrl 128 0.0384 0.0024 4 /rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy_sinktest/nx128_hx16/out > log_ctrl.txt 2>&1
$P invariants.py hxyx4 64 0.0096 0.0096 4 /rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/rb_c0p3_hxy_sinktest/hxy_x4/out > log_hxyx4.txt 2>&1
wait
