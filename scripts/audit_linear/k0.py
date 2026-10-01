import numpy as np
src=open('/rds/project/rds-aSo1XX0UOlw/ir-kenn3/ep_turbulence_paper_letter/tools/verify_endmatter_numerical.py').read()
head=src.split('# ---------------- geometry of kappa')[0]
mr='def min_ratio'+src.split('def min_ratio')[1].split('for name, kap in')[0]
g={}; exec(head+mr,g); min_ratio=g['min_ratio']; Nl=g['Nl']
for a in (0.0,1e-6,1e-3):
    print(a, min_ratio(a*np.ones(Nl))[0])
