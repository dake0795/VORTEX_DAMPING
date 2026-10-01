import numpy as np, os, sys
B='/rds/project/rds-aSo1XX0UOlw/ir-kenn3/DIPOLE_TEST/nl_reducedbox_20260925/'
nx,ny,nz=64,16,32; N=nx*ny*nz
FR=4+8+4+4+N*16+4
def frames(path):
    sz=os.path.getsize(path); n=sz//FR
    return n
def read(path,i):
    with open(path,'rb') as f:
        f.seek(i*FR); b=f.read(FR)
    t=np.frombuffer(b[4:12],'<f8')[0]
    phi=np.frombuffer(b[20:20+N*16],'<c16').reshape((nx,ny,nz),order='F')
    return t,phi
def times(path):
    n=frames(path); ts=[]
    with open(path,'rb') as f:
        for i in range(n):
            f.seek(i*FR+4); ts.append(np.frombuffer(f.read(8),'<f8')[0])
    return np.array(ts)
if __name__=='__main__':
    for p in [B+f'rb_c0p3_hxy/leg_000{l}/out/field.dat' for l in [3,4,5,6]]+[B+f'rb_c0p3_hxy_ampltest/{a}/out/field.dat' for a in ['amp0p3','amp0p03']]:
        ts=times(p); print(p.split('/')[-4:-2], len(ts), ts[0], ts[-1], np.median(np.diff(ts)))
