import revision_studies as r
import simulate as s
from scipy.cluster.vq import kmeans2
from scipy.spatial.distance import cdist
import numpy as np,json

a=np.load(s.DATA/'learning.npz');X=a['X'];Y=a['Y'];itr=a['itr'];iva=a['iva'];ite=a['ite']
centers,_=kmeans2(X[itr],256,iter=30,minit='++',seed=20261005)
dist=cdist(X,centers,metric='sqeuclidean');candidates=[]
for gamma in [.3,.6,1.,2.,4.,8.]:
    Z=np.c_[np.ones(len(X)),X,np.exp(-gamma*dist)]
    for ridge in [1e-8,1e-5,.01]:
        beta=np.linalg.solve(Z[itr].T@Z[itr]+ridge*np.eye(Z.shape[1]),Z[itr].T@Y[itr])
        error=np.mean(((Z[iva]@beta-Y[iva])/r.NET['yscale'])**2);candidates.append((error,gamma,ridge,beta))
loss,gamma,ridge,beta=min(candidates,key=lambda v:v[0]);Z=np.c_[np.ones(len(X)),X,np.exp(-gamma*dist)]
relative=np.linalg.norm((Z[ite]@beta-Y[ite])/r.NET['yscale'])/np.linalg.norm(Y[ite]/r.NET['yscale'])
net={'weights':[centers,np.array([gamma]),beta,np.zeros(4),np.eye(4),np.zeros(4)],'yscale':np.ones(4)}
out={'centers':256,'center_selection':'kmeans++ on training only, 30 iterations, seed 20261005','gamma':gamma,'ridge':ridge,'validation_scaled_mse':loss,'test_scaled_relative_L2':float(relative),'coefficient_and_center_bytes':int(centers.nbytes+beta.nbytes+8),'nominal':{}}
np.savez(s.DATA/'rbf.npz',centers=centers,gamma=gamma,beta=beta,ridge=ridge)
for dr,label in [(1.,'exit'),(-1.,'entry')]:
    met,tr=r.run(8,net=net,direction=dr)
    out['nominal'][label]=met;np.savez_compressed(s.DATA/f'RBF_{label}.npz',trace=tr)
    r.write('rbf_report.json',out)
print(json.dumps(out,indent=2),flush=True)
