"""Horizon, architecture, sensing, delay, actuator and structured stress studies."""
import revision_studies as r
import simulate as s
import numpy as np,json,time,itertools
from scipy.linalg import expm

def train_variant(hidden=(48,32),N=6000,epochs=950,activation=0,scaled=True,seed=20261003):
    a=np.load(s.DATA/'learning.npz');X=a['X'];Y=a['Y'];itr=a['itr'][:int(.7*N)];iva=a['iva'];ite=a['ite']
    rng=np.random.default_rng(seed);dims=[6,*hidden,4]
    ys=np.maximum(Y[itr].std(0),s.LIMITS*.002) if scaled else np.ones(4)
    target=Y/ys;weights=[]
    for n,m in zip(dims[:-1],dims[1:]):weights.extend([rng.normal(0,np.sqrt(2/(n+m)),(n,m)),np.zeros(m)])
    mom=[np.zeros_like(w) for w in weights];var=[np.zeros_like(w) for w in weights];best=np.inf;bw=None;step=0
    def forward(xx):
        z1=xx@weights[0]+weights[1];h1=np.tanh(z1) if activation==0 else np.maximum(0,z1)
        z2=h1@weights[2]+weights[3];h2=np.tanh(z2) if activation==0 else np.maximum(0,z2)
        return h1,h2,h2@weights[4]+weights[5]
    for ep in range(epochs):
        for sub in np.array_split(rng.permutation(itr),max(1,len(itr)//256)):
            h1,h2,pred=forward(X[sub]);dy=(pred-target[sub])/(len(sub)*4)
            grad4=h2.T@dy;grad5=dy.sum(0);dh2=(dy@weights[4].T)*((1-h2*h2) if activation==0 else (h2>0))
            grad2=h1.T@dh2;grad3=dh2.sum(0);dh1=(dh2@weights[2].T)*((1-h1*h1) if activation==0 else (h1>0))
            grad=[X[sub].T@dh1,dh1.sum(0),grad2,grad3,grad4,grad5];step+=1;lr=.0025*(.35+.65*(1-ep/epochs))
            for j in range(6):
                mom[j]=.9*mom[j]+.1*grad[j];var[j]=.999*var[j]+.001*grad[j]**2
                weights[j]-=lr*mom[j]/(1-.9**step)/(np.sqrt(var[j]/(1-.999**step))+1e-8)
        if ep%10==0 or ep==epochs-1:
            val=np.mean((forward(X[iva])[2]-target[iva])**2)
            if val<best:best=val;bw=[x.copy() for x in weights]
    weights=bw;pred=forward(X[ite])[2]*ys
    met={'hidden':list(hidden),'N_total_nested_design':N,'N_training':len(itr),'validation_N':len(iva),'test_N':len(ite),'epochs':epochs,'activation':'tanh' if activation==0 else 'ReLU','output_scaling':scaled,'seed':seed,'test_scaled_relative_L2':float(np.linalg.norm((pred-Y[ite])/r.NET['yscale'])/np.linalg.norm(Y[ite]/r.NET['yscale'])),'validation_mse_in_training_scaling':float(best)}
    return {'weights':weights,'yscale':ys,'activation':activation},met

def polynomial_fit():
    a=np.load(s.DATA/'learning.npz');X=a['X'];Y=a['Y'];itr=a['itr'];iva=a['iva'];ite=a['ite']
    powers=np.array([p for p in itertools.product(range(4),repeat=6) if sum(p)<=3],dtype=np.int64)
    Z=np.prod(X[:,None,:]**powers[None,:,:],axis=2)
    candidates=[]
    for ridge in [1e-8,1e-6,1e-4,.01,1.]:
        coef=np.linalg.solve(Z[itr].T@Z[itr]+ridge*np.eye(len(powers)),Z[itr].T@Y[itr])
        loss=np.mean(((Z[iva]@coef-Y[iva])/r.NET['yscale'])**2);candidates.append((loss,ridge,coef))
    loss,ridge,coef=min(candidates,key=lambda v:v[0]);pred=Z[ite]@coef
    np.savez(s.DATA/'polynomial.npz',coef=coef,powers=powers)
    return (coef,powers),{'degree':3,'basis_terms':len(powers),'ridge':ridge,'selected_validation_loss':loss,'test_scaled_relative_L2':float(np.linalg.norm((pred-Y[ite])/r.NET['yscale'])/np.linalg.norm(Y[ite]/r.NET['yscale']))}

def main():
    start=time.perf_counter();out={'horizons':{},'nominal':{},'ablations':{},'networks':[],'sensing':{},'delays':{},'drivers':{},'structured':{}}
    poly,polyreport=polynomial_fit();out['polynomial_learning']=polyreport
    prepared=r.prepare_plant(s.PARAM)
    def calc(label,kind=4,category='nominal',**kw):
        met,tr=r.run(kind,prepared=prepared,**kw);out[category][label]=met
        np.savez_compressed(s.DATA/f'{category}_{label}.npz',trace=tr)
        r.write('additional_report.json',out);print(category,label,met['cost'],met['failure'],flush=True)
    for dr in [1.,-1.]:
        suffix='exit' if dr>0 else 'entry'
        for H in [0.,6.,12.,24.,36.,48.]:calc(f'H{int(H)}_{suffix}',2 if H==0 else 3,category='horizons',H=H,direction=dr,filter=False)
        for label,kind,filt in [('LQR',2,False),('Oracle24',3,False),('RawNN',4,False),('FilteredNN',4,True),('Polynomial',6,True),('HubPD',0,False),('Passivity',7,False),('Frozen24',3,False)]:
            calc(label+'_'+suffix,kind,direction=dr,filter=filt,frozen=label=='Frozen24',polynomial_model=poly)
    for label,kw in [('NoGate',{'gate':False,'filter':False}),('Gate',{'gate':True,'filter':False}),('GateFilter',{'gate':True,'filter':True}),('NoResidual',{'rmd':False,'filter':True}),('eps1e-10',{'eps':1e-10}),('eps1e-6',{'eps':1e-6}),('StressRaw',{'bias':80.,'filter':False,'end':90.}),('StressFiltered',{'bias':80.,'filter':True,'end':90.})]:calc(label,category='ablations',**kw)
    variants=[((16,8),6000,950,0,True),((32,16),6000,950,0,True),((96,64),6000,950,0,True),((48,32),1500,950,0,True),((48,32),3000,950,0,True),((48,32),6000,200,0,True),((48,32),6000,500,0,True),((48,32),6000,950,1,True),((48,32),6000,950,0,False)]
    for i,args in enumerate(variants):
        net,met=train_variant(*args);met['closed_loop'],tr=r.run(4,net=net,prepared=prepared)
        out['networks'].append(met);np.savez(s.DATA/f'network_variant_{i}.npz',**{f'w{j}':w for j,w in enumerate(net['weights'])},yscale=net['yscale'],activation=net['activation'])
        r.write('additional_report.json',out);print('network',i,met['test_scaled_relative_L2'],flush=True)
    for level in ['low','high']:
        angle=.02 if level=='low' else .1;strain=5e-9 if level=='low' else 25e-9
        noise=[.1 if level=='low' else .5,-.1 if level=='low' else -.5,.05 if level=='low' else .2,.05 if level=='low' else .2,angle/206264.806,-angle/206264.806,strain,-strain,strain,angle/206264.806,angle/206264.806,strain,strain,strain]
        for label,kind,filt in [('LQR',2,False),('FilteredNN',4,True)]:calc(label+'_'+level,kind,category='sensing',noise=noise,observer=True,filter=filt)
    for delay in [0,3,6,13]:
        for label,kind,filt in [('LQR',2,False),('FilteredNN',4,True)]:calc(label+f'_core{delay}',kind,category='delays',core_delay=delay,filter=filt)
    for delay in [1,2,3]:calc(f'residual{delay}',category='delays',residual_delay=delay)
    for bandwidth in [50,200,1000]:calc(f'BW{bandwidth}',category='drivers',bandwidth=bandwidth)
    # Coupled thickness changes: D~h^3, mu~h, transverse G~1/h, all assumptions explicit.
    scenarios=[]
    for scale in [.85,1.15]:
        p=s.PARAM.copy();p['h']*=scale;p['D']*=scale**3;p['mu']*=scale;p['G']/=scale;scenarios.append((f'thickness_{scale}',p,1.,.3,1.,1.))
    for sign in [-1,1]:
        p=s.PARAM.copy()
        for key,factor in [('D',.8 if sign<0 else 1.2),('mu',1.1 if sign<0 else .9),('c',.8 if sign<0 else 1.2),('G',.7 if sign<0 else 1.3),('zeta',.5)]:p[key]*=factor
        scenarios.append((f'corner_{sign}',p,.7,.38 if sign<0 else .22,1.1,1.))
    for label,p,eff,eta,flux,direction in scenarios:
        for label2,kind,filt in [('LQR',2,False),('FilteredNN',4,True)]:
            met,tr=r.run(kind,p=p,eff=eff,eta=eta,flux=flux,direction=direction,filter=filt)
            out['structured'][label+'_'+label2]=met;r.write('additional_report.json',out)
    out['elapsed_s']=time.perf_counter()-start;r.write('additional_report.json',out)

if __name__=='__main__':main()
