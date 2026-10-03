"""Reproducible thermoelastic neural preview benchmark. All claims come from saved runs."""
from pathlib import Path
import os, sys, json, time, argparse
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
ROOT=Path(__file__).resolve().parents[2] if Path(__file__).parent.name=='code' else Path(__file__).resolve().parents[1]
os.environ['MPLCONFIGDIR']=str(Path(__file__).resolve().parent.parent/'.cache'/'matplotlib')
sys.path.insert(0,str(ROOT/'research_work'/'vendor'))
import numpy as np
from scipy.linalg import solve_continuous_are, expm, eigvalsh
from scipy.optimize import root
from scipy.io import savemat
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent.parent
DATA=OUT/'data'; FIG=OUT/'figures'; CODE=OUT/'code'
for p in [OUT,DATA,FIG,CODE]: p.mkdir(exist_ok=True,parents=True)
PARAM=dict(L=9.,width=3.,D=2e4,mu=7.4,Rhub=1.,Jhub=2500.,h=0.025,
    alpha_cte=2e-6,c=1044.,absorptivity=.92,emissivity=.82,sigma=5.670374419e-8,
    G=35.,qIR=120.,S=1350.,eta=.30,beta=np.deg2rad(20),w=3.,te=20.,zeta=.002,
    Ep=63e9,bp=.30,hp=.00030,rhop=7800.,d31=190e-12,patch_length=.50)
LIMITS=np.array([.30,150.,150.,150.])
BETAS=np.array([1.875104068711961,4.694091132974174,7.854757438237612,10.995540734875467,14.13716839104647,17.27875965739948,20.42035225104125,23.56194490180644,26.7035375555183,29.84513020910282,32.98672286269283,36.12831551628262,39.26990816987242,42.41150082346221,45.553093477052,48.69468613064179,51.83627878423159])
BETAS=np.r_[BETAS,np.pi*(np.arange(18,83)-.5)]

def modes(x,n,p=PARAM):
    lam=BETAS[:n]/p['L']; v=x[:,None]*lam
    beta=BETAS[:n];ee=np.exp(-beta);den=1-ee*ee+2*ee*np.sin(beta)
    gam=(1+ee*ee+2*ee*np.cos(beta))/den
    first=np.exp(v-beta)*(np.sin(beta)-np.cos(beta)-ee)/den
    second=np.exp(-v)*(1+ee*(np.sin(beta)+np.cos(beta)))/den
    hyp=first+second;hypp=first-second
    norm=(np.sin(beta)-np.cos(beta)-ee)/den+ee*(1+ee*(np.sin(beta)+np.cos(beta)))/den-np.cos(beta)+gam*np.sin(beta)
    phi=(hyp-np.cos(v)+gam*np.sin(v))/norm
    phip=lam*(hypp+np.sin(v)+gam*np.cos(v))/norm
    phipp=lam**2*(hyp+np.cos(v)-gam*np.sin(v))/norm
    return phi,phip,phipp

def model(n=3,p=PARAM,patches=True):
    # Split integration at every patch edge to avoid indicator quadrature errors.
    starts=np.array([0.,1.8,4.5]); ends=starts+p['patch_length']
    edges=sorted(set([0.,p['L'],*starts,*ends]))
    gx,gw=np.polynomial.legendre.leggauss(100 if n>17 else 60)
    xs=[];ws=[]
    for a,b in zip(edges[:-1],edges[1:]):
        xs.extend((a+b)/2+(b-a)*gx/2);ws.extend(gw*(b-a)/2)
    x=np.array(xs);weight=np.array(ws)
    phi,phip,phipp=modes(x,n,p)
    on=sum(((x>=a)&(x<=b)).astype(float) for a,b in zip(starts,ends)) if patches else np.zeros_like(x)
    zc=(p['h']+p['hp'])/2
    mup=p['mu']+on*2*p['rhop']*p['bp']*p['hp']
    Dp=p['D']+on*2*p['Ep']*p['bp']*(p['hp']**3/12+p['hp']*zc**2)
    M=np.zeros((n+1,n+1));M[0,0]=p['Jhub']+np.sum(weight*mup*(p['Rhub']+x)**2)
    M[0,1:]=M[1:,0]=np.einsum('i,ij->j',weight*mup*(p['Rhub']+x),phi)
    M[1:,1:]=phi.T@((weight*mup)[:,None]*phi)
    K=np.zeros_like(M);K[1:,1:]=phipp.T@((weight*Dp)[:,None]*phipp)
    # Modal damping defined through the bare clamped beam frequencies.
    om=BETAS[:n]**2/p['L']**2*np.sqrt(p['D']/p['mu'])
    C=np.zeros_like(M);C[1:,1:]=np.diag(2*p['zeta']*om*np.diag(M[1:,1:]))
    G=np.zeros(n+1);G[1:]=np.einsum('i,ij->j',weight*p['D']*p['alpha_cte']/p['h'],phipp)
    F=np.zeros((n+1,4));F[0,0]=1.
    chip=2*p['Ep']*p['bp']*p['d31']*zc
    for j,(a,b) in enumerate(zip(starts,ends)):
        _,dp,_=modes(np.array([a,b]),n,p);F[1:,j+1]=chip*(dp[1]-dp[0])
    d=n+1;Mi=np.linalg.inv(M)
    A=np.block([[np.zeros((d,d)),np.eye(d)],[-Mi@K,-Mi@C]])
    B=np.vstack([np.zeros((d,4)),Mi@F]);E=np.r_[np.zeros(d),Mi@G]
    rq=np.r_[0.,np.linalg.solve(K[1:,1:],G[1:])]
    F1=np.r_[np.zeros(d),-Mi@C@rq];F2=np.r_[np.zeros(d),-rq]
    return dict(n=n,M=M,K=K,C=C,F=F,G=G,A=A,B=B,E=E,rq=rq,F1=F1,F2=F2,
                freq=np.sqrt(eigvalsh(K[1:,1:],M[1:,1:]))/(2*np.pi),tip=modes(np.array([p['L']]),n,p)[0][0])

def static_reference(n=5,p=PARAM):
    edges=[0.,.5,1.8,2.3,4.5,5.,p['L']]
    g,w=np.polynomial.legendre.leggauss(100)
    xx=np.concatenate([(a+b)/2+(b-a)*g/2 for a,b in zip(edges[:-1],edges[1:])]);ww=np.concatenate([w*(b-a)/2 for a,b in zip(edges[:-1],edges[1:])])
    Dadd=2*p['Ep']*p['bp']*(p['hp']**3/12+p['hp']*((p['h']+p['hp'])/2)**2)
    wq=np.zeros(len(xx));MT=p['D']*p['alpha_cte']/p['h']
    for a,b in zip(edges[:-1],edges[1:]):
        Dsegment=p['D']+(Dadd if (a,b) in [(0.,.5),(1.8,2.3),(4.5,5.)] else 0.)
        ce=np.minimum(xx,b);mask=xx>a
        wq[mask]+=MT/Dsegment*(xx[mask]*(ce[mask]-a)-(ce[mask]**2-a*a)/2)
    phi=modes(xx,n,p)[0]
    return np.r_[0.,(phi.T@(ww*wq))/(p['L']/4)]

NC=5;ND=NC+1;NX=2*ND
EVAL_N=37;EVAL_DT=1/6400
DES=model(NC)
# Static residual flexibility is calibrated from the exact piecewise-rigidity shape.
DES['rq']=static_reference(NC);DES['G']=DES['K']@DES['rq'];DES['E']=np.r_[np.zeros(ND),np.linalg.solve(DES['M'],DES['G'])]
DES['F1']=np.r_[np.zeros(ND),-np.linalg.solve(DES['M'],DES['C']@DES['rq'])];DES['F2']=np.r_[np.zeros(ND),-DES['rq']]
sc=np.r_[1e-4,np.repeat(2e-3,3),np.repeat(2e-2,NC-3),5e-4,np.repeat(5e-3,3),np.repeat(5e-2,NC-3)]
Q=np.diag(sc**-2);R=np.diag(LIMITS**-2)
# Solve in scaled coordinates to avoid ill-conditioning of SI coordinates.
Dscale=np.diag(sc);Ai=np.linalg.solve(Dscale,DES['A']@Dscale);Bi=np.linalg.solve(Dscale,DES['B'])
Pi=solve_continuous_are(Ai,Bi,np.eye(NX),R);P=np.linalg.solve(Dscale,np.linalg.solve(Dscale,Pi).T).T
GAIN=np.linalg.solve(R,DES['B'].T@P);ACL=DES['A']-DES['B']@GAIN
W=Q+GAIN.T@R@GAIN
DES.update(Q=Q,R=R,P=P,GAIN=GAIN,W=W,ACL=ACL)

def sun(tau,direction,p=PARAM):
    z=np.tanh(np.asarray(tau)/p['w']);chi=(1+direction*z)/2
    dchi=direction*(1-z*z)/(2*p['w'])
    return chi,dchi

def heat(T,tau,direction,eta=.30,flux=1.,p=PARAM,incidence=None):
    T=np.asarray(T);chi,dchi=sun(tau,direction,p)
    ca=p['mu']/p['width']*p['c']/2
    cosbeta=np.cos(p['beta']) if incidence is None else np.maximum(0.,np.cos(incidence))
    solar=(p['absorptivity']-eta)*p['S']*flux*cosbeta*chi
    diff=T[...,0]-T[...,1];rad=p['emissivity']*p['sigma']*(T**4-3.**4)
    td=np.stack([(solar+p['qIR']-p['G']*diff-rad[...,0])/ca,
                 (p['qIR']+p['G']*diff-rad[...,1])/ca],axis=-1)
    j11=(-p['G']-4*p['emissivity']*p['sigma']*T[...,0]**3)/ca
    j22=(-p['G']-4*p['emissivity']*p['sigma']*T[...,1]**3)/ca
    jt=p['G']/ca
    forcing=(p['absorptivity']-eta)*p['S']*flux*cosbeta*dchi/ca
    tdd=np.stack([j11*td[...,0]+jt*td[...,1]+forcing,
                  jt*td[...,0]+j22*td[...,1]],axis=-1)
    return td,tdd

def heat_step(T,tau,direction,eta,flux,h,p=PARAM):
    f=lambda V,s:heat(V,s,direction,eta,flux,p)[0]
    k1=f(T,tau);k2=f(T+h*k1/2,tau+h/2);k3=f(T+h*k2/2,tau+h/2);k4=f(T+h*k3,tau+h)
    return T+h*(k1+2*k2+2*k3+k4)/6

def preview(T,tau,direction,eta,flux,H=24.,h=.12):
    T=np.atleast_2d(T).copy();N=len(T);tau=np.broadcast_to(tau,(N,));direction=np.broadcast_to(direction,(N,))
    eta=np.broadcast_to(eta,(N,));flux=np.broadcast_to(flux,(N,))
    acc=np.zeros((N,4));steps=int(round(H/h))
    for j in range(steps):
        Tmid=heat_step(T,tau+j*h,direction,eta,flux,h/2)
        td,tdd=heat(Tmid,tau+(j+.5)*h,direction,eta,flux)
        dd=td[:,0]-td[:,1];dda=tdd[:,0]-tdd[:,1]
        integ=np.linalg.solve(ACL.T,expm(ACL.T*((j+1)*h))-expm(ACL.T*(j*h)))
        kernel=-np.linalg.solve(R,DES['B'].T@integ@P)
        kval=np.c_[kernel@DES['F1'],kernel@DES['F2']]
        acc+=dd[:,None]*kval[:,0]+dda[:,None]*kval[:,1]
        T=heat_step(T,tau+j*h,direction,eta,flux,h)
    return acc

def features(T,tau,direction,eta,flux):
    T=np.atleast_2d(T);N=len(T)
    return np.stack([(T[:,0]+T[:,1]-560)/110,(T[:,0]-T[:,1])/40,
                     np.tanh(np.broadcast_to(tau,(N,))/18),np.broadcast_to(direction,(N,)),
                     (np.broadcast_to(eta,(N,))-.30)/.08,(np.broadcast_to(flux,(N,))-1)/.2],axis=1)

def train(seed=20261003,N=6000,epochs=950):
    rng=np.random.default_rng(seed)
    mean=rng.uniform(218,370,N);diff=rng.uniform(-30,35,N)
    T=np.c_[mean+diff/2,mean-diff/2];tau=rng.uniform(-55,65,N);direction=rng.choice([-1.,1.],N)
    eta=rng.uniform(.22,.38,N);flux=rng.uniform(.85,1.15,N)
    X=features(T,tau,direction,eta,flux);Y=preview(T,tau,direction,eta,flux)
    ix=rng.permutation(N);a=int(.7*N);b=int(.85*N);itr,iva,ite=ix[:a],ix[a:b],ix[b:]
    yscale=np.std(Y[itr],axis=0);yscale=np.maximum(yscale,LIMITS*.002)
    Ys=Y/yscale
    dims=[6,48,32,4];weights=[]
    for n,m in zip(dims[:-1],dims[1:]):weights.extend([rng.normal(0,np.sqrt(2/(n+m)),(n,m)),np.zeros(m)])
    mom=[np.zeros_like(w) for w in weights];var=[np.zeros_like(w) for w in weights]
    history=[];best=np.inf;bestw=None;step=0
    def forward(xx,w):
        h1=np.tanh(xx@w[0]+w[1]);h2=np.tanh(h1@w[2]+w[3]);return h1,h2,h2@w[4]+w[5]
    for ep in range(epochs):
        for sub in np.array_split(rng.permutation(itr),max(1,len(itr)//256)):
            x=X[sub];y=Ys[sub];h1,h2,pred=forward(x,weights)
            dy=(pred-y)/(len(sub)*4)
            dw4=h2.T@dy;db5=dy.sum(0);dh2=(dy@weights[4].T)*(1-h2*h2)
            dw2=h1.T@dh2;db3=dh2.sum(0);dh1=(dh2@weights[2].T)*(1-h1*h1)
            dw0=x.T@dh1;db1=dh1.sum(0)
            grad=[dw0,db1,dw2,db3,dw4,db5];step+=1
            lr=.0025*(.35+.65*(1-ep/epochs))
            for j in range(6):
                mom[j]=.9*mom[j]+.1*grad[j];var[j]=.999*var[j]+.001*grad[j]**2
                weights[j]-=lr*(mom[j]/(1-.9**step))/(np.sqrt(var[j]/(1-.999**step))+1e-8)
        if ep%10==0 or ep==epochs-1:
            v=np.mean((forward(X[iva],weights)[2]-Ys[iva])**2)
            tr=np.mean((forward(X[itr],weights)[2]-Ys[itr])**2);history.append([ep,tr,v])
            if v<best:best=v;bestw=[w.copy() for w in weights]
            if ep%100==0:print(f'train epoch {ep} validation {v:.6g}',flush=True)
    pred=forward(X[ite],bestw)[2]*yscale
    errors=pred-Y[ite]
    net=dict(weights=bestw,yscale=yscale)
    np.savez(DATA/'network.npz',**{f'w{i}':w for i,w in enumerate(bestw)},yscale=yscale)
    np.savez(DATA/'learning.npz',X=X,Y=Y,itr=itr,iva=iva,ite=ite,pred=pred,history=np.array(history))
    report=dict(N=N,split=[a,b-a,N-b],seed=seed,epochs=epochs,architecture=dims,
                test_rmse=np.sqrt(np.mean(errors**2,0)).tolist(),test_relative_L2=float(np.linalg.norm(errors/yscale)/np.linalg.norm(Y[ite]/yscale)),
                test_relative_L2_in_chosen_SI_units=float(np.linalg.norm(errors)/np.linalg.norm(Y[ite])),
                max_absolute_error=np.max(abs(errors),axis=0).tolist(),best_validation_scaled_mse=float(best))
    return net,report

def load_net():
    a=np.load(DATA/'network.npz');return dict(weights=[a[f'w{i}'] for i in range(6)],yscale=a['yscale'])
def nn(T,tau,direction,eta,flux,net):
    w=net['weights'];x=features(T,tau,direction,eta,flux)
    td,_=heat(T,tau,direction,eta,flux);_,dchi=sun(tau,direction)
    activity=float(np.sum(td**2)+(.3*dchi)**2)
    gate=activity/(activity+1e-8)
    return gate*(np.tanh(np.tanh(x@w[0]+w[1])@w[2]+w[3])@w[4]+w[5])[0]*net['yscale']

def project(u,a,b,baseline=None,Rdiag=np.diag(R)):
    v=np.clip(u,-LIMITS,LIMITS)
    if a@v<=b+1e-10:return v,0.,0.
    if baseline is not None:
        base=np.clip(baseline,-LIMITS,LIMITS)
        if a@base<=b+1e-10:
            denom=a@(v-base)
            alpha=np.clip((b-a@base)/max(denom,1e-30),0.,1.)
            return base+alpha*(v-base),1.,0.
    minval=-np.abs(a)@LIMITS
    if minval>b+1e-10:return -LIMITS*np.sign(a),1.,float(minval-b)
    lo=0.;hi=1.
    for _ in range(60):
        if a@np.clip(u-hi*a/Rdiag,-LIMITS,LIMITS)<=b:break
        hi*=2
    for _ in range(32):
        mid=(lo+hi)/2
        if a@np.clip(u-mid*a/Rdiag,-LIMITS,LIMITS)>b:lo=mid
        else:hi=mid
    return np.clip(u-hi*a/Rdiag,-LIMITS,LIMITS),1.,0.

def initial_temperature(direction,p,eta,flux):
    T0=(p['qIR']/(p['emissivity']*p['sigma'])+3**4)**.25
    if direction>0:return np.array([T0,T0])
    sol=root(lambda T:heat(T,1000,-1,eta,flux,p)[0],np.array([340.,330.]))
    # sunset direction -1 at tau <<0 is full sunlight.
    sol=root(lambda T:heat(T,-1000,-1,eta,flux,p)[0],np.array([340.,330.]))
    assert sol.success
    return sol.x

def simulate(kind,net,direction=1.,p=PARAM,eta=.30,flux=1.,acteff=1.,nplant=EVAL_N,dt=EVAL_DT,end=180.,kick=False,save=False,learn_bias=0.):
    plant=model(nplant,p);d=nplant+1
    A=plant['A'];Bp=plant['B'].copy();Bp[:,1:]*=acteff
    Z=np.block([[A,Bp,plant['E'][:,None]],[np.zeros((5,2*d+5))]])
    ZZ=expm(Z*dt);Ad=ZZ[:2*d,:2*d];Bd=ZZ[:2*d,2*d:2*d+4];Ed=ZZ[:2*d,-1]
    T=initial_temperature(direction,p,eta,flux)
    td,_=heat(T,-p['te'],direction,eta,flux,p)
    y=np.r_[plant['rq']*(T[0]-T[1]),plant['rq']*(td[0]-td[1])]
    if kick:y[0]=1e-5;y[1]=2e-4
    steps=int(round(end/dt));rec=np.empty((steps+1,18));cost=0.;infeasible=0;active=0
    port=np.c_[np.zeros((4,d)),plant['F'].T]
    selector=np.zeros((NX,2*d));selector[:ND,:ND]=np.eye(ND);selector[ND:,d:d+ND]=np.eye(ND)
    port-=np.c_[np.zeros((4,ND)),DES['F'].T]@selector
    rmd_gain=np.diag([0.,5000.,5000.,5000.])
    oracle=None
    if kind=='Oracle':
        # Batch offline nominal preview labels on a coarser grid with interpolation.
        og=np.arange(0,end+.10001,.1);TT=[];tv=T.copy()
        for t in og:TT.append(tv.copy());tv=heat_step(tv,t-p['te'],direction,eta,flux,.1,p)
        oracle=preview(np.array(TT),og-p['te'],direction,eta,flux)
    loop_count=steps+1
    if os.environ.get('SCI_NO_JIT')!='1':
        from fast_kernel import run_core
        def thermal_pack(pp):return np.array([pp['mu']/pp['width']*pp['c']/2,pp['G'],pp['absorptivity'],pp['S'],np.cos(pp['beta']),pp['emissivity']*pp['sigma'],pp['w'],pp['qIR']])
        index={'Passive':0,'Fixed LQR':1,'Thermal LQR':2,'Oracle':3,'Neural raw':4,'Neural filtered':5}[kind]
        rec,cost,active,infeasible=run_core(index,steps,dt,p['te'],float(direction),eta,flux,thermal_pack(p),thermal_pack(PARAM),T,y,
            Ad,Bd,Ed,plant['rq'],plant['tip'],selector,DES['rq'],DES['F1'],DES['F2'],DES['A'],DES['B'],GAIN,P,Q,W,R,sc,LIMITS,rmd_gain@port,
            *net['weights'],net['yscale'],oracle if oracle is not None else np.zeros((2,4)),learn_bias)
        loop_count=0
    for j in range(loop_count):
        t=j*dt;tau=t-p['te']
        tdn,tddn=heat(T,tau,direction,eta,flux,PARAM)
        delta=T[0]-T[1];dd=tdn[0]-tdn[1];dda=tddn[0]-tddn[1]
        rr=np.r_[DES['rq']*delta,DES['rq']*dd]
        yc=np.r_[y[:ND],y[d:d+ND]];x=yc-rr
        f=DES['F1']*dd+DES['F2']*dda
        ulqr=-GAIN@x
        if kind=='Passive':u=np.array([-150.*y[0]-550.*y[d],0.,0.,0.])
        elif kind=='Fixed LQR':u=-GAIN@yc
        elif kind=='Thermal LQR':u=ulqr
        elif kind=='Oracle':u=ulqr+np.array([np.interp(t,og,oracle[:,k]) for k in range(4)])
        else:
            u=ulqr+nn(T,tau,direction,eta,flux,net)
            if learn_bias:u[1:]+=learn_bias*np.sin(.7*t)
        auxiliary=np.zeros(4) if kind=='Passive' else -rmd_gain@port@y
        u+=auxiliary
        u0=u.copy();act=0.;slack=0.
        a=2*DES['B'].T@P@x
        f_eff=f+DES['B']@auxiliary
        b=-.25*(x@W@x)+2*abs(x@P@f_eff)-2*x@P@(DES['A']@x+f)
        if kind=='Neural filtered':u,act,slack=project(u,a,b,ulqr+auxiliary)
        else:u=np.clip(u,-LIMITS,LIMITS)
        xfull=y-np.r_[plant['rq']*delta,plant['rq']*dd]
        vibration=plant['tip']@xfull[1:d]
        static_tip=plant['tip']@plant['rq'][1:]*delta
        actual_tip=plant['tip']@y[1:d]
        Vdot=2*x@P@(DES['A']@x+DES['B']@u+f)
        margin=Vdot+.25*x@W@x-2*abs(x@P@f_eff)
        rec[j]=[t,T[0],T[1],delta,static_tip,actual_tip,vibration,y[0],y[d],*u,float(x@Q@x+u@R@u),act,slack,margin,np.linalg.norm(x/sc)]
        active+=act;infeasible+=slack>1e-8
        if j<steps:
            cost+=dt*(x@Q@x+u@R@u)
            Tnext=heat_step(T,tau,direction,eta,flux,dt,p)
            y=Ad@y+Bd@u+Ed*((delta+Tnext[0]-Tnext[1])/2)
            T=Tnext
    idx=rec[:,0]>=5
    m=dict(controller=kind,direction=float(direction),rms_vibration_um=float(np.sqrt(np.mean(rec[idx,6]**2))*1e6),
       peak_vibration_um=float(np.max(abs(rec[:,6]))*1e6),rms_attitude_arcsec=float(np.sqrt(np.mean(rec[idx,7]**2))*206264.806247),
       peak_attitude_arcsec=float(np.max(abs(rec[:,7]))*206264.806247),cost=float(cost),
       voltage_squared_integral=float(np.trapezoid(np.sum(rec[:,10:13]**2,axis=1),rec[:,0])),
       peak_voltage=float(np.max(abs(rec[:,10:13]))),peak_torque=float(np.max(abs(rec[:,9]))),
       filter_active_percent=float(100*active/(steps+1)),filter_infeasible_count=int(infeasible),
       max_nominal_certificate_violation=float(np.max(rec[:,16])),
       saturation_percent=float(np.mean(np.any(abs(rec[:,9:13])>=LIMITS*(1-1e-8),axis=1))*100),
       static_tip_final_mm=float(rec[-1,4]*1000))
    if save:np.savez_compressed(DATA/(kind.replace(' ','_')+('_sunrise' if direction>0 else '_sunset')+'.npz'),trace=rec[::max(1,int(round(.005/dt)))],native_dt=dt,metrics=json.dumps(m))
    return m,rec

def figures(report):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'axes.titlesize':11,
       'legend.fontsize':8,'lines.linewidth':1.4,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':360})
    colors={'Passive':'#777777','Fixed LQR':'#C04940','Thermal LQR':'#D68B13','Oracle':'#24895C','Neural filtered':'#1867B4','Neural raw':'#8650A5'}
    def finish(fig,name):
        fig.savefig(FIG/(name+'.png'),dpi=360,bbox_inches='tight');fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight');fig.savefig(FIG/(name+'.svg'),bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(2,1,figsize=(7,4.6),sharex=True)
    for direction,ls in [(1,'-'),(-1,'--')]:
        tr=np.load(DATA/('Neural_filtered'+('_sunrise' if direction>0 else '_sunset')+'.npz'))['trace']
        ax[0].plot(tr[:,0],tr[:,1],ls,label=('Exit' if direction>0 else 'Entry')+' front')
        ax[0].plot(tr[:,0],tr[:,2],ls,label=('Exit' if direction>0 else 'Entry')+' back')
        ax[1].plot(tr[:,0],tr[:,3],ls,label='Exit' if direction>0 else 'Entry')
    ax[0].set_ylabel('Temperature (K)');ax[1].set_ylabel('Front minus back (K)');ax[1].set_xlabel('Time (s)');ax[0].legend(ncol=2);ax[1].legend();finish(fig,'fig02_thermal')
    fig,ax=plt.subplots(3,1,figsize=(7,6),sharex=True)
    for kind in ['Passive','Thermal LQR','Oracle','Neural filtered']:
        tr=np.load(DATA/(kind.replace(' ','_')+'_sunrise.npz'))['trace'];c=colors[kind]
        ax[0].plot(tr[:,0],tr[:,6]*1e6,label=kind,color=c)
        ax[1].plot(tr[:,0],tr[:,7]*206264.806247,label=kind,color=c)
        ax[2].plot(tr[:,0],tr[:,10],label=kind,color=c)
    ax[0].set_ylabel('Dynamic tip error (µm)');ax[1].set_ylabel('Hub angle (arcsec)');ax[2].set_ylabel('Patch 1 voltage (V)');ax[2].set_xlabel('Time (s)');ax[0].legend(ncol=2);ax[0].set_xlim(0,90);finish(fig,'fig03_exit_response')
    fig,ax=plt.subplots(2,1,figsize=(7,4.7),sharex=True)
    for kind in ['Passive','Thermal LQR','Oracle','Neural filtered']:
        tr=np.load(DATA/(kind.replace(' ','_')+'_sunset.npz'))['trace'];c=colors[kind]
        ax[0].plot(tr[:,0],tr[:,6]*1e6,label=kind,color=c);ax[1].plot(tr[:,0],tr[:,7]*206264.806247,color=c)
    ax[0].set_ylabel('Dynamic tip error (µm)');ax[1].set_ylabel('Hub angle (arcsec)');ax[1].set_xlabel('Time (s)');ax[0].legend(ncol=2);ax[0].set_xlim(0,90);finish(fig,'fig04_entry_response')
    learning=np.load(DATA/'learning.npz');Y=learning['Y'][learning['ite']];yp=learning['pred']
    fig,ax=plt.subplots(1,2,figsize=(7,3.1))
    for k in range(1,4):ax[0].scatter(Y[:,k],yp[:,k],s=5,alpha=.35,label=f'Patch {k}')
    lim=np.max(abs(Y[:,1:]));ax[0].plot([-lim,lim],[-lim,lim],color='black',lw=.7);ax[0].set_xlabel('Optimal preview voltage (V)');ax[0].set_ylabel('Neural prediction (V)');ax[0].legend()
    hist=learning['history'];ax[1].semilogy(hist[:,0],hist[:,1],label='Training');ax[1].semilogy(hist[:,0],hist[:,2],label='Validation');ax[1].set_xlabel('Epoch');ax[1].set_ylabel('Scaled mean squared error');ax[1].legend();finish(fig,'fig05_learning')
    mc=report['monte_carlo'];fig,ax=plt.subplots(1,2,figsize=(7,3.4))
    for a,key,title in zip(ax,['rms_vibration_um','rms_attitude_arcsec'],['Dynamic tip RMS (µm)','Attitude RMS (arcsec)']):
        vals=[[r[key] for r in mc if r['controller']==k] for k in ['Thermal LQR','Neural filtered']]
        a.boxplot(vals,tick_labels=['Thermal LQR','Neural preview'],showfliers=True);a.set_ylabel(title)
    finish(fig,'fig06_robustness')
    fig,ax=plt.subplots(2,1,figsize=(7,4.3),sharex=True)
    tr=np.load(DATA/'Neural_filtered_sunrise.npz')['trace']
    ax[0].plot(tr[:,0],tr[:,4]*1e3,label='Quasistatic thermal shape');ax[0].plot(tr[:,0],tr[:,5]*1e3,'--',label='Total tip deflection');ax[0].set_ylabel('Tip deflection (mm)');ax[0].legend()
    ax[1].plot(tr[:,0],tr[:,16],label='Nominal inequality residual');ax[1].axhline(0,color='black',lw=.6);ax[1].set_ylabel('CLF residual (1/s)');ax[1].set_xlabel('Time (s)');finish(fig,'fig07_equilibrium_certificate')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--reuse',action='store_true');parser.add_argument('--quick',action='store_true');parser.add_argument('--mc',type=int,default=40);args=parser.parse_args()
    start=time.time()
    if args.reuse:
        net=load_net();learning=json.loads((DATA/'report.json').read_text())['learning']
    else:net,learning=train(N=2400 if args.quick else 6000,epochs=300 if args.quick else 950)
    report=dict(parameters={k:float(v) for k,v in PARAM.items()},learning=learning,controller_modes=NC,evaluation_modes=EVAL_N,simulation_step=EVAL_DT,residual_damping_gain=5000.,nominal=[],monte_carlo=[],verification={})
    for direction in [1.,-1.]:
        for kind in ['Passive','Fixed LQR','Thermal LQR','Oracle','Neural raw','Neural filtered']:
            m,_=simulate(kind,net,direction,save=True);report['nominal'].append(m);print(json.dumps(m),flush=True)
    rng=np.random.default_rng(84563)
    for j in range(args.mc):
        p=PARAM.copy();p['D']*=rng.uniform(.8,1.2);p['mu']*=rng.uniform(.9,1.1);p['c']*=rng.uniform(.8,1.2);p['G']*=rng.uniform(.7,1.3);p['zeta']*=rng.uniform(.5,1.5)
        eta=rng.uniform(.22,.38);flux=rng.uniform(.9,1.1);eff=rng.uniform(.7,1.);direction=rng.choice([-1.,1.])
        ap=model(EVAL_N,p);dd=EVAL_N+1;hh=np.zeros((NX,2*dd));hh[:ND,:ND]=np.eye(ND);hh[ND:,dd:dd+ND]=np.eye(ND)
        bb=ap['B'].copy();bb[:,1:]*=eff
        residual_port=np.c_[np.zeros((4,dd)),ap['F'].T]-np.c_[np.zeros((4,ND)),DES['F'].T]@hh
        spectral=float(np.max(np.linalg.eigvals(ap['A']-bb@GAIN@hh-bb@np.diag([0.,5000.,5000.,5000.])@residual_port).real))
        for kind in ['Thermal LQR','Neural filtered']:
            m,_=simulate(kind,net,direction,p,eta,flux,eff,end=120.);m.update(sample=j,acteff=eff,eta=eta,flux=flux,D=p['D'],mu=p['mu'],c=p['c'],G=p['G'],zeta=p['zeta'],linear_spectral_abscissa=spectral);report['monte_carlo'].append(m)
        if j%5==0:print(f'Uncertainty sample {j+1}/{args.mc}',flush=True)
    care_res=DES['A'].T@P+P@DES['A']-P@DES['B']@np.linalg.solve(R,DES['B'].T)@P+Q
    report['verification']['care_relative_residual']=float(np.linalg.norm(care_res)/np.linalg.norm(Q))
    report['verification']['P_min_eigenvalue']=float(eigvalsh(P)[0]);report['verification']['M_min_eigenvalue']=float(eigvalsh(DES['M'])[0])
    report['verification']['closed_loop_max_real_part']=float(np.max(np.linalg.eigvals(ACL).real))
    pl7=model(EVAL_N);dd=EVAL_N+1;Hsel=np.zeros((NX,2*dd));Hsel[:ND,:ND]=np.eye(ND);Hsel[ND:,dd:dd+ND]=np.eye(ND)
    rp=np.c_[np.zeros((4,dd)),pl7['F'].T]-np.c_[np.zeros((4,ND)),DES['F'].T]@Hsel
    controller=GAIN@Hsel+np.diag([0.,5000.,5000.,5000.])@rp
    report['verification']['high_order_linear_max_real_part']=float(np.max(np.linalg.eigvals(pl7['A']-pl7['B']@controller).real))
    sysaug=expm(np.block([[pl7['A'],pl7['B']],[np.zeros((4,2*dd+4))]])*EVAL_DT)
    report['verification']['high_order_sampled_spectral_radius']=float(np.max(abs(np.linalg.eigvals(sysaug[:2*dd,:2*dd]-sysaug[:2*dd,2*dd:]@controller))))
    low=model(3);scl=np.r_[1e-4,np.repeat(2e-3,3),5e-4,np.repeat(5e-3,3)];dl=np.diag(scl)
    pil=solve_continuous_are(np.linalg.solve(dl,low['A']@dl),np.linalg.solve(dl,low['B']),np.eye(8),R)
    pp=np.linalg.solve(dl,np.linalg.solve(dl,pil).T).T;kl=np.linalg.solve(R,low['B'].T@pp)
    pl5=model(5);hl=np.zeros((8,12));hl[:4,:4]=np.eye(4);hl[4:,6:10]=np.eye(4)
    report['verification']['three_mode_spillover_counterexample_real_part']=float(np.max(np.linalg.eigvals(pl5['A']-pl5['B']@kl@hl).real))
    report['verification']['clamped_patched_frequencies_hz']=model(5)['freq'].tolist()
    bare=model(5,patches=False);exact=BETAS[:5]**2/PARAM['L']**2*np.sqrt(PARAM['D']/PARAM['mu'])/(2*np.pi)
    report['verification']['bare_frequency_max_relative_error']=float(np.max(abs(bare['freq']/exact-1)))
    m1,tr1=simulate('Neural filtered',net,nplant=EVAL_N,dt=EVAL_DT)
    m2,tr2=simulate('Neural filtered',net,nplant=EVAL_N,dt=EVAL_DT/2)
    m3,tr3=simulate('Neural filtered',net,nplant=41,dt=EVAL_DT/2)
    report['verification']['time_refinement_rms_relative_change']=float(abs(m2['rms_vibration_um']/m1['rms_vibration_um']-1))
    report['verification']['modal_refinement_rms_relative_change']=float(abs(m3['rms_vibration_um']/m2['rms_vibration_um']-1))
    report['verification']['modal_refinement_rms_absolute_change_um']=float(abs(m3['rms_vibration_um']-m2['rms_vibration_um']))
    report['verification']['modal_refinement_cost_relative_change']=float(abs(m3['cost']/m2['cost']-1))
    report['verification']['modal_refinement_attitude_relative_change']=float(abs(m3['rms_attitude_arcsec']/m2['rms_attitude_arcsec']-1))
    # Reference derivatives against a centered finite difference away from sharp events.
    T=np.array([300.,285.]);tau=1.7;h=1e-3;td,tdd=heat(T,tau,1.)
    Tm=heat_step(T,tau,1.,.30,1.,-h);Tp=heat_step(T,tau,1.,.30,1.,h)
    num=(heat(Tp,tau+h,1.)[0]-heat(Tm,tau-h,1.)[0])/(2*h)
    report['verification']['thermal_second_derivative_relative_error']=float(np.linalg.norm(num-tdd)/np.linalg.norm(tdd))
    # Costate quadrature refinement independently of NN training.
    sam=np.array([[240.,240.],[330.,318.],[280.,260.]])
    ar=preview(sam,np.array([-3.,0.,10.]),np.array([1.,-1.,1.]),.30,1.,h=.12)
    af=preview(sam,np.array([-3.,0.,10.]),np.array([1.,-1.,1.]),.30,1.,h=.06)
    report['verification']['preview_quadrature_relative_change']=float(np.linalg.norm(ar-af)/np.linalg.norm(af))
    # Plant dissipation and angular momentum without forcing or external wheel torque.
    plant=model(5);y=np.r_[0.,np.repeat(1e-4,5),0.,np.repeat(1e-4,5)];Ap=plant['A'];Ad=expm(Ap*.05)
    momentum0=plant['M'][0]@y[6:];energy=[];mom=[]
    for _ in range(1000):
        energy.append((y[6:]@plant['M']@y[6:]+y[:6]@plant['K']@y[:6])/2);mom.append(plant['M'][0]@y[6:]);y=Ad@y
    report['verification']['unforced_max_energy_increase']=float(np.max(np.diff(energy)))
    report['verification']['unforced_relative_momentum_drift']=float(np.max(abs(np.array(mom)-momentum0))/max(abs(momentum0),1e-12))
    # Deliberately corrupted learned voltages test the safeguard as a separate ablation.
    stress=[]
    for kind in ['Neural raw','Neural filtered']:
        m,_=simulate(kind,net,end=90.,learn_bias=80.);stress.append(m)
    report['verification']['corrupted_network_ablation']=stress
    report['elapsed_seconds']=time.time()-start
    (DATA/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    np.savetxt(DATA/'nominal_metrics.csv',np.array([[r[k] for k in ['direction','rms_vibration_um','peak_vibration_um','rms_attitude_arcsec','cost','peak_voltage','filter_active_percent']] for r in report['nominal']]),delimiter=',',header='direction,rms_vibration_um,peak_vibration_um,rms_attitude_arcsec,cost,peak_voltage,filter_active_percent',comments='')
    grid=np.linspace(0,30,1501);xinitial=sc*np.linspace(.02,.15,NX)
    val=solve_ivp(lambda t,x:ACL@x+DES['F1']*.2*np.sin(.4*t)+DES['F2']*.08*np.cos(.4*t),[0,30],xinitial,t_eval=grid,rtol=1e-11,atol=1e-13,method='DOP853')
    assert val.success
    savemat(DATA/'independent_validation.mat',{k:DES[k] for k in ['M','K','C','F','G','A','B','E','rq','F1','F2','Q','R','P','GAIN','W','ACL']}|{'sc':sc,'beta_roots':BETAS,'parameters_vector':np.array([PARAM['L'],PARAM['D'],PARAM['mu']]),'validation_t':grid,'validation_x':val.y.T,'validation_x0':xinitial},do_compression=True)
    figures(report)
    print('VERIFICATION '+json.dumps(report['verification']),flush=True)

if __name__=='__main__':main()

