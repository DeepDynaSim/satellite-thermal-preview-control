"""Additional paired studies requested in peer review; run commands recorded in README."""
import simulate as s
from pathlib import Path
import json,time,argparse,itertools,platform,concurrent.futures
import numpy as np
from scipy.linalg import expm,solve_discrete_are
from scipy.stats import qmc
from revision_kernel import run_stream,forecast_many,forecast_one,polynomial

DATA=s.DATA;NET=s.load_net();DT=s.EVAL_DT
METRICS=['cost','tip_rms_m','angle_rms_arcsec','filter_fraction','infeasible_count','max_nominal_true_margin','true_margin_violation_fraction','correction_rms_normalized','correction_max_normalized','observer_error_max_scaled','failure','duration_s']

def thermal_parameters(p):
    return np.array([p['mu']*p['c']/(2*p['width']),p['G'],p['absorptivity'],p['S'],np.cos(p['beta']),p['emissivity']*p['sigma'],p['w'],p['qIR']])

def cached_kernels(H=24.,h=.12):
    tau=np.arange(int(round(H/h))+1)*h
    ee=[expm(s.ACL.T*t) for t in tau]
    k=np.array([-np.linalg.solve(s.R,s.DES['B'].T@np.linalg.solve(s.ACL.T,b-a)@s.P)@np.c_[s.DES['F1'],s.DES['F2']] for a,b in zip(ee[:-1],ee[1:])])
    return np.ascontiguousarray(k)

KERNELS={H:cached_kernels(H) for H in [6.,12.,24.,36.,48.]}

def oracle_series(T,end,direction,eta,flux,H=24.,frozen=False,pactual=None):
    grid=np.arange(0,end+.10001,.1);temp=np.empty((len(grid),2));temp[0]=T
    actual=s.PARAM if pactual is None else pactual
    for j in range(len(grid)-1):temp[j+1]=s.heat_step(temp[j],grid[j]-20,direction,eta,flux,.1,actual)
    if H==0:return np.zeros((len(grid),4))
    if frozen:
        td,tdd=s.heat(temp,grid-20,direction,eta,flux)
        return np.c_[td[:,0]-td[:,1],tdd[:,0]-tdd[:,1]]@KERNELS[H].sum(0).T
    return forecast_many(temp,grid-20,np.full(len(grid),direction),np.full(len(grid),eta),np.full(len(grid),flux),thermal_parameters(s.PARAM),KERNELS[H],.12)

def observer_design(dt=DT,process_cov=1e-6):
    nd=s.ND;C=np.zeros((5,s.NX));C[0,0]=1;C[1,nd]=1
    edges=np.array([[0,.5],[1.8,2.3],[4.5,5.]])
    for j,(a,b) in enumerate(edges):
        der=s.modes(np.array([a,b]),s.NC)[1]
        C[j+2,1:nd]=s.PARAM['h']/(2*(b-a))*(der[1]-der[0])
    obsad=expm(s.DES['A']*dt);d=s.NX
    zz=expm(np.block([[s.DES['A'],s.DES['B'],s.DES['E'][:,None]],[np.zeros((5,d+5))]])*dt)
    # Dual discrete Riccati in normalized coordinates, update at 1 ms.
    scales=np.array([.02/206264.806,.02/206264.806,5e-9,5e-9,5e-9])
    Ad1=expm(s.DES['A']*.001);D=np.diag(s.sc);An=np.linalg.solve(D,Ad1@D);Cn=C@D
    cov=solve_discrete_are(An.T,Cn.T,np.eye(d)*process_cov,np.diag(scales**2))
    L=D@cov@Cn.T@np.linalg.inv(Cn@cov@Cn.T+np.diag(scales**2))
    return obsad,zz[:d,d:d+4],zz[:d,-1],L,C

OBS=observer_design()

def prepare_plant(p,n=37,eff=1.,dt=DT):
    plant=s.model(n,p);d=n+1;B=plant['B'].copy();B[:,1:]*=eff
    zz=expm(np.block([[plant['A'],B,plant['E'][:,None]],[np.zeros((5,2*d+5))]])*dt)
    selector=np.zeros((s.NX,2*d));selector[:s.ND,:s.ND]=np.eye(s.ND);selector[s.ND:,d:d+s.ND]=np.eye(s.ND)
    fullport=np.c_[np.zeros((4,d)),plant['F'].T]
    port=fullport-np.c_[np.zeros((4,s.ND)),s.DES['F'].T]@selector
    sensor=np.zeros((5,2*d));sensor[0,0]=1;sensor[1,d]=1
    for j,(a,b) in enumerate([[0,.5],[1.8,2.3],[4.5,5.]]):
        der=s.modes(np.array([a,b]),n,p)[1]
        sensor[j+2,1:d]=p['h']/(2*(b-a))*(der[1]-der[0])
    return plant,(zz[:2*d,:2*d],zz[:2*d,2*d:2*d+4],zz[:2*d,-1],selector,port,fullport,sensor)

def run(kind=4,p=None,n=37,eff=1.,direction=1.,eta=.3,flux=1.,end=120.,H=24.,frozen=False,
        gate=True,filter=True,eps=1e-8,rmd=True,bias=0.,stride=640,core_delay=0,residual_delay=0,bandwidth=None,
        noise=None,observer=False,net=None,prepared=None,seed=7733,polynomial_model=None,dt=DT,obsdesign=None,normal_columns=None):
    p=s.PARAM.copy() if p is None else p.copy();net=NET if net is None else net
    plant,stuff=prepare_plant(p,n,eff,dt) if prepared is None else prepared
    Ad,Bd,Ed,selector,port,fullport,sensor=stuff
    T=s.initial_temperature(direction,p,eta,flux);td=s.heat(T,-20,direction,eta,flux,p)[0]
    y=np.r_[plant['rq']*(T[0]-T[1]),plant['rq']*(td[0]-td[1])]
    oracle=oracle_series(T,end,direction,eta,flux,H,frozen,p) if kind==3 else np.zeros((2,4))
    if not rmd:port=np.zeros_like(port)
    if kind==7:port=fullport
    port=np.diag([0.,5000.,5000.,5000.])@port
    noise=np.zeros(14) if noise is None else np.array(noise)
    columns=(10 if len(noise)>=20 else 7) if normal_columns is None else normal_columns
    normals=np.random.default_rng(seed).normal(size=(int(end/.001)+2,columns)) if observer or np.any(noise) else np.zeros((1,columns))
    coef,powers=(np.zeros((1,4)),np.zeros((1,6),dtype=np.int64)) if polynomial_model is None else polynomial_model
    start=time.perf_counter()
    rec,values,channels=run_stream(kind,int(round(end/dt)),dt,20.,direction,eta,flux,
        thermal_parameters(p),thermal_parameters(s.PARAM),T,y,Ad,Bd,Ed,plant['rq'],plant['tip'],selector,
        s.DES['rq'],s.DES['F1'],s.DES['F2'],s.DES['A'],s.DES['B'],s.GAIN,s.P,s.Q,s.W,s.R,s.sc,s.LIMITS,port,
        *net['weights'],net['yscale'],net.get('activation',0),oracle,bias,eps,gate,filter,coef,powers,stride,
        core_delay,residual_delay,1. if bandwidth is None else 1-np.exp(-2*np.pi*bandwidth*dt),noise,observer,
        *(OBS if obsdesign is None else obsdesign),sensor,normals)
    result=dict(zip(METRICS,map(float,values)));result.update(channels=channels.tolist(),runtime_s=time.perf_counter()-start)
    return result,rec

def write(name,obj):
    (DATA/name).write_text(json.dumps(obj,indent=2,allow_nan=False),encoding='utf-8')

def check():
    T=np.array([[270.,260.],[340.,330.],[300.,307.]])
    tau=np.array([-12.,3.,22.]);dr=np.array([1.,-1.,1.]);eta=np.array([.3,.25,.35]);fl=np.ones(3)
    pred=forecast_many(T,tau,dr,eta,fl,thermal_parameters(s.PARAM),KERNELS[24.],.12)
    exact=s.preview(T,tau,dr,eta,fl)
    test={'cached_preview_max_abs':float(np.max(np.abs(pred-exact)))}
    for name,kind,filt in [('LQR',2,False),('RawNN',4,False),('FilteredNN',4,True)]:
        start=time.perf_counter();a,rec=run(kind,end=30.,filter=filt,stride=1)
        old=s.simulate({'LQR':'Thermal LQR','RawNN':'Neural raw','FilteredNN':'Neural filtered'}[name],NET,end=30.)
        # Original interface returns metrics and decimated traces (inspect its schema).
        test[name]={'new':a,'original_cost':old[0]['cost'],'native_trace_max_abs':float(np.max(np.abs(rec-old[1]))),'elapsed':time.perf_counter()-start}
        assert abs(a['cost']-old[0]['cost'])<1e-8
        assert test[name]['native_trace_max_abs']<1e-5
    write('streaming_validation.json',test);print(json.dumps(test,indent=2),flush=True)

def uncertainty(N=320,workers=4,resume=False):
    seed=20261004;design=qmc.LatinHypercube(9,seed=seed).random(N)
    bounds=np.array([[.8,1.2],[.9,1.1],[.8,1.2],[.7,1.3],[.5,1.5],[.7,1.],[.22,.38],[.9,1.1],[0,1]])
    values=qmc.scale(design,bounds[:,0],bounds[:,1]);values[:,8]=np.where(values[:,8]<.5,-1.,1.)
    np.savez(DATA/'lhs_design.npz',unit=design,parameters=values,seed=seed)
    def case(i):
        a=values[i];p=s.PARAM.copy()
        for key,factor in zip(['D','mu','c','G','zeta'],a[:5]):p[key]*=factor
        prepared=prepare_plant(p,37,a[5]);res={}
        for label,kind,filt in [('LQR',2,False),('Oracle24',3,False),('RawNN',4,False),('FilteredNN',4,True)]:
            r,_=run(kind,p=p,eff=a[5],eta=a[6],flux=a[7],direction=a[8],filter=filt,prepared=prepared)
            res[label]=r
        path=DATA/'lhs_cases';path.mkdir(exist_ok=True);write(f'lhs_cases/case_{i:04d}.json',{'index':i,'parameters':a.tolist(),'results':res})
        return i,res
    # Compile before threading to avoid duplicate LLVM compilation.
    run(2,end=.01)
    results=[None]*N;start=time.perf_counter()
    if resume:
        for i in range(N):
            path=DATA/'lhs_cases'/f'case_{i:04d}.json'
            if path.exists():
                previous=json.loads(path.read_text())
                assert np.array_equal(np.array(previous['parameters']),values[i]),'Resume design differs'
                assert set(previous['results'])=={'LQR','Oracle24','RawNN','FilteredNN'}
                results[i]=previous['results']
        print('Resuming',sum(x is not None for x in results),'verified saved pairs',flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(case,i):i for i in range(N) if results[i] is None}
        for f in concurrent.futures.as_completed(futures):
            i,r=f.result();results[i]=r
            done=sum(x is not None for x in results)
            if done%8==0:print(f'LHS {done}/{N}; elapsed {time.perf_counter()-start:.1f}s',flush=True)
    rng=np.random.default_rng(55041);summary={}
    for label in ['Oracle24','RawNN','FilteredNN']:
        summary[label]={}
        for met in ['cost','tip_rms_m','angle_rms_arcsec']:
            base=np.array([x['LQR'][met] for x in results]);target=np.array([x[label][met] for x in results]);ratio=target/base;diff=target-base
            indices=rng.integers(0,N,(10000,N));bootmedian=np.median(ratio[indices],axis=1);bootmean=np.mean(diff[indices],axis=1)
            summary[label][met]={'median_ratio':float(np.median(ratio)),'q05_q95':np.quantile(ratio,[.05,.95]).tolist(),'min_max':np.array([ratio.min(),ratio.max()]).tolist(),'median_ratio_bootstrap95':np.quantile(bootmedian,[.025,.975]).tolist(),'mean_paired_difference':float(diff.mean()),'mean_difference_bootstrap95':np.quantile(bootmean,[.025,.975]).tolist(),'wins':int(np.sum(ratio<1))}
    write('lhs_report.json',{'N':N,'seed':seed,'bounds':bounds.tolist(),'parameter_order':['D/D0','mu/mu0','c/c0','G/G0','zeta/zeta0','piezo_efficiency','eta','solar_scale','direction'],'duration_s':120,'native_dt':DT,'results':results,'paired_summary':summary,'elapsed_s':time.perf_counter()-start})

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('task',choices=['check','lhs']);parser.add_argument('--N',type=int,default=320);parser.add_argument('--workers',type=int,default=4);parser.add_argument('--resume',action='store_true');args=parser.parse_args()
    if args.task=='check':check()
    else:uncertainty(args.N,args.workers,args.resume)
