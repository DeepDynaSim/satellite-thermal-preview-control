"""Independent cubic Hermite FE assembly and conservative spanwise heat model."""
import revision_studies as r
import simulate as s
import numpy as np,json,time
from scipy.linalg import eigh,expm
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import eigsh
from scipy.integrate import solve_ivp
from numba import njit
from revision_kernel import features
from fast_kernel import thermal,projection

def hermite(v,ell):
    N=np.array([1-3*v*v+2*v**3,ell*(v-2*v*v+v**3),3*v*v-2*v**3,ell*(-v*v+v**3)])
    Nd=np.array([(-6*v+6*v*v)/ell,1-4*v+3*v*v,(6*v-6*v*v)/ell,-2*v+3*v*v])
    Ndd=np.array([(-6+12*v)/ell**2,(-4+6*v)/ell,(6-12*v)/ell**2,(-2+6*v)/ell])
    return N,Nd,Ndd

def assemble(ne=120,nm=37,p=s.PARAM,patches=True,ncells=24):
    edges=np.unique(np.round(np.r_[np.linspace(0,p['L'],ne+1),0,.5,1.8,2.3,4.5,5.],12))
    ne=len(edges)-1;nf=2*(ne+1);M=np.zeros((nf,nf));K=np.zeros_like(M);mh=np.zeros(nf)
    Mhh=p['Jhub'];g=np.zeros((nf,ncells));force=np.zeros((nf,3));probe=np.zeros((3,nf))
    projection=np.zeros((s.NC,nf));gram=np.zeros((s.NC,s.NC))
    gx,gw=np.polynomial.legendre.leggauss(5);zc=(p['h']+p['hp'])/2
    starts=np.array([0,1.8,4.5]);ends=starts+.5
    for e,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
        dofs=np.arange(2*e,2*e+4);ell=b-a;mid=(a+b)/2
        zone=np.flatnonzero((mid>=starts)&(mid<=ends));on=bool(len(zone)) and patches
        mu=p['mu']+(2*p['rhop']*p['bp']*p['hp'] if on else 0)
        D=p['D']+(2*p['Ep']*p['bp']*(p['hp']**3/12+p['hp']*zc**2) if on else 0)
        # Split every element at thermal-cell boundaries as well.
        cuts=np.unique(np.round(np.r_[a,b,np.linspace(0,p['L'],ncells+1)[(np.linspace(0,p['L'],ncells+1)>a+1e-11)&(np.linspace(0,p['L'],ncells+1)<b-1e-11)]],12))
        for aa,bb in zip(cuts[:-1],cuts[1:]):
            for x,w in zip((aa+bb)/2+(bb-aa)*gx/2,gw*(bb-aa)/2):
                N,Nd,Ndd=hermite((x-a)/ell,ell);radius=p['Rhub']+x
                M[np.ix_(dofs,dofs)]+=w*mu*np.outer(N,N);K[np.ix_(dofs,dofs)]+=w*D*np.outer(Ndd,Ndd)
                mh[dofs]+=w*mu*radius*N;Mhh+=w*mu*radius**2
                cell=min(int(x/p['L']*ncells),ncells-1)
                g[dofs,cell]+=w*p['D']*p['alpha_cte']/p['h']*Ndd
                if on:
                    j=zone[0];force[dofs,j]+=w*2*p['Ep']*p['bp']*p['d31']*zc*Ndd
                    probe[j,dofs]+=w*p['h']/(2*.5)*Ndd
                phi=s.modes(np.array([x]),s.NC,p)[0][0]
                projection[:,dofs]+=w*np.outer(phi,N);gram+=w*np.outer(phi,phi)
    free=np.arange(2,nf);Mc=M[np.ix_(free,free)];Kc=K[np.ix_(free,free)]
    # Shift-invert avoids loss of the lowest eigenvalues when the mesh's largest
    # eigenvalue grows as ell^-4. Generalized dense bisection was inadequate on
    # the 480/960-element meshes; those exploratory outputs are superseded.
    eig,V=eigsh(csc_matrix(Kc),k=nm,M=csc_matrix(Mc),sigma=0.,which='LM',tol=1e-11,v0=np.ones(len(free)))
    order=np.argsort(eig);eig=eig[order];V=V[:,order];assert np.all(eig>0);om=np.sqrt(eig)
    mass=np.zeros((nm+1,nm+1));mass[0,0]=Mhh;mass[0,1:]=mass[1:,0]=mh[free]@V;mass[1:,1:]=np.eye(nm)
    stiff=np.diag(np.r_[0,eig]);damp=np.diag(np.r_[0,2*p['zeta']*om])
    act=np.zeros((nm+1,4));act[0,0]=1;act[1:,1:]=V.T@force[free]
    gt=np.vstack((np.zeros(ncells),V.T@g[free]));G=gt.sum(1)
    mi=np.linalg.inv(mass);d=nm+1
    A=np.block([[np.zeros((d,d)),np.eye(d)],[-mi@stiff,-mi@damp]])
    B=np.vstack((np.zeros((d,4)),mi@act));E=np.r_[np.zeros(d),mi@G]
    rq=np.r_[0,np.linalg.solve(stiff[1:,1:],G[1:])]
    hs=np.zeros((s.ND,d));hs[0,0]=1;hs[1:,1:]=np.linalg.solve(gram,projection[:,free]@V)
    selector=np.block([[hs,np.zeros_like(hs)],[np.zeros_like(hs),hs]])
    tip=V[-2] # last displacement in free node order
    sensor=np.zeros((5,2*d));sensor[0,0]=1;sensor[1,d]=1;sensor[2:5,1:d]=probe[:,free]@V
    port=np.c_[np.zeros((4,d)),act.T]-np.c_[np.zeros((4,s.ND)),s.DES['F'].T]@selector
    plant={'M':mass,'K':stiff,'C':damp,'F':act,'G':G,'A':A,'B':B,'E':E,'rq':rq,'tip':tip,'gt':gt,'nm':nm,'edges':edges,'free':free,'V':V,'full_K':Kc,'full_g':g[free],'freq':om/(2*np.pi),'H':selector,'sensor':sensor}
    zz=expm(np.block([[A,B,E[:,None]],[np.zeros((5,2*d+5))]])*r.DT)
    stuff=(zz[:2*d,:2*d],zz[:2*d,2*d:2*d+4],zz[:2*d,-1],selector,port,np.c_[np.zeros((4,d)),act.T],sensor)
    return plant,stuff

def spatial_rhs(t,flat,ncells,p,gradient=.25,shadow_sweep=3.,conductivity=.5):
    T=flat.reshape(2,ncells);dx=p['L']/ncells;x=(np.arange(ncells)+.5)*dx
    chi=(1+np.tanh((t-20-shadow_sweep*x/p['L'])/p['w']))/2
    solar=(p['absorptivity']-p['eta'])*p['S']*np.cos(p['beta'])*chi*(1+gradient*(2*x/p['L']-1))
    lap=np.empty_like(T);lap[:,1:-1]=(T[:,:-2]-2*T[:,1:-1]+T[:,2:])/dx**2
    lap[:,0]=(T[:,1]-T[:,0])/dx**2;lap[:,-1]=(T[:,-2]-T[:,-1])/dx**2
    ca=p['mu']*p['c']/(2*p['width']);delta=T[0]-T[1]
    td=np.array([solar+p['qIR']-p['G']*delta,p['qIR']+p['G']*delta])-p['emissivity']*p['sigma']*(T**4-3**4)+conductivity*lap
    return (td/ca).ravel()

@njit(cache=True)
def spatial_mechanics(kind,dt,steps,Ad,Bd,Eg,Tgrid,delta_grid,static_grid,tip,H,rq,
    A,B,K,P,Q,W,R,limits,portgain,w0,b0,w1,b1,w2,b2,yscale,sc,pnom,F1,F2,refgrid,fgrid,fieldaware,residualmap):
    d=len(tip)+1;y=np.zeros(2*d);y[:d]=static_grid[0]
    tr=np.zeros((steps//64+1,12));nr=0;cost=0.;ss=0.;sa=0.;active=0.;infeasible=0.;peak=np.zeros(4)
    for j in range(steps+1):
        t=j*dt;i=min(int(t/.01),len(Tgrid)-2);v=(t-i*.01)/.01
        T=Tgrid[i]*(1-v)+Tgrid[i+1]*v;delta=delta_grid[i]*(1-v)+delta_grid[i+1]*v
        stat=static_grid[i]*(1-v)+static_grid[i+1]*v
        td,tdd,dc=thermal(T,t-20,1.,.3,1.,pnom);diff=T[0]-T[1];dd=td[0]-td[1]
        reference=np.concatenate((rq*diff,rq*dd));f=F1*dd+F2*(tdd[0]-tdd[1])
        fmean=f.copy()
        if fieldaware:
            reference=refgrid[i]*(1-v)+refgrid[i+1]*v
            f=fgrid[i]*(1-v)+fgrid[i+1]*v
        x=H@y-reference;base=-K@x
        u=base.copy();aux=-portgain@y
        if kind==4:
            xx=features(T,t-20,1.,.3,1.)
            pred=(np.tanh(np.tanh(xx@w0+b0)@w1+b1)@w2+b2)*yscale
            activity=np.sum(td*td)+(.3*dc)**2;u+=activity/(activity+1e-8)*pred
            if fieldaware:u+=residualmap@(f-fmean)
        u+=aux;act=0.;slack=0.
        if kind==4:
            a=2*B.T@P@x;bound=-.25*(x@W@x)+2*abs(x@P@(f+B@aux))-2*x@P@(A@x+f)
            u,act,slack=projection(u,a,bound,base+aux,limits,np.diag(R).copy())
        else:u=np.minimum(np.maximum(u,-limits),limits)
        vib=tip@(y[:d]-stat)[1:];running=x@Q@x+u@R@u
        ss+=vib*vib;sa+=y[0]**2;active+=act;infeasible+=slack>1e-8;peak=np.maximum(peak,np.abs(u))
        if j%64==0:
            tr[nr,0]=t;tr[nr,1:3]=T;tr[nr,3]=vib;tr[nr,4]=y[0];tr[nr,5:9]=u;tr[nr,9]=running;tr[nr,10]=act;tr[nr,11]=slack;nr+=1
        if j<steps:
            cost+=dt*running
            tn=t+dt;ii=min(int(tn/.01),len(Tgrid)-2);vv=(tn-ii*.01)/.01
            dn=delta_grid[ii]*(1-vv)+delta_grid[ii+1]*vv
            y=Ad@y+Bd@u+Eg@((delta+dn)/2)
    return tr[:nr],np.array([cost,np.sqrt(ss/(steps+1)),np.sqrt(sa/(steps+1))*206264.806,active/(steps+1),infeasible]),peak

def main(spatial_only=False):
    start=time.perf_counter();out={'meshes':[],'definition':'Cubic Hermite consistent FE, clamped patched eigenvectors and modal damping 2*zeta*omega; independent of analytic cantilever trial basis.'}
    if spatial_only:out=json.loads((s.DATA/'finite_element_report.json').read_text())
    for ne in ([] if spatial_only else [60,120,240]):
        plant,stuff=assemble(ne);base,_=assemble(ne,6,patches=False)
        exact=s.BETAS[:6]**2/s.PARAM['L']**2*np.sqrt(s.PARAM['D']/s.PARAM['mu'])/(2*np.pi)
        static_full=np.linalg.solve(plant['full_K'],plant['full_g'].sum(1))[-2]
        static_exact=s.model(41)['tip']@s.model(41)['rq'][1:]
        case={'requested_elements':ne,'actual_elements':len(plant['edges'])-1,'bare_frequency_relative_errors':((base['freq']-exact)/exact).tolist(),'static_tip_per_K_m':float(static_full),'static_vs_41Galerkin_relative':float(abs(static_full/static_exact-1)),'dynamic':{}}
        for label,kind,filt in [('LQR',2,False),('FilteredNN',4,True)]:
            met,tr=r.run(kind,filter=filt,end=120.,prepared=(plant,stuff));case['dynamic'][label]=met
            np.savez_compressed(s.DATA/f'FE_{ne}_{label}.npz',trace=tr)
        out['meshes'].append(case);r.write('finite_element_report.json',out);print('FE mesh',ne,case['dynamic'],flush=True)
    p=s.PARAM.copy();T0=s.initial_temperature(1,p,.3,1)[0];spatial={}
    for nc in [12,24,48]:
        grid=np.arange(0,120.0001,.01);sol=solve_ivp(lambda t,y:spatial_rhs(t,y,nc,p),[0,120],np.full(2*nc,T0),rtol=2e-10,atol=2e-9,t_eval=grid)
        assert sol.success
        TT=sol.y.reshape(2,nc,-1).transpose(2,0,1)
        # Conduction sums cancel under insulated end fluxes; face transfer cancels too.
        rng=np.random.default_rng(4403);Tc=rng.uniform(250,350,(2,nc));rhs=spatial_rhs(13,Tc.ravel(),nc,p).reshape(2,nc)
        dx=p['L']/nc;x=(np.arange(nc)+.5)*dx;chi=(1+np.tanh((13-20-3*x/p['L'])/3))/2
        solar=(p['absorptivity']-.3)*p['S']*np.cos(p['beta'])*chi*(1+.25*(2*x/p['L']-1))
        expected=np.sum(solar+2*p['qIR']-p['emissivity']*p['sigma']*(Tc[0]**4+Tc[1]**4-2*3**4))*dx*p['width']
        energy_rate=np.sum(rhs)*p['mu']*p['c']/2*dx
        spatial[nc]={'mean_diff_final':float(np.mean(TT[-1,0]-TT[-1,1])),'span_diff_range_final':np.ptp(TT[-1,0]-TT[-1,1]).item(),'energy_balance_error_W':float(energy_rate-expected),'max_T_K':float(TT.max()),'min_T_K':float(TT.min()),'dynamic':{}}
        np.savez_compressed(s.DATA/f'spatial_thermal_{nc}.npz',time=grid,temperature=TT)
        plant,stuff=assemble(240,ncells=nc);d=38;gt=plant['gt'];E=np.vstack((np.zeros((d,nc)),np.linalg.solve(plant['M'],gt)))
        zz=expm(np.block([[plant['A'],plant['B'],E],[np.zeros((4+nc,2*d+4+nc))]])*r.DT)
        delt=TT[:,0]-TT[:,1];static=np.zeros((len(grid),d));static[:,1:]=delt@np.linalg.solve(plant['K'][1:,1:],gt[1:]).T
        # A spatial reference needs measured field modes, not just two mean temperatures.
        # Analytical heat rates provide reference velocity; central differences of
        # these smooth rates provide acceleration on the validated 0.01 s grid.
        dtd=np.array([spatial_rhs(t,tt.ravel(),nc,p).reshape(2,nc) for t,tt in zip(grid,TT)])
        ddelta=dtd[:,0]-dtd[:,1];acc=np.gradient(ddelta,.01,axis=0,edge_order=2)
        localrq=np.zeros((d,nc));localrq[1:]=np.linalg.solve(plant['K'][1:,1:],gt[1:]);Hq=plant['H'][:s.ND,:d]
        rfield=localrq@delt.T;rfielddot=localrq@ddelta.T;rfieldddot=localrq@acc.T
        refgrid=np.c_[(Hq@rfield).T,(Hq@rfielddot).T]
        fgrid=np.c_[np.zeros((len(grid),s.ND)),(-np.linalg.solve(s.DES['M'],s.DES['C']@(Hq@rfielddot))-(Hq@rfieldddot)).T]
        residualmap=-np.linalg.solve(s.R,s.DES['B'].T@np.linalg.solve(s.ACL.T,expm(s.ACL.T*24)-np.eye(s.NX))@s.P)
        for label,kind,fieldaware in [('LQR',2,False),('FilteredNN',4,False),('FieldLQR',2,True),('FieldNN',4,True)]:
            tr,met,peak=spatial_mechanics(kind,r.DT,int(120/r.DT),np.ascontiguousarray(zz[:2*d,:2*d]),np.ascontiguousarray(zz[:2*d,2*d:2*d+4]),np.ascontiguousarray(zz[:2*d,2*d+4:]),TT.mean(2),delt,static,plant['tip'],plant['H'],s.DES['rq'],s.DES['A'],s.DES['B'],s.GAIN,s.P,s.Q,s.W,s.R,s.LIMITS,np.diag([0,5000,5000,5000])@stuff[4],*r.NET['weights'],r.NET['yscale'],s.sc,r.thermal_parameters(p),s.DES['F1'],s.DES['F2'],refgrid,fgrid,fieldaware,residualmap)
            spatial[nc]['dynamic'][label]=dict(zip(['cost','tip_rms_m','angle_rms_arcsec','filter_fraction','infeasible_count'],map(float,met)))
            np.savez_compressed(s.DATA/f'spatial_FE_{nc}_{label}.npz',trace=tr,peak=peak)
        print('spatial',nc,spatial[nc],flush=True)
    out['spatial_thermal']=spatial;out['elapsed_s']=time.perf_counter()-start;r.write('finite_element_report.json',out)

if __name__=='__main__':
    import sys
    main('--spatial-only' in sys.argv)
