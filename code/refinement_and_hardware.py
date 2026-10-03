import revision_studies as r
import simulate as s
from finite_element_validation import assemble
import numpy as np,json,time
from scipy.linalg import expm

def main():
    report={'modal':[],'fine_FE':[],'algebra':{},'hardware':{}}
    for n in [25,37,41,49,65]:
        plant,stuff=r.prepare_plant(s.PARAM,n);ad,bd,ed,H,port,_,_=stuff
        fb=s.GAIN@H+np.diag([0,5000,5000,5000])@port
        rho=np.max(np.abs(np.linalg.eigvals(ad-bd@fb)))
        met,tr=r.run(4,n=n,prepared=(plant,stuff));report['modal'].append({'n':n,'sampled_linear_radius':float(rho),'metrics':met})
        np.savez_compressed(s.DATA/f'modal_{n}.npz',trace=tr);r.write('refinement_hardware_report.json',report);print('modal',n,rho,met['cost'],met['failure'],flush=True)
    Dadd=2*s.PARAM['Ep']*s.PARAM['bp']*(s.PARAM['hp']**3/12+s.PARAM['hp']*((s.PARAM['h']+s.PARAM['hp'])/2)**2)
    edges=[0,.5,1.8,2.3,4.5,5,9.];tip_exact=0.
    for a,b in zip(edges[:-1],edges[1:]):
        D=s.PARAM['D']+(Dadd if (a,b) in [(0,.5),(1.8,2.3),(4.5,5.)] else 0.)
        tip_exact+=s.PARAM['D']*s.PARAM['alpha_cte']/s.PARAM['h']/D*(9*(b-a)-(b*b-a*a)/2)
    for ne in [480,960]:
        plant,stuff=assemble(ne)
        static=np.linalg.solve(plant['full_K'],plant['full_g'].sum(1))[-2]
        case={'requested_elements':ne,'actual_elements':len(plant['edges'])-1,'static_exact_tip_per_K_m':tip_exact,'static_tip_per_K_m':float(static),'static_relative_error':float(abs(static/tip_exact-1)),'dynamic':{}}
        for label,kind,filt in [('LQR',2,False),('FilteredNN',4,True)]:
            met,tr=r.run(kind,filter=filt,prepared=(plant,stuff));case['dynamic'][label]=met;np.savez_compressed(s.DATA/f'FE_{ne}_{label}.npz',trace=tr)
        report['fine_FE'].append(case);r.write('refinement_hardware_report.json',report);print('fineFE',ne,case['dynamic']['FilteredNN']['tip_rms_m'],flush=True)
    report['modal65_time_refinement']=[]
    for dt in [1/12800,1/25600]:
        plant,stuff=r.prepare_plant(s.PARAM,65,dt=dt);ad,bd,ed,H,port,_,_=stuff
        rho=np.max(np.abs(np.linalg.eigvals(ad-bd@(s.GAIN@H+np.diag([0,5000,5000,5000])@port))))
        met,tr=r.run(4,n=65,prepared=(plant,stuff),dt=dt,stride=int(.1/dt))
        report['modal65_time_refinement'].append({'dt':dt,'sampled_linear_radius':float(rho),'metrics':met});r.write('refinement_hardware_report.json',report)
        np.savez_compressed(s.DATA/f'modal65_dt{int(1/dt)}.npz',trace=tr);print('modal65 finer',dt,rho,met['cost'],flush=True)
    # Independent complete-command versus shifted-core half-space algebra.
    rng=np.random.default_rng(337713);error=[]
    for _ in range(10000):
        x=rng.normal(size=s.NX)*s.sc;f=rng.normal(size=s.NX)*s.sc;uD=rng.normal(size=4)*s.LIMITS*.03;core=rng.normal(size=4)*s.LIMITS*.1
        a=2*s.DES['B'].T@s.P@x;fe=f+s.DES['B']@uD
        bound_full=-.25*x@s.W@x+2*abs(x@s.P@fe)-2*x@s.P@(s.DES['A']@x+f)
        bound_core=-.25*x@s.W@x+2*abs(x@s.P@fe)-2*x@s.P@(s.DES['A']@x+fe)
        error.append(abs((a@(core+uD)-bound_full)-(a@core-bound_core)))
    report['algebra']={'cases':10000,'max_residual_difference':float(max(error))}
    # Conservation of isolated independent FE plant.
    plant,stuff=assemble(120);rng=np.random.default_rng(665);y=rng.normal(size=76)*np.r_[np.full(38,1e-6),np.full(38,1e-5)]
    h0=plant['M'][0]@y[38:];energy=lambda yy:.5*(yy[38:]@plant['M']@yy[38:]+yy[:38]@plant['K']@yy[:38])
    en=[energy(y)];mom=[h0];ad=expm(plant['A']*.001)
    for j in range(10000):y=ad@y;en.append(energy(y));mom.append(plant['M'][0]@y[38:])
    report['FE_conservation']={'maximum_energy_increase_J':float(np.max(np.diff(en))),'momentum_max_abs_drift':float(np.max(np.abs(np.array(mom)-h0)))}
    r.write('refinement_hardware_report.json',report)

if __name__=='__main__':main()
