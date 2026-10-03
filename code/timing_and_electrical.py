"""Warmed CPU timings and illustrative, explicitly assumed electronics accounting."""
import revision_studies as r
import simulate as s
from revision_kernel import features,forecast_one,polynomial
from fast_kernel import thermal
from numba import njit
import numpy as np,time,json,os,platform,winreg

@njit(cache=True)
def surrogate(T,tau,dr,eta,flux,p,w0,b0,w1,b1,w2,b2,ys):
    xx=features(T,tau,dr,eta,flux);td,tdd,dc=thermal(T,tau,dr,eta,flux,p)
    a=np.sum(td*td)+(.3*dc)**2
    return a/(a+1e-8)*(np.tanh(np.tanh(xx@w0+b0)@w1+b1)@w2+b2)*ys

@njit(cache=True)
def poly_surrogate(T,tau,dr,eta,flux,p,coef,powers):
    xx=features(T,tau,dr,eta,flux);td,tdd,dc=thermal(T,tau,dr,eta,flux,p)
    a=np.sum(td*td)+(.3*dc)**2
    return a/(a+1e-8)*polynomial(xx,coef,powers)

@njit(cache=True)
def rbf_surrogate(T,tau,dr,eta,flux,p,centers,gamma,beta):
    xx=features(T,tau,dr,eta,flux);td,tdd,dc=thermal(T,tau,dr,eta,flux,p)
    a=np.sum(td*td)+(.3*dc)**2;basis=np.concatenate((np.ones(1),xx,np.exp(-gamma*np.sum((centers-xx)**2,axis=1))))
    return a/(a+1e-8)*(basis@beta)

def main(wait=True):
    if wait:
        while True:
            try:
                a=json.loads((s.DATA/'additional_report.json').read_text());fe=json.loads((s.DATA/'finite_element_report.json').read_text());ref=json.loads((s.DATA/'refinement_hardware_report.json').read_text())
                ready=(s.DATA/'lhs_report.json').exists() and (s.DATA/'rbf_report.json').exists() and 'elapsed_s' in a and len(a.get('residual_sensor_noise',{}))==4 and len(fe.get('spatial_thermal',{}))==3 and len(fe['spatial_thermal']['48']['dynamic'])==4 and len(ref.get('modal65_time_refinement',[]))==2 and 'FE_conservation' in ref
                if ready:break
            except (FileNotFoundError,json.JSONDecodeError):pass
            time.sleep(5)
    net=r.NET;p=r.thermal_parameters(s.PARAM);poly=np.load(s.DATA/'polynomial.npz');coef=poly['coef'];powers=poly['powers']
    rb=np.load(s.DATA/'rbf.npz');centers=rb['centers'];gamma=float(rb['gamma']);beta=rb['beta']
    rng=np.random.default_rng(76621);N=1000;mean=rng.uniform(218,370,N);diff=rng.uniform(-30,35,N)
    states=np.c_[mean+diff/2,mean-diff/2];tau=rng.uniform(-55,65,N);directions=rng.choice([-1.,1.],N);eta=rng.uniform(.22,.38,N);flux=rng.uniform(.85,1.15,N)
    funcs={'NN':lambda i:surrogate(states[i],tau[i],directions[i],eta[i],flux[i],p,*net['weights'],net['yscale']),
           'CachedOracle24':lambda i:forecast_one(states[i],tau[i],directions[i],eta[i],flux[i],p,r.KERNELS[24.],.12),
           'CubicPolynomial':lambda i:poly_surrogate(states[i],tau[i],directions[i],eta[i],flux[i],p,coef,powers),
           'RBF256':lambda i:rbf_surrogate(states[i],tau[i],directions[i],eta[i],flux[i],p,centers,gamma,beta)}
    for f in funcs.values():
        for i in range(200):f(i)
    timing={}
    for name,fun in funcs.items():
        individual=[];batch=[]
        for repeat in range(15):
            start=time.perf_counter_ns()
            for i in range(N):
                t=time.perf_counter_ns();fun(i);individual.append(time.perf_counter_ns()-t)
            batch.append((time.perf_counter_ns()-start)/N)
        timing[name]={'median_per_call_us':float(np.median(individual)/1000),'p99_per_call_us':float(np.quantile(individual,.99)/1000),'median_batch_amortized_us':float(np.median(batch)/1000),'repeats':15,'states_per_repeat':N}
    with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,r'HARDWARE\DESCRIPTION\System\CentralProcessor\0') as key:cpu=winreg.QueryValueEx(key,'ProcessorNameString')[0].strip()
    weights=net['weights'];memory={
        'NN_weights_and_scales_bytes':int(sum(x.nbytes for x in weights)+net['yscale'].nbytes),
        'NN_dense_multiply_accumulates':int(weights[0].size+weights[2].size+weights[4].size),'NN_tanh_calls':80,
        'Oracle24_cached_thermal_kernel_bytes':int(r.KERNELS[24.].nbytes),'Oracle24_forecast_RK4_stages':200*8,'Oracle24_thermal_acceleration_evaluations':200,
        'Polynomial_coefficients_bytes':int(coef.nbytes),'Polynomial_integer_exponent_bytes':int(powers.nbytes),'Polynomial_terms':len(powers),
        'RBF_centers_coefficients_gamma_bytes':int(centers.nbytes+beta.nbytes+8),'RBF_exponential_calls':len(centers)}
    report={'CPU':cpu,'platform':platform.platform(),'Python':platform.python_version(),'NumPy':np.__version__,'arithmetic':'float64, compiled Numba, fastmath=False, BLAS threads 1','cache':'compiled kernels and mechanical forecast exponentials warmed, no training/file IO in measured intervals','timings':timing,'storage_and_operations':memory,'electrical':{}}
    # Full native trace is required for derivative-based current/loss estimates.
    eps_r=1700.;eps0=8.8541878128e-12;C=eps_r*eps0*s.PARAM['bp']*.5/s.PARAM['hp'];Rdriver=100.;tan_delta=.02
    for label,kind,filt in [('LQR',2,False),('Oracle24',3,False),('RawNN',4,False),('FilteredNN',4,True)]:
        met,tr=r.run(kind,filter=filt,stride=1)
        voltage=tr[:,10:13];dv=np.gradient(voltage,r.DT,axis=0);current=C*dv
        # Two equal ceramics per pair, opposite voltage; resistive loss summed over six ceramics.
        Eres=2*Rdriver*np.trapezoid(np.sum(current**2,axis=1),tr[:,0])
        # Explicit narrow-band dielectric proxy; spectral Parseval energy on mean-removed finite record.
        v=voltage-voltage.mean(0);freq=np.fft.rfftfreq(len(v),r.DT);spec=np.fft.rfft(v,axis=0);omega=2*np.pi*freq
        weights=np.full(len(freq),2.);weights[0]=1.;weights[-1]=1. if len(v)%2==0 else 2.
        Ediel=2*C*tan_delta*r.DT/len(v)*np.sum(weights[:,None]*omega[:,None]*abs(spec)**2)
        H=np.r_[0,-np.cumsum((tr[1:,9]+tr[:-1,9])/2*r.DT)]
        report['electrical'][label]={'assumed_C_each_ceramic_F':C,'assumed_eps_r':eps_r,'assumed_driver_resistance_ohm':Rdriver,'assumed_tan_delta':tan_delta,'rms_current_each_patch_A':np.sqrt(np.mean(current**2,axis=0)).tolist(),'peak_current_A':np.max(abs(current),axis=0).tolist(),'six_ceramic_resistive_loss_J':float(Eres),'six_ceramic_spectral_dielectric_proxy_J':float(Ediel),'maximum_stored_electrical_J':float(np.max(C*np.sum(voltage**2,axis=1))),'wheel_peak_abs_momentum_change_Nms':float(np.max(abs(H))),'wheel_net_momentum_change_Nms':float(H[-1]),'duration_s':120}
        r.write('timing_electrical_report.json',report)
    r.write('timing_electrical_report.json',report);print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
