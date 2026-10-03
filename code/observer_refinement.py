"""Retain failed observer design and test low-bandwidth alternative without cherry picking."""
import revision_studies as r
import simulate as s
import numpy as np,json,time

plant,stuff=r.prepare_plant(s.PARAM);Ad,Bd,Ed,H,port,fullport,Cplant=stuff
port=np.diag([0,5000,5000,5000])@port
def spectral(obs):
    Ao,Bo,Eo,L,C=obs;ny=76;no=12
    F=np.block([[Ad-Bd@port,-Bd@s.GAIN],[-Bo@port,Ao-Bo@s.GAIN]])
    correction=np.block([[np.eye(ny),np.zeros((ny,no))],[L@Cplant,np.eye(no)-L@C]])
    transition=np.eye(ny+no)
    # 32 native updates contain exactly five 1 ms measurement updates.
    for j in range(32):
        t=j*r.DT
        if j==0 or int((t-r.DT)/.001)!=int(t/.001):transition=correction@transition
        transition=F@transition
    return float(np.max(abs(np.linalg.eigvals(transition)))),float(np.log(np.max(abs(np.linalg.eigvals(transition))))/(32*r.DT))

out={'spectrum':[],'selected_cases':{},'selection_rule':'Choose the smallest tested process covariance with spectral radius below one; disclose all candidate spectra and original failed cases.'}
designs={}
for q in [1e-6,1e-8,1e-10,1e-12]:
    obs=r.observer_design(process_cov=q);rho,rate=spectral(obs);out['spectrum'].append({'normalized_process_covariance':q,'floquet_radius_5ms':rho,'equivalent_spectral_rate_s-1':rate});designs[q]=obs
stable=[v['normalized_process_covariance'] for v in out['spectrum'] if v['floquet_radius_5ms']<1]
out['selected_process_covariance']=min(stable) if stable else None
r.write('observer_refinement_report.json',out)
if stable:
    selected=designs[min(stable)]
    for level in ['low','high']:
        angle=.02 if level=='low' else .1;strain=5e-9 if level=='low' else 25e-9
        noise=[.1 if level=='low' else .5,-.1 if level=='low' else -.5,.05 if level=='low' else .2,.05 if level=='low' else .2,angle/206264.806,-angle/206264.806,strain,-strain,strain,angle/206264.806,angle/206264.806,strain,strain,strain]
        for label,kind,filt in [('LQR',2,False),('FilteredNN',4,True)]:
            met,tr=r.run(kind,filter=filt,observer=True,noise=noise,prepared=(plant,stuff),obsdesign=selected,normal_columns=10)
            out['selected_cases'][label+'_'+level]=met;r.write('observer_refinement_report.json',out)
            np.savez_compressed(s.DATA/f'observer_refined_{label}_{level}.npz',trace=tr)
            print(label,level,met['cost'],flush=True)
out['complete']=True;r.write('observer_refinement_report.json',out)
