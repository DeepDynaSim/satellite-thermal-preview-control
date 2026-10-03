import revision_studies as r
import simulate as s
import numpy as np,json,time

while True:
    try:
        report=json.loads((s.DATA/'refinement_hardware_report.json').read_text())
        if 'FE_conservation' in report:break
    except (FileNotFoundError,json.JSONDecodeError):pass
    time.sleep(5)
report['modal65_time_refinement']=[]
for dt in [1/12800,1/25600]:
    plant,stuff=r.prepare_plant(s.PARAM,65,dt=dt);ad,bd,ed,H,port,_,_=stuff
    rho=np.max(np.abs(np.linalg.eigvals(ad-bd@(s.GAIN@H+np.diag([0,5000,5000,5000])@port))))
    met,tr=r.run(4,n=65,prepared=(plant,stuff),dt=dt,stride=int(.1/dt))
    report['modal65_time_refinement'].append({'dt':dt,'sampled_linear_radius':float(rho),'metrics':met})
    r.write('refinement_hardware_report.json',report)
    np.savez_compressed(s.DATA/f'modal65_dt{int(1/dt)}.npz',trace=tr);print('modal65 finer',dt,rho,met['cost'],flush=True)
