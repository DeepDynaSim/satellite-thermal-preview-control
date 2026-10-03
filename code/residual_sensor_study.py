import revision_studies as r
import simulate as s
import json,time,numpy as np

while True:
    try:
        report=json.loads((s.DATA/'additional_report.json').read_text())
        if 'elapsed_s' in report:break
    except (FileNotFoundError,json.JSONDecodeError):pass
    time.sleep(5)
report['residual_sensor_noise']={}
for level,sigma in [('low',1e-6),('high',1e-5)]:
    noise=np.zeros(20);noise[14:17]=5000*sigma;noise[17:20]=5000*sigma/2*np.array([1,-1,1])
    for label,kind,filt in [('LQR',2,False),('FilteredNN',4,True)]:
        met,tr=r.run(kind,noise=noise,filter=filt,seed=99004)
        report['residual_sensor_noise'][label+'_'+level]={'assumed_port_noise_std_Nm_per_Vs':sigma,'noise_hold_s':.001,'metrics':met}
        np.savez_compressed(s.DATA/f'residualnoise_{label}_{level}.npz',trace=tr)
        r.write('additional_report.json',report)
