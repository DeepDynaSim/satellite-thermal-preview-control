import revision_studies as r
import simulate as s
import numpy as np,json
out=json.loads((s.DATA/'observer_refinement_report.json').read_text())
design=r.observer_design(process_cov=out['selected_process_covariance'])
noise=np.zeros(14);noise[2:4]=.005;noise[9:11]=.02/206264.806;noise[11:14]=5e-9
for label,kind,filt in [('LQR',2,False),('FilteredNN',4,True)]:
    met,tr=r.run(kind,noise=noise,observer=True,obsdesign=design,filter=filt,normal_columns=10)
    out['selected_cases'][label+'_calibrated']=met
    np.savez_compressed(s.DATA/f'observer_refined_{label}_calibrated.npz',trace=tr)
out['calibrated_case_assumptions']={'all_biases':0.,'temperature_noise_std_K':.005,'angle_noise_arcsec':.02,'gyro_noise_arcsec_per_s':.02,'strain_noise_std':5e-9}
r.write('observer_refinement_report.json',out)
