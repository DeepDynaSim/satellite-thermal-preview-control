"""Publication figures drawn exclusively from saved computed results."""
import simulate as s
import numpy as np,json,time
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,FancyArrowPatch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':10,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,'lines.linewidth':1.25})
COL={'LQR':'#D68B13','Oracle24':'#24895C','RawNN':'#8650A5','FilteredNN':'#1867B4','Polynomial':'#B04A46','HubPD':'#73777B','Passivity':'#254A57','Frozen24':'#AB6720'}
def save(fig,name):
    for ext in ['png','pdf','svg']:fig.savefig(s.FIG/(name+'.'+ext),dpi=420,bbox_inches='tight')
    plt.close(fig)
def trace(category,label):return np.load(s.DATA/f'{category}_{label}.npz')['trace']
def load(name):return json.loads((s.DATA/name).read_text())

def geometry():
    fig=plt.figure(figsize=(7.4,4.8));ax=fig.add_subplot(211,projection='3d',computed_zorder=False)
    x=np.linspace(0,9,41);yy=np.linspace(-1.5,1.5,9);X,Y=np.meshgrid(x,yy);Z=.002*(X/9)**2*25
    ax.plot_surface(X,Y,Z,cmap='Blues',alpha=1.,linewidth=.3,edgecolor='#C8DDF2',antialiased=True,zorder=2)
    verts=[[(0,-.7,-.7),(0,.7,-.7),(0,.7,.7),(0,-.7,.7)],[(-1.4,-.7,-.7),(-1.4,.7,-.7),(-1.4,.7,.7),(-1.4,-.7,.7)],[(0,-.7,.7),(0,.7,.7),(-1.4,.7,.7),(-1.4,-.7,.7)],[(0,-.7,-.7),(0,-.7,.7),(-1.4,-.7,.7),(-1.4,-.7,-.7)]]
    ax.add_collection3d(Poly3DCollection(verts,facecolor='#D3D7DC',edgecolor='#4B5665',linewidth=.8,zorder=4))
    for a in [0,1.8,4.5]:
        za=.05*(a/9)**2+.02;zb=.05*((a+.5)/9)**2+.02
        ax.add_collection3d(Poly3DCollection([[(a,-.15,za),(a+.5,-.15,zb),(a+.5,.15,zb),(a,.15,za)]],facecolor='#D99C29',edgecolor='#785512',zorder=5))
    for xx in [2,4,6,8]:ax.quiver(xx,.9,1.7,0,0,-1.2,color='#C78419',arrow_length_ratio=.15,linewidth=1.4)
    ax.text2D(.52,.88,'Solar heating',transform=ax.transAxes,fontsize=9)
    ax.text2D(.12,.76,'Hub',transform=ax.transAxes,fontsize=9,zorder=10)
    ax.text2D(.57,.13,'9 m × 3 m appendage',transform=ax.transAxes,fontsize=9)
    ax.text2D(.24,.22,'Gold: PZT pairs',transform=ax.transAxes,fontsize=8)
    ax.set_xlim(-1.5,9);ax.set_ylim(-1.8,1.8);ax.set_zlim(-.8,2);ax.set_box_aspect((10.5,3.6,2.8),zoom=1.0);ax.view_init(20,-65);ax.set_axis_off()
    ax2=fig.add_subplot(212);ax2.set_xlim(0,10);ax2.set_ylim(0,3);ax2.axis('off')
    boxes=[(.05,.95,1.8,'Temperature\nforecast'),(2.35,.95,1.8,'Moving thermal\nreference'),(4.65,.95,2.0,'LQR and learned\npreview'),(7.15,.95,2.6,'Command filter\nand input bounds')]
    for a,b,c,label in boxes:
        ax2.add_patch(Rectangle((a,b),c,1.0,facecolor='#F7F8FA',edgecolor='#48596A'));ax2.text(a+c/2,b+.5,label,ha='center',va='center')
    for a,b in [(1.85,2.35),(4.15,4.65),(6.65,7.15)]:ax2.annotate('',xy=(b,1.45),xytext=(a,1.45),arrowprops={'arrowstyle':'->'})
    ax2.plot([8.45,8.45,5.65,5.65],[.95,.5,.5,.95],color='#48596A');ax2.text(6.9,.15,'Observer and residual strain-rate feedback',ha='center',fontsize=8)
    ax2.text(.1,2.5,'Benchmark schematic; deflection exaggerated for visibility',fontsize=8)
    fig.subplots_adjust(hspace=-.12);save(fig,'fig01_architecture')

def main(wait=True):
    if wait:
        while True:
            try:
                if len(load('timing_electrical_report.json')['electrical'])==4:break
            except (FileNotFoundError,json.JSONDecodeError):pass
            time.sleep(5)
    additional=load('additional_report.json');lhs=load('lhs_report.json');fe=load('finite_element_report.json');ref=load('refinement_hardware_report.json');cpu=load('timing_electrical_report.json')
    geometry()
    fig,ax=plt.subplots(2,1,figsize=(7,4.4),sharex=True)
    for suffix,ls in [('exit','-'),('entry','--')]:
        tr=trace('nominal','FilteredNN_'+suffix)
        for k,label in [(1,'front'),(2,'back')]:ax[0].plot(tr[:,0],tr[:,k],ls,label=suffix.title()+' '+label)
        ax[1].plot(tr[:,0],tr[:,3],ls,label=suffix.title())
    ax[0].set_ylabel('Temperature (K)');ax[1].set_ylabel('Front minus back (K)');ax[1].set_xlabel('Time (s)');ax[0].legend(ncol=2);ax[1].legend();save(fig,'fig02_thermal')
    for suffix,num in [('exit',3),('entry',4)]:
        fig,ax=plt.subplots(3,1,figsize=(7,5.5),sharex=True)
        for label in ['LQR','Oracle24','FilteredNN','Frozen24']:
            tr=trace('nominal',label+'_'+suffix);c=COL[label]
            ax[0].plot(tr[:,0],tr[:,6]*1e6,label=label,color=c);ax[1].plot(tr[:,0],tr[:,7]*206264.806,color=c);ax[2].plot(tr[:,0],tr[:,10],color=c)
        ax[0].legend(ncol=4);ax[0].set_ylabel('Dynamic tip (µm)');ax[1].set_ylabel('Hub angle (arcsec)');ax[2].set_ylabel('Patch 1 voltage (V)');ax[2].set_xlabel('Time (s)');ax[2].set_xlim(0,90);save(fig,f'fig0{num}_{suffix}_response')
    a=np.load(s.DATA/'learning.npz');Y=a['Y'][a['ite']];yp=a['pred'];fig,ax=plt.subplots(2,2,figsize=(7,5.0))
    for k,v in enumerate(ax.flat):
        v.scatter(Y[:,k],yp[:,k],s=7,alpha=.38,color='#1867B4',rasterized=True);low=min(Y[:,k].min(),yp[:,k].min());high=max(Y[:,k].max(),yp[:,k].max());v.plot([low,high],[low,high],color='#2F3D49',lw=.8)
        unit='N m' if k==0 else 'V';v.set_xlabel(f'Analytic label ({unit})');v.set_ylabel(f'NN prediction ({unit})');v.set_title('Hub torque' if k==0 else f'Patch {k}')
    fig.tight_layout();save(fig,'fig05_learning')
    fig,ax=plt.subplots(1,3,figsize=(7.4,3.2))
    for v,key,title in zip(ax,['cost','tip_rms_m','angle_rms_arcsec'],['Cost ratio','Tip RMS ratio','Attitude RMS ratio']):
        base=np.array([case['LQR'][key] for case in lhs['results']]);target=np.array([case['FilteredNN'][key] for case in lhs['results']]);ratio=target/base
        v.hist(ratio,bins=24,color='#1867B4',alpha=.82);v.axvline(1,color='#B04A46',ls='--',lw=1);v.set_xlabel('Filtered NN / LQR');v.set_title(title)
        st=lhs['paired_summary']['FilteredNN'][key];ci=st['median_ratio_bootstrap95'];v.text(.03,.97,f'Median {st["median_ratio"]:.3f}\n95% CI [{ci[0]:.3f}, {ci[1]:.3f}]',transform=v.transAxes,va='top',fontsize=7)
    ax[0].set_ylabel('Paired cases');fig.tight_layout();save(fig,'fig06_robustness')
    tr=trace('nominal','FilteredNN_exit');fig,ax=plt.subplots(2,1,figsize=(7,4.3),sharex=True)
    ax[0].plot(tr[:,0],tr[:,4]*1e3,label='Quasistatic reference');ax[0].plot(tr[:,0],tr[:,5]*1e3,'--',label='Total tip');ax[0].set_ylabel('Tip deflection (mm)');ax[0].legend()
    ax[1].plot(tr[:,0],tr[:,16]);ax[1].axhline(0,color='black',lw=.7);ax[1].set_ylabel('Nominal inequality residual');ax[1].set_xlabel('Time (s)');save(fig,'fig07_equilibrium_certificate')
    fig,ax=plt.subplots(2,2,figsize=(7,5.2),sharex=True)
    for suffix,ls in [('exit','-o'),('entry','--s')]:
        cases=[additional['horizons'][f'H{H}_{suffix}'] for H in [0,6,12,24,36,48]]
        for v,key,factor,ylabel in zip(ax.flat,['cost','tip_rms_m','angle_rms_arcsec','effort'],[1,1e6,1,1],['Quadratic cost','Dynamic tip RMS (µm)','Angle RMS (arcsec)','Weighted command mean square']):
            yy=[(np.sum(np.array(c['channels'])[1]**2/s.LIMITS**2) if key=='effort' else c[key]*factor) for c in cases];v.plot([0,6,12,24,36,48],yy,ls,label=suffix.title());v.set_ylabel(ylabel);v.grid(alpha=.15)
    ax[0,0].legend();ax[1,0].set_xlabel('Preview horizon (s)');ax[1,1].set_xlabel('Preview horizon (s)');fig.tight_layout();save(fig,'fig08_horizon')
    meshes=fe['meshes']+ref['fine_FE'];fig,ax=plt.subplots(1,3,figsize=(7.4,3.2))
    for lab in ['LQR','FilteredNN']:
        ax[0].semilogx([v['actual_elements'] for v in meshes],[v['dynamic'][lab]['tip_rms_m']*1e6 for v in meshes],'-o',label=lab,color=COL[lab])
    ax[0].set_xlabel('Hermite elements');ax[0].set_ylabel('Tip RMS (µm)');ax[0].legend(fontsize=7)
    ax[1].semilogy([v['n'] for v in ref['modal']],[v['metrics']['tip_rms_m']*1e6 for v in ref['modal']],'-o',color=COL['FilteredNN'],label='6400 updates/s')
    for z in ref['modal65_time_refinement']:ax[1].scatter(65,z['metrics']['tip_rms_m']*1e6,marker='^',color='#24895C')
    ax[1].text(.03,.97,'65 modes at 6400/s:\nlinear sampled instability',transform=ax[1].transAxes,va='top',fontsize=7);ax[1].set_xlabel('Galerkin elastic modes');ax[1].set_ylabel('Filtered NN tip RMS (µm)')
    st=np.load(s.DATA/'spatial_thermal_48.npz');d=st['temperature'][:,0]-st['temperature'][:,1];im=ax[2].imshow(d.T,aspect='auto',origin='lower',extent=[0,120,0,9],cmap='inferno');ax[2].set_xlabel('Time (s)');ax[2].set_ylabel('Span coordinate (m)');fig.colorbar(im,ax=ax[2],label='Front minus back (K)',fraction=.045)
    fig.tight_layout();save(fig,'fig09_spatial_convergence')
    obs=load('observer_refinement_report.json');fig,ax=plt.subplots(1,3,figsize=(7.5,3.5))
    for label in ['LQR','FilteredNN']:
        keys=[f'{label}_core{k}' for k in [0,3,6,13]];v=[additional['delays'][k] for k in keys];ax[0].semilogy(np.array([0,3,6,13])*s.EVAL_DT*1000,[c['cost'] for c in v],'-o',color=COL[label],label=label)
    ax[0].set_xlabel('Core command delay (ms)');ax[0].set_ylabel('Quadratic cost');ax[0].legend()
    vv=[additional['nominal']['FilteredNN_exit']['cost']]+[additional['delays'][f'residual{i}']['cost'] for i in [1,2,3]]
    ax[1].semilogy(np.array([0,1,2,3])*s.EVAL_DT*1000,vv,'-o',color=COL['FilteredNN']);ax[1].set_xlabel('Residual delay (ms)');ax[1].set_ylabel('Quadratic cost')
    for j,case in enumerate(['low','high','calibrated']):
        values=[obs['selected_cases'][label+'_'+case]['cost'] for label in ['LQR','FilteredNN']]
        ax[2].bar(np.array([0,1])+j*3,values,color=[COL['LQR'],COL['FilteredNN']])
    ax[2].set_xticks([.5,3.5,6.5],['Low\nbias','High\nbias','Zero bias\n5 mK']);ax[2].set_ylabel('Stable observer cost');ax[2].set_yscale('log');fig.tight_layout();save(fig,'fig10_sensing_delay')
    fig,ax=plt.subplots(1,2,figsize=(7,3.3));names=['CachedOracle24','NN','CubicPolynomial','RBF256'];med=[cpu['timings'][n]['median_per_call_us'] for n in names];p99=[cpu['timings'][n]['p99_per_call_us'] for n in names]
    ax[0].bar(range(4),med,color=['#24895C','#1867B4','#B04A46','#254A57']);ax[0].scatter(range(4),p99,color='black',marker='_',label='99th percentile');ax[0].set_xticks(range(4),['Cached\npreview','Neural','Cubic\npolynomial','RBF']);ax[0].set_ylabel('Warm call time (µs)');ax[0].legend(fontsize=7)
    labels=['NoGate','Gate','GateFilter','NoResidual'];values=[additional['ablations'][v]['cost'] for v in labels];ax[1].bar(range(4),values,color='#1867B4');ax[1].set_xticks(range(4),['NN\nno gate','NN\ngate','NN\nfilter','NN\nno RMD']);ax[1].set_ylabel('Quadratic cost');ax[1].set_yscale('log');fig.tight_layout();save(fig,'fig11_timing_ablation')

if __name__=='__main__':main()
