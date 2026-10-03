% Independent MATLAB validation of Python SciPy matrices and trajectories.
% Uses native MATLAB care and ode113 rather than importing Python solvers.
taskRoot=fileparts(fileparts(mfilename('fullpath')));
if isfolder(fullfile(taskRoot,'data'))
    dataPath=fullfile(taskRoot,'data');
else
    dataPath=fullfile(taskRoot,'SCI_Thermal_Neural_Control','data');
end
v=load(fullfile(dataPath,'independent_validation.mat'));
v.F1=v.F1(:);v.F2=v.F2(:);
Ds=diag(v.sc);As=Ds\v.A*Ds;Bs=Ds\v.B;
[Ps,~,Ks]=care(As,Bs,eye(length(v.sc)),v.R);
Pm=Ds\Ps/Ds;Km=Ks/Ds;
result=struct();
result.MATLAB_version=version;
result.P_relative_difference=norm(Pm-v.P,'fro')/norm(v.P,'fro');
result.K_relative_difference=norm(Km-v.GAIN,'fro')/norm(v.GAIN,'fro');
result.CARE_relative_residual=norm(v.A'*Pm+Pm*v.A-Pm*v.B*(v.R\v.B')*Pm+v.Q,'fro')/norm(v.Q,'fro');
options=odeset('RelTol',1e-11,'AbsTol',1e-13);
[tt,xx]=ode113(@(t,x)(v.A-v.B*Km)*x+v.F1*.2*sin(.4*t)+v.F2*.08*cos(.4*t),v.validation_t,v.validation_x0,options);
result.trajectory_max_scaled_difference=max(abs((xx-v.validation_x)./v.sc),[],'all');
result.trajectory_relative_L2_difference=norm(xx-v.validation_x,'fro')/norm(v.validation_x,'fro');
freq=v.beta_roots(1:5).^2/v.parameters_vector(1)^2*sqrt(v.parameters_vector(2)/v.parameters_vector(3))/(2*pi);
result.analytical_clamped_frequencies_Hz=freq;
% Verify the finite horizon HJB identity with independently evaluated terms.
rng(13579);identityError=zeros(100,1);
for j=1:100
    x=v.sc'.*randn(length(v.sc),1)*.2;s=randn(length(v.sc),1);f=randn(length(v.sc),1)*.001;
    sd=-(v.A-v.B*Km)'*s-Pm*f;
    cd=-2*s'*f+s'*v.B*(v.R\v.B')*s;
    us=-(v.R\v.B')*(Pm*x+s);
    residual=2*sd'*x+cd+x'*v.Q*x+us'*v.R*us+2*(Pm*x+s)'*(v.A*x+v.B*us+f);
    identityError(j)=abs(residual)/(1+abs(x'*v.Q*x)+abs(us'*v.R*us));
end
result.HJB_identity_max_relative_residual=max(identityError);
fid=fopen(fullfile(dataPath,'matlab_validation.json'),'w');fwrite(fid,jsonencode(result,PrettyPrint=true),'char');fclose(fid);
disp(result);
assert(result.P_relative_difference<1e-8);
assert(result.K_relative_difference<1e-8);
assert(result.trajectory_max_scaled_difference<1e-5);
assert(result.HJB_identity_max_relative_residual<1e-8);
