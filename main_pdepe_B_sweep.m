%% build functions

initialize_chem_params

params = inline_script(which('model_params'),{},{'B'});

eval_rhs_str = inline_script(which('eval_Rx'));

rhs_str =['function Rx = rhs_fun_tot(t,u,B)' newline 'u=transpose(u);' newline params newline eval_rhs_str newline 'Rx=transpose(Rx);' newline 'end'];

fid=fopen(strcat(work_dir(),'rhs_fun_tot.m'),'w');
fwrite(fid,rhs_str,'char');
fclose(fid);

clear rhs_fun_tot


%% run the PDE
tic;
rhs = pdepe_fun(3.5);
ic_fun = pdepe_ic();


Ttot=2e5;
Xmax=80;
tspan=linspace(0,Ttot,1e4);
% tspan=logspace(-5,log10(Ttot),100);
xmesh=linspace(0,Xmax,5e1);
options=odeset('RelTol',1e-8, 'AbsTol', 1e-8, 'InitialStep',0.01, 'MaxStep',10);

sol = pdepe(1, rhs, ic_fun, @zeroflux, xmesh, tspan, options);

toc

figure(1);clf();
i_rac=find(strcmp('Rac',chems));

imagesc(tspan, xmesh, sol(:,:,i_rac)');

ylabel('Space (\mum)')
xlabel('Time');
colorbar





