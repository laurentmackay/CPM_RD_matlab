%% build functions

initialize_chem_params

params = inline_script(which('model_params'),{},{'B'});

eval_rhs_str = inline_script(which('eval_Rx'));

rhs_str =['function Rx = rhs_fun_tot(t,u,B)' newline 'u=transpose(u);' newline params newline eval_rhs_str newline 'Rx=transpose(Rx);' newline 'end'];

fid=fopen(strcat(work_dir(),'rhs_fun_tot.m'),'w');
fwrite(fid,rhs_str,'char');
fclose(fid);

clear rhs_fun_tot


%% run the PDE for various values of B
Nrep=5;
tic;
% B_vals=[linspace(2,2.45,4) linspace(2.5,3,10) linspace(3.1,5,6)];
B_vals=linspace(2,0.2,10);

polarized=cell(length(B_vals),1);
t_polarized=cell(length(B_vals),1);
for i=1:length(B_vals)
    B=B_vals(i);
    eval('model_fp');
    eval(params)
    

    eval('model_anon')
    tol=1e-14;
    

    
    fp0=fp;
    while any(abs(rhs_anon(fp))>tol)
        
        [T_vec,Y_vec] = ode15s(@(t,u) rhs_anon(u)',[0 1e4],fp,odeset('NonNegative',1:N_species));
        fp=Y_vec(end,:);
        
        
    end

rhs = pdepe_fun(B);
ic_fun = @(x) fp(1:N_slow)'*(1+0.05*randn());


Ttot=2e6;
Xmax=30;
tspan=linspace(0,Ttot,1e4);
% tspan=logspace(-5,log10(Ttot),100);
xmesh=linspace(0,Xmax,1e3);
options=odeset('RelTol',1e-5, 'AbsTol', 1e-5, 'InitialStep',1e-6);
i_rac=find(strcmp(chems,'Rac'));

for j=1:Nrep
    tic;
    sol = pdepe(0, rhs, ic_fun, @zeroflux, xmesh, tspan, options);
    toc
    delta=max(sol(:,:,i_rac),[],2)-min(sol(:,:,i_rac),[],2);
    polar = delta(end)>0.05;
    polarized{i}{end+1}=polar
    if polar
        t_polarized{i}{end+1}=tspan(find(delta>0.05,1));
    end
end
end
% 
% figure(1);clf();
% i_rac=find(strcmp('Rac',chems));
% 
% imagesc(tspan, xmesh, sol(:,:,i_rac)');
% 
% ylabel('Space (\mum)')
% xlabel('Time');
% colorbar




%%

figure(33);
subplot(1,2,1)
p=cellfun(@(x) mean([x{:}]), polarized)
e=cellfun(@(x) std([x{:}]), polarized)/sqrt(Nrep);
% h=errorbar(B_vals,p,e );
h=plot(B_vals,p);
set(h,'linewidth',2)
xlabel('B (s^{-1})');
ylabel('Prob. Self-Polarization')
set(gca,'FontSize',14)


for i=find(cellfun(@isempty, t_polarized))'
t_polarized{i}={};
end

subplot(1,2,2)
t=cellfun(@(x) mean([x{:}]), t_polarized)
e=cellfun(@(x) std([x{:}]), t_polarized)/sqrt(Nrep);

h=errorbar(B_vals,t/3600,e/3600 );
set(h,'linewidth',2)
xlabel('B (s^{-1})');
ylabel('Self-Polarization Time (hours)')
set(gca,'FontSize',14)


