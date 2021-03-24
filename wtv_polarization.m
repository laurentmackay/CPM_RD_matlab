model_name='chem_Rx_Pax_Kathy';
deploy_model(model_name,1);
mk_fun2('main_FVM_polarization',{},{'B','lam_p_0','dt','copyNum','cpmstep0','model_name'})

set_experiment('lam_p_0_Inf_ext');
save_dir=results_dir();


B_vals=linspace(3,10,30);



N_reps=3;

% p=parpool(3);

cpmstep0=2e1;
dt=cpmstep0/5;

parfor i0=1:length(B_vals)*N_reps
    i=ceil(i0/N_reps);
    j=mod(i0-1,N_reps)+1;
    main_FVM_polarization_fun(B_vals(i), Inf, save_dir, dt,j, cpmstep0,model_name);
end
