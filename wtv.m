model_name='chem_Rx_Pax_Kathy';
deploy_model(model_name,1);
mk_fun2('main_FVM',{},{'B','lam_p_0','dt','copyNum','cpmstep0','model_name'})

B_vals=linspace(2,3,30);

N_reps=3;

% p=parpool(1);

parfor i0=1:length(B_vals)*N_reps
    i=ceil(i0/N_reps);
    j=mod(i0-1,N_reps)+1;
    main_FVM_fun(B_vals(i), Inf,10,j,2e2,model_name);
end
