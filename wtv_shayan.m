model_name='chem_Rx_Pax_Asheesh';
deploy_model(model_name,1);
mk_fun2('main_FVM',{},{'B','copyNum','model_name','Ttot'});

set_experiment('random_walk_B_sweep_200k_large_range');

save_dir=results_dir();

B_vals=linspace(2.5,25,10);

N_reps=5;

p=parpool(2);

Ttot=2e5;

parfor i0=1:length(B_vals)*N_reps
    i=ceil(i0/N_reps);
    j=mod(i0-1,N_reps)+1;
    main_FVM_fun(B_vals(i), save_dir, j, model_name, Ttot);
end
