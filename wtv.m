% mk_rxn_files('chem_Rx_Pax_Kathy');
mk_fun2('main_FVM',{},{'B','lam_p_0','dt','copyNum'})

B_vals=linspace(2,3,1);

N_reps=2;

% p=parpool(1);

parfor i0=1:length(B_vals)*N_reps
    i=ceil(i0/N_reps);
    j=mod(i0-1,N_reps)+1;
    main_FVM_fun(B_vals(i), Inf,10,j);
end
