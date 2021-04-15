mk_rxn_files('chem_Rx_Pax2');

params = inline_script('model_params');
rhs_str = inline_script('eval_Rx');

str =['function Rx = rhs_fun(t,u)' newline 'u=transpose(u);' newline params newline rhs_str newline 'Rx=transpose(Rx);' newline 'end'];

fid=fopen('rhs_fun.m','w');
fwrite(fid,str,'char');
fclose(fid);

[T,Y] = ode15s(@ rhs_fun,[0 1e4],N0(1,:));

figure(1);clf();
plot(T,Y)

figure(2);clf();
plot(T,Y(:,2))