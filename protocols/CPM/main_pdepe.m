initialize_chem_params

rhs = pdepe_fun();
ic_fun = pdepe_ic();

Ttot=1e4;
Xmax=20;
tspan=linspace(0,Ttot,100);
xmesh=linspace(0,Xmax,1e3);

sol = pdepe(0, rhs, ic_fun, @zeroflux, xmesh, tspan);


figure(1);clf();
i_rac=find(strcmp('Rac',chems));

imagesc(tspan, xmesh, sol(:,:,i_rac)');
ylabel('Space ($\mu$m)')
xlabel('Time');
colorbar





