iter=0;
Nsteps=floor(Ttot/min(picstep,cpmstep))+1;



center=zeros(2,Nsteps);
Results=zeros([shape,N_species,Nsteps]);
Times=zeros(1,Nsteps);

areas=zeros(1,Nsteps);
perims=zeros(1,Nsteps);

Ham0=zeros(1,Nsteps);
Ham=zeros(1,Nsteps);


%save the initial conditions
save_results