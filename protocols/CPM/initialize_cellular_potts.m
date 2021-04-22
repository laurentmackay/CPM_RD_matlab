
%parameter for the CPM scaled by the size of a lattice cell

%these parameters have been determined for h=1, and must be rescaled as h changes.
%their h-dependent scaling factors are taken from:
%
%Magno R, Grieneisen VA, Marée AF. 
%The biophysical nature of cells: potential cell behaviours revealed by analytical and computational studies of cell surface mechanics. 
%BMC Biophys. 2015;8:8. Published 2015 May 12. doi:10.1186/s13628-015-0022-x

lam_a=0.3*h^4; %energy cost of area change
lam_p_0=0.3;
lam_p=lam_p_0*h^2; %energy cost of permiter change
J=0.1*h; %energy cost of change in medium contact

B_0=0.7;
B_rho=(B_0/0.3)*h^2;%chemical potential rho
B_R=(B_0/0.3)*(.18/.13)*h^2; %chemical potential rac
%(defined such that they have no net effect at the saddle)


a=A; %ideal area      values from abira
per=Per*(1 + (sqrt(2)-1)/2); %ideal permiter       values from abira 128 for perfect circle data 295
Hb=0; %membranes resistance to movement
T=0.5; %"temperture" strength of noise






H0=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian
dH_chem=0; %initial chemotactic contribution to the hamiltonian

grow_count=0;
shrink_count=0;