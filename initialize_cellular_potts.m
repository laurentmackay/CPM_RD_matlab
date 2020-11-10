
%parameter for the CPM scaled by hgth of a lattice point

%these parameters have been determined for h=1, and must be rescaled as h changes.
%their h-dependent scaling factors are taken from:
%
%Magno R, Grieneisen VA, Marée AF. 
%The biophysical nature of cells: potential cell behaviours revealed by analytical and computational studies of cell surface mechanics. 
%BMC Biophys. 2015;8:8. Published 2015 May 12. doi:10.1186/s13628-015-0022-x

lam_a=3*h^4; %energy cost of area change
lam_p=60*h^2; %energy cost of permiter change
J=0*h; %energy cost of change in medium contact

B_rho=2e3*h^2;%chemical potential rho
B_R=2e3*(.18/.13)*h^2; %chemical potential rac
%(defined such that they have no net effect at the saddle)

% a=1308/h^2; %ideal area      values from abira
% per=128/h; %ideal permiter       values from abira 128 for perfect circle data 295
a=A; %ideal area      values from abira
per=Per; %ideal permiter       values from abira 128 for perfect circle data 295
Hb=0; %membranes resistance to movement
T=900; %"temperture" strength of noise






H0=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian

grow_count=0;
shrink_count=0;