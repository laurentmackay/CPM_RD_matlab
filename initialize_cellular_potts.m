
%parameter for the CPM scaled by hgth of a lattice point

lam_a=3*h^4; %energy cost of area change
lam_p=14*h^2; %energy cost of permiter change
J=0*h; %energy cost of change in medium contact

B_rho=2e3*h^2;%chemical potential rho
B_R=2e3*(.3/.13)*h^2; %chemical potential rac
%(defined such that they have no net effect at the saddle)

a=1308/h^2; %ideal area      values from abira
per=128/h; %ideal permiter       values from abira 128 for perfect circle data 295
Hb=0; %membranes resistance to movement
T=100; %"temperture" strength of noise






H0=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian

grow_count=0;
shrink_count=0;