
%parameter for the CPM scaled by length of a lattice point

lam_a=3*len^4; %energy cost of area change
lam_p=9*len^2; %energy cost of permiter change
J=0*len; %energy cost of change in medium contact

B_rho=2e4;%chemical potential rho
B_R=2e4*(.3/.13); %chemical potential rac
%(defined such that they have no net effect at the saddle)

a=1308/len^2; %ideal area      values from abira
per=128/len; %ideal permiter       values from abira 128 for perfect circle data 295
Hb=0; %membranes resistance to movement
T=100; %"temperture" strength of noise
vmax=3/60; %max speed of the cell





H0=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian