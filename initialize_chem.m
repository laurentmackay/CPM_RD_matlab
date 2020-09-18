initialize_chem_params

%~~~~~~~~setting up the chemical state of the cell~~~~~~~~~~~~~~~~~~~~~~
time=0;
reactions=0;

D_1=0.43;                  %inactive rho/rac
D_2=0.02;                  %active rho/rac
D_3=0.03;                  %pax
D = [D_1 D_1 D_2 D_2 D_3 D_3];

D=[D_1 D_2 D_1 D_2];

N_instantaneous=50;

%Parameters from (Tang et al., 2018) 
B_1=5;

I_rho=0.016;
I_R=0.020;

L_rho=0.34; delta_rho=0.016;
L_R=0.34;  delta_R=0.025; alpha_R=15; Rtot=7.5;
delta_P=0.0004; I_K=0.009;
L_K=5.77;

k_X=41.7; k_G=5.71; k_C=5;
GIT=0.11; PIX=0.069; Paxtot=2.3;
n=4; m=4; gamma=0.3;
PAKtot = gamma*Rtot;

alpha=alpha_R/Rtot;



%%%

%to find total number of molecules took concentrations from paper and
% assumed sperical cell r=5um gives N=3e5*[X]
totalRho = 2250000/SF;
totalRac = 2250000/SF;
totalPax = 690000/SF;
Rho_Square = totalRho/(A);    %Average number of Rho per square
Rac_Square = totalRac/(A);    %Average number of Rac per square
Pax_Square = totalPax/(A);    %Average number of Pax per square

R_eq=0; % i seet the equillbirum values to 0 so the Rac always causes exansion Rho always causes retraction
rho_eq=0;

% rough estimate of uninduced state
RhoRatio_u = 0.44;
RacRatio_u = 0.085;%0.085
% PaxRatio_u = 0.22;%0.22


% rough estimate of the induced state 
RhoRatio_i = 0.2;
RacRatio_i = 0.4; %0.215;

% RhoRatio_i = RhoRatio_u;
% RacRatio_i = RacRatio_u;



%0.215;
% PaxRatio_i = 0.33; %0.33

% RhoRatio_i = 0.35;
% RacRatio_i = 0.12; %0.215;
% PaxRatio_i = 0.27; %0.33
% 
% RhoRatio_i = 0.2;
% RacRatio_i = 0.2; %0.215;
% PaxRatio_i = 0.31; %0.33

RhoRatio=[RhoRatio_u; RhoRatio_i];
RacRatio=[RacRatio_u; RacRatio_i];
% PaxRatio=[PaxRatio_u; PaxRatio_i];


Rho0 = Rho_Square*RhoRatio;           %active Rho
Rhoi0 = Rho_Square - Rho0;               %inactive Rho that's ready to convert to active Rho

Rac0 = Rac_Square*RacRatio;           %active Rac
Raci0 = Rac_Square - Rac0;        %inactive Rac that's ready to convert to active Rac



%Setting up initial state of the cell
N0=[Raci0 Rac0 Rhoi0 Rho0];

%setting up intial chemcial states 
% x=reshape(round(N0.*cell_mask(:)),[N,N,N_species]);  %puts molecules in thier places 
mask=induced_mask&cell_mask;
[tmp,tmp2]=meshgrid((0:N_species-1)*sz,find(mask));
i_induced=tmp+tmp2;

x=zeros([shape ,N_species]); % where the chemical information is stored

x(i_induced)=repmat(N0(2,:),nnz(mask),1);

mask=~induced_mask&cell_mask;
[tmp,tmp2]=meshgrid((0:N_species-1)*sz,find(mask));
i_induced=tmp+tmp2;
x(i_induced)=repmat(N0(1,:),nnz(mask),1);

x=round(x);

% reaction=zeros(N,N,9);
% 
% reaction(:,:,3) = delta_rho;                                             %From active rho to inactive rho
% reaction(:,:,4) = delta_R;                                               %From active Rac to inactive Rac
% reaction(:,:,6) = delta_P;                                               %From phosphorylated Pax to unphosphorylated Pax


alpha_chem=zeros([shape N_rx]);
alpha_rx=zeros(1,N_rx);
alpha_diff=zeros(6,1); %the 6 here is not the same as N_rx...need to figure out a system
ir0=((1:N_rx)-1)*sz;

RhoRatio=zeros(shape);
RacRatio=zeros(shape);
RbarRatio=zeros(shape);
PaxRatio=zeros(shape);
K=zeros(shape);
K_is=zeros(shape);
I_Ks=zeros(shape);

vox=cell_inds(1:A);

% [x,sz,alpha_rx,...
%     alpha_chem,PaxRatio,RhoRatio,K_is,K,...
%     RacRatio,RbarRatio,I_Ks,reaction,ir0,...
%     k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
%     PAKtot,i]=update_alpha_chem(x,sz,alpha_rx,...
%     alpha_chem,PaxRatio,RhoRatio,K_is,K,...
%     RacRatio,RbarRatio,I_Ks,reaction,ir0,...
%     k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
%     PAKtot,i);

% update_all=true;
update_alpha_chem

RhoRatio0=RhoRatio;
RacRatio0=RacRatio;
