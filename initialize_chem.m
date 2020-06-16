%~~~~~~~~setting up the chemical state of the cell~~~~~~~~~~~~~~~~~~~~~~
time=0;
reactions=0;

D_1=0.43;                  %inactive rho/rac
D_2=0.02;                  %active rho/rac
D_3=0.03;                  %pax
D = [D_1 D_1 D_2 D_2 D_3 D_3];
N_instantaneous=50;

%Parameters from (Tang et al., 2018) 
B_1=5;

I_rho=0.016; L_rho=0.34; delta_rho=0.016;
L_R=0.34; I_R=0.003; delta_R=0.025; alpha_R=15; Rtot=7.5;
delta_P=0.0004; I_K=0.009;
L_K=5.77;

k_X=41.7; k_G=5.71; k_C=5;
GIT=0.11; PIX=0.069; Paxtot=2.3;
n=4; m=4; gamma=0.3;
PAKtot = gamma*Rtot;

alpha=alpha_R/Rtot;

h=len;

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
RacRatio_u = 0.085;
PaxRatio_u = 0.22;


% rough estimate of the induced state 
RhoRatio_i = 0.21;
RacRatio_i = 0.4; %0.215;
PaxRatio_i = 0.33;

RhoRatio=[RhoRatio_u; RhoRatio_i];
RacRatio=[RacRatio_u; RacRatio_i];
PaxRatio=[PaxRatio_u; PaxRatio_i];


numberofC = Rho_Square*RhoRatio;           %active Rho
numberofD = Rac_Square*RacRatio;           %active Rac
numberofF = Pax_Square*PaxRatio;           %phosphorylated Pax

K_is=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio).*(1+alpha_R*RacRatio)+k_G*k_X*GIT*PIX); %intial value of some ratios
K=alpha_R*RacRatio.*K_is.*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio);
P_i=1-PaxRatio.*(1+k_G*k_X*k_C*GIT*PIX*PAKtot*K_is.*(1+alpha_R*RacRatio));

if P_i<0
    error('No Initial Inactive Pax')
end

numberofA = Rho_Square - numberofC;               %inactive Rho that's ready to convert to active Rho
numberofB = Rac_Square*(1-RacRatio-gamma*K);      %inactive Rac that's ready to convert to active Rac
numberofE = Pax_Square*P_i;                       %unphosphorylated Pax that's ready to convert to phosphorylated Pax


numberofG = Rac_Square-numberofB-numberofD;       %RacGTP that is in a complex
numberofH = Pax_Square-numberofF-numberofE;       %phosphorylated Paxilin that is in a complex

%Setting up initial state of the cell
N0=[numberofA numberofB numberofC numberofD numberofE numberofF numberofG numberofH];

%setting up intial chemcial states 
% x=reshape(round(N0.*cell_mask(:)),[N,N,N_species]);  %puts molecules in thier places 
mask=induced_mask&cell_mask;
[tmp,tmp2]=meshgrid((0:N_species-1)*sz,find(mask));
i_induced=tmp+tmp2;
x(i_induced)=repmat(N0(2,:),nnz(mask),1);

mask=~induced_mask&cell_mask;
[tmp,tmp2]=meshgrid((0:N_species-1)*sz,find(mask));
i_induced=tmp+tmp2;
x(i_induced)=repmat(N0(1,:),nnz(mask),1);

x=round(x);

reaction=zeros(N,N,9);

reaction(:,:,3) = delta_rho;                                             %From active rho to inactive rho
reaction(:,:,4) = delta_R;                                               %From active Rac to inactive Rac
reaction(:,:,6) = delta_P;                                               %From phosphorylated Pax to unphosphorylated Pax


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
%% this should exactly match what is at the end of CPM_chem_func.m
i=cell_inds(1:A);

[tmp,tmp2]=meshgrid(ir0,i);
i_chem=tmp+tmp2;

a_c_0=alpha_chem(i_chem);

%update Ratios
RacRatio(i)=x(i+(4-1)*sz)./(x(i+(4-1)*sz)+x(i+(2-1)*sz)+x(i+(7-1)*sz));
RbarRatio(i)=x(i+(7-1)*sz)./(x(i+(4-1)*sz)+x(i+(2-1)*sz)+x(i+(7-1)*sz));
RhoRatio(i)=x(i+(3-1)*sz)./(x(i+(3-1)*sz)+x(i+(1-1)*sz));
PaxRatio(i)=x(i+(6-1)*sz)./(x(i+(6-1)*sz)+x(i+(5-1)*sz)+x(i+(8-1)*sz));

RacRatio(isnan(RacRatio))=0;
RbarRatio(isnan(RbarRatio))=0;
RhoRatio(isnan(RhoRatio))=0;
PaxRatio(isnan(PaxRatio))=0;

%update other parameters    
K_is(i)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i)).*(1+alpha_R*RacRatio(i))+k_G*k_X*GIT*PIX);
K(i)=alpha_R*RacRatio(i).*K_is(i).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(i));%RbarRatio(I)/gamma;         %changed from paper
I_Ks(i)=I_K*(1-K_is(i).*(1+alpha_R*RacRatio(i)));

reaction(i+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(i)+RbarRatio(i)).^m));            %From inactive rho to active rho changed from model
reaction(i+(2-1)*sz) = (I_R+I_Ks(i)).*(L_rho^m./(L_rho^m+RhoRatio(i).^m));                %From inactive Rac to active Rac
reaction(i+(5-1)*sz) = B_1*(K(i).^m./(L_K^m+K(i).^m));


alpha_chem(i_chem) = reaction(i_chem).*x(i_chem);
alpha_rx=sum(alpha_chem(i_chem));

