%~~~~~~~~setting up the chemical state of the cell~~~~~~~~~~~~~~~~~~~~~~
time=0;
reactions=0;

D_1=0.0143;%0.43;                  %inactive rho/rac
D_2=0.0007;%0.02;                  %active rho/rac
D_3=0.001%0.03;                  %pax
D = [D_1 D_1 D_2 D_2 D_3 D_3];
N_instantaneous=50;

%Parameters from (Tang et al., 2018) 
B_1=4.5;

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
RhoRatio = 0.3;
RacRatio = 0.13;
PaxRatio = 0.05;

% rough estimate of the induced state 
RhoRatio = 0.13;
RacRatio = 0.63;
PaxRatio = 0.35;


numberofC = Rho_Square*RhoRatio;           %active Rho
numberofD = Rac_Square*RacRatio;           %active Rac
numberofF = Pax_Square*PaxRatio;           %phosphorylated Pax

K_is=1/((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio)*(1+alpha_R*RacRatio)+k_G*k_X*GIT*PIX); %intial value of some ratios
K=alpha_R*RacRatio*K_is*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio);
P_i=1-PaxRatio*(1+k_G*k_X*k_C*GIT*PIX*PAKtot*K_is*(1+alpha_R*RacRatio));

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
x=reshape(round(N0.*cell_mask(:)),[N,N,N_species]);  %puts molecules in thier places 

reaction=zeros(N,N,9);

reaction(:,:,3) = delta_rho;                                             %From active rho to inactive rho
reaction(:,:,4) = delta_R;                                               %From active Rac to inactive Rac
reaction(:,:,6) = delta_P;                                               %From phosphorylated Pax to unphosphorylated Pax


alpha_chem=zeros([shape 6]);
alpha_rx=zeros(1,N_rx);
alpha_diff=zeros(6,1);
ir0=((1:N_rx)-1)*sz;

RhoRatio=zeros(shape);
RacRatio=zeros(shape);
RbarRatio=zeros(shape);
PaxRatio=zeros(shape);
K=zeros(shape);
K_is=zeros(shape);
I_Ks=zeros(shape);
%%
I=cell_inds(1:A);

[tmp,tmp2]=meshgrid(ir0,I);
I_chem=tmp+tmp2;

a_c_0=alpha_chem(I_chem);

%update Ratios
RacRatio(I)=x(I+(4-1)*sz)./(x(I+(4-1)*sz)+x(I+(2-1)*sz)+x(I+(7-1)*sz));
RbarRatio(I)=x(I+(7-1)*sz)./(x(I+(4-1)*sz)+x(I+(2-1)*sz)+x(I+(7-1)*sz));
RhoRatio(I)=x(I+(3-1)*sz)./(x(I+(3-1)*sz)+x(I+(1-1)*sz));
PaxRatio(I)=x(I+(6-1)*sz)./(x(I+(6-1)*sz)+x(I+(5-1)*sz)+x(I+(8-1)*sz));

RacRatio(isnan(RacRatio))=0;
RbarRatio(isnan(RbarRatio))=0;
RhoRatio(isnan(RhoRatio))=0;
PaxRatio(isnan(PaxRatio))=0;

%update other parameters    
K_is(I)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(I)).*(1+alpha_R*RacRatio(I))+k_G*k_X*GIT*PIX);
K(I)=alpha_R*RacRatio(I).*K_is(I).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(I));%RbarRatio(I)/gamma;         %changed from paper
I_Ks(I)=I_K*(1-K_is(I).*(1+alpha_R*RacRatio(I)));

reaction(I+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(I)+RbarRatio(I)).^m));            %From inactive rho to active rho changed from model
reaction(I+(2-1)*sz) = (I_R+I_Ks(I)).*(L_rho^m./(L_rho^m+RhoRatio(I).^m));                %From inactive Rac to active Rac
reaction(I+(5-1)*sz) = B_1*(K(I).^m./(L_K^m+K(I).^m));


alpha_chem(I_chem) = reaction(I_chem).*x(I_chem);
alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));


