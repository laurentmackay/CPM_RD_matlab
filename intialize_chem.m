%~~~~~~~~setting up the chemical state of the cell~~~~~~~~~~~~~~~~~~~~~~
time=0;
reactions=0;

D_1=0.43;                  %inactive rho/rac
D_2=0.02;                  %active rho/rac
D_3=0.03;                  %pax
D = [D_1 D_1 D_2 D_2 D_3 D_3];
N_instantaneous=50;

%Parameters from (Tang et al., 2018) 
B_1=4.26;  
gamma = 0.3;
L_R=0.34;
delta_rho=0.016;
L_rho=0.34;
m=4;                                    %More commonly known as 'n' in literature
L_K=5.77;
delta_P=0.00041;
k_X=41.7;
k_G=5.71;
k_C=5;
GIT=0.11;
PIX=0.069;
Paxtot=2.3;
alpha_R=15;
Rtot=7.5; % from Suplemental 1
alpha=alpha_R/Rtot;
delta_R = 0.025;
PAKtot = gamma*Rtot;

I_R = 0.003;
I_K = 0.009;
I_rho = 0.016;
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

RhoRatio = 0.3;
RacRatio = 0.13;
PaxRatio = 0.15;


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
reaction(:,:,3) = delta_rho*ones(N,N);                                             %From active rho to inactive rho
reaction(:,:,4) = delta_R*ones(N,N);                                               %From active Rac to inactive Rac
reaction(:,:,6) = delta_P*ones(N,N);                                               %From phosphorylated Pax to unphosphorylated Pax

alpha_chem=zeros(N,N,6);
alpha_diff=zeros(6,1);


% determining fractional expression per grid point
RhoRatio=x(:,:,3)./(x(:,:,3)+x(:,:,1));
RhoRatio(isnan(RhoRatio))=0;
RacRatio=x(:,:,4)./(x(:,:,4)+x(:,:,2)+x(:,:,7));
RacRatio(isnan(RacRatio))=0;
PaxRatio=x(:,:,6)./(x(:,:,6)+x(:,:,5)+x(:,:,8));
PaxRatio(isnan(PaxRatio))=0;
RbarRatio=x(:,:,7)./(x(:,:,4)+x(:,:,2)+x(:,:,7));  % this is gamma*K    
RbarRatio(isnan(RbarRatio))=0;

%----reactions propensites that vary lattice ot lattice 
K_is=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio).*(1+alpha_R*RacRatio)+k_G*k_X*GIT*PIX);
K=RbarRatio/gamma;         %changed from paper
I_Ks=I_K*(1-K_is.*(1+alpha_R*RacRatio));
reaction(:,:,1) = I_rho*(L_R^m./(L_R^m +(RacRatio+RbarRatio).^m));            %From inactive rho to active rho changed from model
reaction(:,:,2) = (I_R+I_Ks).*(L_rho^m./(L_rho^m+RhoRatio.^m));
reaction(:,:,5) = B_1*(K.^m./(L_K^m+K.^m));
alpha_chem=reaction(:,:,1:6).*x(:,:,1:6);