%~~~~~~~~setting up the chemical state of the cell~~~~~~~~~~~~~~~~~~~~~~
time=0;
reactions=0;
D_1=0.43/10;                  %inactive rho/rac
D_2=0.02/10;                  %active rho/rac
D_3=0.03/10;                  %pax
D = [D_1 D_1 D_2 D_2 D_3 D_3];
N_instantaneous=50;
h=len/N;

%Parameters
B_1=1;   %orginal 4.26
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
Rtot=7.5; % from S1
alpha=alpha_R/Rtot;
delta_R = 0.025;
PAKtot = gamma*Rtot;

I_R = 0.003;
I_K = 0.009;
I_rho = 0.016;

%%%

%to find total number of molecules took concentrations from paper and
% assumed sperical cell r=5um gives N=3e5*[X]
totalRho = 2250000/10;
totalRac = 2250000/10;
totalPax = 690000/10;
Rho_Square = totalRho/(A);    %Average number of Rho per square
Rac_Square = totalRac/(A);    %Average number of Rac per square
Pax_Square = totalPax/(A);    %Average number of Pax per square

R_eq=0.13;
rho_eq=0.3;

RhoRatio = 0.3;
RacRatio = 0.13;
PaxRatio = 0.15;
Pax_reactions = 0;
Rac_reactions = 0;
Rho_reactions = 0;
diff_reactions= 0;

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
% x(:,:,(2-1)) = (numberofA)*cell_mask;     %inactive Rho
% x(:,:,(3-1)) = (numberofB)*cell_mask;     %inactive Rac
% x(:,:,(4-1)) = (numberofC)*cell_mask;     %active Rho
% x(:,:,(5-1)) = (numberofD)*cell_mask;     %active Rac
% x(:,:,(6-1)) = (numberofE)*cell_mask;     %unphosphorylated Pax
% x(:,:,(7-1)) = (numberofF)*cell_mask;     %phosphorylated Pax
% x(:,:,(8-1)) = (numberofG)*cell_mask;     %complexed Rac (Rbar)
% x(:,:,(9-1)) = (numberofH)*cell_mask;     %complexed Pax (Pbar)
% x=round(x);
x=reshape(round(N0.*cell_mask(:)),[N,N,N_species]);

reaction=zeros(N,N,9);
reaction(:,:,3) = delta_rho*ones(N,N);                                             %From active rho to inactive rho
reaction(:,:,4) = delta_R*ones(N,N);                                               %From active Rac to inactive Rac
reaction(:,:,6) = delta_P*ones(N,N);                                               %From phosphorylated Pax to unphosphorylated Pax

alpha_chem=zeros(N,N,6);
alpha_diff=zeros(6,1);

%intial values of some dynamic values

% determining fractional expression per grid point
RhoRatio=x(:,:,(4-1))./(x(:,:,(4-1))+x(:,:,(2-1)));
RhoRatio(isnan(RhoRatio))=0;
RacRatio=x(:,:,(5-1))./(x(:,:,(5-1))+x(:,:,(3-1))+x(:,:,(8-1)));
RacRatio(isnan(RacRatio))=0;
PaxRatio=x(:,:,(7-1))./(x(:,:,(7-1))+x(:,:,(6-1))+x(:,:,(9-1)));
PaxRatio(isnan(PaxRatio))=0;
RbarRatio=x(:,:,(8-1))./(x(:,:,(5-1))+x(:,:,(3-1))+x(:,:,(8-1)));  % this is gamma*K    
RbarRatio(isnan(RbarRatio))=0;

%----reactions that vary lattice ot lattice 
K_is=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio).*(1+alpha_R*RacRatio)+k_G*k_X*GIT*PIX);
K=RbarRatio/gamma;         %changed from paper
I_Ks=I_K*(1-K_is.*(1+alpha_R*RacRatio));
reaction(:,:,1) = I_rho*(L_R^m./(L_R^m +(RacRatio+RbarRatio).^m));            %From inactive rho to active rho changed from model
reaction(:,:,2) = (I_R+I_Ks).*(L_rho^m./(L_rho^m+RhoRatio.^m));
reaction(:,:,5) = B_1*(K.^m./(L_K^m+K.^m));
alpha_chem=reaction(:,:,1:6).*x(:,:,1:6);