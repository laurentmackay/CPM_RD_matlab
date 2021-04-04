initialize_chem_params

i_chem_0 = ((1:N_species)-1)*sz;

%~~~~~~~~setting up the chemical state of the cell~~~~~~~~~~~~~~~~~~~~~~


% D_1=0.43;                  %inactive rho/rac
% D_2=0.02;                  %active rho/rac
% D_3=0.03;                  %pax
% D_3=1;                  %pax
% D = [D_1 D_1 D_2 D_2 D_3 D_3];
%
% D=[D_1 D_2 D_1 D_2 D_3 D_3 0 0];

% D=[D_1 D_2 D_1 D_2];

% D=[D_1 D_1 D_1 D_1];

%to find total number of molecules took concentrations from paper and
% assumed sperical cell r=5um gives N=3e5*[X]
totalRho = 2250000/SF;
totalRac = 2250000/SF;
totalPax = 690000/SF;
Rho_Square = totalRho/(A);    %Average number of Rho per square
Rac_Square = totalRac/(A);    %Average number of Rac per square
Pax_Square = totalPax/(A);    %Average number of Pax per square

Rho_Square = 1;    %Average number of Rho per square
Rac_Square = 1;    %Average number of Rac per square
Pax_Square = 1;    %Average number of Pax per square

N_instantaneous=50;

model_params

VolCell=(0.5*10^-6)*(h*10^-6)^2; %volume of a single lattice cell in m^3
muM = 6.02214*10^23*VolCell*10^3*10^-6; %Rescales from \muM to molecules

%%%


if length(D)==9
    Rho_Square=muM*3.1;
    Rac_Square=muM*7.5;
    C_Square=muM*2.4;
end


if length(D)==9
    Rho_Square=muM*3.1;
    Rac_Square=muM*7.5;
    C_Square=muM*2.4;
end



% rough estimate of uninduced state

RhoRatio_u = 0.8;
RacRatio_u = 0.045;%0.085

PaxRatio_u = 0.082;
% PaxRatio_u = 0.002;

% rough estimate of the induced state
RhoRatio_i = 0.02;
RacRatio_i = 0.35; %0.215;

if length(D)==9
    RacRatio_u = 0.35; %0.215;
    RacRatio_i = 0.85; %0.215;
    CRatio=[0.6; 0.5];
    
    RhoRatio_i = 0.05;
    RhoRatio_u = 0.4;
end

if length(D)==9
    RacRatio_u = 0.35; %0.215;
    RacRatio_i = 0.85; %0.215;
    CRatio=[0.6; 0.5];
    
    RhoRatio_i = 0.05;
    RhoRatio_u = 0.4;
end


% RhoRatio_i = RhoRatio_u;
% RacRatio_i = RacRatio_u;



%0.215;
PaxRatio_i = 0.2;

RhoRatio=[RhoRatio_u; RhoRatio_i];
RacRatio=[RacRatio_u; RacRatio_i];
PaxRatio=[PaxRatio_u; PaxRatio_i];


K_is=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio).*(1+alpha_R*RacRatio)+k_G*k_X*GIT*PIX); %intial value of some ratios
K=alpha_R*RacRatio.*K_is.*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio);
P_i=1-PaxRatio.*(1+k_G*k_X*k_C*GIT*PIX*PAKtot*K_is.*(1+alpha_R*RacRatio));

Rho0 = Rho_Square*RhoRatio;           %active Rho
Rhoi0 = Rho_Square - Rho0;

Rac0 = Rac_Square*RacRatio;           %active Rac


    fp=0; %trick for inlining
    eval('model_fp');
    rhs_anon=0; %trick for inlining
    eval('model_anon')
    tol=1e-14;
    
    relax=any(abs(rhs_anon(fp))>tol);
    
    fp0=fp;
    while any(abs(rhs_anon(fp))>tol)
        
        [T_vec,Y_vec] = ode15s(@(t,u) rhs_anon(u)',[0 1e4],fp,odeset('NonNegative',1:N_species));
        fp=Y_vec(end,:);
        
        
    end
    
    if relax
        disp('Relaxed to a new fixed point:')
        disp(strjoin(strcat(chems,'=',string(fp))),', ')
        fid=fopen(which('model_fp'),'w');
        fwrite(fid,['fp = [' num2str(fp,12) '];'],'char');
        fclose(fid);
    end

    
    %     N0(1,2)/Rac_Square
    %     N0(1,4)/Rho_Square
    %     N0(1,6)/Pax_Square
    if plotting && relax
        figure(3);clf();
        plot(T_vec,Y_vec);
        legend(chems)
        drawnow;
        
        yl=ylim();
        ylim([0 yl(end)]);
    end



% N0=[Raci0 Rac0 Rhoi0 Rho0];


%setting up intial chemcial states
% x=reshape(round(N0.*cell_mask(:)),[N,N,N_species]);  %puts molecules in their places

induced_ic

mask=induced_mask&cell_mask;
[tmp,tmp2]=meshgrid((0:N_species-1)*sz,find(mask));
i_induced=tmp+tmp2;

x=zeros([shape ,N_species]); % where the chemical information is stored

if ~isempty(ic)
    x(i_induced)=repmat(ic,nnz(mask),1);
    mask=~induced_mask&cell_mask;
else
    mask=cell_mask;
end

[tmp,tmp2]=meshgrid((0:N_species-1)*sz,find(mask));
x(tmp+tmp2)=repmat(fp,nnz(mask),1);

% x=round(x);
x=(1+noise*rand(size(x))).*x;
x(x<0)=0;
% reaction=zeros(N,N,9);
%
% reaction(:,:,3) = delta_rho;                                             %From active rho to inactive rho
% reaction(:,:,4) = delta_R;                                               %From active Rac to inactive Rac
% reaction(:,:,6) = delta_P;                                               %From phosphorylated Pax to unphosphorylated Pax


alpha_chem=zeros([shape N_rx]);
alpha_rx=zeros(1,N_rx);
alpha_diff=zeros(6,1); %the 6 here is not the same as N_rx...need to figure out a system
ir0=((1:N_rx)-1)*sz;

% RacRatio0=zeros(shape);
% RhoRatio=zeros(shape);
% RacRatio=zeros(shape);
% RbarRatio=zeros(shape);
% PaxRatio=zeros(shape);
% K=zeros(shape);
% K_is=zeros(shape);
% I_Ks=zeros(shape);
% 
% feedback=zeros(shape);
% Q_R=zeros(shape);
% Q_C=zeros(shape);
% Q_rho=zeros(shape);
% Q_P=zeros(shape);
vox=cell_inds(1:A);

eval_model


