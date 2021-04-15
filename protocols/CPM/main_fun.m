function [B_1,copyNum] = main_fun(B_1,copyNum)



if isempty(getCurrentTask())
    % clc
    % clearvars -except var
    % close all;
    set(0,'DefaultFigureVisible','on');
    %restarting the environment


    % this saves videos to avi change to whatever's convenient
    % only avi works for Beluga. If want to save space, use MPEG-4
    vid = VideoWriter(['results'],'Motion JPEG AVI');
    %frames=[]; % store the video frames in the mat file *take too much space
    %so don't use it for long-time simulation

    open(vid);
end

N_species=8; %number of chemical species
finaltime= 4*3600;
nrx=3e4; %number of times reactions are carried out in a chem_func loop
SF=2; % speed factor that reduces molecule number for speed
%%%CAUTION: don't divide too much so that each lattice has enough molecules%%%
Gsize=100; %length of the grid in micrometers
N=30; % number of points used to discretize the grid
sz=N^2;
h=Gsize/N; %length of a lattice square
[j, i,] = meshgrid(1:N,1:N); %the i and j need to be reversed for some reason (\_(:0)_/)

x=zeros(N,N,N_species); % where the chemical information is stored

% setting the initial shape of the cell to a circle of radius r in mircons
R=20;
cell_mask=(i-N/2).^2 +(j-N/2).^2 < (R/h)^2;
cell_maskp=cell_mask; % intializing a variable for CPM_step.m

% vectorizing the diffusion direction for grid points (periodic)
up = sub2ind([N,N],circshift(i,1,1),j);
down = sub2ind([N,N],circshift(i,-1,1),j);
left = sub2ind([N,N],i,circshift(j,1,2));
right = sub2ind([N,N],i,circshift(j,-1,2));

%combing diffusion directions into jump array
jump = zeros(sz,4);
jump(:,1) = up(:);
jump(:,2) = down(:);
jump(:,3) = left(:);
jump(:,4) = right(:);

%functions that measure perimeter and center of mass
perim = @(x) nnz(x&~x(up)) + nnz(x&~x(down)) + nnz(x&~x(left)) + nnz(x&~x(right));
com = @(x) [sum(sum(i.*x)),sum(sum(j.*x))]/nnz(x);


%parameter for the CPM scaled by length of a lattice point

lam_a=3*h^4; %energy cost of area change
lam_p=14*h^2; %energy cost of permiter change
J=0*h; %energy cost of change in medium contact

B_rho=2e4;%chemical potential rho
B_R=2e4*(.3/.13); %chemical potential rac
%(defined such that they have no net effect at the saddle)

a=1308/h^2; %ideal area      values from abira
per=128/h; %ideal permiter       values from abira 128 for perfect circle data 295
Hb=0; %membranes resistance to movement
T=100; %"temperture" strength of noise
vmax=3/60; %max speed of the cell



Per=perim(cell_mask); %current lattice permiter
A=nnz(cell_mask); %current lattice area

cell_inds=zeros(N*N,1);
cell_inds(1:A)=find(cell_mask);

H0=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian

time=0;
reactions=0;

D_1=0.43;                  %inactive rho/rac
D_2=0.02;                  %active rho/rac
D_3=0.03;                  %pax
D = [D_1 D_1 D_2 D_2 D_3 D_3];
N_instantaneous=50;

gamma = 0.3;
L_R=0.34;
delta_rho=0.016;
L_rho=0.34;
m=4;             %More commonly known as 'n' in literature
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



totalRho = 2250000/SF;
totalRac = 2250000/SF;
totalPax = 690000/SF;
Rho_Square = totalRho/(A);    %Average number of Rho per square
Rac_Square = totalRac/(A);    %Average number of Rac per square
Pax_Square = totalPax/(A);    %Average number of Pax per square




RhoRatio=0.35;
RacRatio=0.16;
PaxRatio=0.27;

numberofC = Rho_Square*RhoRatio;           %active Rho
numberofD = Rac_Square*RacRatio;           %active Rac
numberofF = Pax_Square*PaxRatio;           %phosphorylated Pax

K_is=1/((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio)*(1+alpha_R*RacRatio)+k_G*k_X*GIT*PIX);
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

N0=[numberofA numberofB numberofC numberofD numberofE numberofF numberofG numberofH];

x=reshape(round(N0.*cell_mask(:)),[N,N,N_species]);  %puts molecules in thier places

reaction=zeros(N,N,9);
reaction(:,:,3) = delta_rho*ones(N,N);   %From active rho to inactive rho
reaction(:,:,4) = delta_R*ones(N,N);     %From active Rac to inactive Rac
reaction(:,:,6) = delta_P*ones(N,N);     %From phosphorylated Pax to unphosphorylated Pax

alpha_chem=zeros(N,N,6);

RhoRatio=x(:,:,3)./(x(:,:,3)+x(:,:,1));
RhoRatio(isnan(RhoRatio))=0;
RacRatio=x(:,:,4)./(x(:,:,4)+x(:,:,2)+x(:,:,7));
RacRatio(isnan(RacRatio))=0;
PaxRatio=x(:,:,6)./(x(:,:,6)+x(:,:,5)+x(:,:,8));
PaxRatio(isnan(PaxRatio))=0;

K_is=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio).*(1+alpha_R*RacRatio)+k_G*k_X*GIT*PIX);
K=alpha_R*RacRatio.*K_is.*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio);
I_Ks=I_K*(1-K_is.*(1+alpha_R*RacRatio));
reaction(:,:,1) = I_rho*(L_R^m./(L_R^m +(RacRatio+gamma.*K).^m));
reaction(:,:,2) = (I_R+I_Ks).*(L_rho^m./(L_rho^m+RhoRatio.^m));
reaction(:,:,5) = B_1*(K.^m./(L_K^m+K.^m));
alpha_chem=reaction(:,:,1:6).*x(:,:,1:6);

timecheck=0;
tic
picstep=0.5*h/vmax; %timepoints where we take a frame for the video
z=1;


N_result=floor(finaltime/picstep)+1;
Results=zeros(N,N,N_species+1,N_result); %an array where we store results
Timeseries=zeros(1,N_result);
TRac=zeros(1,N_result);
TRho=zeros(1,N_result);
% %included in pic.m. Use these two lines when don't plot at the same time
Results(:,:,1,z)=cell_mask;
Results(:,:,2:(N_species+1),z)=x;

center=zeros(N_result,2); %an array where we store the COM
center(z,:)=com(cell_mask);


if usejava('desktop')
    close all
    % to make a video all frames must be the same size setting its position just
    % stops some bugs
    figure('Position', [200 75 1000 900])
    fs=14; %axis font size

    subplot(2,2,1)
    imagesc(cell_mask,[0 1]);
    % colorbar;
    center(z,:)=com(cell_mask);
    hold on
    plot(center(1:z,2),center(1:z,1),'r')
    hold off
    ax = gca;
    ax.FontSize = fs;
    % title('Cell', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')


    subplot(2,2,2)
    imagesc(RhoRatio,[0.1 0.6]);
    colorbar;
    ax = gca;
    ax.FontSize = fs;
    title('Rho', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')

    subplot(2,2,3)
    imagesc(RacRatio,[0 0.3]);
    colorbar
    ax = gca;
    ax.FontSize = fs;
    title('Rac', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')

    subplot(2,2,4)
    imagesc(PaxRatio,[0 0.4]);
    colorbar
    ax = gca;
    ax.FontSize = fs;
    title('Pax', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')

    %colormap jet
    %saveas(gcf,'CPM.png') %if you want a an image of a frame


    drawnow
    %adding videos the frame

    %frames=[frames getframe(gcf)];
    frame=getframe(gcf);
end

reactions=0; %intializing a reaction counter

%diffusion timestep
eps=0.00005;
pmax=0.004;%0.5;%max allowed pT for one cell


%intializing variables for enumerate_diffusion.m making sure their size is
%constant
diffuse_mask=false(size(jump,2),sz);
dt=pmax*(h^2)/(max(D)*size(jump,1));%auto-determine timestep
num_vox_diff=zeros(1,sz);
pT0 = zeros(sz,length(D));
pi = zeros(size(jump,2),sz);

diffusing_species=1:6; %only the first 6 species diffuse


for drx=1:size(jump,2) %itterating over all possible directions
    diffuse_mask(drx,:)=cell_mask(jump(:,drx))&cell_mask(:); %finding cites where diffusion is possible

end

for vox=1:size(diffuse_mask,2)
    num_vox_diff(vox)=nnz(diffuse_mask(:,vox)); %total number of pssoible diffusion reactions in a voxel
end

for i=1:A
    vox=cell_inds(i);
    pT0(vox,:) = num_vox_diff(vox)*D/(h^2);
    pi(:,vox)=diffuse_mask(:,vox)'./sum(diffuse_mask(:,vox));
end

%vectorized indices for the reaction and diffusion propensites
id0=(diffusing_species-1)*sz;
ir0=((1:size(alpha_chem,3))-1)*sz;

%intializing total chemical reaction propensities
alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));

%arrays recording Ratio changes
%after run, plot TRac/TRho against Timeseries
% Timeseries=[];
% TRac=[];
% TRho=[];

last_time=time; %used to time the CPM_step
rx_count=zeros([N,N]); %record number of reactions in each lattice square
dt_diff=zeros(size(D)); % timing for each species in each square since last diffusion
rx_speedup=2;
P_diff=0.5;

tic
while time<finaltime
    A=nnz(cell_mask); %current area
    cell_inds(1:A)=find(cell_mask);

    while (time-last_time)<picstep

        alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
        % can be a function to be put under a c wrapper (slow part)

        [x,alpha_rx,alpha_chem,time,PaxRatio,RhoRatio,K_is,K,RacRatio,...
            I_Ks,reaction,jump,ir0,id0,cell_inds,k_X,PIX,k_G,k_C,GIT,...
            Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,alpha,...
            PAKtot,rx_count,dt_diff] =  CPM_chem_func(x,alpha_rx,...
            alpha_chem,time,PaxRatio,RhoRatio,K_is,K,rx_speedup,RacRatio,...
            I_Ks,reaction,jump,ir0,id0,cell_inds,k_X,PIX,k_G,k_C,GIT,...
            Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,alpha,...
            PAKtot,nrx,A,pi,pT0,dt_diff,rx_count,P_diff,Rac_Square,...
            Pax_Square,Rho_Square);

        reactions=reactions+nrx; %reaction counter

        if time>=timecheck+picstep % takes video frames
            pic
            % %included in pic.m. Use these two lines when don't plot at the same time
            z=z+1;
            if z<=N_result
                Results(:,:,1,z) = cell_mask;
                Results(:,:,2:(N_species+1),z) = x;
                Timeseries(z) = time;
                TRac(z) = mean(RacRatio(cell_mask));
                TRho(z) = mean(RhoRatio(cell_mask));
                center(z,:)=com(cell_mask);
            end
            timecheck=timecheck+picstep;
            time;
            reactions;
            toc;
        end

    end

    last_time=time;

    for kk=1:Per %itterates CPM step Per times
        CPM_step
    end

    diffusing_species=1:6; %only the first 6 species diffuse
    
    
    for drx=1:size(jump,2) %itterating over all possible directions
        diffuse_mask(drx,:)=cell_mask(jump(:,drx))&cell_mask(:); %finding cites where diffusion is possible
    
    end
    
    for vox=1:size(diffuse_mask,2)
        num_vox_diff(vox)=nnz(diffuse_mask(:,vox)); %total number of pssoible diffusion reactions in a voxel
    end
    
    for i=1:A
        vox=cell_inds(i);
        pT0(vox,:) = num_vox_diff(vox)*D/(h^2);
        pi(:,vox)=diffuse_mask(:,vox)'./sum(diffuse_mask(:,vox));
    end
end



%saving final results
toc;
if isempty(getCurrentTask())
    close(vid);
    fn=['results/final_B_' num2str(B_1) '.mat'];
    ls results
    disp(['saving to: ' fn]);
    save(fn);
else
    fn=['results/final_B_' num2str(B_1) '_copy' int2str(copyNum) '.mat'];
    disp(['saving to: ' fn]);
    ls results
    save(fn);

end

end
