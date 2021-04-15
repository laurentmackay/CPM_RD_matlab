


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

intialize_chem %all reaction-diffusion parameter are getting initialized

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


pic %takes a frame for the video

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

enumerate_diffusion %determines the possible diffusion reactions

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

    enumerate_diffusion %recaucluates diffusable sites
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

