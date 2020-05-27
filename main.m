clc
clearvars
close all;
set(0,'DefaultFIgureVisible','on');
%restarting the environment

mkdir 'results'
vid = VideoWriter(['results'],'MPEG-4'); % this saves videos to mp4 change to whatever's convenient

open(vid);

N_species=8; %number of chemical species
finaltime=100; %time the simulation end
SF=500; % speed factor I divide molecule number by this for speed
Gsize=100; %length of the grid in um
N=50; % number of points used to discretize the grid
sz=N^2;
len=Gsize/N; %length of a latice square
[j, i,] = meshgrid(1:N,1:N); %the i and j need to be reversed for some reason (\_(:0)_/)

x=zeros(N,N,N_species); % where the chemical information is stored

% setting the initial shape of the cell to a circle of radius r
R=20;
cell_mask=(i-N/2).^2 +(j-N/2).^2 < (R/len)^2;


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

lam_a=1*len^4; %energy cost of area change
lam_p=3*len^2; %energy cost of permiter change
J=0*len; %energy cost of change in medium contact

B_rho=2e4;%chemical potential rho
B_R=2e4*(.3/.13); %chemical potential rac
%(defined such that they have no net effect at the saddle)

a=1308/len^2; %ideal area      values from abira
per=128/len; %ideal permiter       values from abira 128 for perfect circle data 295
Hb=0; %membranes resistance to movement
T=100; %"temperture" strength of noise
vmax=3/60; %max speed of the cell



Per=perim(cell_mask); %current lattice permiter
A=nnz(cell_mask); %current lattice area

H0=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian


intialize_chem %all reaction-diffusion parameter are getting initialized

timecheck=0;
tic
picstep=(len)/vmax; %timepoints where we take a frame for the video
z=1;

center=zeros(floor(finaltime/picstep)+1,2); %an array where we store the COM
center(z,:)=com(cell_mask);

Results=zeros(N,N,N_species+1,floor(finaltime/picstep)+1); %an array where we store results

pic %takes a frame for the video

nrx=2e4; %number of times reactions are carried out in a chem_func loop
reactions=0; %intializing a reaction counter

%intializing variables for enumerate_diffusion.m making sure there size is
%constant
ij0=(1:(sz))';
diffuse_mask=false(size(jump,2),sz);
num_diffuse=zeros(1,size(jump,2));
ij_diffuse=zeros(4,(N)*(N));
diffusing_species_sum=zeros(size(jump,2),6);

enumerate_diffusion %determines the possible diffusion reactions in a way that be convereted to c

%vectorized indcies for the reaction and diffusion propensites
id0=(diffusing_species-1)*sz;
ir0=((1:size(alpha_chem,3))-1)*sz;

%intializing total chemical reaction and diffusion-reaction propensities
cell_inds=zeros(N*N,1);
cell_inds(1:A)=find(cell_mask);
alpha_diff=sum(diffusing_species_sum).*D/(h*h);
alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
alpha_rx2=alpha_rx;
a_total=sum(alpha_diff)+sum(alpha_rx(:));


if (1/a_total)*nrx>(len)/(4*vmax) %makes sure that you don't stay in the CPM__chem func loop for to long
    error('cell moving to fast consider lowering nrx')
end



last_time=time; %used to time the CMP_step
tic
while time<finaltime
    A=nnz(cell_mask); %current area
    cell_inds(1:A)=find(cell_mask); %all cell sites padded with 0s
    
    while (time-last_time)<(len/vmax)
        
        % has to be a function to be put under a c wrapper (slow part)
        [x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
            alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
            RacRatio,RbarRatio,I_Ks,reaction,ij_diffuse,jump,ir0,id0,...
            k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
            alpha,PAKtot] =  CPM_chem_func(x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
            alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
            RacRatio,RbarRatio,I_Ks,reaction,ij_diffuse,jump,ir0,id0,cell_inds,...
            k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
            alpha,PAKtot,nrx,A);
        
        reactions=reactions+nrx; %reaction counter
        
        if time>=timecheck+picstep % takes video frames
            pic
            z=z+1;
            center(z,:)=com(cell_mask);
            timecheck=timecheck+picstep;
            time
            reactions
            toc
        end
        
    end
    
    last_time=time;
    
    for kk=1:Per %itterates CPM step Per times
        CPM_step
    end
    
    enumerate_diffusion %recaucluates diffusable cites
end

%calculates speed by the distance the COM moved every 120s 
%thats (how they do it experimentally)
dx=zeros(floor(finaltime/120),1);
dx(1)=sqrt(sum((center(120/picstep,:)-center(1,:)).^2));
for i=2:length(dx-1)
    dx(i)=sqrt(sum((center(i*120/picstep,:)-center((i-1)*120/picstep,:)).^2));
end
v=dx/120*3600;

%saving final results 
toc
close(vid);
cur=pwd;
cd results
save(['final'])
cd(cur)
