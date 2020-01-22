
clc
clearvars
% close all;
% set(0,'DefaultFIgureVisible','on');
count=0;

folder = [pwd '/results'];
mkdir 'results'
vid = VideoWriter('results/vid','Indexed AVI');
open(vid);

N_species=8;
N=100;
sz=100*100;
len=100;
[j, i,] = meshgrid(1:N,1:N); %the i and j need to be reversed for some reason (�\_(:0)_/�)

x=zeros(N,N,N_species);
R=20;
% x(:,:,(1-1)) = ; %intial location of the cell
cell_mask=(i-N/2).^2 +(j-N/2).^2 < (N*R/len)^2;


up = sub2ind([N,N],circshift(i,1,1),j);
down = sub2ind([N,N],circshift(i,-1,1),j); 

left = sub2ind([N,N],i,circshift(j,1,2));
right = sub2ind([N,N],i,circshift(j,-1,2));

jump = zeros(sz,4);
jump(:,1) = up(:);
jump(:,2) = down(:);
jump(:,3) = left(:);
jump(:,4) = right(:);


perim = @(x) nnz(x&~x(up)) + nnz(x&~x(down)) + nnz(x&~x(left)) + nnz(x&~x(right));

com = @(x) [sum(sum(i.*x)),sum(sum(j.*x))]/nnz(x);


lam_a=3.1; %energy cost of area change
lam_p=8.8; %energy cost of permiter change
J=0; %energy cost of change in medium contact
a=1308/len^2; %ideal area      values from abira
per=295/len; %ideal permiter       values from abira
Hb=0; %membranes resistance to movement
T=30; %"temperture" strength of noise
B_R=10000;
B_rho=10000;


Per=perim(cell_mask); %current lattice
A=nnz(cell_mask); %current area
H0=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian

figure(1);
clf(1);
intialize_chem

xp=x;
N_update=1;
iplot=50;

timecheck=0;
tic
stepratio=100;
finaltime=0.1;
picstep=finaltime/2;
z=1;
center=zeros(floor(finaltime/picstep),2);
center(z,:)=com(cell_mask);
Results=zeros(N_species+1,N,N,floor(finaltime/picstep));
pic
vmax=50/3600;
last_time=time;

ij0=(1:(sz))';
diffuse_mask=false(size(jump,2),sz);
enumerate_diffusion

diffusing_species_sum=zeros(size(jump,2),6);
diffusing_species_sum2=zeros(size(jump,2),6);

id0=(diffusing_species-1)*sz;
ir0=((1:size(alpha_chem,3))-1)*sz;

for drx=1:length(ij_diffuse)
    diffusing_species_sum(drx,:)=sum(x(id0 + ij_diffuse{drx}));
end



tic
while time<finaltime
    cell_inds=find(cell_mask);
    alpha_rx=sum(alpha_chem(ir0 + cell_inds));
    alpha_rx2=alpha_rx;
    j=1:size(alpha_chem,1);
    nrx=1e2;
	
    while (time-last_time)<vmax/(h)
	
        [x,diffusing_species_sum,D,h,alpha_rx,...
alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
RacRatio,RbarRatio,I_Ks,reaction,ij_diffuse,jump,ir0,id0,cell_inds,...
k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
alpha,PAKtot] = CPM_chem_func(x,diffusing_species_sum,D,h,alpha_rx,...
alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
RacRatio,RbarRatio,I_Ks,reaction,ij_diffuse,jump,ir0,id0,cell_inds,...
k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
alpha,PAKtot,nrx);
        reactions=reactions+nrx;
        
        pic
    end
    disp(['t = ' num2str(time) ', dt = '  num2str(time-last_time) ])
    
    last_time=time;
    
%     for kk=1:Per
%         CPM_step
%     end
    
    
    

enumerate_diffusion

    reactions
end
toc
close(vid);
% cd results 
% save('final')