% clc
% clearvars
% close all;
% set(0,'DefaultFIgureVisible','on');
%restarting the environment

% mkdir 'results'
%vid = VideoWriter(['results'],'MPEG-4'); % this saves videos to mp4 change to whatever's convenient

figure(1);clf();
set(gcf,'defaultaxesfontsize',14);
d=[0.04, 0.04];
panelA=subplot(2,2,1); annotatePlot('A',22,d);
panelB=subplot(2,2,2); annotatePlot('B',22,d);
panelC=subplot(2,2,3); annotatePlot('C',22,d);
panelD=subplot(2,2,4); annotatePlot('D',22,d);
Results=[];
Times=[];

%open(vid);

N_species=8; %number of chemical species
N_rx=6; %number of reactions (this should determined automatically)
Ttot=3.2e3; %time the simulation end
SF=1; % speed factor I divide molecule number by this for speed
Gsize=20; %length of the grid in um
N=30; % number of points used to discretize the grid
shape=[N,N];
sz=prod(shape);
len=Gsize/N; %length of a latice square
[j, i] = meshgrid(1:shape(2),1:shape(1)); %the i and j need to be reversed for some reason (\_(:0)_/)

div=0.1;

initialize_cell_geometry
initialize_cellular_potts
initialize_chem %all reaction-diffusion parameter are getting initialized

lastplot=0;
tic
picstep=0.05;(len)/vmax; %timepoints where we take a frame for the video
z=1;

center=zeros(floor(Ttot/picstep)+1,2); %an array where we store the COM
center(z,:)=com(cell_mask);

% Results=zeros(N,N,N_species+1,floor(Ttot/picstep)+1); %an array where we store results
pic_WP %takes a frame for the video

nrx=4e5; %number of times reactions are carried out in a chem_func loop
reactions=0; %intializing a reaction counter


%diffusion timestep
eps=0.00005;
pmax=0.004;%0.5;%max allowed pT for one cell

dt=pmax*(h^2)/(max(D)*size(jump,1));%auto-determine timestep

%intializing variables for enumerate_diffusion.m making sure their size is
%constant
ij0=(1:(sz))';
diffuse_mask=false(size(jump,2),sz);
num_diffuse=zeros(1,size(jump,2));
ij_diffuse=zeros(4,(N)*(N));
diffusing_species_sum=zeros(size(jump,2),6);
num_vox_diff=zeros(1,sz);
pT0 = zeros(sz,length(D));
pi = zeros(size(jump,2),sz);
dt_diff=zeros(size(D));

enumerate_diffusion %determines the possible diffusion reactions in a way that be convereted to c

%vectorized indcies for the reaction and diffusion propensites
id0=(diffusing_species-1)*sz;


%initializing total chemical reaction and diffusion-reaction propensities
alpha_diff=sum(diffusing_species_sum).*D/(h*h);
alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
alpha_rx2=alpha_rx;
a_total=sum(alpha_diff)+sum(alpha_rx(:));


if (1/a_total)*nrx>(len)/(4*vmax) %makes sure that you don't stay in the CPM__chem func loop for to long
    error('cell moving to fast consider lowering nrx')
end

numDiff=0;
numReac=0;
%arrays recording Ratio change
%after run, plot TRac/TRho over Timeseries
Timeseries=[];
TRac=[];
TRho=[];
TPax=[];

last_time=time; %used to time the CMP_step
tic
rx_speedup=2;
rx_count=zeros(shape);
dt_diff=zeros(size(D));
P_diff=0.5;

mk_fun("SSA0")

while time<Ttot
    A=nnz(cell_mask); %current area
    cell_inds(1:A)=find(cell_mask); %all cell sites padded with 0s
    
    while (time-last_time)<Ttot
        
        alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
        
        % has to be a function to be put under a c wrapper (slow part)
%         [x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
%             alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
%             RacRatio,RbarRatio,I_Ks,reaction,jump,ir0,id0,cell_inds,...
%             k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
%             alpha,PAKtot,rx_count,dt_diff] =  SSA_mex(x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
%             alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
%             RacRatio,RbarRatio,I_Ks,reaction,jump,ir0,id0,cell_inds,...
%             k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
%             alpha,PAKtot,nrx,A,num_vox_diff,pi,pT0,dt_diff,rx_count,0.5,rx_speedup);

        [A,B_1,D,GIT,I_K,I_Ks,I_R,I_rho,K,K_is,L_K,L_R,L_rho,PAKtot,PIX,P_diff,...
        PaxRatio,Pax_Square,Paxtot,RacRatio,Rac_Square,RhoRatio,Rho_Square,a_total,alpha,...
        alpha_R,alpha_chem,alpha_rx,cell_inds,diffuse_mask,diffusing_species_sum,dt_diff,...
        gamma,h,id0,ir0,jump,k_C,k_G,k_X,m,nrx,pT0,pi,reaction,rx_count,rx_speedup,time,...
        x] = SSA0_fun(A,B_1,D,GIT,I_K,I_Ks,I_R,I_rho,K,K_is,L_K,L_R,L_rho,PAKtot,PIX,P_diff,...
        PaxRatio,Pax_Square,Paxtot,RacRatio,Rac_Square,RhoRatio,Rho_Square,a_total,alpha,alpha_R,...
        alpha_chem,alpha_rx,cell_inds,diffuse_mask,diffusing_species_sum,dt_diff,gamma,h,id0,...
        ir0,jump,k_C,k_G,k_X,m,nrx,pT0,pi,reaction,rx_count,rx_speedup,time,x);
        
        reactions=reactions+nrx; %reaction counter
        
        if time>=lastplot+picstep % takes video frames
            pic_WP
            
%             Timeseries=[Timeseries time];
%             
%             TRac=[TRac sum(sum(x(:,:,4)))/sum(sum(sum(x(:,:,[4 2 7]))))];
%             TRho=[TRho sum(sum(x(:,:,3)))/sum(sum(sum(x(:,:,[1 3]))))];
%             TPax=[TPax sum(sum(x(:,:,6)))/sum(sum(sum(x(:,:,[6 5 8]))))];
% 
%             z=z+1;
%             center(z,:)=com(cell_mask);
            lastplot=time;
            time
%             reactions
            toc
        end
        
    end
    
    last_time=time;
    
%     for kk=1:Per %itterates CPM step Per times
%         CPM_step
%     end
    
    enumerate_diffusion %recaucluates diffusable sites
end



%calculates speed by the distance the COM moved every 120s 
%thats (how they do it experimentally)
% dx=zeros(floor(finaltime/120),1);
% dx(1)=sqrt(sum((center(120/picstep,:)-center(1,:)).^2));
% for i=2:length(dx-1)
%     dx(i)=sqrt(sum((center(i*120/picstep,:)-center((i-1)*120/picstep,:)).^2));
% end
% v=dx/120*3600;
% 
% %saving final results 
toc


% figure(1);
% plot(Timeseries,TRho,Timeseries,TRac,Timeseries,TPax);
% legend('Rho','Rac','Pax','Location','Best');
% xlabel('Time');
% title(['B = ' num2str(B_1)]);


%close(vid);
% cur=pwd;
% cd results
% save(['final'])
% cd(cur)
