% clc
% clearvars
% close all;
% set(0,'DefaultFIgureVisible','on');
%restarting the environment

mkdir 'results'
%vid = VideoWriter(['results'],'MPEG-4'); % this saves videos to mp4 change to whatever's convenient

%open(vid);

N_species=8; %number of chemical species
finaltime=2e4; %time the simulation end
SF=10; % speed factor I divide molecule number by this for speed
Gsize=100; %length of the grid in um
N=20; % number of points used to discretize the grid
sz=N^2;
len=Gsize/N; %length of a latice square
[j, i,] = meshgrid(1:N,1:N); %the i and j need to be reversed for some reason (\_(:0)_/)

x=zeros(N,N,N_species); % where the chemical information is stored

initialize_cell_geometry
initialize_cellular_potts
initialize_chem %all reaction-diffusion parameter are getting initialized

timecheck=0;
tic
picstep=(len)/vmax; %timepoints where we take a frame for the video
z=1;

center=zeros(floor(finaltime/picstep)+1,2); %an array where we store the COM
center(z,:)=com(cell_mask);

Results=zeros(N,N,N_species+1,floor(finaltime/picstep)+1); %an array where we store results

%pic %takes a frame for the video

nrx=3e4; %number of times reactions are carried out in a chem_func loop
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
while time<finaltime
    A=nnz(cell_mask); %current area
    cell_inds(1:A)=find(cell_mask); %all cell sites padded with 0s
    
    while (time-last_time)<(len/vmax)
        
        alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
        
        % has to be a function to be put under a c wrapper (slow part)
        [x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
            alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
            RacRatio,RbarRatio,I_Ks,reaction,ij_diffuse,jump,ir0,id0,cell_inds,...
            k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
            alpha,PAKtot,numDiff,numReac] =  CPM_chem_func_mex(x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
            alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
            RacRatio,RbarRatio,I_Ks,reaction,ij_diffuse,jump,ir0,id0,cell_inds,...
            k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
            alpha,PAKtot,nrx,A,numDiff,numReac);
        
        reactions=reactions+nrx; %reaction counter
        
        if time>=timecheck+picstep % takes video frames
            %pic
            Timeseries=[Timeseries time];
            
            TRac=[TRac sum(sum(x(:,:,4)))/sum(sum(sum(x(:,:,[4 2 7]))))];
            TRho=[TRho sum(sum(x(:,:,3)))/sum(sum(sum(x(:,:,[1 3]))))];
            TPax=[TPax sum(sum(x(:,:,6)))/sum(sum(sum(x(:,:,[6 5 8]))))];

            z=z+1;
            center(z,:)=com(cell_mask);
            timecheck=timecheck+picstep;
            time
            reactions
            toc
        end
        
    end
    
    last_time=time;
    
%     for kk=1:Per %itterates CPM step Per times
%         CPM_step
%     end
    
    enumerate_diffusion %recaucluates diffusable cites
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


figure(1);
plot(Timeseries,TRho,Timeseries,TRac,Timeseries,TPax);
legend('Rho','Rac','Pax','Location','Best');
xlabel('Time');
title(['B = ' num2str(B_1)]);


%close(vid);
% cur=pwd;
% cd results
% save(['final'])
% cd(cur)
