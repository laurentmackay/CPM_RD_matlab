
plotting=usejava('desktop') && isempty(getCurrentTask())

if plotting %do not display pictures when running in parallel...i.e., on the cluster
    
    figure(1);clf();
    set(gcf,'defaultaxesfontsize',14);
    d=[0.04, 0.04];
    h_vec={};
    for i=1:9
        h_vec{i}=subplot(3,3,i);
    end
%     panelA=subplot(2,2,1); annotatePlot('A',22,d);
%     panelB=subplot(2,2,2); annotatePlot('B',22,d);
%     panelC=subplot(2,2,3); annotatePlot('C',22,d);
%     panelD=subplot(2,2,4); annotatePlot('D',22,d);
    
end



%open(vid);


nrx=3e4; %number of times reactions are carried out in a chem_func loop


Ttot=3*3.6e3; %Total simulation time
SF=2; % speed factor I divide molecule number by this for speed
Gsize=80; %length of the grid in um
N=80; % number of points used to discretize the grid
shape=[N,N];
sz=prod(shape);
h=Gsize/(N-1); %length of a latice square

vmax=3/60; %max speed of the cell
picstep=5;
cpmsteps=15;

cpmstep0=h/vmax;
cpmstep=cpmstep0/cpmsteps;



[j, i] = meshgrid(1:shape(2),1:shape(1)); %the i and j need to be reversed for some reason (\_(:0)_/)

div=0.1;

%prepare some .m files to model the chemical reactions from the reactions specified in `chem_Rx` file
mk_rxn_files('chem_Rx_Maree2012')

restart=false;

if ~restart
    initialize_cell_geometry
    initialize_cellular_potts
    initialize_chem_params
    initialize_chem %all reaction-diffusion parameter are getting initialized
end

lastplot=0;
lastcpm=0;
tic
%timepoints where we take a frame for the video

initialize_results


% Results=zeros(N,N,N_species+1,floor(Ttot/picstep)+1); %an array where we store results

if plotting
    pic_gen %takes a frame for the video
    gif('test.gif','frame',h_vec{2})
end
reactions=0; %intializing a reaction counter


%diffusion timestep
eps=0.00005;
pmax=0.004;%0.5;%max allowed pT for one cell

dt=pmax*(h^2)/(max(D)*size(jump,1));%auto-determine timestep

%intializing variables for enumerate_diffusion.m making sure their size is
%constant
N_dim=size(jump,2);
ij0=(1:(sz))';
diffuse_mask=false(N_dim,sz);
num_diffuse=zeros(1,size(jump,2));
ij_diffuse=zeros(4,(N)*(N));
diffusing_species_sum=zeros(N_dim,length(D));
num_vox_diff=zeros(1,sz);
pT0 = zeros(sz,length(D));
pi = zeros(N_dim,sz);
dt_diff=zeros(size(D));

enumerate_diffusion %determines the possible diffusion reactions in a way that be convereted to c

%vectorized indcies for the reaction and diffusion propensites
id0=(diffusing_species-1)*sz;


%initializing total chemical reaction and diffusion-reaction propensities
alpha_diff=sum(diffusing_species_sum).*D/(h*h);
alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
alpha_rx2=alpha_rx;
a_total=sum(alpha_diff)+sum(alpha_rx(:));


if (1/a_total)*nrx>h/(4*vmax) %makes sure that you don't stay in the CPM__chem func loop for to long
    error('cell moving to fast consider lowering nrx')
end

numDiff=0;
numReac=0;
cpmcounter=0;
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

SSA='SSA02';
SSA_fn=mk_fun(SSA,'gamma','alpha','pi','jump','beta');

d0=sum(sum(sum(x(:,:,:),3)))-(totalRac+totalRho);
SSA_call=[getFunctionHeader(SSA_fn) ';'];

disp(SSA_call)

while time<Ttot
    A=nnz(cell_mask); %current area
    cell_inds(1:A)=find(cell_mask); %all cell sites padded with 0s (initially)
    
    while (time-last_time)<Ttot
        
        alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
%         disp('gonna do SSA')
%         if sum(sum(sum(x(:,:,:),3)))-(totalRac+totalRho)~=d0
%             disp('main: molecules changed')
%         end
        
        x0=x;
        
%         disp('tryng SSA')    
        %run the SSA
        eval(['try' newline SSA_call newline 'catch err' newline 'disp(err.stack(1));' newline 'end']);
    
        try
%         [A,D,I_R,I_rho,L_R,L_rho,P_diff,RacRatio,Rac_Square,RhoRatio,Rho_Square,alpha_chem,alpha_rx,cell_inds,delta_R,delta_rho,diffuse_mask,diffusing_species_sum,dt_diff,h,id0,ir0,jump,m,nrx,pT0,pi,rx_count,rx_speedup,time,x] = SSA02_fun(A,D,I_R,I_rho,L_R,L_rho,P_diff,RacRatio,Rac_Square,RhoRatio,Rho_Square,alpha_chem,alpha_rx,cell_inds,delta_R,delta_rho,diffuse_mask,diffusing_species_sum,dt_diff,h,id0,ir0,jump,m,nrx,pT0,pi,rx_count,rx_speedup,time,x);
        catch err
        disp(err.stack.file);
        disp(err.stack.line);
        end
        
%        disp('just did SSA')
%         if sum(sum(sum(x(:,:,:),3)))-(totalRac+totalRho)~=d0
%             disp('main: molecules changed')
%         end
%         
        reactions=reactions+nrx; %reaction counter
        
        
        
        
        if time>=lastcpm+cpmstep
            
            for kk=1:(2*Per)/cpmsteps %itterates CPM step Per times
%                 disp('doing CPM')
                try
                    CPM_step
                catch
                    time=Ttot;
                    break;
                end
                
%                 if sum(sum(sum(x(:,:,:),3)))-(totalRac+totalRho)~=d0
%                     disp('main: molecules changed')
%                 end
            end
            
            enumerate_diffusion %recalcluates diffusable sites
            lastcpm=time;
            cpmcounter=cpmcounter+1;
        end
        
        if time>=lastplot+picstep || time==lastcpm % takes video frames
            pic_gen
            
            lastplot=time;      
            if cpmcounter==cpmsteps
                if plotting
                    gif
                end
                if ~isempty(getCurrentTask())
                    disp([num2str(copyNum) ': B=' num2str(B_1) ', t=' num2str(time)])
                end
                save_results
                cpmcounter=0;
            end
            
            
        end
        
    end
    
    last_time=time;
    
    
    
    enumerate_diffusion %recaucluates diffusable sites
end

toc

if isempty(getCurrentTask())  
%     close(vid);
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
