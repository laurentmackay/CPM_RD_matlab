function [B_1,copyNum] = main3_fun(B_1,copyNum)

if isempty(getCurrentTask()) %do not display pictures when running in parallel...i.e., on the cluster

    figure(1);clf();
    set(gcf,'defaultaxesfontsize',14);
    d=[0.04, 0.04];
    panelA=subplot(2,2,1); annotatePlot('A',22,d);
    panelB=subplot(2,2,2); annotatePlot('B',22,d);
    panelC=subplot(2,2,3); annotatePlot('C',22,d);
    panelD=subplot(2,2,4); annotatePlot('D',22,d);

end



%open(vid);


nrx=3e4; %number of times reactions are carried out in a chem_func loop


Ttot=3*3.6e3; %time the simulation end
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
mk_rxn_files('chem_Rx')

restart=false;

if ~restart
    
    
    
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
    
    R=0.2*N/2;
    cell_mask=(i-N/2).^2 +(j-N/2).^2 < R^2;
    induced_mask=cell_mask & (i-min(i(cell_mask)))<=2*div*(max(i(cell_mask))-min(i(cell_mask)));
    induced_mask(:)=0;
    
    i0=i;
    j0=j;
    
    Per=perim(cell_mask); %current lattice permiter
    A=nnz(cell_mask); %current lattice area
    
    
    cell_maskp=cell_mask; % intializing a variable for CPM_step.m
    cell_inds=zeros(N*N,1);
    cell_inds(1:A)=find(cell_mask);
    
    
    
    
    
    
    lam_a=3*h^4; %energy cost of area change
    lam_p=60*h^2; %energy cost of permiter change
    J=0*h; %energy cost of change in medium contact
    
    B_rho=2e3*h^2;%chemical potential rho
    B_R=2e3*(.18/.13)*h^2; %chemical potential rac
    
    a=A; %ideal area      values from abira
    per=Per; %ideal permiter       values from abira 128 for perfect circle data 295
    Hb=0; %membranes resistance to movement
    T=900; %"temperture" strength of noise
    
    
    
    
    
    
    H0=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian
    dH_chem=0; %initial chemotactic contribution to the hamiltonian
    
    grow_count=0;
    shrink_count=0;
    N_species = 8;
    N_rx = 6;
    initialize_chem_params
    
    time=0;
    reactions=0;
    
    D_1=0.43;                  %inactive rho/rac
    D_2=0.02;                  %active rho/rac
    D_3=0.03;                  %pax
    D = [D_1 D_1 D_2 D_2 D_3 D_3];
    
    D=[D_1 D_2 D_1 D_2 D_3 D_3 0 0];
    
    
    
    N_instantaneous=50;
    
    
    
    I_rho=0.016;
    I_R=0.003;
    I_K=0.009;
    delta_R=0.025;
    delta_rho=0.016;
    delta_P=0.0004;
    
    L_rho=0.34;
    L_R=0.34;
    L_K=5.77;
    
    
    alpha_R=15; Rtot=7.5;
    
    
    
    k_X=41.7; k_G=5.71; k_C=5;
    GIT=0.11; PIX=0.069; Paxtot=2.3;
    n=3; m=4; gamma=0.3;
    PAKtot = gamma*Rtot;
    
    alpha=alpha_R/Rtot;
    
    
    
    
    totalRho = 2250000/SF;
    totalRac = 2250000/SF;
    totalPax = 690000/SF;
    Rho_Square = totalRho/(A);    %Average number of Rho per square
    Rac_Square = totalRac/(A);    %Average number of Rac per square
    Pax_Square = totalPax/(A);    %Average number of Pax per square
    
    
    
    
    RhoRatio_u = 0.55;
    RacRatio_u = 0.12;%0.085
    
    PaxRatio_u = 0.22;
    
    
    RhoRatio_i = 0.2;
    RacRatio_i = 0.5; %0.215;
    
    
    
    
    
    PaxRatio_i = 0.33;
    
    
    RhoRatio=[RhoRatio_u; RhoRatio_i];
    RacRatio=[RacRatio_u; RacRatio_i];
    PaxRatio=[PaxRatio_u; PaxRatio_i];
    
    K_is=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio).*(1+alpha_R*RacRatio)+k_G*k_X*GIT*PIX); %intial value of some ratios
    K=alpha_R*RacRatio.*K_is.*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio);
    P_i=1-PaxRatio.*(1+k_G*k_X*k_C*GIT*PIX*PAKtot*K_is.*(1+alpha_R*RacRatio));
    
    
    Rho0 = Rho_Square*RhoRatio;           %active Rho
    Rhoi0 = Rho_Square - Rho0;               %inactive Rho that's ready to convert to active Rho
    
    Rac0 = Rac_Square*RacRatio;           %active Rac
    if length(D)==8
        Raci0 = Rac_Square*(1-RacRatio-gamma*K);        %inactive Rac that's ready to convert to active Rac
    else
        Raci0 = Rac_Square*(1-RacRatio);
    end
    
    Pax0 = Pax_Square*PaxRatio;           %active Rac
    Paxi0 = Pax_Square*P_i;        %inactive Rac that's ready to convert to active Rac
    
    
    RacPak0 = Rac_Square - Rac0 - Raci0;
    GPP0 = Pax_Square - Pax0 - Paxi0;
    
    N0=[Raci0 Rac0 Rhoi0 Rho0 Paxi0 Pax0];
    
    N0=[Raci0 Rac0 Rhoi0 Rho0 Paxi0 Pax0 RacPak0 GPP0];
    
    
    mask=induced_mask&cell_mask;
    [tmp,tmp2]=meshgrid((0:N_species-1)*sz,find(mask));
    i_induced=tmp+tmp2;
    
    x=zeros([shape ,N_species]); % where the chemical information is stored
    
    x(i_induced)=repmat(N0(2,:),nnz(mask),1);
    
    mask=~induced_mask&cell_mask;
    [tmp,tmp2]=meshgrid((0:N_species-1)*sz,find(mask));
    i_induced=tmp+tmp2;
    x(i_induced)=repmat(N0(1,:),nnz(mask),1);
    
    x=round(x);
    
    
    
    alpha_chem=zeros([shape N_rx]);
    alpha_rx=zeros(1,N_rx);
    alpha_diff=zeros(6,1); %the 6 here is not the same as N_rx...need to figure out a system
    ir0=((1:N_rx)-1)*sz;
    
    RhoRatio=zeros(shape);
    RacRatio=zeros(shape);
    RbarRatio=zeros(shape);
    PaxRatio=zeros(shape);
    K=zeros(shape);
    K_is=zeros(shape);
    I_Ks=zeros(shape);
    
    vox=cell_inds(1:A);
    
    
    update_alpha_chem
    
    RhoRatio0=RhoRatio;
    RacRatio0=RacRatio;
end

lastplot=0;
lastcpm=0;
tic
%timepoints where we take a frame for the video

iter=0;


Nsteps=floor(Ttot/min(cpmstep0))+1;


center=zeros(2,Nsteps);
Results=zeros([shape,N_species,Nsteps]);
Times=zeros(1,Nsteps);

areas=zeros(1,Nsteps);
perims=zeros(1,Nsteps);

Ham0=zeros(1,Nsteps);
Hchem=zeros(1,Nsteps);


save_results


% Results=zeros(N,N,N_species+1,floor(Ttot/picstep)+1); %an array where we store results
if usejava('desktop') && isempty(getCurrentTask())
    % to make a video all frames must be the same size setting its position just
    % stops some bugs
    fs=14; %axis font size

    subplot(2,2,1)
    plotCellIm(panelA,double(cell_mask),cell_mask,i0,j0)
    colorbar;

    hold on
    try
    plot(center(2,1:iter),center(1,1:iter),'r')
    catch e
        disp(e)
    end
    hold off
    ax = gca;
    ax.FontSize = fs;
    % title('Cell', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')


    subplot(2,2,2)
    plotCellIm(panelB,RhoRatio,cell_mask,i0,j0)
    caxis([0 1])
    colorbar;
    ax = gca;
    ax.FontSize = fs;
    title('Rho', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')

    subplot(2,2,3)
    plotCellIm(panelC,RacRatio,cell_mask,i0,j0)
    caxis([0 0.4])
    colorbar
    ax = gca;
    ax.FontSize = fs;
    set(gca,'Color',[1 1 1]*1)
    title('Rac', 'Fontsize', 24)
    axis tight



     subplot(2,2,4)
    plotCellIm(panelD,PaxRatio,cell_mask,i0,j0)
    caxis([0 1])
    colorbar
    ax = gca;
    ax.FontSize = fs;
    set(gca,'Color',[1 1 1]*1)
    title('Pax', 'Fontsize', 24)
    axis tight




    % xlabel('X')
    % ylabel('Y')
    % xlabel('X')
    % ylabel('Y')

    %colormap jet
    %saveas(gcf,'CPM.png') %if you want a an image of a frame


    drawnow
    %adding videos the frame

    %frames=[frames getframe(gcf)];
    frame=getframe(gcf);
end
if usejava('desktop') && isempty(getCurrentTask())
    gif('test.gif','frame',panelC)
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

diffusing_species=1:N_species; %only the first 6 species diffuse


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
SSA_fn=mk_fun(SSA,'gamma','alpha','pi');
SSA_call=[getFunctionHeader(SSA_fn) ';'];
d0=sum(sum(sum(x(:,:,:),3)))-(totalRac+totalRho);
while time<Ttot
    A=nnz(cell_mask); %current area
    cell_inds(1:A)=find(cell_mask); %all cell sites padded with 0s (initially)

    while (time-last_time)<Ttot

        alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
%         disp('gonna do SSA')
        if sum(sum(sum(x(:,:,:),3)))-(totalRac+totalRho)~=d0
            disp('molecules changed')
        end

        x0=x;

        %run the SSA
        eval(['try' newline SSA_call newline 'catch err' newline 'disp(err);' newline 'end']);
    disp('tryng SSA')
        try
%         [A,D,I_R,I_rho,L_R,L_rho,P_diff,RacRatio,Rac_Square,RhoRatio,Rho_Square,alpha_chem,alpha_rx,cell_inds,delta_R,delta_rho,diffuse_mask,diffusing_species_sum,dt_diff,h,id0,ir0,jump,m,nrx,pT0,pi,rx_count,rx_speedup,time,x] = SSA02_fun(A,D,I_R,I_rho,L_R,L_rho,P_diff,RacRatio,Rac_Square,RhoRatio,Rho_Square,alpha_chem,alpha_rx,cell_inds,delta_R,delta_rho,diffuse_mask,diffusing_species_sum,dt_diff,h,id0,ir0,jump,m,nrx,pT0,pi,rx_count,rx_speedup,time,x);
        catch err
        disp(err.stack.file);
        disp(err.stack.line);
        end

%        disp('just did SSA')
        if sum(sum(sum(x(:,:,:),3)))-(totalRac+totalRho)~=d0
            disp('molecules changed')
        end

        reactions=reactions+nrx; %reaction counter




        if time>=lastcpm+cpmstep

            for kk=1:(2*Per)/cpmsteps %itterates CPM step Per times
                disp('doing CPM')
                try
                    adj_empty = ~cell_mask(up) | ~cell_mask(down) | ~cell_mask(left) | ~cell_mask(right); %cells at the boundry
                    adj_full = cell_mask(up) | cell_mask(down) | cell_mask(left) | cell_mask(right); %cells at the boundry
                    
                    bndry_cell = cell_mask & adj_empty;
                    bndry_empty = ~cell_mask & adj_full;
                    bndry = find( bndry_cell | bndry_empty );
                    
                    bndry_up=cell_mask  & ~cell_mask(up);
                    bndry_down=cell_mask  & ~cell_mask(down);
                    bndry_l=cell_mask  & ~cell_mask(left);
                    bndry_r=cell_mask  & ~cell_mask(right);
                    
                    bndry_ud= bndry_up | bndry_down;
                    bndry_lr= bndry_l | bndry_r;
                    
                    bndrys=[bndry_up(:) bndry_down(:) bndry_l(:) bndry_r(:)];
                    
                    
                    if any(cell_maskp~=cell_mask)
                        error('not reseting')
                    end
                    
                    
                    rho_eq=mean(RhoRatio(find(cell_mask)));
                    R_eq=mean(RacRatio(find(cell_mask)));
                    Ncell_mask=squeeze(sum(sum(x))); %for a sanity check
                    A0=A;
                    
                    no_holes=false;
                    
                    while ~no_holes
                        vox_trial = bndry(randi(length(bndry)));
                    
                        r=randi(size(jump,2));
                        vox_ref=jump(sub2ind([sz,4],vox_trial,r));
                        cell_maskp(vox_trial) = cell_mask(vox_ref);% make a new trial configuration
                    
                        Per=perim(cell_maskp); % perimeter
                        A=nnz(cell_maskp); % area
                        HA=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian after the possible change
                        dH=HA-H0;
                        no_holes = getfield(bwconncomp(cell_maskp,4),'NumObjects')==1 && getfield(bwconncomp(~cell_maskp,4),'NumObjects')==1 ;%makes sure the cell stays connected and no hole
                        if ~no_holes
                            cell_maskp(vox_trial)=cell_mask(vox_trial);%revert change
                        end
                    end
                    
                    
                    reacted = 0;
                    if  no_holes
                        %check if growing or shrinking
                        grow= cell_maskp(vox_trial) & ~cell_mask(vox_trial);
                        shrink= ~cell_maskp(vox_trial) & cell_mask(vox_trial);
                    
                    
                        if grow
                            f=1;
                            dH_chem=B_rho*(RhoRatio(vox_ref)-rho_eq)-B_R*(RacRatio(vox_ref)-R_eq);
                    
                        elseif shrink
                            f=-1;
                            dH_chem=-B_rho*(RhoRatio(vox_trial)-rho_eq)+B_R*(RacRatio(vox_trial)-R_eq);
                    
                        end
                    
                    
                        if (grow || shrink) && rand<exp(-(dH+dH_chem+Hb)/T)
                            reacted=1;
                                cm0=cell_mask;
                            cell_mask=cell_maskp; %changing cell shape
                    
                            if shrink
                                bndry_up=cell_mask  & ~cell_mask(up);
                                bndry_down=cell_mask  & ~cell_mask(down);
                                bndry_l=cell_mask  & ~cell_mask(left);
                                bndry_r=cell_mask  & ~cell_mask(right);
                    
                                bndry_ud= bndry_up | bndry_down;
                                bndry_lr= bndry_l | bndry_r;
                    
                            end
                    
                                %recalculate parameters
                            Per=perim(cell_mask);
                            A=nnz(cell_mask);
                    
                    
                    
                            if grow
                                inds=cell_inds(1:A-1);
                                cell_inds(1:A)=find(cell_mask);
                            else
                                cell_inds(1:A)=find(cell_mask);
                                inds=cell_inds(1:A);
                            end
                    
                    
                    
                                if grow
                                    dist=max(abs(i0(vox_ref)-i0(inds)),abs(j0(vox_ref)-j0(inds)));
                                else
                                    dist=max(abs(i0(vox_trial)-i0(inds)),abs(j0(vox_trial)-j0(inds)));
                                end
                    
                                min_dist=10;
                                transport_mask=((D~=0).*D/min(D(D~=0))+(D==0).*prod(shape))*min_dist>dist;
                    
                                transport_mask(find(vox_trial==inds),:)=false;
                    
                            x0=x;
                    
                    
                                for i=1:length(D)
                    
                                    inds2=inds(transport_mask(:,i))+(i-1)*sz;
                                    i_trial=vox_trial+(i-1)*sz;
                    
                    
                                    sum0=sum(x(inds+(i-1)*sz));
                    
                    
                                    if grow
                                        us=x(vox_ref+(i-1)*sz);
                                        Ts=sum(x(inds2));
                                        f=Ts/(Ts+us);
                                        x(i_trial)=us;
                                        inds2=[inds2; i_trial];
                                    else
                                        ut=x(i_trial);
                                        f=1+(ut/sum(x(inds2)));
                                        x(i_trial)=0;
                                    end
                                    if ~isfinite(f)
                                        error("NAN");
                                    end
                                    x(inds2)=floor(x(inds2)*f)+[0; diff(floor(cumsum(rem(x(inds2)*f,1.0))+1e-5))]; %the 1e-5 is a fudge-factor to prevent underflow erros, they are typically of the order 1e-10 so the 1e-5 dominates
                    
                                    if ~grow
                                        sum1=sum(x(inds+(i-1)*sz));
                                        if sum1-sum0~=ut
                                            disp('check this out boss')
                                            ifix=jump(vox_trial,find(cell_mask(jump(vox_trial,:)),1))+(i-1)*sz;
                                            x(ifix)=x(ifix)-(sum1-sum0-ut);
                                        end
                                    else
                                        sum1=sum(x([inds; vox_trial]+(i-1)*sz));
                                        if sum1~=sum0
                                            error('we got a wild ass over here')
                                        end
                                    end
                    
                                end
                    
                    
                            I=[vox_trial vox_ref]; %places where molecule number has changed
                            H0=HA; %changing the hamiltonn to the new one
                    
                            %recalculate parameters
                            Per=perim(cell_mask);
                            A=nnz(cell_mask);
                            cell_inds(1:A)=find(cell_mask);
                    
                    
                            if grow
                                vox=cell_inds(1:A);
                            else
                                vox=[cell_inds(1:A); vox_trial];
                            end
                    
                            update_alpha_chem
                    
                            alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
                            if grow
                                disp('grow');
                                grow_count=grow_count+1;
                            else
                                disp('shrink');
                                shrink_count=shrink_count+1;
                            end
                    
                        end
                    end
                    
                    if ~reacted
                        %no move
                        cell_maskp=cell_mask;
                        Per=perim(cell_mask); % perimter
                        A=nnz(cell_mask); % area
                        cell_inds(1:A)=find(cell_mask);
                    end
                    
                    
                    
                    Ncell_maskp=squeeze(sum(sum(x)));
                    if any(Ncell_mask~=Ncell_maskp)
                        error('molecule loss')
                    end
                    
                    if min(cell_mask(:))<0
                        error('Oh no! D: (negtive numbers)')
                    end
                catch
                    time=Ttot;
                    break;
                end

                if sum(sum(sum(x(:,:,:),3)))-(totalRac+totalRho)~=d0
                    disp('molecules changed')
                end
            end

            diffusing_species=1:N_species; %only the first 6 species diffuse
            
            
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
            lastcpm=time;
            cpmcounter=cpmcounter+1;
        end

        if time>=lastplot+picstep || time==lastcpm % takes video frames
            if usejava('desktop') && isempty(getCurrentTask())
                % to make a video all frames must be the same size setting its position just
                % stops some bugs
                fs=14; %axis font size
            
                subplot(2,2,1)
                plotCellIm(panelA,double(cell_mask),cell_mask,i0,j0)
                colorbar;
            
                hold on
                try
                plot(center(2,1:iter),center(1,1:iter),'r')
                catch e
                    disp(e)
                end
                hold off
                ax = gca;
                ax.FontSize = fs;
                % title('Cell', 'Fontsize', 24)
                % xlabel('X')
                % ylabel('Y')
            
            
                subplot(2,2,2)
                plotCellIm(panelB,RhoRatio,cell_mask,i0,j0)
                caxis([0 1])
                colorbar;
                ax = gca;
                ax.FontSize = fs;
                title('Rho', 'Fontsize', 24)
                % xlabel('X')
                % ylabel('Y')
            
                subplot(2,2,3)
                plotCellIm(panelC,RacRatio,cell_mask,i0,j0)
                caxis([0 0.4])
                colorbar
                ax = gca;
                ax.FontSize = fs;
                set(gca,'Color',[1 1 1]*1)
                title('Rac', 'Fontsize', 24)
                axis tight
            
            
            
                 subplot(2,2,4)
                plotCellIm(panelD,PaxRatio,cell_mask,i0,j0)
                caxis([0 1])
                colorbar
                ax = gca;
                ax.FontSize = fs;
                set(gca,'Color',[1 1 1]*1)
                title('Pax', 'Fontsize', 24)
                axis tight
            
            
            
            
                % xlabel('X')
                % ylabel('Y')
                % xlabel('X')
                % ylabel('Y')
            
                %colormap jet
                %saveas(gcf,'CPM.png') %if you want a an image of a frame
            
            
                drawnow
                %adding videos the frame
            
                %frames=[frames getframe(gcf)];
                frame=getframe(gcf);
            end
            gif
            lastplot=time;
            if cpmcounter==cpmsteps
                iter=iter+1;
                
                center(:,iter)=com(cell_mask);
                Results(:,:,1,iter)=cell_mask;
                Results(:,:,2:(N_species+1),iter)=x; %storing the results
                Times(iter)=time;
                
                areas(iter)=A;
                perims(iter)=Per;
                
                Ham0(iter)=H0;
                Hchem(iter)=dH_chem;
                cpmcounter=0;
            end


        end

    end

    last_time=time;



    diffusing_species=1:N_species; %only the first 6 species diffuse
    
    
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
end
