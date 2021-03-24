function [B,lam_p_0,save_dir,dt,copyNum,cpmstep0,model_name] = ...
main_FVM_polarization_fun(B,lam_p_0,save_dir,dt,copyNum,cpmstep0,model_name)


plotting=usejava('desktop') && isempty(getCurrentTask());
try
    inputname(1);
catch
    deploy_model(model_name);
end


if plotting 

    pic_fig=figure(1);clf();
panel1=subplot(2,2,1);
panel2=subplot(2,2,2);
panel3=subplot(2,2,3);
panel4=subplot(2,2,4);
end





nrx=1e5; 

noise=0.0005;

Ttot=2e5; 

SF=2; 
Gsize=80; 
N=150; 
shape=[N,N];
sz=prod(shape);
h=Gsize/(N-1); 
cpm_wait=5; 

vmax=3/60; 
picstep=5;
cpmsteps=5;


cpmstep=cpmstep0/cpmsteps;



[j, i] = meshgrid(1:shape(2),1:shape(1)); 

div=0.1;



restart=false;
tic

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

Per=perim(cell_mask); 
A=nnz(cell_mask); 


cell_maskp=cell_mask; 
cell_inds=zeros(N*N,1);
cell_inds(1:A)=find(cell_mask);


adj_empty = ~cell_mask(up) | ~cell_mask(down) | ~cell_mask(left) | ~cell_mask(right); 
adj_full = cell_mask(up) | cell_mask(down) | cell_mask(left) | cell_mask(right); 

bndry_cell = cell_mask & adj_empty;
bndry_empty = ~cell_mask & adj_full;
bndry_mask = bndry_cell | bndry_empty;
bndry = find( bndry_mask );

bndry_up=cell_mask  & ~cell_mask(up);
bndry_down=cell_mask  & ~cell_mask(down);
bndry_l=cell_mask  & ~cell_mask(left);
bndry_r=cell_mask  & ~cell_mask(right);

bndry_ud= bndry_up | bndry_down;
bndry_lr= bndry_l | bndry_r;

bndrys=[bndry_up(:) bndry_down(:) bndry_l(:) bndry_r(:)];

N_species = 6;
N_rx = 6;
D = [0.43         0.02         0.43         0.02          0.1           20];
N_slow = 6;
chems={'Raci','Rac','Rhoi','Rho','Paxi','Pax'};





N_species = 6;
N_rx = 6;
D = [0.43         0.02         0.43         0.02          0.1           20];
N_slow = 6;
chems={'Raci','Rac','Rhoi','Rho','Paxi','Pax'};





i_chem_0 = ((1:N_species)-1)*sz;






totalRho = 2250000/SF;
totalRac = 2250000/SF;
totalPax = 690000/SF;
Rho_Square = totalRho/(A);    
Rac_Square = totalRac/(A);    
Pax_Square = totalPax/(A);    

Rho_Square = 1;    
Rac_Square = 1;    
Pax_Square = 1;    

N_instantaneous=50;


I_rho=0.016000000;
L_rho=0.340000000;
delta_rho=0.016000000;
L_R=0.340000000;
I_R=0.003000000;
delta_R=0.025000000;
alpha_R=15.000000000;
delta_P=0.000400000;
I_K=0.009000000;
L_K=5.770000000;
k_X=41.700000000;
k_G=5.710000000;
k_C=5.000000000;
GIT=0.110000000;
PIX=0.069000000;
Paxtot=2.300000000;
n=4.000000000;
m=4.000000000;
alpha_PAK=0.300000000;
Rho_Square=1.000000000;
Rac_Square=1.000000000;
Pax_Square=1.000000000;
PAKtot=2.250000000;
Rtot=7.5;
cnsrv_1=Rac_Square;
cnsrv_2=Rac_Square;
cnsrv_3=Pax_Square;
VolCell=(0.5*10^-6)*(h*10^-6)^2; 
muM = 6.02214*10^23*VolCell*10^3*10^-6; 



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




RhoRatio_u = 0.8;
RacRatio_u = 0.045;

PaxRatio_u = 0.082;

RhoRatio_i = 0.02;
RacRatio_i = 0.35; 

if length(D)==9
    RacRatio_u = 0.35; 
    RacRatio_i = 0.85; 
    CRatio=[0.6; 0.5];
    
    RhoRatio_i = 0.05;
    RhoRatio_u = 0.4;
end

if length(D)==9
    RacRatio_u = 0.35; 
    RacRatio_i = 0.85; 
    CRatio=[0.6; 0.5];
    
    RhoRatio_i = 0.05;
    RhoRatio_u = 0.4;
end





PaxRatio_i = 0.2;

RhoRatio=[RhoRatio_u; RhoRatio_i];
RacRatio=[RacRatio_u; RacRatio_i];
PaxRatio=[PaxRatio_u; PaxRatio_i];


K_is=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio).*(1+alpha_R*RacRatio)+k_G*k_X*GIT*PIX); 
K=alpha_R*RacRatio.*K_is.*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio);
P_i=1-PaxRatio.*(1+k_G*k_X*k_C*GIT*PIX*PAKtot*K_is.*(1+alpha_R*RacRatio));

Rho0 = Rho_Square*RhoRatio;           
Rhoi0 = Rho_Square - Rho0;

Rac0 = Rac_Square*RacRatio;           


if length(D)~=9
    
    
    
    
    
    Pax0 = Pax_Square*PaxRatio;           
    if length(D)==6
        Paxi0=Pax_Square-Pax0;
    else
        Paxi0 = Pax_Square*P_i;        
    end
    
    K_PAK = (1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*Pax0/Pax_Square).*alpha_PAK*PAKtot.*K_is;
    
    if length(D)==8
        Raci0 = Rac_Square*(1-(1+K_PAK).*RacRatio);        
    else
        Raci0 = Rac_Square*(1-RacRatio);
    end
    
    RacPak0 = K_PAK .* Rac0;
    GPP0 = (k_G*k_X*k_C*GIT*PIX*K_is*PAKtot.*(1+alpha_R*Rac0/Rac_Square)) .*Pax0;
    
    RacPak0 =  Rac_Square - Rac0 - Raci0;
    GPP0 = Pax_Square - Pax0 - Paxi0;
    
    
    N0=[Raci0 Rac0 Rhoi0 Rho0 Paxi0 Pax0];
    
    N0=[Raci0 Rac0 Rhoi0 Rho0 Paxi0 Pax0 RacPak0 GPP0];
    
    
    
    
    
    
    
    
    fp=0; 
    eval('model_fp');
    rhs_anon=0; 
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
    if  nnz(induced_mask)==0
        N0(1,1:N_species)=fp;
        
        
        
        
    end
    
    
    
    
    if plotting && relax
        figure(3);clf();
        plot(T_vec,Y_vec);
        legend(chems)
        drawnow;
        
        yl=ylim();
        ylim([0 yl(end)]);
    end
else
    
    
    unscale=5e0;
    
    I_C=2.95*muM/unscale;
    I_R=0.5*muM/unscale;
    I_rho=3.3*muM/unscale;
    
    a1=1.15*muM;
    a2=1*muM;
    n=3;
    alpha=4.5/unscale;
    beta=0.3/unscale;
    d_C=1/unscale;
    d_R=1/unscale;
    d_rho=1/unscale;
    
    ff=0.4;
    
    R_b=3*muM;
    rho_b=1.25*muM;
    
    rescale=35;
    
    I_P1=10.5*muM/rescale;
    delta_P1=0.21;
    kPI5K=0.084;
    k21=0.021;
    
    kPI3K=0.00072*5;
    kPTEN=0.432/5;
    
    P3b=.2*muM;
    
    C0=CRatio*C_Square;
    Ci0=C_Square-C0;
    
    Raci0 = Rac_Square*(1-RacRatio);
    
    P10=[1; 1]*I_P1/delta_P1;
    
    P20=P10.*(kPI5K*(1+Rac0/R_b)/2)/k21;
    P30=P20.*(kPI3K*(1+Rac0/R_b)/2)./( kPTEN*(1+Rho0/rho_b)/2);
    
    N0=[Raci0 Rac0 Ci0 C0 Rhoi0 Rho0 P10 P20 P30];
end




mask=induced_mask&cell_mask;
[tmp,tmp2]=meshgrid((0:N_species-1)*sz,find(mask));
i_induced=tmp+tmp2;

x=zeros([shape ,N_species]); 

x(i_induced)=repmat(N0(2,1:N_species),nnz(mask),1);

mask=~induced_mask&cell_mask;
[tmp,tmp2]=meshgrid((0:N_species-1)*sz,find(mask));
x(tmp+tmp2)=repmat(N0(1,1:N_species),nnz(mask),1);

x=(1+noise*rand(size(x))).*x;
x(x<0)=0;


alpha_chem=zeros([shape N_rx]);
alpha_rx=zeros(1,N_rx);
alpha_diff=zeros(6,1); 
ir0=((1:N_rx)-1)*sz;

vox=cell_inds(1:A);

RacRatio0 = x(:,:,2) ./ Rac_Square;
RhoRatio = x(:,:,4) ./ Rho_Square;
PaxRatio = x(:,:,6) ./ Pax_Square;
K_is=1.0./((1.0+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
I_Ks=I_K.*(1.0-K_is.*(1+alpha_R.*RacRatio0));
P_i=1.0-PaxRatio.*(1+k_G.*k_X.*k_C.*GIT.*PIX.*PAKtot.*K_is.*(1+alpha_R.*RacRatio0));
RacRatio = (x(:,:,2) + alpha_PAK.*K) ./ Rac_Square;
Q_R = ((1-RacRatio)./(1-RacRatio0)).*(I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio).^m));
Q_P = (P_i./(1-PaxRatio)).*B.*(K.^n./(L_K.^n+K.^n));




lam_a=1*h^4; 

lam_p=lam_p_0*h^2; 
J=0*h; 

B_0=1.5;
B_rho=(B_0/0.3)*h^2;
B_R=(B_0/0.3)*(.18/.13)*h^2; 


a=A; 
per=Per; 
Hb=0; 
T=0.3; 






H0=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; 
dH_chem=0; 

grow_count=0;
shrink_count=0;
end



lastplot=0;
lastcpm=0;


iter=0;

time=0;
reactions=0;

Nsteps=floor(Ttot/min(cpmstep0*cpm_wait))+1;


center=zeros(2,Nsteps);
Results=zeros([shape,N_species+1,Nsteps]);
Times=zeros(1,Nsteps);

areas=zeros(1,Nsteps);
perims=zeros(1,Nsteps);

Ham0=zeros(1,Nsteps);
Hchem=zeros(1,Nsteps);


iter=iter+1;

center(:,iter)=com(cell_mask);
Results(:,:,1,iter)=cell_mask;
Results(:,:,2:end,iter)=x; 
Times(iter)=time;

areas(iter)=A;
perims(iter)=Per;

Ham0(iter)=H0;
Hchem(iter)=dH_chem;


u = reshape(x,[sz ,size(x,3)]);
if plotting

tp__0=tic;

plotCellIm(panel1,reshape(RacRatio0,shape),cell_mask,i0,j0);
caxis(panel1,'auto');
colorbar(panel1);
title(panel1,'RacRatio0', 'Fontsize', 24);

plotCellIm(panel2,reshape(RacRatio,shape),cell_mask,i0,j0);
caxis(panel2,'auto');
colorbar(panel2);
title(panel2,'RacRatio', 'Fontsize', 24);

plotCellIm(panel3,reshape(RhoRatio,shape),cell_mask,i0,j0);
caxis(panel3,'auto');
colorbar(panel3);
title(panel3,'RhoRatio', 'Fontsize', 24);

plotCellIm(panel4,reshape(PaxRatio,shape),cell_mask,i0,j0);
caxis(panel4,'auto');
colorbar(panel4);
title(panel4,'PaxRatio', 'Fontsize', 24);

sgtitle(pic_fig,['t=' num2str(time) ', t_{plot}=' num2str(double(tic-tp__0)*1e-6), ', t_{sim}=' num2str(toc)], 'Fontsize', 10,'FontWeight','bold')

end
if plotting && usejava('desktop') && isempty(getCurrentTask())
    delete test.gif
    gif('test.gif','frame',panel1)
end

time=0;






N_dim=size(jump,2);
ij0=(1:(sz))';
diffuse_mask=false(N_dim,sz);
num_diffuse=zeros(1,size(jump,2));
ij_diffuse=zeros(4,(N)*(N));
diffusing_species_sum=zeros(N_dim,length(D));
num_vox_diff=zeros(1,sz);
pT0 = zeros(sz,length(D));
pi = zeros(N_dim,sz);

diffusing_species=1:N_species; 


for drx=1:size(jump,2) 
    diffuse_mask(drx,:)=cell_mask(jump(:,drx))&cell_mask(:); 
    
end

for vox=1:size(diffuse_mask,2)
    num_vox_diff(vox)=nnz(diffuse_mask(:,vox)); 
end

for i=1:A
    vox=cell_inds(i);
    pT0(vox,:) = num_vox_diff(vox)*D/(h^2);
    pi(:,vox)=diffuse_mask(:,vox)'./sum(diffuse_mask(:,vox));
end
id0=(diffusing_species-1)*sz;


alpha_diff=sum(diffusing_species_sum).*D/(h*h);
alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
alpha_rx2=alpha_rx;
a_total=sum(alpha_diff)+sum(alpha_rx(:));


numDiff=0;
numReac=0;
cpmcounter=0;
Timeseries=[];
TRac=[];
TRho=[];
TPax=[];

last_time=time; 
rx_speedup=2;
rx_count=zeros(shape);
dt_diff=zeros(size(D));
P_diff=0.5;


d0=sum(x(:));



if isempty(getCurrentTask());  end



T_integration = cpmstep;
keep_running=true;

while time<Ttot && keep_running
    A=nnz(cell_mask); 
    cell_inds(1:A)=find(cell_mask); 
    
    while (time-last_time)<Ttot

       
        


u = x(cell_inds(1:A) + i_chem_0);

interior=~bndrys&cell_mask(:);

vox=repmat((1:sz)',1,N_dim);
vox=vox(interior)';




row = find(interior);
N_ind = length(row);
i = [row'; row'];
j=[ vox(:)';  jump(row)';];

Delta=repmat([-1/h; +1/h],N_ind,1);
u_x = sparse(i,j,Delta,numel(interior),sz);

[i2,j2,v ] = find(u_x);
i2=mod(i2-1,sz)+1;
u_xx = sparse(i2,j2,v/h,sz,sz);

u_xx=u_xx(cell_inds(1:A),cell_inds(1:A));




RacRatio0 = u(:,2) ./ Rac_Square;
RhoRatio = u(:,4) ./ Rho_Square;
PaxRatio = u(:,6) ./ Pax_Square;
K_is=1.0./((1.0+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
I_Ks=I_K.*(1.0-K_is.*(1+alpha_R.*RacRatio0));
P_i=1.0-PaxRatio.*(1+k_G.*k_X.*k_C.*GIT.*PIX.*PAKtot.*K_is.*(1+alpha_R.*RacRatio0));
RacRatio = (u(:,2) + alpha_PAK.*K) ./ Rac_Square;
Q_R = ((1-RacRatio)./(1-RacRatio0)).*(I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio).^m));
Q_P = (P_i./(1-PaxRatio)).*B.*(K.^n./(L_K.^n+K.^n));
f_Raci = -(Q_R.*u(:,1))+ (delta_R.*u(:,2));
f_Rhoi = -(Q_rho.*u(:,3))+ (delta_rho.*u(:,4));
f_Paxi = -(Q_P.*u(:,5))+ (delta_P.*u(:,6));

Rx = [f_Raci,...
-f_Raci,...
f_Rhoi,...
-f_Rhoi,...
f_Paxi,...
-f_Paxi];
u_prev = u;

t0=time;

eye = speye(A);
MAT_list = arrayfun(@(Di)( eye*3/(2*dt)-Di*u_xx),D,'UniformOutput',0); 
if any(u(:)<0)
    disp('negatory pig pen')
end

while time-t0<T_integration
    

    
    Rx_prev=Rx;
    RacRatio0 = u(:,2) ./ Rac_Square;
RhoRatio = u(:,4) ./ Rho_Square;
PaxRatio = u(:,6) ./ Pax_Square;
K_is=1.0./((1.0+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
I_Ks=I_K.*(1.0-K_is.*(1+alpha_R.*RacRatio0));
P_i=1.0-PaxRatio.*(1+k_G.*k_X.*k_C.*GIT.*PIX.*PAKtot.*K_is.*(1+alpha_R.*RacRatio0));
RacRatio = (u(:,2) + alpha_PAK.*K) ./ Rac_Square;
Q_R = ((1-RacRatio)./(1-RacRatio0)).*(I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio).^m));
Q_P = (P_i./(1-PaxRatio)).*B.*(K.^n./(L_K.^n+K.^n));
f_Raci = -(Q_R.*u(:,1))+ (delta_R.*u(:,2));
f_Rhoi = -(Q_rho.*u(:,3))+ (delta_rho.*u(:,4));
f_Paxi = -(Q_P.*u(:,5))+ (delta_P.*u(:,6));

Rx = [f_Raci,...
-f_Raci,...
f_Rhoi,...
-f_Rhoi,...
f_Paxi,...
-f_Paxi];
b=(2*u-u_prev/2)/dt + (2*Rx-Rx_prev); 
    u_prev=u;

    for i = 1:N_species
        u(:,i) = MAT_list{i}\b(:,i);
        
        if any(u(:)<0)
            disp('wild ass over here')
        end
    
    end
    

    time=time+dt;


end
x(cell_inds(1:A) + i_chem_0) = u(:);
u = reshape(x,[sz ,size(x,3)]);
RacRatio0 = u(:,2) ./ Rac_Square;
RhoRatio = u(:,4) ./ Rho_Square;
PaxRatio = u(:,6) ./ Pax_Square;
K_is=1.0./((1.0+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
I_Ks=I_K.*(1.0-K_is.*(1+alpha_R.*RacRatio0));
P_i=1.0-PaxRatio.*(1+k_G.*k_X.*k_C.*GIT.*PIX.*PAKtot.*K_is.*(1+alpha_R.*RacRatio0));
RacRatio = (u(:,2) + alpha_PAK.*K) ./ Rac_Square;
Q_R = ((1-RacRatio)./(1-RacRatio0)).*(I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio).^m));
Q_P = (P_i./(1-PaxRatio)).*B.*(K.^n./(L_K.^n+K.^n));
f_Raci = -(Q_R.*u(:,1))+ (delta_R.*u(:,2));
f_Rhoi = -(Q_rho.*u(:,3))+ (delta_rho.*u(:,4));
f_Paxi = -(Q_P.*u(:,5))+ (delta_P.*u(:,6));

Rx = [f_Raci,...
-f_Raci,...
f_Rhoi,...
-f_Rhoi,...
f_Paxi,...
-f_Paxi];

if time>=lastcpm+cpmstep
            
            for kk=1:(2*Per)/cpmsteps 
                try
                    if all(isfinite([lam_a,lam_p]))
    
    adj_empty = ~cell_mask(up) | ~cell_mask(down) | ~cell_mask(left) | ~cell_mask(right); 
adj_full = cell_mask(up) | cell_mask(down) | cell_mask(left) | cell_mask(right); 

bndry_cell = cell_mask & adj_empty;
bndry_empty = ~cell_mask & adj_full;
bndry_mask = bndry_cell | bndry_empty;
bndry = find( bndry_mask );

bndry_up=cell_mask  & ~cell_mask(up);
bndry_down=cell_mask  & ~cell_mask(down);
bndry_l=cell_mask  & ~cell_mask(left);
bndry_r=cell_mask  & ~cell_mask(right);

bndry_ud= bndry_up | bndry_down;
bndry_lr= bndry_l | bndry_r;

bndrys=[bndry_up(:) bndry_down(:) bndry_l(:) bndry_r(:)];
is_discrete = all(mod(x(cell_inds(1:A)),1)==0);
    
    
    if any(cell_maskp~=cell_mask)
        error('not reseting')
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    rho_eq=mean(RhoRatio(find(cell_mask)));
    R_eq=mean(RacRatio(find(cell_mask)));
    Ncell_mask=squeeze(sum(sum(x))); 
    A0=A;
    
    no_holes=false;
    
    while ~no_holes
        vox_trial = bndry(randi(length(bndry)));
        
        r=randi(size(jump,2));
        vox_ref=jump(sub2ind([sz,4],vox_trial,r));
        cell_maskp(vox_trial) = cell_mask(vox_ref);
        
        Per=perim(cell_maskp); 
        A=nnz(cell_maskp); 
        HA=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; 
        dH=HA-H0;
        no_holes = getfield(bwconncomp(cell_maskp,4),'NumObjects')==1 && getfield(bwconncomp(~cell_maskp,4),'NumObjects')==1 ;
        if ~no_holes
            cell_maskp(vox_trial)=cell_mask(vox_trial);
        end
    end
    
    
    reacted = 0;
    if  no_holes
        
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
            
            
            
            cell_mask=cell_maskp; 
            
            if shrink
                bndry_up=cell_mask  & ~cell_mask(up);
                bndry_down=cell_mask  & ~cell_mask(down);
                bndry_l=cell_mask  & ~cell_mask(left);
                bndry_r=cell_mask  & ~cell_mask(right);
                
                bndry_ud= bndry_up | bndry_down;
                bndry_lr= bndry_l | bndry_r;
                
            end
            
            
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
            
            min_dist=5600;
            
            
            
            
            
            
            x0=x;
            inds2=inds+i_chem_0;
            i_trial=vox_trial+i_chem_0;
            
            
            if grow
                us=x(vox_ref+i_chem_0);
                Ts=sum(x(inds2));
                f=Ts./(Ts+us);
                x(i_trial)=us;
                inds2=[inds2; i_trial];
            else
                ut=x(i_trial);
                f=1+(ut./sum(x(inds2)));
                x(i_trial)=0;
            end
            
            if is_discrete
                x(inds2)=floor(x(inds2).*f)+[zeros(1,N_species); diff(floor(cumsum(rem(x(inds2).*f,1.0))+1e-5))]; 
            else
                x(inds2)=x(inds2).*f;
            end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            I=[vox_trial vox_ref]; 
            H0=HA; 
            
            
            Per=perim(cell_mask);
            A=nnz(cell_mask);
            cell_inds(1:A)=find(cell_mask);
            
            
            if grow
                vox=cell_inds(1:A);
            else
                vox=[cell_inds(1:A); vox_trial];
            end
            
            
            
            alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
            if grow
                
                grow_count=grow_count+1;
            else
                
                shrink_count=shrink_count+1;
            end
            
        end
    end
    
    if ~reacted
        
        cell_maskp=cell_mask;
        Per=perim(cell_mask); 
        A=nnz(cell_mask); 
        cell_inds(1:A)=find(cell_mask);
    else
        RacRatio0 = x(:,:,2) ./ Rac_Square;
RhoRatio = x(:,:,4) ./ Rho_Square;
PaxRatio = x(:,:,6) ./ Pax_Square;
K_is=1.0./((1.0+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
I_Ks=I_K.*(1.0-K_is.*(1+alpha_R.*RacRatio0));
P_i=1.0-PaxRatio.*(1+k_G.*k_X.*k_C.*GIT.*PIX.*PAKtot.*K_is.*(1+alpha_R.*RacRatio0));
RacRatio = (x(:,:,2) + alpha_PAK.*K) ./ Rac_Square;
Q_R = ((1-RacRatio)./(1-RacRatio0)).*(I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio).^m));
Q_P = (P_i./(1-PaxRatio)).*B.*(K.^n./(L_K.^n+K.^n));
end
    
    
    
    Ncell_maskp=squeeze(sum(sum(x)));
    
    if (is_discrete & any(Ncell_mask~=Ncell_maskp)) | (~is_discrete & any(abs(Ncell_mask-Ncell_maskp)>1e-5))
        error('molecule loss')
    end
    
    if min(cell_mask(:))<0
        error('Oh no! D: (negtive numbers)')
    end
end
catch err
                    rethrow(err)
                    break;
                end
                
                adj_empty = ~cell_mask(up) | ~cell_mask(down) | ~cell_mask(left) | ~cell_mask(right); 
adj_full = cell_mask(up) | cell_mask(down) | cell_mask(left) | cell_mask(right); 

bndry_cell = cell_mask & adj_empty;
bndry_empty = ~cell_mask & adj_full;
bndry_mask = bndry_cell | bndry_empty;
bndry = find( bndry_mask );

bndry_up=cell_mask  & ~cell_mask(up);
bndry_down=cell_mask  & ~cell_mask(down);
bndry_l=cell_mask  & ~cell_mask(left);
bndry_r=cell_mask  & ~cell_mask(right);

bndry_ud= bndry_up | bndry_down;
bndry_lr= bndry_l | bndry_r;

bndrys=[bndry_up(:) bndry_down(:) bndry_l(:) bndry_r(:)];
end
            
            diffusing_species=1:N_species; 


for drx=1:size(jump,2) 
    diffuse_mask(drx,:)=cell_mask(jump(:,drx))&cell_mask(:); 
    
end

for vox=1:size(diffuse_mask,2)
    num_vox_diff(vox)=nnz(diffuse_mask(:,vox)); 
end

for i=1:A
    vox=cell_inds(i);
    pT0(vox,:) = num_vox_diff(vox)*D/(h^2);
    pi(:,vox)=diffuse_mask(:,vox)'./sum(diffuse_mask(:,vox));
end
lastcpm=time;
            cpmcounter=cpmcounter+1;
        end
                
        if time>=lastplot+picstep || time==lastcpm 
 
            if cpmcounter==cpmsteps*cpm_wait
              u = reshape(x,[sz ,size(x,3)]);
              if plotting

tp__0=tic;

plotCellIm(panel1,reshape(RacRatio0,shape),cell_mask,i0,j0);
caxis(panel1,'auto');
colorbar(panel1);
title(panel1,'RacRatio0', 'Fontsize', 24);

plotCellIm(panel2,reshape(RacRatio,shape),cell_mask,i0,j0);
caxis(panel2,'auto');
colorbar(panel2);
title(panel2,'RacRatio', 'Fontsize', 24);

plotCellIm(panel3,reshape(RhoRatio,shape),cell_mask,i0,j0);
caxis(panel3,'auto');
colorbar(panel3);
title(panel3,'RhoRatio', 'Fontsize', 24);

plotCellIm(panel4,reshape(PaxRatio,shape),cell_mask,i0,j0);
caxis(panel4,'auto');
colorbar(panel4);
title(panel4,'PaxRatio', 'Fontsize', 24);

sgtitle(pic_fig,['t=' num2str(time) ', t_{plot}=' num2str(double(tic-tp__0)*1e-6), ', t_{sim}=' num2str(toc)], 'Fontsize', 10,'FontWeight','bold')

end
lastplot=time; 
            
                i_rac = find(strcmp(chems,'Rac')); 
                inds=cell_inds(1:A)+sz*(i_rac-1);
                d_Rac=max(max(x(inds)))-min(min(x(inds)));
                if plotting
                    gif
                end
                if ~isempty(getCurrentTask())
                    disp([num2str(copyNum) ': B=' num2str(B) ', t=' num2str(time) ', delta_Rac=' num2str(d_Rac)])
                end
                
                iter=iter+1;

center(:,iter)=com(cell_mask);
Results(:,:,1,iter)=cell_mask;
Results(:,:,2:end,iter)=x; 
Times(iter)=time;

areas(iter)=A;
perims(iter)=Per;

Ham0(iter)=H0;
Hchem(iter)=dH_chem;

cpmcounter=0;
                
                i_rac = find(strcmp(chems,'Rac')); 
                inds=cell_inds(1:A)+sz*(i_rac-1);
                if d_Rac>0.05
                   keep_running=false;
                   break;
                end
                
            end
            
            
        end
        
    end
    
    last_time=time;
    
    
    
    diffusing_species=1:N_species; 


for drx=1:size(jump,2) 
    diffuse_mask(drx,:)=cell_mask(jump(:,drx))&cell_mask(:); 
    
end

for vox=1:size(diffuse_mask,2)
    num_vox_diff(vox)=nnz(diffuse_mask(:,vox)); 
end

for i=1:A
    vox=cell_inds(i);
    pT0(vox,:) = num_vox_diff(vox)*D/(h^2);
    pi(:,vox)=diffuse_mask(:,vox)'./sum(diffuse_mask(:,vox));
end
end

toc

    fn=strcat(save_dir, 'final_B_', num2str(B), '_copy', int2str(copyNum), '.mat');
    disp(['saving to: ' fn]);
    ls results
    save(fn,'-v7.3');


end
