
N_species=8; %number of chemical species
N_rx=6; %number of reactions (this should determined automatically)
Ttot=2e1; %time the simulation end
SF=1; % speed factor I divide molecule number by this for speed
Gsize=10; %length of the grid in um
N=30; % number of points used to discretize the grid
shape=[N,N];
sz=prod(shape);
len=Gsize/N; %length of a lattice square
dx=len;
[j, i] = meshgrid(1:shape(2),1:shape(1)); %the i and j need to be reversed for some reason (\_(:0)_/)

<<<<<<< HEAD
div=0.2;%fraction of the cell in the induced state
=======
div=0.2;
>>>>>>> f97891902ca8a3bc3ad20acaf34593c39421e13f

initialize_cell_geometry
initialize_cellular_potts
initialize_chem 

D=[D 0 0];

ij0=(1:(sz))';
diffuse_mask=false(size(jump,2),sz);
num_diffuse=zeros(1,size(jump,2));
ij_diffuse=zeros(4,(N)*(N));

enumerate_diffusion

<<<<<<< HEAD
pmax=0.5;%max allowed pT for one cell

dt=pmax*(h^2)/(max(D)*size(diffuse_mask,1));


=======
pmax=0.5;

dt=pmax*(h^2)/(max(D)*size(diffuse_mask,1));

>>>>>>> f97891902ca8a3bc3ad20acaf34593c39421e13f
% dt=0.01;
imagesc(x(:,:,2))

jumpp=jump';
<<<<<<< HEAD
%pre-compute single molecule diffusion probabilities for the timestep of
%length dt
=======
>>>>>>> f97891902ca8a3bc3ad20acaf34593c39421e13f
pT = zeros(size(num_vox_diff,2),size(x,3));
pi = zeros(size(diffuse_mask,1),sz);
for i=1:A
    vox=cell_inds(i);
    pT(vox,:) = num_vox_diff(vox)*D*dt/(h^2);
    pi(:,vox)=diffuse_mask(:,vox)'./sum(diffuse_mask(:,vox));
end
max(max(pT))
<<<<<<< HEAD

=======
% pi=diffuse_mask/sum(diffuse_mask)
>>>>>>> f97891902ca8a3bc3ad20acaf34593c39421e13f
while true
tic;
x=Alg3_mex(x,dt,D,dx,jumpp,diffuse_mask,pT,pi,cell_inds,A,1);
toc

imagesc(x(:,:,2));
% sum(sum(x(:,:,2)))
colorbar
drawnow;
end