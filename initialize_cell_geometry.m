


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

% setting the initial shape of the cell to a circle of radius r
R=0.2*N/2;
cell_mask=(i-N/2).^2 +(j-N/2).^2 < R^2;
% cell_mask=true(shape);
% cell_mask=(i>1&i<N)&(j>1&j<N);
induced_mask=cell_mask & (i-min(i(cell_mask)))<=2*div*(max(i(cell_mask))-min(i(cell_mask)));
% induced_mask(:)=0;
% induced_mask=(i<=div*N&j<=N/2)|(i>(1-div)*N&j>N/2);
% induced_mask=induced_mask&cell_mask;

i0=i;
j0=j;

Per=perim(cell_mask); %current lattice permiter
A=nnz(cell_mask); %current lattice area


cell_maskp=cell_mask; % intializing a variable for CPM_step.m
cell_inds=zeros(N*N,1);
cell_inds(1:A)=find(cell_mask);


detect_bndrys
