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
R=20;
cell_mask=(i-N/2).^2 +(j-N/2).^2 < (R/len)^2;

Per=perim(cell_mask); %current lattice permiter
A=nnz(cell_mask); %current lattice area


cell_maskp=cell_mask; % intializing a variable for CPM_step.m
cell_inds=zeros(N*N,1);
cell_inds(1:A)=find(cell_mask);



