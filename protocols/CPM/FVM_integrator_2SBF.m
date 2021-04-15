T_final=1e1;
dt=1e0;
t_plot=2;


% figure(3);clf();

% u = reshape(x,[sz,size(x,3)]);
% u = u()
u = x(cell_inds(1:A) + ((1:N_species)-1)*sz);

% u(i0(:)>40 & u(:,2)>0,2) = u(i0(:)>40 & u(:,2)>0,2)/2;
%
% u(i0(:)<40 & u(:,2)>0,2) = u(i0(:)<40 & u(:,2)>0,2) + sum( u(i0(:)>40 & u(:,2)>0,2))/nnz(u(i0(:)<40 & u(:,2)>0));

interior=~bndrys&cell_mask(:);

vox=repmat((1:sz)',1,N_dim);
vox=vox(interior)';

% ind_0(:) = [ind_0(interior(:,1:2)) ind_0(interior(:,3:4))

%%

dir=1:N_dim;

row = find(interior);
N_ind = length(row);
i = [row'; row'];
j=[ vox(:)';  jump(row)';];

Delta=repmat([-1/h; +1/h],N_ind,1);
u_x = sparse(i,j,Delta,numel(interior),sz);

Delta2 = repelem(repmat([1/h 1/h],1,N_dim/2),sum(interior));
v=nonzeros(Delta2.*u_x(row,:)');
i2 = repelem(vox,2);
j2 = mod(find(u_x(row,:)')-1,sz)+1;
u_xx = sparse(i2,j2,v,sz,sz);

u_xx=u_xx(cell_inds(1:A),cell_inds(1:A));

%%

%     u_xx = (u_x(jump(:,1),:)-u_x(jump(:,2),:)+...
%         u_x(jump(:,3),:)-u_x(jump(:,4),:))/h;

t=0;
t_last=0;
im=nan(shape);
eval_Rx
%  dt * u_xx*(D(1:N_slow).*u(:,1:N_slow))
tic;
eye = speye(A);
A_mat = arrayfun(@(Di)( eye*3/(2*dt)-Di*u_xx),D,'UniformOutput',0);
u_prev = u(:,1:N_slow);
while t<=T_final
    
    Rx_prev=Rx;
    eval_Rx
    b=(2*u(:,1:N_slow)-u_prev/2)/dt + 2*Rx-Rx_prev;
    u_prev=u(:,1:N_slow);
    for i = 1:N_slow
    
        u(:,i) = A_mat{i}\b(:,i);
    
    end
    project_fast
    
    
    t=t+dt;
    %
%     if t-t_last>=t_plot || t>=T_final
% %         t_plot/toc
%         im(cell_inds(1:A))=u(:,2)/Rac_Square;
%         imagesc(im); colorbar();
%         title(['time = ' num2str(t)])
%         drawnow;
%         tic
%         t_last=t;
%     end
end
x(cell_inds(1:A) + ((1:N_species)-1)*sz) = u(:);
u = reshape(x,[sz,size(x,3)]);
eval_Rx
