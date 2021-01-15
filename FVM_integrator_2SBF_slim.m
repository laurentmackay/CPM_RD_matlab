


u = x(cell_inds(1:A) + i_chem_0);

interior=~bndrys&cell_mask(:);

vox=repmat((1:sz)',1,N_dim);
vox=vox(interior)';


%%


row = find(interior);
N_ind = length(row);
i = [row'; row'];
j=[ vox(:)';  jump(row)';];

Delta=repmat([-1/h; +1/h],N_ind,1);
u_x = sparse(i,j,Delta,numel(interior),sz);

[i2,j2,v ] = find(u_x);
i2=mod(i2-1,sz)+1;
u_xx = sparse(i2,j2,v/h,sz,sz);

%system reduction
u_xx=u_xx(cell_inds(1:A),cell_inds(1:A));




eval_Rx
u_prev = u;

t0=time;

eye = speye(A);
MAT_list = arrayfun(@(Di)( eye*3/(2*dt)-Di*u_xx),D,'UniformOutput',0);



while time-t0<T_integration
    

    
    Rx_prev=Rx;
    
    eval_Rx
    b=(2*u-u_prev/2)/dt + 2*Rx-Rx_prev;
    u_prev=u;

    for i = 1:N_species
        u(:,i) = MAT_list{i}\b(:,i);
    end
  
    time=time+dt;


end
%expanding back in to the full system
x(cell_inds(1:A) + i_chem_0) = u(:);
u = reshape(x,[sz ,size(x,3)]);
eval_Rx %re-evaulate on full matrix
