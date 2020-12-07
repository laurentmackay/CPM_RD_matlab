
u = reshape(x,[sz,size(u,3)]);

interior=~bndrys&cell_mask(:);

inds=repmat((1:z)',N_dim);

v=cell2mat(arrayfun(@(i) 1j^(2*(i-1))*ones(1,nnz(interior(:,i))),1:N_dim,'UniformOutput',0));

u_x = sparse(inds(interior),jump(interior),v,sz,sz);


u_xx = (u_x(:,jump(:,1))-u_x(:,jump(:,2))+...
        u_x(:,jump(:,3))-u_x(:,jump(:,4)));
    
    t=0;
    
eval_Rx

Rx_prev=Rx;

while t<=T_final
    u = u + dt * neg_div_J*(D.*u) + 3*dt/2*Rx - dt/2*Rx_prev;
    t=t+dt;
    
    eval_Rx
    Rx_prev=Rx;
end
