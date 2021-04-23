    N=1e1;
   
    i=(1:N)';
    L=1
    [j,i]=meshgrid(1:N,1:N);
     sz=numel(i);
     N_dim=ndims(i);
     shape=size(i);
     
    h=L/(N-1);
    if N_dim==2
        up=sub2ind(shape,circshift(i,-1),j);
        down=sub2ind(shape,circshift(i,1),j);
        left=sub2ind(shape,i,circshift(j,-1,2));
        right= sub2ind(shape,i,circshift(j,1,2));
        jump=[up(:) down(:)  left(:) right(:)];
    else

        left=sub2ind(shape,circshift(i,-1));
        right= sub2ind(shape,circshift(j,1));
        jump=[left(:) right(:)];
    end
    
    v=zeros([shape N_dim]);
    v(:,:,1)=double(i>0.2*N & i<=0.7*N & j>0.2*N & j<=0.7*N);
    v(:,:,2)=-double(i>0.2*N & i<=0.7*N & j>0.2*N & j<=0.7*N);
    v0=v;
    v=reshape(v,[sz N_dim]);
    
    jump_for = jump(:,1:2:2*N_dim-1);
    jump_back = jump(:,2:2:2*N_dim);
    
%     
%     lambda_tilde = sign(v).*sqrt(abs(v(jump_for)).*abs(v))
%     u_for=sparse(repmat((1:sz)',1,N_dim),jump_for+[0:(N_dim-1)]*sz,1)
    
    
    lambda_tilde = arrayfun(@(i) sign(v(:,i)).*sqrt(abs(v(jump_for(:,i),i))).*abs(v(:,i)),1:N_dim,'UniformOutput',false) ;
    u_for = arrayfun(@(i) sparse(1:sz,jump_for(:,i),1),1:N_dim,'UniformOutput',false);
    f_for = arrayfun(@(i) u_for{i}.*v(jump_for(:,i),i), 1:N_dim,'UniformOutput',false);
    f=arrayfun(@(i) spdiags(v((1:sz),i),0,sz,sz),1:N_dim,'UniformOutput',false);
    
    f_star = arrayfun(@(i) (f_for{i}+f{i})/2 - (sign(lambda_tilde{i}).*(f_for{i}-f{i}))/2 ,1:N_dim,'UniformOutput',0);
%     f_star=sparse(sz,sz);
    
    for k_=1:N_dim
        tmp=f_star{k_};
        tmp(v(jump_for(:,k_))==0,:)=0;
        tmp(jump_back(v(jump_back(:,k_))==0,k_),:)=0;
        f_star{k_}=tmp;
    end

%     f_star_minus = (f_minus+f)/2-(spdiags(lambda_minus,0,sz,sz)*D_plus)/2
    
    J_adv_2 = arrayfun(@(i) -(f_star{i} - f_star{i}(jump_back(:,i),:))/(2*h),1:N_dim,'UniformOutput',false);
    J_adv = J_adv_2{1}+J_adv_2{2}
    
    
    u0 = (1+0.1*(rand(size(i))-0.5)).*sign(abs(sum(abs(v0),3)));
    
    
    
    sum(J_adv*u0(:))
    
    
    u=u0(:);
    dt=1e-7;
    figure(1);clf();

    s=pcolor(u0);
    set(s,'EdgeColor','none');
    colorbar();
    
%     
    while true
        u=u+J_adv*u*dt;
            s=pcolor(reshape(u(:),shape));
            set(s,'EdgeColor','none');
            colorbar();
%         plot(x,u);
        drawnow
        
    end
%     