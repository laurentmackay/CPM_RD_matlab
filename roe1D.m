    N=1e1;
    sz=N;
    i=(1:N)';
    x=linspace(0,1,N)';
    h=(x(end)-x(1))/(N-1);
    jump=[circshift(i,-1) circshift(i,1)];
    v=double(x>0.2 & x<0.7);
    
    lambda_plus=(v(:)+abs(v(:)))/2;
    lambda_minus=(v(:)-abs(v(:)))/2;

    
    u_plus = sparse(1:sz , jump(:,1)',1);

    u_minus = sparse(1:sz , jump(:,2)',1);

    D_plus =  (u_plus-speye(sz));
    D_minus =  (u_minus-speye(sz));


    
    f_plus = u_plus.*v(jump(:,1));
    f_plus(v==0,:)=0;
    
    f_minus = u_minus.*v(jump(:,2));
    f_minus(v==0,:)=0;
    full(f_minus)

    f=spdiags(v(:),0,sz,sz);
    
    f_star_plus = (f_plus+f)/2-(D_plus.*sign(lambda_plus))/2;

    
    f_star_minus = (f_minus+f)/2+(D_minus.*sign(lambda_minus))/2;
%     f_star_minus = (f_minus+f)/2-(spdiags(lambda_minus,0,sz,sz)*D_plus)/2
    
    J_adv_2 = -(f_star_plus - f_star_minus)/(2*h);
    
    
    u0 = rand(size(x)).*v;
    
    sum(J_adv_2*u0)
    