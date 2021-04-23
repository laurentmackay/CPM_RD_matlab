    N=1+1e3;
    sz=N;
    i=(1:N)';
    x=linspace(0,1,N)';
    h=(x(end)-x(1))/(N-1);
    jump=[circshift(i,-1) circshift(i,1)];
    v=-double(x>=0.3 & x<=0.7);
    
    lambda_tilde = sign(v).*sqrt(v(jump(:,1)).*v);%ok this is not really Roe's method, but that will require a dynamic stencil
%     lambda_plus=(lambda_tilde+abs(lambda_tilde))/2;
%     lambda_minus=(lambda_tilde-abs(lambda_tilde))/2;
    
    

    
    u_for = sparse(1:sz , jump(:,1)',1);


    
    f_for = u_for.*v(jump(:,1));
    f=spdiags(v(:),0,sz,sz);
    
    f_star = (f_for+f)/2-(sign(lambda_tilde).*(f_for-f))/2;
    f_star(v(jump(:,1))==0,:)=0;
    f_star(jump(v(jump(:,2))==0,2),:)=0;
    
    
    J_adv_2 = -(f_star - f_star(jump(:,2),:))/(2*h);
    
    
    u0 = (1+0.1*(rand(size(x))-0.5)).*abs(v);
    
    
    
    sum(J_adv_2*u0)
    
    
    u=u0;
    dt=0.001;
    figure(1);clf();
    t=0
    while t<100
        u=u+J_adv_2*u*dt;
        plot(x,u,x,u0);
        drawnow
        t=t+dt;
    end
    