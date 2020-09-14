function [x,A,sz,...
    alpha_chem,PaxRatio,RhoRatio,K_is,K,RacRatio,I_Ks,reaction,ir0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    I,Rac_Square,Pax_Square,Rho_Square]=update(x,A,sz,...
    alpha_chem,PaxRatio,RhoRatio,K_is,K,RacRatio,I_Ks,reaction,ir0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    I,Rac_Square,Pax_Square,Rho_Square,gamma)


Rac_Square=sum(sum(sum(x(:,:,[4 2 7]))))/A;
Rho_Square=sum(sum(sum(x(:,:,[3 1]))))/A;
Pax_Square=sum(sum(sum(x(:,:,[5 6 8]))))/A;
for ind=[cell_inds(1:A)' I(1)]
    RacRatio(ind)=x(ind+(4-1)*sz)./Rac_Square;
    RhoRatio(ind)=x(ind+(3-1)*sz)./Rho_Square;
    PaxRatio(ind)=x(ind+(6-1)*sz)./Pax_Square;
    K_is(ind)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(ind)).*(1+alpha_R*RacRatio(ind))+k_G*k_X*GIT*PIX);
    K(ind)=alpha_R*RacRatio(ind).*K_is(ind).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(ind));%RbarRatio(I)/gamma;         %changed from paper
    I_Ks(ind)=I_K*(1-K_is(ind).*(1+alpha_R*RacRatio(ind)));
    reaction(ind+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(ind)+gamma*K(ind)).^m));            %From inactive rho to active rho changed from model
    reaction(ind+(2-1)*sz) = (I_R+I_Ks(ind)).*(L_rho^m./(L_rho^m+RhoRatio(ind).^m));                %From inactive Rac to active Rac
    reaction(ind+(5-1)*sz) = B_1*(K(ind).^m./(L_K^m+K(ind).^m));
    [tmp,tmp2]=meshgrid(ir0,ind);
    tmp3=tmp+tmp2;
    alpha_chem(tmp3) = reaction(tmp3).*x(tmp3);
end
end