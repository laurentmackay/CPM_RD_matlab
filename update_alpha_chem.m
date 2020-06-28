function [x,sz,alpha_rx,...
    alpha_chem,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,reaction,ir0,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    PAKtot,i]=update_alpha_chem(x,sz,alpha_rx,...
    alpha_chem,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,reaction,ir0,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    PAKtot,i,A)

%to properly locate in alpha_chem(ir0+I)
if length(i)>1
[tmp,tmp2]=meshgrid(ir0,i);
    I_rx=tmp+tmp2;
else
    I_rx=i+ir0;
end

a_c_0=alpha_chem(I_rx);

%update Ratios
RacRatio(i)=nan2zero(x(i+(4-1)*sz)./(x(i+(4-1)*sz)+x(i+(2-1)*sz)+x(i+(7-1)*sz)));
RbarRatio(i)=nan2zero(x(i+(7-1)*sz)./(x(i+(4-1)*sz)+x(i+(2-1)*sz)+x(i+(7-1)*sz)));
RhoRatio(i)=nan2zero(x(i+(3-1)*sz)./(x(i+(3-1)*sz)+x(i+(1-1)*sz)));
PaxRatio(i)=nan2zero(x(i+(6-1)*sz)./(x(i+(6-1)*sz)+x(i+(5-1)*sz)+x(i+(8-1)*sz)));



gamma=0.3;

% RacRatio(i)=nan2zero(x(i+(4-1)*sz)./1.434948979591837e+03);
% RbarRatio(i)=nan2zero(x(i+(7-1)*sz)./1.434948979591837e+03);
% RhoRatio(i)=nan2zero(x(i+(3-1)*sz)./1.434948979591837e+03);
% PaxRatio(i)=nan2zero(x(i+(6-1)*sz)./1.077806122448980e+03);

if sum([nnz(isnan(RhoRatio)), nnz(isnan(RacRatio)), nnz(isnan(PaxRatio))])~=0
    disp("woah")
end

% RacRatio(isnan(RacRatio))=0;
% RbarRatio(isnan(RbarRatio))=0;
% RhoRatio(isnan(RhoRatio))=0;
% PaxRatio(isnan(PaxRatio))=0;

%update other parameters    
K_is(i)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i)).*(1+alpha_R*RacRatio(i))+k_G*k_X*GIT*PIX);
K(i)=alpha_R*RacRatio(i).*K_is(i).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(i));%RbarRatio(I)/gamma;         %changed from paper
I_Ks(i)=I_K*(1-K_is(i).*(1+alpha_R*RacRatio(i)));

reaction(i+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(i)+gamma*K(i)).^m));            %From inactive rho to active rho changed from model
reaction(i+(2-1)*sz) = (I_R+I_Ks(i)).*(L_rho^m./(L_rho^m+RhoRatio(i).^m));                %From inactive Rac to active Rac
reaction(i+(5-1)*sz) = B_1*(K(i).^m./(L_K^m+K(i).^m));


% alpha_chem(I_rx) = reaction(I_rx).*x(I_rx);
% % alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
% alpha_rx=squeeze(sum(sum(alpha_chem)))';


%         ai20=alpha_chem(ir0+i2(1));

% 
alpha_chem(I_rx) = reaction(I_rx).*x(I_rx); %chemical reaction
alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);
