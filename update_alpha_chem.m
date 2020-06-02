% function [x,A,sz,diffusing_species_sum,alpha_rx,...
%     alpha_chem,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
%     RacRatio,RbarRatio,I_Ks,N_instantaneous,reaction,ir0,id0,cell_inds,...
%     k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
%     alpha,PAKtot,rx,diffused,i,I,xi0,neg]=update(x,A,sz,diffusing_species_sum,alpha_rx,...
%     alpha_chem,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
%     RacRatio,RbarRatio,I_Ks,N_instantaneous,reaction,ir0,id0,cell_inds,...
%     k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
%     alpha,PAKtot,rx,diffused,i,I,xi0,neg)

%to properly locate in alpha_chem(ir0+I)
[tmp,tmp2]=meshgrid(ir0,I);
I_rx=tmp(:)+tmp2(:);

a_c_0=alpha_chem(I_rx);

%update Ratios
RacRatio(I)=nan2zero(x(I+(4-1)*sz)./(x(I+(4-1)*sz)+x(I+(2-1)*sz)+x(I+(7-1)*sz)));
RbarRatio(I)=nan2zero(x(I+(7-1)*sz)./(x(I+(4-1)*sz)+x(I+(2-1)*sz)+x(I+(7-1)*sz)));
RhoRatio(I)=nan2zero(x(I+(3-1)*sz)./(x(I+(3-1)*sz)+x(I+(1-1)*sz)));
PaxRatio(I)=nan2zero(x(I+(6-1)*sz)./(x(I+(6-1)*sz)+x(I+(5-1)*sz)+x(I+(8-1)*sz)));

if sum([nnz(isnan(RhoRatio)), nnz(isnan(RacRatio)), nnz(isnan(PaxRatio))])~=0
    disp("woah")
end

% RacRatio(isnan(RacRatio))=0;
% RbarRatio(isnan(RbarRatio))=0;
% RhoRatio(isnan(RhoRatio))=0;
% PaxRatio(isnan(PaxRatio))=0;

%update other parameters    
K_is(I)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(I)).*(1+alpha_R*RacRatio(I))+k_G*k_X*GIT*PIX);
K(I)=alpha_R*RacRatio(I).*K_is(I).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(I));%RbarRatio(I)/gamma;         %changed from paper
I_Ks(I)=I_K*(1-K_is(I).*(1+alpha_R*RacRatio(I)));

reaction(I+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(I)+RbarRatio(I)).^m));            %From inactive rho to active rho changed from model
reaction(I+(2-1)*sz) = (I_R+I_Ks(I)).*(L_rho^m./(L_rho^m+RhoRatio(I).^m));                %From inactive Rac to active Rac
reaction(I+(5-1)*sz) = B_1*(K(I).^m./(L_K^m+K(I).^m));


% alpha_chem(I_rx) = reaction(I_rx).*x(I_rx);
% alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));


%         ai20=alpha_chem(ir0+i2(1));

% 
alpha_chem(I_rx) = reaction(I_rx).*x(I_rx); %chemical reaction
alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0);
