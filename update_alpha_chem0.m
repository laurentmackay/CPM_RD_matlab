%to properly locate in alpha_chem(ir0+I)
if length(vox)>1
[tmp,tmp2]=meshgrid(ir0,cell_inds(1:A));
    I_rx=tmp+tmp2;
else
    I_rx=vox+ir0;
end

a_c_0=alpha_chem(I_rx);

%update Ratios
% RacRatio(vox)=nan2zero(x(vox+(4-1)*sz)./(x(vox+(4-1)*sz)+x(vox+(2-1)*sz)+x(vox+(7-1)*sz)));

% RhoRatio(vox)=nan2zero(x(vox+(3-1)*sz)./(x(vox+(3-1)*sz)+x(vox+(1-1)*sz)));
% PaxRatio(vox)=nan2zero(x(vox+(6-1)*sz)./(x(vox+(6-1)*sz)+x(vox+(5-1)*sz)+x(vox+(8-1)*sz)));

% RacRatio(vox)=x(vox+(4-1)*sz)./(x(vox+(4-1)*sz)+x(vox+(2-1)*sz)+x(vox+(7-1)*sz));
% RhoRatio(vox)=x(vox+(3-1)*sz)./(x(vox+(3-1)*sz)+x(vox+(1-1)*sz));
% PaxRatio(vox)=x(vox+(6-1)*sz)./(x(vox+(6-1)*sz)+x(vox+(5-1)*sz)+x(vox+(8-1)*sz));


RacRatio(vox)=x(vox+(4-1)*sz)./Rac_Square;
RhoRatio(vox)=x(vox+(3-1)*sz)./Rho_Square;
PaxRatio(vox)=x(vox+(6-1)*sz)./Pax_Square;


gamma=0.3;

% RacRatio(vox)=nan2zero(x(vox+(4-1)*sz)./1.434948979591837e+03);
% RbarRatio(vox)=nan2zero(x(vox+(7-1)*sz)./1.434948979591837e+03);
% RhoRatio(vox)=nan2zero(x(vox+(3-1)*sz)./1.434948979591837e+03);
% PaxRatio(vox)=nan2zero(x(vox+(6-1)*sz)./1.077806122448980e+03);

% if sum([nnz(isnan(RhoRatio)), nnz(isnan(RacRatio)), nnz(isnan(PaxRatio))])~=0
%     disp("woah")
% end

% RacRatio(isnan(RacRatio))=0;
% RbarRatio(isnan(RbarRatio))=0;
% RhoRatio(isnan(RhoRatio))=0;
% PaxRatio(isnan(PaxRatio))=0;

%update other parameters    
K_is(vox)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(vox)).*(1+alpha_R*RacRatio(vox))+k_G*k_X*GIT*PIX);
K(vox)=alpha_R*RacRatio(vox).*K_is(vox).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(vox));         %changed from paper
I_Ks(vox)=I_K*(1-K_is(vox).*(1+alpha_R*RacRatio(vox)));

reaction(vox+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(vox)+gamma*K(vox)).^m));            %From inactive rho to active rho changed from model
reaction(vox+(2-1)*sz) = (I_R+I_Ks(vox)).*(L_rho^m./(L_rho^m+RhoRatio(vox).^m));                %From inactive Rac to active Rac
reaction(vox+(5-1)*sz) = B_1*(K(vox).^m./(L_K^m+K(vox).^m));


% alpha_chem(I_rx) = reaction(I_rx).*x(I_rx);
% % alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
% alpha_rx=squeeze(sum(sum(alpha_chem)))';


%         ai20=alpha_chem(ir0+i2(1));

% 
alpha_chem(I_rx) = reaction(I_rx).*x(I_rx); %chemical reaction
alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);
