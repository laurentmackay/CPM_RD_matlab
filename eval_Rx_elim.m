f_E = -(a.*(cnsrv_1 - u(:,2)).*u(:,1))+ (d.*(((cnsrv_1 - u(:,2)).*u(:,1).*a)./(d + k)))+ (k.*(((cnsrv_1 - u(:,2)).*u(:,1).*a)./(d + k)));
f_S = -(a.*(cnsrv_1 - u(:,2)).*u(:,1))+ (d.*(((cnsrv_1 - u(:,2)).*u(:,1).*a)./(d + k)));
f_P = (k.*(((cnsrv_1 - u(:,2)).*u(:,1).*a)./(d + k)));
subs__0 = 1./(d + k);
subs__1 = (cnsrv_1 - u(:,2)).*a.*subs__0;
J_gamma_1=(u(:,1).*a)./(d + k);
J_gamma_2=((cnsrv_1 - u(:,2)).*a)./(d + k);
subs2__0 = (J_gamma_2.*f_S)./(J_gamma_1 + 1);
subs2__1 = -subs2__0;

Rx = [(u(:,1).*a.*k.*((((cnsrv_1 - u(:,2)).*u(:,1).*a)./(d + k)) - cnsrv_1))./(d + k),...
-(u(:,1).*a.*k.*((((cnsrv_1 - u(:,2)).*u(:,1).*a)./(d + k)) - cnsrv_1).*(d + k + u(:,1).*a - a.*((((cnsrv_1 - u(:,2)).*u(:,1).*a)./(d + k)) - cnsrv_1)))./((d + k).*(d + k + u(:,1).*a))];