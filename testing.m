function [RhoRatio,RacRatio,PaxRatio,B,ss1,ss2] = testing(RhoRatio,RacRatio,PaxRatio,B)
    %define parameter values
    I_R=0.003;
    I_K=0.009;
    I_rho=0.016;L_R=0.34;n=4;gamma=0.3;delta_rho=0.016;
    L_rho=0.34;delta_R=0.025;m=4;L_K=5.77;delta_P=0.0004;
    k_X=41.7;k_G=5.71;k_C=5;GIT=0.11;PIX=0.069;Paxtot=2.3;alpha_R=15;Rtot=7.5;
    PAKtot = gamma*Rtot;
    D_1=0.43;%inactive Rho/Rac
    D_2=0.02;%active
    D_3=0.03;%pax

    %auxiliary functions
    K_is=@(R,P) 1/((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*P)*(1+alpha_R*R)+k_G*k_X*GIT*PIX);
    K=@(R,P) alpha_R*R*K_is(R,P)*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*P);
    I_Ks=@(R,P) I_K*(1-K_is(R,P)*(1+alpha_R*R));
    P_i=@(R,P) 1-P*(1+k_G*k_X*k_C*GIT*PIX*PAKtot*K_is(R,P)*(1+alpha_R*R));

    %ODEs
    RhoDot = @(Rho,R,P) I_rho*(L_R^n/(L_R^n + (R+gamma*K(R,P))^n))*(1-Rho) - delta_rho*Rho;
    Rdot = @(Rho,R,P) (I_R+I_Ks(R,P))*(L_rho^n/(L_rho^n+Rho^n))*(1-R-gamma*K(R,P)) - delta_R*R;
    Pdot = @(Rho,R,P) B*(K(R,P)^m/(L_K^m+K(R,P)^m))*P_i(R,P)-delta_P*P;

    %system
    %Kathy =  @(T,x) [RhoDot(x(1),x(2),x(3)) Rdot(x(1),x(2),x(3)) Pdot(x(1),x(2),x(3))]';
    Kathy =  @(x) [RhoDot(x(1),x(2),x(3)) Rdot(x(1),x(2),x(3)) Pdot(x(1),x(2),x(3))];

    %initial guesses for optimization
    uninduced_0=[0.99,0.01,0.01];
    trytry=[0.3991 0.1 0.3410];
    try4=[0.3732 0.12, 0.2];
    induced_0=[0.01,0.99,0.99];
    %{
    [T,Y]=ode15s(Kathy,linspace(0,18000,100),try4);
    figure();
    plot(T,Y(:,1));
    %}
    %run optimization from the two initial guesses, presumably it will find the steady states
    [ss1,~,exitflag1]=fsolve(Kathy,uninduced_0,optimset('TolFun',1e-16));
    [ss2,~,exitflag2]=fsolve(Kathy,induced_0,optimset('TolFun',1e-16));
    %you will probably want a smart criterion to decide if these two steady
    %states are "different". fsolve gives an exit flag which is a good start,
    %and should maybe be combined with some citerion on distance between the two
    %steady states for added robustness.

    %use the optimized values to define a search range for the saddle
    prec=0.001; %this value should be small but not 1e-12 small
    ss10=ss1+prec*(ss2-ss1);
    ss20=ss2-prec*(ss2-ss1);
    bounds=sort([ss10;ss20]);

    %solve for the saddle
    [ss3,~,~,exitflag3]=lsqnonlin(@(x) (Kathy(x)-[0 0 0]),mean([ss10;ss20]),bounds(1,:),bounds(2,:),optimset('TolFun',1e-16,'TolPCG',0.001,'TolX',1e-12,'MaxFunEvals',1e3));
    RhoRatio = ss3(1);
    RacRatio = ss3(2);
    PaxRatio = ss3(3);
end