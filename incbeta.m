function  y = incbeta(x,a,b)
    %numerical approximation to the incomplete beta function
    %MATLAB port of code found in:
    %https://github.com/codeplea/incbeta/blob/master/incbeta.c
    
    flip=false;
    %/*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
    if (x > (a+1.0)/(a+b+2.0)) 
        flip=true;
        a0=a;
        a=b;
        b=a0;
        x=1.0-x;
        %/*Use the fact that beta is symmetrical.*/
    end
    
    STOP=1.0e-3;
    TINY=1.0e-30;
    
    %*Find the first part before the continued fraction.*
     lbeta_ab = gammaln(a)+gammaln(b)-gammaln(a+b);
     front = exp(log(x)*a+log(1.0-x)*b-lbeta_ab) / a;

    %/*Use Lentz's algorithm to evaluate the continued fraction.*/
    f = 1.0; c = 1.0; d = 0.0;

    
    for i=0:200
%     for (i = 0; i <= 200; ++i) {
        m = i/2;


        if (i == 0) 
            numerator = 1.0; %/*First numerator is 1.0.*/
        elseif (mod(i,2) == 0) 
            numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m)); %/*Even term.*/
        else
            numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1)); %/*Odd term.*/
        end

        %/*Do an iteration of Lentz's algorithm.*/
        d = 1.0 + numerator * d;
        if (abs(d) < TINY) 
            d = TINY;
        end
        d = 1.0 / d;

        c = 1.0 + numerator / c;
        if (abs(c) < TINY)
            c = TINY; 
        end

         cd = c*d;
        f = f*cd;

        %/*Check for stop.*/
        if (abs(1.0-cd) < STOP) 
            y = front * (f-1.0);
            if flip
                y=1-y;
            end
            return;
        end
    end

    y=NaN;
end

