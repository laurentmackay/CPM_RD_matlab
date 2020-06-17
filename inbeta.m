function [f B]=inbeta(z,a,b,N)
%F=INBETA(X,P,Q) In complete beta function for complex x
% there are branch cuts [-infinity,0] and [1 infinity] on the real axis
% f = int(y.^(a-1).*(1-y).^(b-1),0..x);
% B is the complete integral
% The calculation is based in Muir's continued fraction expansion
% Biometrika Vol. 22 284--297 (1930-31)
% discussed in
% L.A. Aroian Annals of Mathematical Statistics Vol. 12 218-223 (1941) 
% and correection vol 30 1265 (1959)
%
% for large $z$ we use a Taylor expansion in powers of 1/z 
%
% The number of terms N can be adjusted
%
% Copyright Jim McElwaine 2010
%
% Call without arguments to run a test with random a and b in the
% complex plane
% The maximum error in the derivative is returned as a test

B = beta(a,b);
if B~=0
    f=0;
    for k=N:-1:1
      f = k*(b-k)*z./((a+2*k-1).*(a+2*k).*(1+f));
      j=k-1;
      f = -(a+j)*(a+b+j)*z./((a+2*j).*(a+2*j+1).*(1+f));
    end
    f=z.^a.*(1-z).^b./(a.*(1+f));

    if z>(a+1)/(a+b+2)
        f = B-f;
    end
    f=f/B;
else
    f=0;
end
return;
