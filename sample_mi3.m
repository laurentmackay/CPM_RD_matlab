function Y = sample_mi(n,p)
if nargin < 3,
    if nargin == 2,
        m = 1;
    else
        error('You need to input at least two arguments.');
        return,
    end;
end;
if (length(n)~=1) | (fix(n) ~= n) | (n < 0),
   error('n must be a positive integer.');
   return,
end;
P = sum(p);
if P ~= 1,
    error('The sum of the input probabilities must be equal 1.')
    return,
end;
 
    o = ones(1,n);
    s = cumsum(p/P);
    r = rand(1,n);
    for j = 1:length(p);
        o = o + (r > s(j)); 
    end;
    for j = 1:length(p);
        X(1,j) = sum(o == j); 
    end;

Y = X./n;
return,
    
end

