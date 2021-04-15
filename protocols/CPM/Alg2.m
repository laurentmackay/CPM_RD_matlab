function i = Alg2(p,xi,p0)
i=0;
b=p0;

while xi>=b
    i=i+1;
    try
    b=b+p(i);
    catch e
        disp(e)
    end
end

end