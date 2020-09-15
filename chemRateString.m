function rate = chemRateString(chem,r)
if r>0
    if r~=1
        rate = ['.*' chem '.^' num2str(r)];
    else
        rate=['.*' chem];
    end
else
    rate='';
end
end

