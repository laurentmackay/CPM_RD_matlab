function mi = Alg4(m,mi)

mp=sum(mi);
while mp>m
i=Alg2(mi/mp);
mi(i)=mi(i)-1;
mp=mp-1;
end

end

