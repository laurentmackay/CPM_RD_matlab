function y = zero2nan(x)
y=x;
y(x==0)=NaN;
end

