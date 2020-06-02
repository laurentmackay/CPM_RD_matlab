function y = nan2zero(x)
y=x;
y(isnan(x))=0;
end

