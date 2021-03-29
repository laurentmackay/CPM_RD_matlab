function y=find_crossings(x,thresh)

    y=find([0 diff(x>thresh)]);
    
end

