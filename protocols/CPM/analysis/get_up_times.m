function dt = up_times(x,thresh)
    inds = find_crossings(x,thresh);
    
    if x(1)>thresh
        dt=diff([0 inds]);
    else
        dt=diff(inds);
    end
end

