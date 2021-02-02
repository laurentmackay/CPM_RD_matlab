function out = breaklines(str,n,cont)
    function l=breakline(l)
        breaks = regexp(l,'[\+\-\*\/ ]');
        breaker = find(breaks<n_eff-1,1,'last');
        breaker = breaks(breaker);
        while ~isempty(breaker) && breaker < breaks(end)
            l = [l(1:breaker) cont newline l(breaker+1:end)];
            breaks = regexp(l,'[\+\-\*\/ ]');
            breaker = breaker + length(breaker)+length(newline);
            breaker = find((breaks-breaker)<n_eff-1,1,'last');
            breaker = breaks(breaker);
        end
    end
if nargin<2
    n=80;
end
if nargin<3
    cont='...';
end

n_eff=n-length(cont);

lines = textscan([str newline], '%s', 'Delimiter', newline );
lines=lines{1};
too_long = cellfun(@length,lines)' > n_eff;
for i=find(too_long)
    lines{i} = breakline(lines{i});
end

out = strjoin(lines, newline);
end

