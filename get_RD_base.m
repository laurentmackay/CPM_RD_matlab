function base = get_RD_base()
global RD_base

if isempty(RD_base)
    p=mfilename('fullpath');
    RD_base = regexprep(p,[ '[^\' filesep ']+$'],'');
end

base=RD_base;

end

