function dep = getScriptDep(f)
if isstring(f) || ischar(f)
    vars=getVars(f);
    local=getInitialized(f);
    [~,i]=setdiff(vars,local);
    dep=vars(sort(i));
else
    dep={};
    for i=1:length(f)
    dep=union(dep,getScriptDep(string(f{i})));
    end
end
end

