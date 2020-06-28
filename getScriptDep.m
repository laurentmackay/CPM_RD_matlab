function dep = getScriptDep(f)
if isstring(f) || ischar(f)
    vars=getVars(f);
    local=getInitialized(f);
    dep=setdiff(vars,local);
else
    dep={};
    for i=1:length(f)
    dep=union(dep,getScriptDep(string(f{i})));
    end
end
end

