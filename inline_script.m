function [str, deps, init] = inline_script(script, override)
    function bool = isScript(f)
        fid = fopen(f);
        try
            l1=fgetl(fid);
        catch e
            disp(e);
        end
        while isempty(l1)&&~isnumeric(l1)
            l1=fgetl(fid);
        end
        fclose(fid);
        if ~isnumeric(l1)
            
            m=regexp(l1,"[\s]?(function)[\s]+",'ONCE');
            bool=isempty(m);
        else
            bool=true;
        end
    end

if nargin==1
    override={};
end

if ~isScript(strcat(script,".m"))
    error('the file provided is not a script') ;
end



listing=dir('*.m');
mfiles= {listing.name}';
inds=cellfun(@(f) isScript(f), mfiles);
scripts=strrep(mfiles(inds),".m","");



source = fileread(strcat(script,".m"));
source=regexprep(source,"\n(%[^\n]*\n)+","\n"); %remove block comments
source=regexprep(source,"%[^\n]*\n","\n"); %remove end of line comments

scripts_match = strcat("(?<![\'",'\"])',scripts,'([^a-zA-Z_$0-9]|$)');%dont match script names inside of strings

[matched_scripts, starts ] = regexp(source,scripts_match,'match','start');
inds=~cellfun('isempty',matched_scripts);
matched_scripts = scripts(inds);
[~,foo]=sort(cellfun(@(x) x(1),starts(inds)));
matched_scripts = matched_scripts(foo);


script_match=strcat("[^\']?(",matched_scripts, ")(?:[ \f\t\r\n$]+|$)");
script_reps = {};

deps0=getVars(strcat(script,".m"));
deps={};
init0=getInitialized(strcat(script,".m"));
init=init0;
if ~isempty(matched_scripts)
    for s=matched_scripts'
        [rep,new_deps,new_init] = inline_script(s,override);
        script_reps{end+1}=[rep newline];
        new_deps=setdiff(new_deps,init);
        init=union(init,new_init);
%         [script ' : ' script_name ]
        deps=union(deps, new_deps);

    end
    special = {'\\a','\\b','\\f','\\r','\\t','\\v'}';
    special_rep = cellfun(@(x) ['\\' x],special,'UniformOutput',false);
    deps=unique([deps{:} setdiff(deps0,init) override{:}]);
    str=regexprep(source,script_match,regexprep(script_reps,special,special_rep));
else
    deps=setdiff(deps0,init);
    str=source;
end


end

