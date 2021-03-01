function [str, deps, init] = inline_script(script, manual_deps, override)
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

if nargin<2
    manual_deps={};
end

if nargin<3
    override={};
elseif ~iscell(override)
    override=cellstr(override);
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


if ~isempty(override) %remove overrides
    override_refs=strcat('(?<=(?:^|\n|\;)[ \t\f]*)',override,'(?=[ \t\f]*\=)[^\n\;]+[\;]?');
    source=regexprep(source,override_refs,'');
end

scripts_match = strcat("(?<![\'",'\"_A-Za-z0-9])',scripts,'([^a-zA-Z_$0-9]|$)');%dont match script names inside of strings

[matched_scripts, starts ] = regexp(source,scripts_match,'match','start');
inds=~cellfun('isempty',matched_scripts);
matched_scripts = scripts(inds);
[~,foo]=sort(cellfun(@(x) x(1),starts(inds)));
matched_scripts = matched_scripts(foo);


script_match=strcat("(?<![\'])(",matched_scripts, ")(?:[ \f\t\r\n$]+|$)");
script_reps = {};

deps0=getVars(strcat(script,".m"));
deps={};
init0=getInitialized(strcat(script,".m"));
init=init0;
if ~isempty(matched_scripts)
    for s=matched_scripts'
        [rep,new_deps,new_init] = inline_script(s,manual_deps, override);
        script_reps{end+1}=[rep newline];
        new_deps=setdiff(new_deps,init,'stable');
        init=union(init,new_init);
%         [script ' : ' script_name ]
        deps=union(deps, new_deps);

    end
    special = {'\\a','\\b','\\f','\\r','\\t','\\v'}';
    special_rep = cellfun(@(x) ['\\' x],special,'UniformOutput',false);
    deps=unique([deps{:} setdiff(deps0,init,'stable') manual_deps{:} override{:}],'stable');
    str=regexprep(source,script_match,regexprep(script_reps,special,special_rep));
else
    deps=setdiff(deps0,init);
    str=source;
end




end

