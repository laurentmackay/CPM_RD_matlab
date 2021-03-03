function [str, deps, init] = inline_script(script, manual_deps, override, supp_path)
    function bool = isScript(f)
        fid = fopen(f);
        try
            l1=fgetl(fid);
        catch e
%             rethrow(e)
            disp(e);
        end
        if exist('l1','var')
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
        else
            bool=false;
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

if ~strcmp(script(end-1:end),'.m')
    script=strcat(script,'.m');
end

try

if ~isScript(script)
    error('the file provided is not a script') ;
end

catch err
    disp(err)
    end






listing=dir('*.m');
mfiles= {listing.name}';
inds=cellfun(@(f) isScript(f), mfiles);
scripts=strrep(mfiles(inds),".m","");

% p_scripts_match = strcat("(?<![\'",'\"_A-Za-z0-9])','([a-zA-Z][a-zA-Z_0-9]*)','([^a-zA-Z_$0-9]|$)');%dont match script names inside of strings

%append the file extenstion if needed

source = fileread(script);
source=regexprep(source,"\n(%[^\n]*\n)+","\n"); %remove block comments
source=regexprep(source,"%[^\n]*\n","\n"); %remove end of line comments


if ~isempty(override) %remove overrides
    override_refs=strcat('(?<=(?:^|\n|\;)[ \t\f]*)',override,'(?=[ \t\f]*\=)[^\n\;]+[\;]?');
    source=regexprep(source,override_refs,'');
end

% scripts_match = strcat("(?<![\'",'\"_A-Za-z0-9])',scripts,'([^a-zA-Z_$0-9]|$)');%dont match script names inside of strings

% [matched_scripts, starts ] = regexp(source,scripts_match,'match','start');


% p_scripts_match = strcat("(?<[^\'",'\"_A-Za-z0-9\*\-\+\/\^\(][ \t\f]*)','[a-zA-Z][a-zA-Z_0-9]*','(?=(?:[^a-zA-Z_$0-9\=\(\+\-\/\*\^]|$|\;|\r|\n))');%dont match script names inside of strings

% p_scripts_match = strcat("(?<=(?:[^\'",'\"_A-Za-z0-9\*\-\+\/\^\(\[]|^))','[ \t\f]*[a-zA-Z][a-zA-Z_0-9]*[ \t\f]*','(?=(?:$|\;|\r|\n))');%dont match script names inside of strings

scripts_match = strcat('(?<=(?:[\n\;]|^))','[ \t\f]*[a-zA-Z][a-zA-Z_0-9]*[ \t\f]*','(?=(?:$|\;|\r|\n))');%dont match script names inside of strings

[matched_scripts, starts ] = regexp(source,scripts_match,'match','start');
matched_scripts=strtrim(matched_scripts);
inds = cellfun(@(x) exist(x,'builtin'),matched_scripts)==0;
matched_scripts = matched_scripts(inds);
starts=starts(inds);

inds = cellfun(@(x) exist(x,'file'),matched_scripts)~=0;
matched_scripts = matched_scripts(inds);
starts=starts(inds);



ms0=matched_scripts;
matched_scripts = cellfun(@which,matched_scripts,'UniformOutput',0);

inds = cellfun(@isScript,matched_scripts)~=0;
matched_scripts = matched_scripts(inds);
ms0=ms0(inds);
starts=starts(inds);

% inds=~cellfun('isempty',matched_scripts);
% matched_scripts = scripts(inds);
if ~isempty(starts)
    [~,foo]=sort(starts);
    matched_scripts = matched_scripts(foo);
end

script_match=strcat("(?<![\'])(",ms0, ")(?:[ \f\t\r\n$]+|$)");
script_reps = {};

deps0=getVars(script);
deps={};
init0=getInitialized(script);
init=init0;
if ~isempty(matched_scripts)
    for s=matched_scripts
        [rep,new_deps,new_init] = inline_script(s{1},manual_deps, override);
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

