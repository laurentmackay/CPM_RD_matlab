function fn=mk_fun(script,varargin)
%varargin should be strings denoting any function names in that are in the matlab path
%which should be considered as local variables
%
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

suffix='_fun';
fn=[script suffix '.m'];

if ~isScript(strcat(script,".m"))
   error('the file provided is not a script') ;
end
if nargin==2 %remove any assignments to the rep variable
    clean = @(x) regexprep(x,{'%','\\','\t',[ varargin{1} '[ ]?=[^\n]+']},{'%%','\\\\','\\t',''});
else
    clean = @(x) regexprep(x,{'%','\\','\t'},{'%%','\\\\','\\t'});
end
listing=dir('*.m');
mfiles= {listing.name}';
inds=cellfun(@(f) isScript(f), mfiles);
scripts=strrep(mfiles(inds),".m","");
matches=strjoin(scripts,"|");
reg=strcat("[\s]?(" ,matches, ")[\s]?");

source = fileread(strcat(script,".m"));
source=regexprep(source,"%[^\n]*\n","\n"); %remove block comments

[matched_scripts, starts ] = regexp(source,strcat(scripts,'[^a-zA-Z_$0-9]'),'match','start');
inds=~cellfun('isempty',matched_scripts);
matched_scripts = scripts(inds);
[~,foo]=sort(cellfun(@(x) x(1),starts(inds)));
matched_scripts = matched_scripts(foo);

deps0=getVars(strcat(script,".m"));
deps={};
init0=getInitialized(strcat(script,".m"));
init=init0;
for s=matched_scripts'
    script_name=strcat(s,".m");
    new_deps=getScriptDep(script_name);
    new_deps=setdiff(new_deps,init);
    init=union(init,getInitialized(script_name));
    deps=union(deps, new_deps);
    
end
deps=unique([deps{:} setdiff(deps0,init) varargin{:}]);

maxwidth=80;

header="function [";
line_length=strlength(header);
for d=deps
    if line_length+length(d)+1<=maxwidth
        header=strcat(header,d,",");
        line_length=line_length+strlength(d)+1;
    else
        header=strcat(header,"...\n",d,",");
        line_length=strlength(d)+1;
    end
end

header=strcat(extractBetween(header,1,strlength(header)-1),"] = ");
line_length=line_length+3;
add=strcat(script,suffix,"(");

if line_length+length(add)<=maxwidth+3
    header=strcat(header,add);
    line_length=line_length+strlength(add);
else
    header=strcat(header,"...\n",add);
    line_length=strlength(add);
end

for d=deps
    if line_length+length(d)+1<=maxwidth+3
        header=strcat(header,d,",");
        line_length=line_length+strlength(d)+1;
    else
        header=strcat(header,"...\n",d,",");
        line_length=strlength(d)+1;
    end
end
header=strcat(extractBetween(header,1,strlength(header)-1),")\n");

in = fopen(strcat(script,".m"));



out = fopen(fn,'wt');

fprintf(out,header);

line=fgetl(in);
script_match=strcat("[^\']?(",matched_scripts, ")(?:[ \f\t\r\n$]+|$)");

while ~isnumeric(line)
    %
    [match_script,match_start]=regexp(line,script_match,'tokens','start');
    i_match=~cellfun('isempty',match_script);
    match_script=match_script(i_match);
    match_start=match_start(i_match);
    
    if isempty(match_script)
        fprintf(out,strcat(clean(line),'\n'));
    elseif length(match_script)==1
        
        offset=regexp(line,matched_scripts(i_match));
        indent=string(line(1:offset-1));
        try
        mid=fopen(strcat(string(match_script{1}),'.m'));

            
        line=fgetl(mid);
        catch err
            disp(err)
        end
        %inline script calls
        while ~isnumeric(line)
            if isempty(line) || ~(first(char(line))=='%')
                fprintf(out,strcat(indent,clean(line),'\n'));
            end
            line=fgetl(mid);
        end
        
         fclose(mid);
    end
    
    line=fgetl(in);
end

fclose(in);
fprintf(out,'end\n');
fclose(out);

end

