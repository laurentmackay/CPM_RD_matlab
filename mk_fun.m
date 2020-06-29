function out = mk_fun(script)
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

append="_fun";


if ~isScript(strcat(script,".m"))
   error('the file provided is not a script') ;
end

clean = @(x) strrep(x,"%","%%");

listing=dir('*.m');
mfiles= {listing.name}';
inds=cellfun(@(f) isScript(f), mfiles);
scripts=strrep(mfiles(inds),".m","");
matches=strjoin(scripts,"|");
reg=strcat("[\s]?(" ,matches, ")[\s]?");

source = fileread(strcat(script,".m"));
matched_scripts = scripts(~cellfun('isempty',regexp(source,scripts,'match')));

deps=getScriptDep(strcat(script,".m"));
init=getInitialized(strcat(script,".m"));
for s=matched_scripts'
    new_deps=setdiff(getScriptDep(strcat(s,".m")),init);
    deps=[deps, new_deps];
end
deps=unique(deps);

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
add=strcat(script,append,"(");

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



out = fopen(strcat(script,append,".m"),'wt');

fprintf(out,header);

line=fgetl(in);
while ~isnumeric(line)
    
    match=regexp(line,scripts,'match');
    i_match=~cellfun('isempty',match);
    match=match(i_match);
    
    
    if isempty(match)
        fprintf(out,strcat(clean(line),'\n'));
    elseif length(match)==1
        
        offset=regexp(line,scripts(i_match));
        indent=string(line(1:offset-1));
        
        mid=fopen(strcat(string(match{1}),'.m'));
        
        line=fgetl(mid);
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

