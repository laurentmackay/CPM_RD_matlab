function out = mk_fun(script)
    function bool = isScript(f)
        fid = fopen(f);
        
        l1=fgetl(fid);
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




if ~isScript(strcat(script,".m"))
   error('the file provided is not a script') ;
end

clean = @(x) strrep(x,"%","%%");


mfiles=strtrim(string(ls("*.m")));
inds=arrayfun(@(f) isScript(f), mfiles);
scripts=strrep(mfiles(inds),".m","");
matches=strjoin(scripts,"|");
reg=strcat("[\s]?(" ,matches, ")[\s]?");

in = fopen(strcat(script,".m"));

out = fopen(strcat(script,"_fun.m"),'wt');


line=fgets(in);
while ~isnumeric(line)
    
    match=regexp(line,scripts,'match');
    match=match(~cellfun('isempty',match));
    
%     disp(line);
    
    if isempty(match)
        fprintf(out,clean(line));
    elseif length(match)==1
        mid=fopen(strcat(string(match{1}),'.m'));

        line=fgetl(mid);
        while ~isnumeric(line)
            if isempty(line) || ~(first(char(line))=='%')
                fprintf(out,strcat(clean(line),'\n'));
            end
            line=fgetl(mid);
        end
        
         fclose(mid);
    end
    
    line=fgets(in);
end

fclose(in);
fclose(out);

end

