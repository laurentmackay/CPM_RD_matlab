function init = getInitialConditions(f,chems)


str=fileread(f);

str=regexprep(str,"%[ \%\f\w\=\(\)\+\;\:\.\*\,\]\[\-\/\'\^\?]+",""); %remove comments
str=regexprep(str,"\'[^\'\n\r]+\'",""); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',""); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',""); %remove function definition

code='A-Za-z0-9_ \t\f\*\+\-\/\,\=\[\]\>\<\&\~\:\|\{\}\^\.\(\)';


dbl='[0-9]+(?:(?:[edf\.][\+\-]?)?[0-9])?';
    init_default=regexp(str,['\*\(0\)[ \t\f]*=[^\n\;]*?([' code '])+[ \t\f]*'],'tokens');

% regexp(str,init_default,'tokens');



name='[a-zA-Z_$][a-zA-Z_$0-9\-]*';


str=strjoin(regexp(str,['[^\n]*' name '\(0\)[^\n\;]*(\n|\;|$)'],"match"),newline);


N=length(chems);

if ~isempty(init_default)
    init = repmat(string(init_default{1}{1}),size(chems));
%     init = repmat(str2num(init_default{1}{1}),[1,N]);
    
else
    init=string(nan(1,N));


end


for i=1:N
    chem=chems{i};
    init_str=regexp(str,[chem '\(0\)[ \t\f]*=?[^\r\n\;]*=[ \t\f]*([' code ']+)'],'tokens');
    if ~isempty(init_str)
        init(i)=init_str{1}{1};
    end
end
end
