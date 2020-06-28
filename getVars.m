function vars = getVars(f)
str=fileread(f);

str=regexprep(str,"%[ \%\f\w\=\(\)\+\;\:\.\*\,\]\[\-\/\']+",""); %remove comments
str=regexprep(str,"\'[^\'\n\r]+\'",""); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',""); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',""); %remove function definition


vars=regexp(str,'[a-zA-Z_$][a-zA-Z_$0-9]*',"match"); %get variables
vars=unique(vars);

vars_bare=regexp(str,'([a-zA-Z_$][a-zA-Z_$0-9]*)[\*\+\-\/\,]',"tokens"); %get variables
vars_bare=cellfun(@(x) x(1),vars_bare);
vars_bare=unique(vars_bare);



mask = cellfun(@(x) exist(x,'builtin')~=5,vars);
mask = mask & cellfun(@(x) exist(x,'file')~=2,vars);

vars=unique([vars(mask) vars_bare]);
end

