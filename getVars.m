function vars = getVars(f)
str=fileread(f);

str=regexprep(str,"%[ \f\w\=\(\)\+\;\:\_\.\*\,\]\[\-]+[\n\r\f]",""); %remove comments
str=regexprep(str,"\'[^\'\n\r]+\'",""); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',""); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',""); %remove function definition


vars=regexp(str,'[a-zA-Z_$][a-zA-Z_$0-9]*',"match"); %get variables
vars=unique(vars);


mask = cellfun(@(x) exist(x,'builtin')~=5,vars);
mask = mask & cellfun(@(x) exist(x,'file')~=2,vars);

vars=vars(mask);
end

