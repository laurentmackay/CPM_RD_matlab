function fh = getFunctionHeader(f)
str=fileread(f);
str=regexprep(str,"\.\.\.[\r\n]+",""); %remove elipses
name='[a-zA-Z_$][a-zA-Z_$0-9]*';
str=regexp(str,['function([^\=]+\=[ \t\f]*' name '\([^\)]+\))'],"tokens"); %remove function definition
fh=strtrim(str{1}{1});
end

