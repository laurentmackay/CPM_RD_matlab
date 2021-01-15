function fh = getFunctionHeader(f)
str0=fileread(f);
str=regexprep(str0,"\.\.\.[\r\n]+",""); %remove elipses
name='[a-zA-Z_$][a-zA-Z_$0-9]*';
str=regexp(str,['function([^\(]+' name '\([^\)]*\))'],"tokens"); %remove function definition
fh=strtrim(str{1}{1});
end

