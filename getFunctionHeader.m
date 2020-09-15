function fh = getFunctionHeader(f)
str=fileread(f);
str=regexprep(str,"\.\.\.[\r\n]+",""); %remove elipses
str=regexp(str,'function([^\=]+\=[^\=]+\))',"tokens"); %remove function definition
fh=strtrim(str{1}{1});
end

