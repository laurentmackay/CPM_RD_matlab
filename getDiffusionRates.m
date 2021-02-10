function D = getDiffusionRates(f,chems)


str=fileread(f);

str=regexprep(str,"%[ \%\f\w\=\(\)\+\;\:\.\*\,\]\[\-\/\'\^\?]+",""); %remove comments
str=regexprep(str,"\'[^\'\n\r]+\'",""); %remove hardcoded strings with single quotes
str=regexprep(str,'\"[^\"\n\r]+\"',""); %remove hardcoded strings with double quotes
str=regexprep(str,'function[^\=]+\=[^\=]+\)',""); %remove function definition

name='[a-zA-Z_$][a-zA-Z_$0-9\-]*';


str=strjoin(regexp(str,['[^\n]*D\(' name '\)[^\n\;]*(\n|\;)'],"match"),newline);


N=length(chems);
D=zeros(1,N);
dbl='0-9edf\.';
for i=1:N
    chem=chems{i};
    Dstr=regexp(str,['D\(' chem '\)[ \t\f]*=[^\n\;]*?([' dbl '])+[ \t\f]*[\r\n\;]'],'tokens');
    if ~isempty(Dstr)
        D(i)=str2double(Dstr{1}{1});
    end
end
end

